"""Seeded bounded local search with independent original-reference validation."""

import hashlib
import math
import random
from collections.abc import Callable
from typing import Annotated, Literal, Self, TypedDict

from Bio.Data import CodonTable
from Bio.Seq import Seq
from pydantic import Field, model_validator

from .execution import ExecutionControl
from .metrics import codon_metrics, gc_fraction
from .models import CDSMode, Construct, FrozenModel, MotifMode, Mutability, Region, Strand
from .motifs import find_motif
from .permissions import effective_permissions
from .profiles import CodonProfile
from .translation import extract_feature
from .validation import (
    CheckStatus,
    Finding,
    ValidationReport,
    check_feasibility,
    validate_candidate,
)

Nonnegative = Annotated[float, Field(ge=0, allow_inf_nan=False)]


class ObjectiveWeights(FrozenModel):
    codon_adaptation: Nonnegative = 1.0
    global_gc: Nonnegative = 0.0
    rare_codons: Nonnegative = 0.0
    custom: Nonnegative = 0.0


class AdditionalLock(Region):
    sequence: str
    reason: str

    @model_validator(mode="after")
    def valid_sequence(self) -> Self:
        if len(self.sequence) != self.end - self.start or set(self.sequence) - set("ACGT"):
            raise ValueError("Additional lock must contain exact uppercase DNA for its interval")
        return self


class SearchConfig(FrozenModel):
    seed: Annotated[int, Field(strict=True)]
    iterations: Annotated[int, Field(ge=0, strict=True)] = 1000
    alternatives: Annotated[int, Field(ge=1, strict=True)] = 1
    algorithm: Literal["hill_climb", "simulated_annealing"] = "simulated_annealing"
    initial_temperature: Nonnegative = 0.05
    n_terminal_lock_nt: Annotated[int, Field(ge=0, strict=True)] = 30
    weights: ObjectiveWeights = Field(default_factory=ObjectiveWeights)
    gc_target: Annotated[float, Field(ge=0, le=1, allow_inf_nan=False)] | None = None
    rare_threshold: Annotated[float, Field(ge=0, le=1, allow_inf_nan=False)] | None = None
    allowed_positions: tuple[Annotated[int, Field(ge=0, strict=True)], ...] | None = None
    additional_locks: tuple[AdditionalLock, ...] = ()

    @model_validator(mode="after")
    def objective_parameters(self) -> Self:
        if self.weights.global_gc and self.gc_target is None:
            raise ValueError("GC objective requires an explicit target")
        if self.weights.rare_codons and self.rare_threshold is None:
            raise ValueError("Rare-codon objective requires an explicit threshold")
        return self


class MutationSite(FrozenModel):
    positions: tuple[int, ...]
    alternatives: tuple[str, ...]
    reverse: bool = False
    feature_id: str | None = None


class NucleotideChange(FrozenModel):
    position: int
    original: str
    designed: str
    reason: str = "constrained local search"


class Design(FrozenModel):
    sequence: str
    objective_loss: float
    components: dict[str, float]
    validation: ValidationReport
    changes: tuple[NucleotideChange, ...]


class SearchResult(FrozenModel):
    status: Literal["success", "infeasible", "budget_exhausted", "invalid_reference"]
    designs: tuple[Design, ...]
    config: SearchConfig
    reference_sha256: str | None
    profile_fingerprint: str | None
    proposals_evaluated: int
    validation: ValidationReport
    warnings: tuple[str, ...]


def _protected_positions(reference: Construct, config: SearchConfig) -> set[int]:
    protected = set()
    if config.allowed_positions is not None:
        allowed = set(config.allowed_positions)
        if any(p >= len(reference.sequence) for p in allowed):
            raise ValueError("Allowed mutation position exceeds the original construct")
        protected.update(set(range(len(reference.sequence))) - allowed)
    for feature in reference.features:
        if feature.is_cds:
            oriented = tuple(range(feature.start, feature.end))
            if feature.strand == Strand.REVERSE:
                oriented = tuple(reversed(oriented))
            count = max(3 if feature.cds_mode == CDSMode.COMPLETE else 0, config.n_terminal_lock_nt)
            protected.update(oriented[feature.frame_offset : feature.frame_offset + count])
    return protected


def mutation_sites(
    reference: Construct, config: SearchConfig, baseline_sequence: str | None = None
) -> tuple[MutationSite, ...]:
    """Propose substitutions; every proposal still needs all-CDS final validation."""
    permissions = effective_permissions(reference)
    baseline_sequence = reference.sequence if baseline_sequence is None else baseline_sequence
    allowed = set(range(len(reference.sequence))) - _protected_positions(reference, config)
    for lock in config.additional_locks:
        if lock.end > len(reference.sequence):
            raise ValueError("Additional lock exceeds construct length")
        allowed.difference_update(range(lock.start, lock.end))
    sites: list[MutationSite] = []
    for feature in reference.features:
        if not feature.is_cds:
            continue
        table = CodonTable.unambiguous_dna_by_id[feature.translation_table]
        dna = extract_feature(baseline_sequence, feature)[feature.frame_offset :]
        oriented = tuple(range(feature.start, feature.end))
        if feature.strand == Strand.REVERSE:
            oriented = tuple(reversed(oriented))
        oriented = oriented[feature.frame_offset :]
        for offset in range(0, len(dna), 3):
            codon = dna[offset : offset + 3]
            if offset == 0 and feature.cds_mode == CDSMode.COMPLETE:
                continue
            if codon in table.stop_codons:
                choices = tuple(sorted(table.stop_codons))
            else:
                aa = table.forward_table[codon]
                choices = tuple(
                    sorted(c for c, amino in table.forward_table.items() if amino == aa)
                )
            positions = oriented[offset : offset + 3]
            choices = tuple(
                choice
                for choice in choices
                if all(
                    base == current or (p in allowed and permissions[p].state != Mutability.LOCKED)
                    for p, base, current in zip(positions, choice, codon, strict=True)
                )
            )
            if len(choices) > 1:
                sites.append(
                    MutationSite(
                        positions=positions,
                        alternatives=choices,
                        reverse=feature.strand == Strand.REVERSE,
                        feature_id=feature.id,
                    )
                )
    for permission in permissions:
        if (
            permission.position in allowed
            and permission.state in {Mutability.FREE_MUTATION, Mutability.DESIGNABLE_RBS}
            and not permission.cds_ids
        ):
            sites.append(MutationSite(positions=(permission.position,), alternatives=tuple("ACGT")))
    return tuple(sites)


def _propose(sequence: str, site: MutationSite, rng: random.Random) -> str:
    current = "".join(sequence[p] for p in site.positions)
    if site.reverse:
        current = str(Seq(current).complement())
    choices = [choice for choice in site.alternatives if choice != current]
    if not choices:
        return sequence
    replacement = rng.choice(choices)
    if site.reverse:
        replacement = str(Seq(replacement).complement())
    bases = list(sequence)
    for position, base in zip(site.positions, replacement, strict=True):
        bases[position] = base
    return "".join(bases)


def validate_search_candidate(
    reference: Construct,
    sequence: str,
    config: SearchConfig,
    *,
    baseline_sequence: str | None = None,
) -> ValidationReport:
    """Fresh original-reference checks AND any explicitly carried design locks."""
    report = validate_candidate(reference, sequence)
    checks = list(report.checks)
    if any(c.code == "REFERENCE_SCHEMA" for c in checks):
        return report
    baseline = reference.sequence if baseline_sequence is None else baseline_sequence
    changed = tuple(
        p
        for p in sorted(_protected_positions(reference, config))
        if sequence[p : p + 1] != baseline[p : p + 1]
    )
    checks.append(
        Finding(
            status=CheckStatus.ERROR if changed else CheckStatus.PASS,
            code="SEARCH_REGION_PROTECTION",
            message="Preserve configured initiation and mutation-region boundaries",
            positions=changed,
        )
    )
    for lock in config.additional_locks:
        same = sequence[lock.start : lock.end] == lock.sequence
        checks.append(
            Finding(
                status=CheckStatus.PASS if same else CheckStatus.ERROR,
                code="ADDITIONAL_LOCK",
                message=lock.reason,
                positions=tuple(range(lock.start, lock.end)),
            )
        )
    return ValidationReport(checks=tuple(checks))


def _motif_penalty(reference: Construct, sequence: str) -> int:
    penalty = 0
    for constraint in reference.constraints:
        matches = find_motif(sequence, constraint, reference.topology)
        if constraint.mode == MotifMode.AVOID:
            penalty += len(matches)
        elif constraint.mode == MotifMode.REQUIRED and not matches:
            penalty += 1
    return penalty


def optimize(
    reference: Construct,
    config: SearchConfig,
    profile: CodonProfile | None = None,
    *,
    initial_sequence: str | None = None,
    custom_objective: Callable[[str], float] | None = None,
    execution: ExecutionControl | None = None,
) -> SearchResult:
    """Repair motifs before scoring valid designs; search exhaustion is not infeasibility.

    Intermediate motif-repair states may violate unresolved motifs, but never locks,
    masks or proteins. Once feasible, every accepted proposal satisfies all constraints.
    No invalid sequence is returned as a Design. Final results are validated afresh.
    """
    config = SearchConfig.model_validate(config.model_dump())
    execution = execution or ExecutionControl()
    execution.checkpoint(0, config.iterations, "Preparing constrained search")
    if profile is not None:
        profile = CodonProfile.model_validate(profile.model_dump())
    if (config.weights.codon_adaptation or config.weights.rare_codons) and profile is None:
        raise ValueError("Codon-adaptation objective requires a sourced codon profile")
    if config.weights.custom and custom_objective is None:
        raise ValueError("Custom objective weight requires a scorer")
    baseline = check_feasibility(reference)
    fingerprint = (
        None
        if any(c.code == "REFERENCE_SCHEMA" for c in baseline.checks)
        else hashlib.sha256(reference.model_dump_json().encode()).hexdigest()
    )

    class ResultContext(TypedDict):
        config: SearchConfig
        reference_sha256: str | None
        profile_fingerprint: str | None

    common: ResultContext = dict(
        config=config,
        reference_sha256=fingerprint,
        profile_fingerprint=profile.fingerprint if profile else None,
    )
    if not baseline.valid:
        return SearchResult(
            **common,
            status="infeasible" if baseline.infeasible else "invalid_reference",
            designs=(),
            proposals_evaluated=0,
            validation=baseline,
            warnings=(),
        )
    reference = Construct.model_validate(reference.model_dump())
    starting_sequence = reference.sequence if initial_sequence is None else initial_sequence
    current = starting_sequence
    report = validate_search_candidate(
        reference, current, config, baseline_sequence=starting_sequence
    )

    def fundamental_error(report: ValidationReport) -> bool:
        return any(
            c.status == CheckStatus.ERROR and c.code not in {"MOTIF_AVOID", "MOTIF_REQUIRED"}
            for c in report.checks
        )

    if fundamental_error(report):
        return SearchResult(
            **common,
            status="invalid_reference",
            designs=(),
            proposals_evaluated=0,
            validation=report,
            warnings=("Initial sequence violates original controls",),
        )
    sites = mutation_sites(reference, config, current)

    def score(sequence: str) -> tuple[float, dict[str, float]]:
        components: dict[str, float] = {}
        if profile and (config.weights.codon_adaptation or config.weights.rare_codons):
            metrics = codon_metrics(
                reference, sequence, profile, rare_threshold=config.rare_threshold
            )
            if config.weights.codon_adaptation:
                components["codon_adaptation"] = (
                    sum(1 - m.adaptation_score for m in metrics)
                    / max(1, len(metrics))
                    * config.weights.codon_adaptation
                )
            if config.weights.rare_codons:
                components["rare_codons"] = (
                    sum(m.rare_fraction or 0 for m in metrics)
                    / max(1, len(metrics))
                    * config.weights.rare_codons
                )
        if config.weights.global_gc:
            assert config.gc_target is not None
            components["global_gc"] = (
                abs(gc_fraction(sequence) - config.gc_target) * config.weights.global_gc
            )
        if config.weights.custom:
            assert custom_objective is not None
            components["custom"] = custom_objective(sequence) * config.weights.custom
        if any(not math.isfinite(value) or value < 0 for value in components.values()):
            raise ValueError("Objective losses must be finite and nonnegative")
        return sum(components.values()), components

    feasible: dict[str, tuple[float, dict[str, float]]] = {}

    def retain(sequence: str, scored: tuple[float, dict[str, float]]) -> None:
        feasible[sequence] = scored
        if len(feasible) > config.alternatives:
            worst = max(feasible, key=lambda seq: (feasible[seq][0], seq))
            del feasible[worst]

    current_penalty = _motif_penalty(reference, current)
    # Scorers only see biologically valid designs, never unresolved repair states.
    current_loss = score(current)[0] if report.valid else 0.0
    if report.valid:
        retain(current, score(current))
    rng = random.Random(config.seed)
    evaluated = 0
    for iteration in range(config.iterations if sites else 0):
        execution.checkpoint(iteration, config.iterations, "Searching valid designs")
        proposal = _propose(current, rng.choice(sites), rng)
        candidate_report = validate_search_candidate(
            reference, proposal, config, baseline_sequence=starting_sequence
        )
        evaluated += 1
        if fundamental_error(candidate_report) or (
            current_penalty == 0 and not candidate_report.valid
        ):
            continue
        penalty = _motif_penalty(reference, proposal)
        loss = 0.0
        if candidate_report.valid:
            loss, components = score(proposal)
            retain(proposal, (loss, components))
        delta = penalty - current_penalty if current_penalty else loss - current_loss
        temperature = config.initial_temperature * (1 - iteration / max(1, config.iterations))
        accept = delta <= 0 or (
            config.algorithm == "simulated_annealing"
            and temperature > 0
            and rng.random() < math.exp(-delta / temperature)
        )
        if accept:
            current, current_penalty, current_loss = proposal, penalty, loss
    designs: list[Design] = []
    execution.checkpoint(config.iterations, config.iterations, "Independently validating designs")
    for sequence in sorted(feasible, key=lambda seq: (feasible[seq][0], seq))[
        : config.alternatives
    ]:
        final_report = validate_search_candidate(
            reference, sequence, config, baseline_sequence=starting_sequence
        )
        final_report.require_valid()
        final_loss, components = score(sequence)
        designs.append(
            Design(
                sequence=sequence,
                objective_loss=final_loss,
                components=components,
                validation=final_report,
                changes=tuple(
                    NucleotideChange(position=i, original=a, designed=b)
                    for i, (a, b) in enumerate(zip(reference.sequence, sequence, strict=True))
                    if a != b
                ),
            )
        )
    warnings = ["Objective scores are engineering preferences, not absolute expression predictions"]
    if len(designs) < config.alternatives:
        warnings.append(
            "Search budget did not produce the requested number of distinct valid designs"
        )
    if profile and profile.reference_kind != "high_expression":
        warnings.append("Genome/custom codon adaptation is not high-expression-reference CAI")
    return SearchResult(
        **common,
        status="success" if len(designs) == config.alternatives else "budget_exhausted",
        designs=tuple(designs),
        proposals_evaluated=evaluated,
        validation=designs[0].validation
        if designs
        else validate_search_candidate(
            reference, current, config, baseline_sequence=starting_sequence
        ),
        warnings=tuple(warnings),
    )
