"""Transparent synthesis-risk heuristics with explicit, versioned policies."""

from collections import defaultdict
from itertools import groupby
from typing import Annotated, Self

from Bio.Seq import Seq
from pydantic import Field, model_validator

from .errors import InvalidDesignError
from .execution import ExecutionControl
from .metrics import gc_fraction, sequence_entropy
from .models import Construct, FrozenModel, Identifier, Topology
from .optimization import SearchConfig, SearchResult, optimize
from .profiles import CodonProfile
from .validation import CheckStatus, ValidationReport, validate_candidate

PositiveInt = Annotated[int, Field(ge=1, strict=True)]


class FractionRange(FrozenModel):
    minimum: Annotated[float, Field(ge=0, le=1, allow_inf_nan=False)]
    maximum: Annotated[float, Field(ge=0, le=1, allow_inf_nan=False)]

    @model_validator(mode="after")
    def ordered(self) -> Self:
        if self.minimum > self.maximum:
            raise ValueError("Minimum must not exceed maximum")
        return self


class SynthesisProfile(FrozenModel):
    id: Identifier = "generic"
    version: Identifier = "explicit-policy-v1"
    source: Identifier = "User-selected engineering policy; PRODUCT_SPEC.md O"
    accessed_date: Annotated[str, Field(pattern=r"^\d{4}-\d{2}-\d{2}$")] = "2026-10-08"
    approximation: Identifier = "Transparent heuristics; no prediction of vendor acceptance"
    global_gc_range: FractionRange | None = None
    gc_window: PositiveInt | None = None
    local_gc_range: FractionRange | None = None
    maximum_local_gc_spread: Annotated[float, Field(ge=0, le=1, allow_inf_nan=False)] | None = None
    homopolymer_limit: PositiveInt | None = None
    repeat_min_length: Annotated[int, Field(ge=2, strict=True)] | None = None
    inverted_repeat_min_length: Annotated[int, Field(ge=2, strict=True)] | None = None
    complexity_window: PositiveInt | None = None
    minimum_entropy: Annotated[float, Field(ge=0, le=2, allow_inf_nan=False)] | None = None

    @model_validator(mode="after")
    def complete_policies(self) -> Self:
        if (self.gc_window is None) != (
            self.local_gc_range is None and self.maximum_local_gc_spread is None
        ):
            raise ValueError("Local GC assessment requires a window and explicit range or spread")
        if (self.complexity_window is None) != (self.minimum_entropy is None):
            raise ValueError(
                "Complexity assessment requires both window and explicit entropy limit"
            )
        return self


class SynthesisFinding(FrozenModel):
    status: CheckStatus
    metric: str
    message: str
    observed: float | None = None
    positions: tuple[int, ...] = ()
    penalty: Annotated[float, Field(ge=0, allow_inf_nan=False)] = 0


class SynthesisReport(FrozenModel):
    profile: SynthesisProfile
    findings: tuple[SynthesisFinding, ...]
    biological_validation: ValidationReport
    global_gc: Annotated[float, Field(ge=0, le=1)]

    @property
    def risk_loss(self) -> float:
        return sum(f.penalty for f in self.findings)

    def require_policy(self) -> None:
        self.biological_validation.require_valid()
        if self.risk_loss:
            raise InvalidDesignError("Design violates the explicitly selected synthesis policy")
        if any(f.observed is None for f in self.findings):
            raise InvalidDesignError(
                "Unconfigured synthesis criteria cannot be certified as passing"
            )


def _windows(sequence: str, width: int, topology: Topology) -> tuple[tuple[int, str], ...]:
    if width > len(sequence):
        return ()
    circular = topology == Topology.CIRCULAR
    scan = sequence + sequence[: width - 1] if circular else sequence
    return tuple(
        (start, scan[start : start + width])
        for start in range(len(sequence) if circular else len(sequence) - width + 1)
    )


def _run_positions(sequence: str, topology: Topology) -> list[tuple[int, ...]]:
    runs = [
        tuple(i for i, _ in group)
        for _, group in groupby(enumerate(sequence), key=lambda pair: pair[1])
    ]
    if topology == Topology.CIRCULAR and len(runs) > 1 and sequence[0] == sequence[-1]:
        runs = [runs[-1] + runs[0], *runs[1:-1]]
    return runs


def scan_synthesis(
    reference: Construct, sequence: str, profile: SynthesisProfile
) -> SynthesisReport:
    """Assess selected criteria; unspecified thresholds remain visibly unassessed."""
    profile = SynthesisProfile.model_validate(profile.model_dump())
    validation = validate_candidate(reference, sequence)
    if not sequence or set(sequence) - set("ACGT") or len(sequence) != len(reference.sequence):
        raise ValueError("Synthesis analysis requires equal-length concrete DNA")
    findings: list[SynthesisFinding] = []

    def emit(metric: str, value: float, penalty: float, positions: tuple[int, ...] = ()) -> None:
        findings.append(
            SynthesisFinding(
                status=CheckStatus.WARNING if penalty else CheckStatus.PASS,
                metric=metric,
                observed=value,
                penalty=penalty,
                positions=positions,
                message=f"{metric}: "
                f"{'outside selected policy' if penalty else 'within selected policy'}",
            )
        )

    def unassessed(metric: str) -> None:
        findings.append(
            SynthesisFinding(
                status=CheckStatus.WARNING,
                metric=metric,
                message="No compatible window or explicit threshold configured; "
                "no vendor acceptance claim",
            )
        )

    gc = gc_fraction(sequence)
    if profile.global_gc_range:
        bounds = profile.global_gc_range
        emit("global_gc", gc, max(0, bounds.minimum - gc, gc - bounds.maximum))
    else:
        unassessed("global_gc")
    if profile.gc_window:
        windows = _windows(sequence, profile.gc_window, reference.topology)
        local_bounds = profile.local_gc_range
        if not windows:
            unassessed("local_gc")
        for start, window in windows:
            local = gc_fraction(window)
            if local_bounds:
                emit(
                    "local_gc",
                    local,
                    max(0, local_bounds.minimum - local, local - local_bounds.maximum),
                    tuple((start + i) % len(sequence) for i in range(len(window))),
                )
        if windows and profile.maximum_local_gc_spread is not None:
            fractions = [gc_fraction(window) for _, window in windows]
            spread = max(fractions) - min(fractions)
            emit("local_gc_spread", spread, max(0, spread - profile.maximum_local_gc_spread))
    else:
        unassessed("local_gc")
    if profile.homopolymer_limit:
        for positions in _run_positions(sequence, reference.topology):
            emit(
                "homopolymer",
                len(positions),
                max(0, len(positions) - profile.homopolymer_limit),
                positions,
            )
    else:
        unassessed("homopolymer")
    if profile.repeat_min_length:
        words: dict[str, list[int]] = defaultdict(list)
        for start, word in _windows(sequence, profile.repeat_min_length, reference.topology):
            words[word].append(start)
        repeats = [(word, starts) for word, starts in words.items() if len(starts) > 1]
        if not words:
            unassessed("direct_repeat")
        elif not repeats:
            emit("direct_repeat", 0, 0)
        for word, starts in repeats:
            emit(
                "direct_repeat",
                len(starts),
                len(starts) - 1,
                tuple(sorted({(s + i) % len(sequence) for s in starts for i in range(len(word))})),
            )
    else:
        unassessed("direct_repeat")
    if profile.inverted_repeat_min_length:
        width = profile.inverted_repeat_min_length
        words = defaultdict(list)
        for start, word in _windows(sequence, width, reference.topology):
            words[word].append(start)
        seen: set[str] = set()
        pairs = 0
        for word, starts in words.items():
            reverse = str(Seq(word).reverse_complement())
            if word in seen:
                continue
            seen.update((word, reverse))
            matched = False
            for a in starts:
                first = {(a + i) % len(sequence) for i in range(width)}
                for b in words.get(reverse, ()):
                    second = {(b + i) % len(sequence) for i in range(width)}
                    if not first.intersection(second):
                        emit("inverted_repeat", width, 1, tuple(sorted(first | second)))
                        pairs += 1
                        matched = True
                        break
                if matched:
                    break
        if not words:
            unassessed("inverted_repeat")
        elif not pairs:
            emit("inverted_repeat", 0, 0)
    else:
        unassessed("inverted_repeat")
    if profile.complexity_window and profile.minimum_entropy is not None:
        windows = _windows(sequence, profile.complexity_window, reference.topology)
        if not windows:
            unassessed("complexity")
        for start, window in windows:
            entropy = sequence_entropy(window)
            emit(
                "complexity",
                entropy,
                max(0, profile.minimum_entropy - entropy),
                tuple((start + i) % len(sequence) for i in range(len(window))),
            )
    else:
        unassessed("complexity")
    return SynthesisReport(
        profile=profile, findings=tuple(findings), biological_validation=validation, global_gc=gc
    )


def load_synthesis_profile(identifier: str) -> SynthesisProfile:
    """Public guidance snapshots, deliberately incomplete approximations."""
    if identifier == "generic":
        return SynthesisProfile()
    if identifier == "genscript-public":
        return SynthesisProfile(
            id=identifier,
            version="gene-synthesis-guidance-2026-10-08",
            source="https://www.genscript.com/gene_synthesis.html",
            approximation="Public GC and repeat warnings; other thresholds unpublished; "
            "not GenScript's proprietary score or an acceptance prediction",
            global_gc_range=FractionRange(minimum=0.2, maximum=0.8),
            repeat_min_length=20,
        )
    if identifier == "twist-public":
        return SynthesisProfile(
            id=identifier,
            version="DOC-001081-REV4-2026-10-08",
            source="https://www.twistbioscience.com/content/dam/twistbioscience/resources/"
            "2023-06/DOC-001081_TechNote_TwistTipVectorDesign-REV4-singles.pdf",
            approximation="Clonal gene design guidance: GC, repeats and local-GC spread only; "
            "repeat Tm and proprietary combined score are not modeled; no acceptance prediction",
            global_gc_range=FractionRange(minimum=0.25, maximum=0.65),
            repeat_min_length=20,
            gc_window=50,
            maximum_local_gc_spread=0.5,
        )
    raise ValueError(f"Unknown synthesis profile: {identifier}; create an explicit custom policy")


class SynthesisRepairResult(FrozenModel):
    search: SearchResult
    before: SynthesisReport
    after: tuple[SynthesisReport, ...]


def repair_synthesis(
    reference: Construct,
    profile: SynthesisProfile,
    config: SearchConfig,
    *,
    initial_sequence: str | None = None,
    codon_profile: CodonProfile | None = None,
    execution: ExecutionControl | None = None,
) -> SynthesisRepairResult:
    """Repair advisory risks; residual warnings stay visible and cannot be certified."""
    initial = reference.sequence if initial_sequence is None else initial_sequence
    before = scan_synthesis(reference, initial, profile)
    weights = config.weights.model_dump()
    weights["custom"] = config.weights.custom or 1.0
    config = SearchConfig.model_validate({**config.model_dump(), "weights": weights})
    if codon_profile is None and config.weights.codon_adaptation:
        raise ValueError("Supply a codon profile or explicitly disable codon-adaptation weight")
    search = optimize(
        reference,
        config,
        codon_profile,
        initial_sequence=initial,
        custom_objective=lambda sequence: scan_synthesis(reference, sequence, profile).risk_loss,
        execution=execution,
    )
    return SynthesisRepairResult(
        search=search,
        before=before,
        after=tuple(scan_synthesis(reference, d.sequence, profile) for d in search.designs),
    )
