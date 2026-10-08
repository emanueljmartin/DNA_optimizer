"""Independent reference-based validation and conservative feasibility proofs."""

from enum import StrEnum

from Bio.Data import IUPACData
from Bio.Seq import Seq
from pydantic import ValidationError
from pydantic_core import PydanticSerializationError

from .errors import InvalidDesignError, TranslationError
from .models import (
    Construct,
    FrozenModel,
    MotifConstraint,
    MotifMode,
    MotifStrand,
    Mutability,
    Topology,
)
from .motifs import find_motif
from .permissions import effective_permissions
from .translation import extract_feature, translate_feature


class CheckStatus(StrEnum):
    PASS = "PASS"
    WARNING = "WARNING"
    ERROR = "ERROR"


class Finding(FrozenModel):
    status: CheckStatus
    code: str
    message: str
    feature_id: str | None = None
    constraint_id: str | None = None
    positions: tuple[int, ...] = ()
    proves_infeasible: bool = False


class ValidationReport(FrozenModel):
    checks: tuple[Finding, ...]

    @property
    def valid(self) -> bool:
        return not any(check.status == CheckStatus.ERROR for check in self.checks)

    @property
    def infeasible(self) -> bool:
        return any(check.proves_infeasible for check in self.checks)

    def require_valid(self) -> None:
        if not self.valid:
            raise InvalidDesignError(
                "; ".join(c.message for c in self.checks if c.status == CheckStatus.ERROR)
            )


def _oriented_motifs(constraint: MotifConstraint) -> tuple[str, ...]:
    motifs: list[str] = []
    if constraint.strand != MotifStrand.REVERSE:
        motifs.append(constraint.motif)
    if constraint.strand != MotifStrand.FORWARD:
        motifs.append(str(Seq(constraint.motif).reverse_complement()))
    return tuple(motifs)


def _required_is_forbidden(required: MotifConstraint, avoid: MotifConstraint, length: int) -> bool:
    """Prove subset recognition and containing scopes, without enumerating designs."""
    if len(required.motif) != len(avoid.motif):
        return False
    rstart, rend = (required.region.start, required.region.end) if required.region else (0, length)
    astart, aend = (avoid.region.start, avoid.region.end) if avoid.region else (0, length)
    if not (astart <= rstart and rend <= aend):
        return False
    if avoid.position is not None and avoid.position != required.position:
        return False
    values = IUPACData.ambiguous_dna_values
    return all(
        any(
            all(set(values[r]) <= set(values[a]) for r, a in zip(rmotif, amotif, strict=True))
            for amotif in _oriented_motifs(avoid)
        )
        for rmotif in _oriented_motifs(required)
    )


def _feasibility_checks(reference: Construct) -> list[Finding]:
    checks: list[Finding] = []
    permissions = effective_permissions(reference)
    avoids = [c for c in reference.constraints if c.mode == MotifMode.AVOID]
    for constraint in reference.constraints:
        if constraint.mode == MotifMode.PRESERVE and not find_motif(
            reference.sequence, constraint, reference.topology
        ):
            checks.append(
                Finding(
                    status=CheckStatus.ERROR,
                    code="PRESERVE_ABSENT",
                    message=f"PRESERVE {constraint.id} has no reference occurrence",
                    constraint_id=constraint.id,
                    proves_infeasible=True,
                )
            )
        if constraint.mode == MotifMode.REQUIRED:
            # Overapproximate mutable bases as arbitrary A/C/G/T. If no possible
            # scoped window can match even under that relaxation, impossibility
            # is proven without guessing synonymous feasibility.
            length = len(reference.sequence)
            starts = (
                range(length)
                if reference.topology == Topology.CIRCULAR
                else range(length - len(constraint.motif) + 1)
            )
            possible = False
            for start in starts:
                if constraint.position is not None and start != constraint.position:
                    continue
                positions = tuple((start + i) % length for i in range(len(constraint.motif)))
                if constraint.region and not all(
                    constraint.region.start <= p < constraint.region.end for p in positions
                ):
                    continue
                if any(
                    all(
                        permissions[p].state != Mutability.LOCKED
                        or reference.sequence[p] in IUPACData.ambiguous_dna_values[symbol]
                        for p, symbol in zip(positions, motif, strict=True)
                    )
                    for motif in _oriented_motifs(constraint)
                ):
                    possible = True
                    break
            if not possible:
                checks.append(
                    Finding(
                        status=CheckStatus.ERROR,
                        code="REQUIRED_IMPOSSIBLE",
                        message=f"REQUIRED {constraint.id}: no compatible location "
                        "under locks/scope",
                        constraint_id=constraint.id,
                        proves_infeasible=True,
                    )
                )
            for avoid in avoids:
                if _required_is_forbidden(constraint, avoid, len(reference.sequence)):
                    checks.append(
                        Finding(
                            status=CheckStatus.ERROR,
                            code="MOTIF_CONFLICT",
                            message=f"REQUIRED {constraint.id} is forbidden by AVOID {avoid.id}",
                            constraint_id=constraint.id,
                            proves_infeasible=True,
                        )
                    )
    for constraint in avoids:
        for occurrence in find_motif(reference.sequence, constraint, reference.topology):
            positions = occurrence.positions(len(reference.sequence))
            if all(permissions[p].state == Mutability.LOCKED for p in positions):
                checks.append(
                    Finding(
                        status=CheckStatus.ERROR,
                        code="IMMUTABLE_FORBIDDEN_MOTIF",
                        message=f"AVOID {constraint.id} is trapped in immutable DNA",
                        constraint_id=constraint.id,
                        positions=positions,
                        proves_infeasible=True,
                    )
                )
    return checks


def _reference_checks(reference: Construct) -> list[Finding]:
    checks: list[Finding] = []
    for feature in reference.features:
        if feature.is_cds:
            try:
                translate_feature(reference.sequence, feature)
            except TranslationError as exc:
                checks.append(
                    Finding(
                        status=CheckStatus.ERROR,
                        code="REFERENCE_CDS_INVALID",
                        message=str(exc),
                        feature_id=feature.id,
                    )
                )
        if feature.rbs_mask is not None:
            dna = extract_feature(reference.sequence, feature)
            if any(m != "N" and m != base for m, base in zip(feature.rbs_mask, dna, strict=True)):
                checks.append(
                    Finding(
                        status=CheckStatus.ERROR,
                        code="REFERENCE_RBS_MASK_INVALID",
                        message=f"Original DNA disagrees with fixed mask for {feature.id}",
                        feature_id=feature.id,
                    )
                )
    return checks


def check_feasibility(reference: Construct) -> ValidationReport:
    """Valid means no proven contradiction; it does NOT prove a solution exists."""
    try:
        reference = Construct.model_validate(reference.model_dump(warnings=False))
    except (ValidationError, PydanticSerializationError) as exc:
        return ValidationReport(
            checks=(Finding(status=CheckStatus.ERROR, code="REFERENCE_SCHEMA", message=str(exc)),)
        )
    checks = _reference_checks(reference) + _feasibility_checks(reference)
    if not checks:
        checks.append(
            Finding(
                status=CheckStatus.WARNING,
                code="FEASIBILITY_UNKNOWN",
                message="No contradiction proven; general synonymous feasibility was not solved",
            )
        )
    return ValidationReport(checks=tuple(checks))


def validate_candidate(reference: Construct, candidate: Construct | str) -> ValidationReport:
    """Recompute ALL invariants using original controls, never candidate controls.

    An invalid original cannot be used as a biological preservation baseline. A
    raw candidate string is supported so illegal DNA can produce a structured report.
    """
    checks: list[Finding] = []
    try:
        reference = Construct.model_validate(reference.model_dump(warnings=False))
    except (ValidationError, PydanticSerializationError) as exc:
        return ValidationReport(
            checks=(Finding(status=CheckStatus.ERROR, code="REFERENCE_SCHEMA", message=str(exc)),)
        )
    if isinstance(candidate, Construct):
        # Round-trip validation also rejects bypasses through model_copy/model_construct.
        try:
            candidate = Construct.model_validate(candidate.model_dump(warnings=False))
        except (ValidationError, PydanticSerializationError) as exc:
            return ValidationReport(
                checks=(
                    Finding(status=CheckStatus.ERROR, code="CANDIDATE_SCHEMA", message=str(exc)),
                )
            )
        controls = (
            "features",
            "constraints",
            "topology",
            "source_sequence",
            "host_profile",
            "synthesis_profile",
            "rbs_settings",
            "fragmentation_settings",
            "model_versions",
        )
        changed = [
            name for name in controls if getattr(candidate, name) != getattr(reference, name)
        ]
        checks.append(
            Finding(
                status=CheckStatus.ERROR if changed else CheckStatus.PASS,
                code="REFERENCE_CONTROLS",
                message=f"Candidate changed reference controls: {', '.join(changed)}"
                if changed
                else "Candidate uses the original biological controls",
            )
        )
        sequence = candidate.sequence
    else:
        sequence = candidate
    if not isinstance(sequence, str) or not sequence or set(sequence) - set("ACGT"):
        checks.append(
            Finding(
                status=CheckStatus.ERROR,
                code="DNA_ALPHABET",
                message="Candidate must be nonempty uppercase concrete A/C/G/T DNA",
            )
        )
        return ValidationReport(checks=tuple(checks))
    checks.append(
        Finding(status=CheckStatus.PASS, code="DNA_ALPHABET", message="DNA alphabet valid")
    )
    if len(sequence) != len(reference.sequence):
        checks.append(
            Finding(
                status=CheckStatus.ERROR,
                code="SEQUENCE_LENGTH",
                message="Only substitutions are supported; candidate length differs from original",
            )
        )
        return ValidationReport(checks=tuple(checks))
    checks.append(
        Finding(
            status=CheckStatus.PASS, code="SEQUENCE_LENGTH", message="Construct length preserved"
        )
    )
    baseline_checks = _reference_checks(reference)
    checks.extend(baseline_checks)
    checks.extend(_feasibility_checks(reference))
    changed_locks = tuple(
        p.position
        for p in effective_permissions(reference)
        if p.state == Mutability.LOCKED and reference.sequence[p.position] != sequence[p.position]
    )
    checks.append(
        Finding(
            status=CheckStatus.ERROR if changed_locks else CheckStatus.PASS,
            code="LOCKED_INTEGRITY",
            message="Immutable original DNA changed"
            if changed_locks
            else "Locked DNA byte-identical",
            positions=changed_locks,
        )
    )
    for feature in reference.features:
        if feature.is_cds:
            try:
                original = translate_feature(reference.sequence, feature)
                translated = translate_feature(sequence, feature)
                same = original == translated
                checks.append(
                    Finding(
                        status=CheckStatus.PASS if same else CheckStatus.ERROR,
                        code="PROTEIN_INTEGRITY",
                        message=f"CDS {feature.id}: protein {'preserved' if same else 'changed'}",
                        feature_id=feature.id,
                    )
                )
            except TranslationError as exc:
                checks.append(
                    Finding(
                        status=CheckStatus.ERROR,
                        code="CDS_TRANSLATION",
                        message=str(exc),
                        feature_id=feature.id,
                    )
                )
        if feature.rbs_mask is not None:
            dna = extract_feature(sequence, feature)
            same = all(m == "N" or m == base for m, base in zip(feature.rbs_mask, dna, strict=True))
            checks.append(
                Finding(
                    status=CheckStatus.PASS if same else CheckStatus.ERROR,
                    code="RBS_MASK_INTEGRITY",
                    message=f"Fixed RBS mask check for {feature.id}",
                    feature_id=feature.id,
                )
            )
    for constraint in reference.constraints:
        occurrences = find_motif(sequence, constraint, reference.topology)
        if constraint.mode == MotifMode.AVOID:
            satisfied = not occurrences
        elif constraint.mode == MotifMode.REQUIRED:
            satisfied = bool(occurrences)
        else:
            original_occurrences = find_motif(reference.sequence, constraint, reference.topology)
            satisfied = bool(original_occurrences) and all(
                occurrence in occurrences
                and all(
                    sequence[p] == reference.sequence[p]
                    for p in occurrence.positions(len(sequence))
                )
                for occurrence in original_occurrences
            )
        checks.append(
            Finding(
                status=CheckStatus.PASS if satisfied else CheckStatus.ERROR,
                code=f"MOTIF_{constraint.mode.value}",
                message=f"{constraint.mode.value} motif {constraint.id}: "
                f"{'satisfied' if satisfied else 'violated'}",
                constraint_id=constraint.id,
            )
        )
    checks.extend(
        Finding(status=CheckStatus.WARNING, code="IMPORT_WARNING", message=message)
        for message in reference.import_warnings
    )
    return ValidationReport(checks=tuple(checks))
