"""Deterministic Gibson tiling with exact independent reconstruction and DNA/DNA Tm."""

import hashlib
import math
from collections import Counter, deque
from typing import Annotated, Literal, Self, TypedDict

import Bio
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from pydantic import Field, model_validator

from .execution import ExecutionControl
from .metrics import gc_fraction, sequence_entropy
from .models import Construct, FrozenModel, Mutability, Topology
from .permissions import effective_permissions
from .synthesis import FractionRange, SynthesisProfile, SynthesisReport, scan_synthesis
from .validation import CheckStatus, Finding, ValidationReport, validate_candidate

Positive = Annotated[float, Field(gt=0, allow_inf_nan=False)]
Nonnegative = Annotated[float, Field(ge=0, allow_inf_nan=False)]


class TmConditions(FrozenModel):
    """Biopython Tm_NN defaults; explicit solution assumptions, not Gibson mix claims."""

    strand1_nm: Positive = 25
    strand2_nm: Positive = 25
    sodium_mm: Positive = 50
    potassium_mm: Nonnegative = 0
    tris_mm: Nonnegative = 0
    magnesium_mm: Nonnegative = 0
    dntps_mm: Nonnegative = 0
    salt_correction: Literal[5, 7] = 5


class GibsonSettings(FrozenModel):
    maximum_ordered_length: Annotated[int, Field(ge=2, strict=True)]
    overlap_length: Annotated[int, Field(ge=2, strict=True)] = 30
    minimum_ordered_length: Annotated[int, Field(ge=1, strict=True)] = 60
    search_window: Annotated[int, Field(ge=1, strict=True)] = 100
    maximum_fragments: Annotated[int, Field(ge=1, strict=True)] = 128
    overlap_gc_range: FractionRange | None = None
    preferred_tm_c: Annotated[float, Field(allow_inf_nan=False)] | None = None
    require_unique_overlaps: bool = True
    avoid_protected_overlaps: bool = False
    conditions: TmConditions = Field(default_factory=TmConditions)
    synthesis_profile: SynthesisProfile | None = None

    @model_validator(mode="after")
    def limits(self) -> Self:
        if self.maximum_ordered_length <= self.overlap_length:
            raise ValueError("Maximum ordered length must exceed overlap length")
        if self.minimum_ordered_length > self.maximum_ordered_length:
            raise ValueError("Minimum ordered length exceeds maximum")
        if self.minimum_ordered_length <= self.overlap_length:
            raise ValueError("Minimum ordered length must exceed overlap for progressing tiles")
        return self


class OverlapMetrics(FrozenModel):
    start: int
    sequence: str
    gc: float
    tm_c: float
    entropy_bits: float
    forward_occurrences: int
    reverse_occurrences: int
    protected_bases: int
    unique: bool
    penalty_components: dict[str, float]


class GibsonFragment(FrozenModel):
    number: int
    start: int
    end: int
    sequence: str
    unique_length: int
    ordered_length: int
    gc: float
    left_overlap: str
    right_overlap: str
    synthesis: SynthesisReport | None = None


class GibsonPlan(FrozenModel):
    status: Literal["success", "infeasible", "budget_exhausted"]
    settings: GibsonSettings
    topology: Topology
    source_sha256: str
    minimum_geometric_fragments: int
    fragments: tuple[GibsonFragment, ...]
    overlaps: tuple[OverlapMetrics, ...]
    validation: ValidationReport
    model_versions: dict[str, str]
    warnings: tuple[str, ...]


def minimum_fragment_count(
    length: int, maximum: int, overlap: int, *, circular: bool = False
) -> int:
    if any(type(value) is not int for value in (length, maximum, overlap)):
        raise ValueError("Fragment geometry requires integer lengths")
    if length < 1 or overlap < 1 or maximum <= overlap:
        raise ValueError("Positive construct/overlap and maximum > overlap are required")
    if circular:
        return max(1, math.ceil(length / (maximum - overlap)))
    return max(1, math.ceil((length - overlap) / (maximum - overlap)))


def overlap_tm(sequence: str, conditions: TmConditions) -> float:
    gc_fraction(sequence)
    if len(sequence) < 2:
        raise ValueError("Nearest-neighbor Tm requires at least two bases")
    conditions = TmConditions.model_validate(conditions.model_dump())
    return float(
        mt.Tm_NN(
            sequence,
            check=False,
            strict=True,
            nn_table=mt.DNA_NN3,
            dnac1=conditions.strand1_nm,
            dnac2=conditions.strand2_nm,
            selfcomp=sequence == str(Seq(sequence).reverse_complement()),
            Na=conditions.sodium_mm,
            K=conditions.potassium_mm,
            Tris=conditions.tris_mm,
            Mg=conditions.magnesium_mm,
            dNTPs=conditions.dntps_mm,
            saltcorr=conditions.salt_correction,
        )
    )


def reconstruct(plan: GibsonPlan) -> str:
    if not plan.fragments:
        raise ValueError("Plan contains no fragments")
    overlap = plan.settings.overlap_length
    dna = plan.fragments[0].sequence
    for fragment in plan.fragments[1:]:
        if dna[-overlap:] != fragment.sequence[:overlap]:
            raise ValueError("Adjacent fragments do not share the declared overlap")
        dna += fragment.sequence[overlap:]
    if plan.topology == Topology.CIRCULAR:
        if dna[-overlap:] != dna[:overlap]:
            raise ValueError("Circular closing overlap does not match the origin")
        dna = dna[:-overlap]
    return dna


def validate_gibson(reference: Construct, sequence: str, plan: GibsonPlan) -> ValidationReport:
    """Fresh biological checks plus independent physical tiling/order verification."""
    checks = list(validate_candidate(reference, sequence).checks)

    def add(code: str, passed: bool, message: str) -> None:
        checks.append(
            Finding(
                status=CheckStatus.PASS if passed else CheckStatus.ERROR, code=code, message=message
            )
        )

    add(
        "GIBSON_SOURCE",
        plan.source_sha256 == hashlib.sha256(sequence.encode()).hexdigest()
        and plan.topology == reference.topology,
        "Fragment plan belongs to the final sequence/topology",
    )
    add("GIBSON_STATUS", plan.status == "success", "Fragment plan completed successfully")
    maximum, overlap = plan.settings.maximum_ordered_length, plan.settings.overlap_length
    add(
        "GIBSON_COUNT",
        len(plan.fragments) <= plan.settings.maximum_fragments
        and plan.minimum_geometric_fragments
        == minimum_fragment_count(
            len(sequence), maximum, overlap, circular=reference.topology == Topology.CIRCULAR
        ),
        "Fragment-count metadata agrees with geometry and selected budget",
    )
    scan = sequence + (sequence[:overlap] if reference.topology == Topology.CIRCULAR else "")
    try:
        rebuilt = reconstruct(plan)
    except (ValueError, IndexError):
        rebuilt = ""
    add("GIBSON_RECONSTRUCTION", rebuilt == sequence, "Ordered fragments reconstruct DNA exactly")
    expected_start = 0
    valid_tiles = bool(plan.fragments)
    circular = reference.topology == Topology.CIRCULAR
    for number, fragment in enumerate(plan.fragments, start=1):
        expected_unique = len(fragment.sequence) - (overlap if number > 1 else 0)
        if reference.topology == Topology.CIRCULAR and number == len(plan.fragments):
            expected_unique -= overlap
        valid_tiles &= (
            fragment.number == number
            and fragment.start == expected_start
            and 0 <= fragment.start < fragment.end <= len(scan)
            and fragment.sequence == scan[fragment.start : fragment.end]
            and fragment.ordered_length == len(fragment.sequence) <= maximum
            and fragment.unique_length == expected_unique
            and fragment.left_overlap
            == (fragment.sequence[:overlap] if number > 1 or circular else "")
            and fragment.right_overlap
            == (fragment.sequence[-overlap:] if number < len(plan.fragments) or circular else "")
            and math.isclose(fragment.gc, gc_fraction(fragment.sequence), abs_tol=1e-12)
            and (
                len(plan.fragments) == 1
                or len(fragment.sequence) >= plan.settings.minimum_ordered_length
            )
        )
        expected_start = fragment.end - overlap
    valid_tiles &= bool(plan.fragments) and plan.fragments[-1].end == len(scan)
    add(
        "GIBSON_GEOMETRY",
        valid_tiles,
        "Fragment order, coordinates and ordered lengths are correct",
    )
    actual_starts = tuple(f.start for f in plan.fragments[1:])
    if reference.topology == Topology.CIRCULAR:
        actual_starts = (0,) + actual_starts
    add(
        "GIBSON_OVERLAPS",
        tuple(o.start for o in plan.overlaps) == actual_starts
        and all(o.sequence == scan[o.start : o.start + overlap] for o in plan.overlaps),
        "Overlap metadata agrees with physical boundaries",
    )
    occurrences = Counter(
        scan[i : i + overlap]
        for i in range(
            len(sequence)
            if reference.topology == Topology.CIRCULAR
            else max(0, len(sequence) - overlap + 1)
        )
    )
    unique = all(
        occurrences[o.sequence] == 1 and occurrences[str(Seq(o.sequence).reverse_complement())] == 0
        for o in plan.overlaps
    )
    locked = [
        permission.state == Mutability.LOCKED for permission in effective_permissions(reference)
    ]
    add(
        "GIBSON_PROTECTED_OVERLAPS",
        not plan.settings.avoid_protected_overlaps
        or all(
            not any(locked[(o.start + i) % len(sequence)] for i in range(overlap))
            for o in plan.overlaps
        ),
        "Selected protected-overlap policy is independently checked",
    )
    add(
        "GIBSON_METRICS",
        all(
            math.isclose(o.gc, gc_fraction(o.sequence), abs_tol=1e-12)
            and math.isclose(o.entropy_bits, sequence_entropy(o.sequence), abs_tol=1e-12)
            and o.protected_bases
            == sum(locked[(o.start + i) % len(sequence)] for i in range(overlap))
            and math.isclose(o.tm_c, overlap_tm(o.sequence, plan.settings.conditions), abs_tol=1e-9)
            and o.forward_occurrences == occurrences[o.sequence]
            and o.reverse_occurrences == occurrences[str(Seq(o.sequence).reverse_complement())]
            and o.unique
            == (
                occurrences[o.sequence] == 1
                and occurrences[str(Seq(o.sequence).reverse_complement())] == 0
            )
            for o in plan.overlaps
        ),
        "Overlap composition, Tm and occurrences are freshly verified",
    )
    add(
        "GIBSON_ORDER_UNIQUENESS",
        unique,
        "Full overlaps identify one position and orientation each",
    )
    return ValidationReport(checks=tuple(checks))


def plan_gibson(
    reference: Construct,
    sequence: str,
    settings: GibsonSettings,
    *,
    execution: ExecutionControl | None = None,
) -> GibsonPlan:
    execution = execution or ExecutionControl()
    execution.checkpoint(0, 1, "Preparing assembly geometry")
    validate_candidate(reference, sequence).require_valid()
    settings = GibsonSettings.model_validate(settings.model_dump())
    length, width = len(sequence), settings.overlap_length
    circular = reference.topology == Topology.CIRCULAR
    minimum = minimum_fragment_count(
        length, settings.maximum_ordered_length, width, circular=circular
    )
    warnings = [
        "DNA/DNA Tm uses declared solution assumptions; not an empirical Gibson success prediction"
    ]

    class PlanMetadata(TypedDict):
        settings: GibsonSettings
        topology: Topology
        source_sha256: str
        minimum_geometric_fragments: int
        model_versions: dict[str, str]

    metadata: PlanMetadata = dict(
        settings=settings,
        topology=reference.topology,
        source_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
        minimum_geometric_fragments=minimum,
        model_versions={"Biopython": Bio.__version__, "tm_table": "DNA_NN3", "planner": "dp-v1"},
    )

    def failure(status: Literal["infeasible", "budget_exhausted"], message: str) -> GibsonPlan:
        report = ValidationReport(
            checks=(
                Finding(
                    status=CheckStatus.ERROR,
                    code="GIBSON_NO_TILING",
                    message=message,
                    proves_infeasible=status == "infeasible",
                ),
            )
        )
        return GibsonPlan(
            **metadata,
            status=status,
            fragments=(),
            overlaps=(),
            validation=report,
            warnings=tuple(warnings),
        )

    if circular and width >= length:
        return failure("infeasible", "Circular overlap must be shorter than one construct turn")
    total = length + width if circular else length
    scan = sequence + (sequence[:width] if circular else "")
    locked = [p.state == Mutability.LOCKED for p in effective_permissions(reference)]
    counts = Counter(
        scan[i : i + width] for i in range(length if circular else max(0, length - width + 1))
    )
    metrics: dict[int, OverlapMetrics] = {}

    def metric(start: int) -> OverlapMetrics:
        if start not in metrics:
            word = scan[start : start + width]
            reverse = str(Seq(word).reverse_complement())
            gc, entropy = gc_fraction(word), sequence_entropy(word)
            tm_c = overlap_tm(word, settings.conditions)
            protected = sum(locked[(start + i) % length] for i in range(width))
            components = {"protected_fraction": protected / width, "composition": (2 - entropy) / 2}
            if settings.overlap_gc_range:
                bounds = settings.overlap_gc_range
                components["gc"] = max(0, bounds.minimum - gc, gc - bounds.maximum)
            if settings.preferred_tm_c is not None:
                components["tm"] = abs(tm_c - settings.preferred_tm_c) / 100
            unique = counts[word] == 1 and counts[reverse] == 0
            components["nonunique"] = 0 if unique else 1
            if settings.synthesis_profile is not None:
                components["synthesis"] = scan_synthesis(
                    Construct(sequence=word), word, settings.synthesis_profile
                ).risk_loss
            metrics[start] = OverlapMetrics(
                start=start,
                sequence=word,
                gc=gc,
                tm_c=tm_c,
                entropy_bits=entropy,
                forward_occurrences=counts[word],
                reverse_occurrences=counts[reverse],
                protected_bases=protected,
                unique=unique,
                penalty_components=components,
            )
        return metrics[start]

    def eligible(start: int) -> bool:
        value = metric(start)
        return (not settings.require_unique_overlaps or value.unique) and (
            not settings.avoid_protected_overlaps or not value.protected_bases
        )

    if circular and not eligible(0):
        return failure(
            "infeasible",
            "The fixed circular origin violates overlap requirements; "
            "choose a new origin explicitly",
        )
    lower = settings.minimum_ordered_length - width
    upper = settings.maximum_ordered_length - width
    upper_count = max(1, (total - width) // lower)
    cap = min(settings.maximum_fragments, upper_count)
    path: tuple[int, ...] | None = None
    for count in range(minimum, cap + 1):
        execution.checkpoint(
            count - minimum, max(1, cap - minimum + 1), "Searching fragment boundaries"
        )
        if count == 1:
            if total <= settings.maximum_ordered_length:
                path = ()
                break
            continue
        if total < count * settings.minimum_ordered_length - (count - 1) * width:
            break
        previous: dict[int, float] = {0: 0.0}
        parents: list[dict[int, int]] = []
        for layer in range(1, count):
            current: dict[int, float] = {}
            parent: dict[int, int] = {}
            queue: deque[int] = deque()
            ordered = sorted(previous)
            cursor = 0
            remaining = count - layer
            earliest = max(layer * lower, total - width - remaining * upper)
            latest = min(layer * upper, total - width - remaining * lower)
            for start in range(earliest, latest + 1):
                if start % 100 == 0:
                    execution.checkpoint(layer, count, "Scoring overlap candidates")
                while cursor < len(ordered) and ordered[cursor] <= start - lower:
                    prior = ordered[cursor]
                    while queue and previous[queue[-1]] > previous[prior]:
                        queue.pop()
                    queue.append(prior)
                    cursor += 1
                while queue and queue[0] < start - upper:
                    queue.popleft()
                if not queue or not eligible(start):
                    continue
                ideal = layer * (total - width) / count
                balance = max(0, abs(start - ideal) - settings.search_window) / max(1, length)
                prior = queue[0]
                current[start] = (
                    previous[prior] + sum(metric(start).penalty_components.values()) + balance
                )
                parent[start] = prior
            previous = current
            parents.append(parent)
            if not previous:
                break
        if previous and len(parents) == count - 1:
            end = min(previous, key=lambda position: (previous[position], position))
            reverse_path = []
            for parent in reversed(parents):
                reverse_path.append(end)
                end = parent[end]
            path = tuple(reversed(reverse_path))
            break
    if path is None:
        status: Literal["budget_exhausted", "infeasible"] = (
            "budget_exhausted" if settings.maximum_fragments < upper_count else "infeasible"
        )
        return failure(
            status,
            "No tiling meets ordered-length and overlap requirements within the fragment limit",
        )
    starts, ends = (0,) + path, tuple(p + width for p in path) + (total,)
    fragments = []
    for number, (start, end) in enumerate(zip(starts, ends, strict=True), start=1):
        dna = scan[start:end]
        unique_length = (
            len(dna) - (width if number > 1 else 0) - (width if circular and end == total else 0)
        )
        assessment = (
            scan_synthesis(Construct(sequence=dna), dna, settings.synthesis_profile)
            if settings.synthesis_profile
            else None
        )
        fragments.append(
            GibsonFragment(
                number=number,
                start=start,
                end=end,
                sequence=dna,
                unique_length=unique_length,
                ordered_length=len(dna),
                gc=gc_fraction(dna),
                left_overlap=dna[:width] if number > 1 or circular else "",
                right_overlap=dna[-width:] if end < total or circular else "",
                synthesis=assessment,
            )
        )
    overlap_starts = ((0,) if circular else ()) + path
    plan = GibsonPlan(
        **metadata,
        status="success",
        fragments=tuple(fragments),
        overlaps=tuple(metric(p) for p in overlap_starts),
        validation=ValidationReport(checks=()),
        warnings=tuple(warnings),
    )
    report = validate_gibson(reference, sequence, plan)
    execution.checkpoint(1, 1, "Reconstruction validated")
    # Nonunique overlaps can be explored explicitly, but cannot be exported as unambiguous.
    return plan.model_copy(update={"validation": report})
