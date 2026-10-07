import pytest
from pydantic import ValidationError

from dna_design import (
    Construct,
    Feature,
    MotifConstraint,
    MotifMode,
    Region,
    Topology,
    check_feasibility,
    find_motif,
    restriction_constraint,
    validate_candidate,
)


def free_construct(
    sequence: str, *constraints: MotifConstraint, topology: str = "linear"
) -> Construct:
    feature = Feature(
        id="f",
        name="f",
        feature_type="spacer",
        start=0,
        end=len(sequence),
        mutability="FREE_MUTATION",
    )
    return Construct(
        sequence=sequence, features=(feature,), constraints=constraints, topology=topology
    )


def test_avoid_synonymous_repair() -> None:
    feature = Feature(id="c", name="c", feature_type="CDS", start=0, end=12, strand=1)
    avoid = restriction_constraint("EcoRI", MotifMode.AVOID)
    reference = Construct(sequence="ATGGAATTCTAA", features=(feature,), constraints=(avoid,))
    assert not validate_candidate(reference, reference.sequence).valid
    assert not check_feasibility(reference).infeasible
    assert validate_candidate(reference, "ATGGAGTTCTAA").valid


def test_required_and_multiple_motifs() -> None:
    required = MotifConstraint(id="r", motif="AGGAGG", mode="REQUIRED")
    avoid = MotifConstraint(id="a", motif="AAAA", mode="AVOID")
    reference = free_construct("AGGAGG", required, avoid)
    assert validate_candidate(reference, "AGGAGG").valid
    assert not validate_candidate(reference, "ACGAGG").valid
    assert not validate_candidate(reference, "AAAAAA").valid


def test_preserve_exact_nucleotides_not_just_degenerate_match() -> None:
    preserve = MotifConstraint(id="p", motif="GAN", mode="PRESERVE", strand="forward", position=0)
    reference = free_construct("GAACCC", preserve)
    assert validate_candidate(reference, "GAATTT").valid
    assert not validate_candidate(reference, "GAGCCC").valid


def test_preserve_all_reference_occurrences() -> None:
    preserve = MotifConstraint(id="p", motif="GAA", mode="PRESERVE", strand="forward")
    reference = free_construct("GAAGAACC", preserve)
    assert not validate_candidate(reference, "GAAGAGCC").valid


def test_iupac_overlapping_matches() -> None:
    constraint = MotifConstraint(id="m", motif="AN", mode="AVOID", strand="forward")
    assert [o.start for o in find_motif("AAAC", constraint)] == [0, 1, 2]
    degenerate = MotifConstraint(id="m", motif="ARY", mode="REQUIRED", strand="forward")
    assert len(find_motif("AGC", degenerate)) == 1


def test_reverse_strand_motif_and_palindrome() -> None:
    reverse = MotifConstraint(id="r", motif="AGG", mode="AVOID", strand="reverse")
    assert [o.start for o in find_motif("CCATCCT", reverse)] == [4]
    palindromic = restriction_constraint("EcoRI", MotifMode.REQUIRED)
    assert {o.strand for o in find_motif("GAATTC", palindromic)} == {-1, 1}


def test_motif_region_requires_full_occurrence_and_pinned_position() -> None:
    scoped = MotifConstraint(
        id="r", motif="AAA", mode="REQUIRED", strand="forward", region=Region(start=1, end=4)
    )
    assert [o.start for o in find_motif("AAAAA", scoped)] == [1]
    pinned = scoped.model_copy(update={"position": 2})
    assert not find_motif("AAAAA", pinned)


def test_circular_origin_spanning_motifs() -> None:
    constraint = restriction_constraint("EcoRI", MotifMode.REQUIRED)
    reference = free_construct("ATTCCCGA", constraint, topology="circular")
    occurrences = find_motif(reference.sequence, constraint, reference.topology)
    assert {o.start for o in occurrences} == {6}
    assert validate_candidate(reference, reference.sequence).valid
    assert not find_motif(reference.sequence, constraint, Topology.LINEAR)
    preserve = restriction_constraint("EcoRI", MotifMode.PRESERVE)
    reference = free_construct("ATTCCCGA", preserve, topology="circular")
    assert not validate_candidate(reference, "GTTCCCGA").valid


def test_circular_one_base_motif_and_full_length() -> None:
    motif = MotifConstraint(id="m", motif="A", mode="REQUIRED", strand="forward")
    assert [o.start for o in find_motif("AAA", motif, Topology.CIRCULAR)] == [0, 1, 2]
    full = MotifConstraint(id="m", motif="GAA", mode="REQUIRED", strand="forward")
    assert [o.start for o in find_motif("AAG", full, Topology.CIRCULAR)] == [2]


def test_trapped_forbidden_motif_reports_proven_infeasibility() -> None:
    avoid = restriction_constraint("EcoRI", MotifMode.AVOID)
    reference = Construct(sequence="GAATTC", constraints=(avoid,))
    report = check_feasibility(reference)
    assert report.infeasible
    assert any(c.code == "IMMUTABLE_FORBIDDEN_MOTIF" for c in report.checks)
    assert validate_candidate(reference, reference.sequence).infeasible


def test_required_forbidden_conflict_and_degenerate_subset() -> None:
    required = MotifConstraint(id="r", motif="GAA", mode="REQUIRED")
    avoid = MotifConstraint(id="a", motif="GAN", mode="AVOID")
    reference = free_construct("ACGTTT", required, avoid)
    report = check_feasibility(reference)
    assert report.infeasible
    assert any(c.code == "MOTIF_CONFLICT" for c in report.checks)


def test_distinct_scopes_and_strands_are_not_falsely_infeasible() -> None:
    required = MotifConstraint(
        id="r", motif="GAA", mode="REQUIRED", strand="forward", region=Region(start=0, end=3)
    )
    avoid = MotifConstraint(
        id="a", motif="GAA", mode="AVOID", strand="forward", region=Region(start=3, end=6)
    )
    assert not check_feasibility(free_construct("GAACCC", required, avoid)).infeasible
    reverse = avoid.model_copy(update={"region": None, "strand": "reverse"})
    assert not check_feasibility(free_construct("GAACCC", required, reverse)).infeasible


def test_preserve_absent_and_preserve_avoid_conflict() -> None:
    preserve = MotifConstraint(id="p", motif="GAA", mode="PRESERVE", strand="forward")
    assert check_feasibility(free_construct("CCCCCC", preserve)).infeasible
    avoid = MotifConstraint(id="a", motif="GAA", mode="AVOID", strand="forward")
    assert check_feasibility(free_construct("GAACCC", preserve, avoid)).infeasible


def test_bad_motifs_and_unknown_enzymes() -> None:
    for motif in ("", "ax", "A?", "AU"):
        with pytest.raises(ValidationError):
            MotifConstraint(id="m", motif=motif, mode="AVOID")
    with pytest.raises(ValueError, match="Unknown restriction"):
        restriction_constraint("DoesNotExist", MotifMode.AVOID)


def test_invalid_scope_and_duplicate_constraints() -> None:
    constraint = MotifConstraint(id="m", motif="AA", mode="REQUIRED", position=5)
    with pytest.raises(ValidationError, match="position"):
        free_construct("AAAA", constraint)
    with pytest.raises(ValidationError, match="length"):
        free_construct("A", MotifConstraint(id="m", motif="AA", mode="REQUIRED"))
    constraint = MotifConstraint(id="m", motif="A", mode="REQUIRED")
    with pytest.raises(ValidationError, match="unique"):
        free_construct("AAAA", constraint, constraint)
