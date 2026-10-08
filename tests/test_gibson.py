"""Synthetic exact reconstruction and assembly geometry, not wet-lab yield tests."""

import random

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from dna_design import Construct, InvalidDesignError
from dna_design.gibson import (
    GibsonSettings,
    TmConditions,
    minimum_fragment_count,
    overlap_tm,
    plan_gibson,
    reconstruct,
    validate_gibson,
)


def synthetic_dna(length: int, seed: int = 1) -> str:
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


@pytest.mark.parametrize(
    "length,maximum,overlap,expected",
    [(100, 100, 30, 1), (101, 100, 30, 2), (240, 100, 30, 3), (241, 100, 30, 4)],
)
def test_minimum_count_includes_overlaps(
    length: int, maximum: int, overlap: int, expected: int
) -> None:
    assert minimum_fragment_count(length, maximum, overlap) == expected
    assert expected * maximum - (expected - 1) * overlap >= length
    assert expected == 1 or (expected - 1) * maximum - (expected - 2) * overlap < length


@pytest.mark.parametrize("topology", ["linear", "circular"])
def test_plan_reconstructs_exactly_with_overlap_lengths_and_order(topology: str) -> None:
    dna = synthetic_dna(650)
    original = Construct(sequence=dna, topology=topology)
    config = GibsonSettings(
        maximum_ordered_length=200, overlap_length=20, minimum_ordered_length=60
    )
    result = plan_gibson(original, dna, config)
    assert result.status == "success"
    assert result.validation.valid
    assert reconstruct(result) == dna
    assert all(f.ordered_length == len(f.sequence) <= 200 for f in result.fragments)
    assert sum(f.unique_length for f in result.fragments) == len(dna)
    assert len(result.fragments) == minimum_fragment_count(
        650, 200, 20, circular=topology == "circular"
    )
    assert all(overlap.unique for overlap in result.overlaps)
    assert result == plan_gibson(original, dna, config)
    assert validate_gibson(original, dna, result).valid


def test_single_linear_and_circular_tiles() -> None:
    dna = synthetic_dna(75, 3)
    settings = GibsonSettings(maximum_ordered_length=100, overlap_length=20)
    for topology in ("linear", "circular"):
        original = Construct(sequence=dna, topology=topology)
        plan = plan_gibson(original, dna, settings)
        assert plan.validation.valid
        assert len(plan.fragments) == 1
        assert plan.fragments[0].ordered_length == (95 if topology == "circular" else 75)
        assert reconstruct(plan) == dna


def test_short_linear_tile_and_infeasible_minimum_lengths() -> None:
    original = Construct(sequence="ACGT")
    plan = plan_gibson(original, original.sequence, GibsonSettings(maximum_ordered_length=100))
    assert plan.validation.valid
    dna = synthetic_dna(101)
    original = Construct(sequence=dna)
    plan = plan_gibson(
        original, dna, GibsonSettings(maximum_ordered_length=100, minimum_ordered_length=90)
    )
    assert plan.status == "infeasible"
    assert not plan.validation.valid


def test_ambiguous_overlap_is_reported_and_cannot_pass_validation() -> None:
    original = Construct(sequence="A" * 250)
    config = GibsonSettings(maximum_ordered_length=100)
    impossible = plan_gibson(original, original.sequence, config)
    assert impossible.status == "infeasible"
    diagnostic = plan_gibson(
        original, original.sequence, config.model_copy(update={"require_unique_overlaps": False})
    )
    assert diagnostic.fragments
    assert not diagnostic.validation.valid
    with pytest.raises(InvalidDesignError):
        diagnostic.validation.require_valid()


def test_fresh_reconstruction_rejects_forged_fragment_metadata_and_dna() -> None:
    dna = synthetic_dna(250)
    original = Construct(sequence=dna)
    plan = plan_gibson(original, dna, GibsonSettings(maximum_ordered_length=100))
    first = plan.fragments[0]
    altered = first.model_copy(update={"sequence": "C" + first.sequence[1:]})
    if altered.sequence == first.sequence:
        altered = first.model_copy(update={"sequence": "A" + first.sequence[1:]})
    tampered = plan.model_copy(update={"fragments": (altered,) + plan.fragments[1:]})
    assert not validate_gibson(original, dna, tampered).valid
    altered = first.model_copy(update={"ordered_length": 1, "unique_length": 1})
    assert not validate_gibson(
        original, dna, plan.model_copy(update={"fragments": (altered,) + plan.fragments[1:]})
    ).valid
    reversed_plan = plan.model_copy(update={"fragments": tuple(reversed(plan.fragments))})
    assert not validate_gibson(original, dna, reversed_plan).valid


def test_protected_boundaries_and_fragment_budget_are_explicit() -> None:
    dna = synthetic_dna(250)
    original = Construct(sequence=dna)  # Unannotated DNA is locked.
    config = GibsonSettings(maximum_ordered_length=100, avoid_protected_overlaps=True)
    assert plan_gibson(original, dna, config).status == "infeasible"
    assert (
        plan_gibson(original, dna, config.model_copy(update={"maximum_fragments": 1})).status
        == "budget_exhausted"
    )
    with pytest.raises(InvalidDesignError):
        plan_gibson(original, ("C" if dna[0] == "A" else "A") + dna[1:], config)


def test_export_metadata_and_protected_overlap_policy_are_freshly_checked() -> None:
    dna = synthetic_dna(250)
    original = Construct(sequence=dna)
    plan = plan_gibson(original, dna, GibsonSettings(maximum_ordered_length=100))
    first = plan.fragments[0].model_copy(update={"right_overlap": "A" * 30, "gc": 0})
    assert not validate_gibson(
        original, dna, plan.model_copy(update={"fragments": (first,) + plan.fragments[1:]})
    ).valid
    overlap = plan.overlaps[0].model_copy(update={"tm_c": 0})
    assert not validate_gibson(
        original, dna, plan.model_copy(update={"overlaps": (overlap,) + plan.overlaps[1:]})
    ).valid
    forged = plan.model_copy(
        update={"settings": plan.settings.model_copy(update={"avoid_protected_overlaps": True})}
    )
    assert not validate_gibson(original, dna, forged).valid
    assert not validate_gibson(
        original, dna, plan.model_copy(update={"minimum_geometric_fragments": 100})
    ).valid


def test_tm_matches_public_biopython_example_and_shows_conditions() -> None:
    word = "CGTTCCAAAGATGTGGGCATGAGCTTAC"
    assert overlap_tm(word, TmConditions()) == pytest.approx(60.32, abs=0.01)
    assert overlap_tm(word, TmConditions(sodium_mm=100)) > overlap_tm(word, TmConditions())
    with pytest.raises(ValueError):
        overlap_tm("A", TmConditions())


@settings(max_examples=25, deadline=None)
@given(st.integers(120, 600), st.integers(0, 65535))
def test_randomized_tiling_has_no_missing_or_duplicated_bases(length: int, seed: int) -> None:
    dna = synthetic_dna(length, seed)
    original = Construct(sequence=dna)
    plan = plan_gibson(
        original,
        dna,
        GibsonSettings(maximum_ordered_length=100, overlap_length=12, minimum_ordered_length=30),
    )
    assert plan.status == "success"
    assert plan.validation.valid
    assert reconstruct(plan) == dna
    assert sum(f.unique_length for f in plan.fragments) == length
