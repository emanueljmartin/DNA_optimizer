import pytest
from hypothesis import given
from hypothesis import strategies as st
from pydantic import ValidationError

from dna_design import (
    Construct,
    Feature,
    Mutability,
    Region,
    internal_to_ui,
    ui_to_internal,
)


@pytest.mark.parametrize("dna", ["A", "ACGT", "ATGGAATAA"])
def test_legal_dna(dna: str) -> None:
    assert Construct(sequence=dna).source_sequence == dna


@pytest.mark.parametrize("dna", ["", "N", "acgt", "AC GU", "ACTU", "ACT-", "ACT1"])
def test_illegal_dna(dna: str) -> None:
    with pytest.raises(ValidationError):
        Construct(sequence=dna)


def test_default_mutability_and_frozen_controls() -> None:
    cds = Feature(id="c", name="c", feature_type="CDS", start=0, end=9, strand=1)
    assert cds.mutability == Mutability.SYNONYMOUS_ONLY
    promoter = Feature(id="p", name="p", feature_type="promoter", start=0, end=3)
    assert promoter.mutability == Mutability.LOCKED
    with pytest.raises(ValidationError, match="frozen"):
        cds.start = 1


@pytest.mark.parametrize("feature_type", ["promoter", "operator", "terminator", "RBS", "custom"])
def test_functional_noncoding_defaults_locked(feature_type: str) -> None:
    assert (
        Feature(id="f", name="f", feature_type=feature_type, start=0, end=3).mutability
        == Mutability.LOCKED
    )


@pytest.mark.parametrize("start,end", [(-1, 3), (3, 3), (3, 2), (True, 3), (1.2, 3)])
def test_bad_intervals(start: int, end: int) -> None:
    with pytest.raises(ValidationError):
        Region(start=start, end=end)


def test_controls_bounds_and_identifiers(cds: Feature) -> None:
    with pytest.raises(ValidationError, match="length"):
        Construct(sequence="A", features=(cds,))
    with pytest.raises(ValidationError, match="unique"):
        Construct(sequence="ATGGAATAA", features=(cds, cds))
    with pytest.raises(ValidationError, match="strand"):
        Feature(id="c", name="c", feature_type="CDS", start=0, end=9)
    with pytest.raises(ValidationError, match="multiple"):
        Feature(id="c", name="c", feature_type="CDS", start=0, end=8, strand=1)
    with pytest.raises(ValidationError, match="genetic code"):
        Feature(
            id="c", name="c", feature_type="CDS", start=0, end=9, strand=1, translation_table=999
        )


def test_dual_stop_sense_code_is_explicitly_unsupported() -> None:
    with pytest.raises(ValidationError, match="Context-dependent"):
        Feature(
            id="c", name="c", feature_type="CDS", start=0, end=9, strand=1, translation_table=27
        )


def test_coordinate_boundaries() -> None:
    assert ui_to_internal(1, 1, 1) == (0, 1)
    assert internal_to_ui(0, 1, 1) == (1, 1)
    assert ui_to_internal(1, 10, 10) == (0, 10)
    assert ui_to_internal(4, 9, 9) == (3, 9)
    assert internal_to_ui(3, 9, 9) == (4, 9)


@pytest.mark.parametrize("start,end", [(0, 1), (2, 1), (1, 11), (True, 1), (1.0, 2)])
def test_ui_coordinate_errors(start: int, end: int) -> None:
    with pytest.raises(ValueError):
        ui_to_internal(start, end, 10)


@pytest.mark.parametrize("start,end", [(-1, 1), (0, 0), (2, 1), (0, 11), (False, 2)])
def test_internal_coordinate_errors(start: int, end: int) -> None:
    with pytest.raises(ValueError):
        internal_to_ui(start, end, 10)


@given(st.integers(min_value=1, max_value=100000), st.integers(min_value=1, max_value=10000))
def test_coordinate_round_trip(start: int, span: int) -> None:
    end = start + span - 1
    assert internal_to_ui(*ui_to_internal(start, end, end), end) == (start, end)


@given(st.text(alphabet="ACGT", min_size=1, max_size=500))
def test_json_round_trip(dna: str) -> None:
    construct = Construct(sequence=dna, random_seed=42, model_versions={"test": "1"})
    assert Construct.model_validate_json(construct.model_dump_json()) == construct
