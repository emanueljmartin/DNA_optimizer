import pytest
from Bio.Seq import Seq
from hypothesis import given
from hypothesis import strategies as st
from pydantic import ValidationError

from dna_design import (
    CDSMode,
    Construct,
    Feature,
    Mutability,
    TranslationError,
    effective_permissions,
    translate_feature,
    validate_candidate,
)


def test_synonymous_passes_nonsynonymous_fails(coding_construct: Construct) -> None:
    assert validate_candidate(coding_construct, "ATGGAGTAA").valid
    report = validate_candidate(coding_construct, "ATGGACTAA")
    assert not report.valid
    assert any(c.code == "PROTEIN_INTEGRITY" for c in report.checks)


def test_complete_cds_start_stop_integrity(coding_construct: Construct) -> None:
    assert translate_feature(coding_construct.sequence, coding_construct.features[0]) == "ME"
    assert not validate_candidate(coding_construct, "ATGGAACAA").valid
    assert not validate_candidate(coding_construct, "ATGTAATAA").valid
    assert not validate_candidate(coding_construct, "ACGGAATAA").valid
    # Bacterial alternative starts are initiator methionine, not elongation valine.
    assert validate_candidate(coding_construct, "GTGGAATAA").valid


def test_reverse_strand_cds() -> None:
    feature = Feature(id="r", name="r", feature_type="CDS", start=0, end=9, strand=-1)
    sequence = str(Seq("ATGGAATAA").reverse_complement())
    reference = Construct(sequence=sequence, features=(feature,))
    assert translate_feature(sequence, feature) == "ME"
    assert validate_candidate(reference, str(Seq("ATGGAGTAA").reverse_complement())).valid
    assert not validate_candidate(reference, str(Seq("ATGGACTAA").reverse_complement())).valid


def test_alternative_translation_tables() -> None:
    # Table 4 translates TGA as W; table 11 treats it as an internal stop.
    feature = Feature(
        id="c", name="c", feature_type="CDS", start=0, end=9, strand=1, translation_table=4
    )
    reference = Construct(sequence="ATGTGATAA", features=(feature,))
    assert translate_feature(reference.sequence, feature) == "MW"
    assert validate_candidate(reference, "ATGTGGTAA").valid
    feature11 = feature.model_copy(update={"translation_table": 11})
    with pytest.raises(TranslationError, match="stop"):
        translate_feature(reference.sequence, feature11)


def test_overlapping_shifted_cds_preserve_all() -> None:
    # GAA -> GAG preserves E in the main frame but changes AAT -> AGT (N -> S)
    # in the shifted CDS. Fragment semantics are explicit for this synthetic case.
    first = Feature(id="a", name="a", feature_type="CDS", start=0, end=9, strand=1)
    second = Feature(
        id="b", name="b", feature_type="CDS", start=1, end=7, strand=1, cds_mode=CDSMode.FRAGMENT
    )
    reference = Construct(sequence="ATGGAATAA", features=(first, second))
    report = validate_candidate(reference, "ATGGAGTAA")
    assert not report.valid
    failures = [
        c.feature_id for c in report.checks if c.code == "PROTEIN_INTEGRITY" and c.status == "ERROR"
    ]
    assert failures == ["b"]
    # A degenerate change can be synonymous for two overlapping forward CDSs.
    a = Feature(id="a", name="a", feature_type="CDS", start=0, end=6, strand=1, cds_mode="fragment")
    b = Feature(id="b", name="b", feature_type="CDS", start=3, end=9, strand=1, cds_mode="fragment")
    overlap = Construct(sequence="GAAGAAGAA", features=(a, b))
    assert validate_candidate(overlap, "GAAGAGGAA").valid


def test_overlapping_opposite_strands() -> None:
    a = Feature(id="a", name="a", feature_type="CDS", start=0, end=9, strand=1)
    b = Feature(
        id="b", name="b", feature_type="CDS", start=0, end=9, strand=-1, cds_mode="fragment"
    )
    reference = Construct(sequence="ATGGAATAA", features=(a, b))
    assert not validate_candidate(reference, "ATGGAGTAA").valid


def test_locked_dominates_and_annotation_order_does_not_matter(cds: Feature) -> None:
    lock = Feature(id="l", name="lock", feature_type="protected region", start=3, end=6)
    free = Feature(
        id="f",
        name="free",
        feature_type="spacer",
        start=0,
        end=9,
        mutability=Mutability.FREE_MUTATION,
    )
    for features in ((cds, lock, free), (free, lock, cds)):
        reference = Construct(sequence="ATGGAATAA", features=features)
        assert effective_permissions(reference)[4].state == Mutability.LOCKED
        assert effective_permissions(reference)[0].state == Mutability.SYNONYMOUS_ONLY
        assert not validate_candidate(reference, "ATGGAGTAA").valid
        assert validate_candidate(reference, reference.sequence).valid


def test_free_cannot_weaken_functional_annotation() -> None:
    promoter = Feature(id="p", name="p", feature_type="promoter", start=0, end=3)
    free = Feature(
        id="f", name="f", feature_type="spacer", start=0, end=6, mutability="FREE_MUTATION"
    )
    reference = Construct(sequence="ACGTAC", features=(promoter, free))
    assert not validate_candidate(reference, "TCGTAC").valid
    assert validate_candidate(reference, "ACGAAA").valid


def test_unannotated_and_uninterpretable_synonymous_positions_lock() -> None:
    reference = Construct(sequence="ACGT")
    assert all(p.state == Mutability.LOCKED for p in effective_permissions(reference))
    assert not validate_candidate(reference, "ACGA").valid
    feature = Feature(
        id="f", name="f", feature_type="custom", start=0, end=4, mutability="SYNONYMOUS_ONLY"
    )
    reference = Construct(sequence="ACGT", features=(feature,))
    assert not validate_candidate(reference, "ACGA").valid


def test_explicit_cds_free_mutation_still_preserves_protein(cds: Feature) -> None:
    feature = cds.model_copy(update={"mutability": Mutability.FREE_MUTATION})
    reference = Construct(sequence="ATGGAATAA", features=(feature,))
    assert not validate_candidate(reference, "ATGGACTAA").valid


@pytest.mark.parametrize("strand", [1, -1])
def test_fragment_frame_offset_and_leading_base_lock(strand: int) -> None:
    feature = Feature(
        id="f",
        name="f",
        feature_type="CDS",
        start=0,
        end=7,
        strand=strand,
        frame_offset=1,
        cds_mode="fragment",
    )
    sequence = "AGAAGAA"
    candidate = "AGAAGAG"
    if strand == -1:
        sequence = str(Seq(sequence).reverse_complement())
        candidate = str(Seq(candidate).reverse_complement())
    reference = Construct(sequence=sequence, features=(feature,))
    assert translate_feature(sequence, feature) == "EE"
    assert validate_candidate(reference, candidate).valid
    position = 0 if strand == 1 else 6
    bases = list(sequence)
    bases[position] = "C" if bases[position] != "C" else "A"
    assert not validate_candidate(reference, "".join(bases)).valid


def test_complete_cds_cannot_use_frame_offset() -> None:
    with pytest.raises(ValidationError, match="frame offset"):
        Feature(id="f", name="f", feature_type="CDS", start=0, end=7, strand=1, frame_offset=1)


@pytest.mark.parametrize("strand", [1, -1])
def test_rbs_mask_fixed_and_mutable_positions(strand: int) -> None:
    feature = Feature(
        id="r",
        name="r",
        feature_type="RBS",
        start=0,
        end=6,
        strand=strand,
        mutability="DESIGNABLE_RBS",
        rbs_mask="NNAGGN",
    )
    sequence = "AAAGGT"
    mutable = "CCAGGA"
    fixed = "AAAAGT"
    if strand == -1:
        sequence, mutable, fixed = (
            str(Seq(s).reverse_complement()) for s in (sequence, mutable, fixed)
        )
    reference = Construct(sequence=sequence, features=(feature,))
    assert validate_candidate(reference, mutable).valid
    assert not validate_candidate(reference, fixed).valid


@pytest.mark.parametrize("mask", [None, "N", "NNXNNN"])
def test_invalid_rbs_masks(mask: str | None) -> None:
    with pytest.raises(ValidationError):
        Feature(
            id="r",
            name="r",
            feature_type="RBS",
            start=0,
            end=6,
            strand=1,
            mutability="DESIGNABLE_RBS",
            rbs_mask=mask,
        )


@given(st.lists(st.sampled_from(("GAA", "GAG")), min_size=1, max_size=40))
def test_property_synonymous_preservation(codons: list[str]) -> None:
    dna = "ATG" + "".join(codons) + "TAA"
    feature = Feature(id="c", name="c", feature_type="CDS", start=0, end=len(dna), strand=1)
    reference = Construct(sequence=dna, features=(feature,))
    alternative = "ATG" + "".join("GAG" if c == "GAA" else "GAA" for c in codons) + "TAG"
    assert validate_candidate(reference, alternative).valid


@given(st.text(alphabet="ACGT", min_size=1, max_size=200), st.data())
def test_property_any_lock_mutation_fails(dna: str, data: st.DataObject) -> None:
    position = data.draw(st.integers(0, len(dna) - 1))
    alternative = data.draw(
        st.sampled_from(tuple(base for base in "ACGT" if base != dna[position]))
    )
    feature = Feature(id="l", name="l", feature_type="protected region", start=0, end=len(dna))
    reference = Construct(sequence=dna, features=(feature,))
    candidate = dna[:position] + alternative + dna[position + 1 :]
    assert not validate_candidate(reference, candidate).valid
