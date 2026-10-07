import pytest
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from hypothesis import given
from hypothesis import strategies as st

from dna_design import (
    CheckStatus,
    Construct,
    Feature,
    InvalidDesignError,
    MotifConstraint,
    ParseError,
    UnsupportedLocationError,
    check_feasibility,
    parse_fasta,
    parse_genbank,
    parse_raw,
    translate_feature,
    validate_candidate,
)

from .conftest import genbank_text


def test_raw_and_fasta_import() -> None:
    assert parse_raw(" acg\n t \t").sequence == "ACGT"
    fasta = parse_fasta(">synthetic description\nacgt\nACGT\n")
    assert fasta.name == "synthetic description"
    assert fasta.sequence == "ACGTACGT"
    assert fasta.source_sequence == fasta.sequence


@pytest.mark.parametrize("text", ["", "ATGU", "ACGN", ">header\nACGT"])
def test_bad_raw(text: str) -> None:
    with pytest.raises(ParseError):
        parse_raw(text)


@pytest.mark.parametrize(
    "text", ["ACGT", ">\nACGT", ">one\nACGT\n>two\nACGT", ">one\nACGN", ">one\n"]
)
def test_bad_fasta(text: str) -> None:
    with pytest.raises(ParseError):
        parse_fasta(text)


@pytest.mark.parametrize("strand", [1, -1])
def test_genbank_annotations_and_strands(strand: int) -> None:
    promoter = SeqFeature(
        SimpleLocation(0, 2, strand=1), type="promoter", qualifiers={"note": ["locked"]}
    )
    text = genbank_text(strand=strand, extra_features=(promoter,), topology="circular")
    construct = parse_genbank(text)
    assert construct.topology == "circular"
    cds = construct.features[0]
    assert (cds.start, cds.end, cds.strand) == (0, 9, strand)
    assert cds.qualifiers["translation"] == ("ME",)
    assert translate_feature(construct.sequence, cds) == "ME"
    assert construct.features[1].mutability == "LOCKED"
    assert construct.record_annotations["source_features"][0]["qualifiers"]["organism"] == [
        "synthetic"
    ]
    assert construct.original_record_text == text
    assert Construct.model_validate_json(construct.model_dump_json()) == construct


def test_genbank_import_then_synonymous_validation() -> None:
    construct = parse_genbank(genbank_text())
    assert validate_candidate(construct, "ATGGAGTAA").valid
    assert not validate_candidate(construct, "ATGGACTAA").valid


def test_unsupported_compound_location_explicitly_reported() -> None:
    compound = CompoundLocation([SimpleLocation(0, 3, strand=1), SimpleLocation(6, 9, strand=1)])
    with pytest.raises(UnsupportedLocationError, match="compound"):
        parse_genbank(genbank_text(location=compound))


def test_fuzzy_location_explicitly_reported() -> None:
    location = SimpleLocation(BeforePosition(0), 9, strand=1)
    with pytest.raises(UnsupportedLocationError, match="fuzzy"):
        parse_genbank(genbank_text(location=location))


def test_remote_location_explicitly_reported() -> None:
    location = SimpleLocation(0, 9, strand=1, ref="other")
    with pytest.raises(UnsupportedLocationError, match="remote"):
        parse_genbank(genbank_text(location=location))


def test_genbank_translation_disagreement_rejected() -> None:
    with pytest.raises(ParseError, match="disagrees"):
        parse_genbank(genbank_text(qualifiers={"translation": ["MA"]}))


def test_genbank_frame_offset_and_genetic_code() -> None:
    construct = parse_genbank(
        genbank_text(sequence="AGAAGAA", qualifiers={"codon_start": ["2"], "translation": ["EE"]})
    )
    assert construct.features[0].frame_offset == 1
    assert construct.features[0].cds_mode == "fragment"
    construct = parse_genbank(
        genbank_text(
            sequence="ATGTGATAA", qualifiers={"transl_table": ["4"], "translation": ["MW"]}
        )
    )
    assert construct.features[0].translation_table == 4


@pytest.mark.parametrize(
    "key", ["transl_except", "exception", "ribosomal_slippage", "pseudo", "pseudogene"]
)
def test_genbank_unsupported_recoding_exception(key: str) -> None:
    with pytest.raises(ParseError, match="unsupported"):
        parse_genbank(genbank_text(qualifiers={key: ["unsupported"]}))


def test_multiple_genbank_and_invalid_records_rejected() -> None:
    with pytest.raises(ParseError, match="one GenBank"):
        parse_genbank(genbank_text() + genbank_text())
    with pytest.raises(ParseError):
        parse_genbank("not genbank")


def test_candidate_cannot_redefine_original_protein(coding_construct: Construct) -> None:
    # Even removing the CDS control cannot bypass translation of the original CDS.
    malicious = Construct(sequence="ATGGACTAA", source_sequence=coding_construct.sequence)
    report = validate_candidate(coding_construct, malicious)
    assert not report.valid
    assert any(c.code == "REFERENCE_CONTROLS" and c.status == "ERROR" for c in report.checks)
    assert any(c.code == "PROTEIN_INTEGRITY" and c.status == "ERROR" for c in report.checks)


def test_candidate_cannot_redefine_protected_dna_or_source() -> None:
    lock = Feature(id="l", name="l", feature_type="protected region", start=0, end=4)
    reference = Construct(sequence="ACGT", features=(lock,))
    malicious = Construct(sequence="ACGA", source_sequence="ACGA", features=())
    assert not validate_candidate(reference, malicious).valid
    # Independently passing a string still uses reference.sequence, not source_sequence.
    reference = Construct(sequence="ACGT", source_sequence="ACGA", features=(lock,))
    assert not validate_candidate(reference, "ACGA").valid


def test_candidate_cannot_drop_motif_or_alter_topology() -> None:
    required = MotifConstraint(id="r", motif="GAA", mode="REQUIRED", strand="forward")
    free = Feature(
        id="f", name="f", feature_type="spacer", start=0, end=3, mutability="FREE_MUTATION"
    )
    reference = Construct(sequence="GAA", features=(free,), constraints=(required,))
    malicious = reference.model_copy(
        update={"sequence": "CCC", "constraints": (), "topology": "circular"}
    )
    report = validate_candidate(reference, malicious)
    assert not report.valid
    assert any(c.code == "MOTIF_REQUIRED" and c.status == "ERROR" for c in report.checks)


def test_reference_controls_checked_even_for_unchanged_candidate(
    coding_construct: Construct,
) -> None:
    malicious = coding_construct.model_copy(update={"features": ()})
    assert not validate_candidate(coding_construct, malicious).valid


def test_bypassed_model_schema_still_checked(coding_construct: Construct) -> None:
    malicious = coding_construct.model_copy(update={"sequence": "NNNNNNNNN"})
    assert not validate_candidate(coding_construct, malicious).valid
    invalid_reference = coding_construct.model_copy(update={"sequence": "NNNNNNNNN"})
    assert not validate_candidate(invalid_reference, "ATGGAATAA").valid
    assert not check_feasibility(invalid_reference).valid


def test_reference_protein_invalid_is_not_validated_against_itself(cds: Feature) -> None:
    invalid_reference = Construct(sequence="ATGTAATAA", features=(cds,))
    report = validate_candidate(invalid_reference, invalid_reference.sequence)
    assert not report.valid
    assert any(c.code == "REFERENCE_CDS_INVALID" for c in report.checks)


def test_original_rbs_mask_mismatch_reported() -> None:
    feature = Feature(
        id="r",
        name="r",
        feature_type="RBS",
        start=0,
        end=3,
        strand=1,
        mutability="DESIGNABLE_RBS",
        rbs_mask="AGN",
    )
    reference = Construct(sequence="AAA", features=(feature,))
    assert not validate_candidate(reference, reference.sequence).valid
    assert not check_feasibility(reference).valid


def test_simultaneous_constraints(cds: Feature) -> None:
    lock = Feature(id="l", name="l", feature_type="protected region", start=0, end=3)
    preserve = MotifConstraint(id="p", motif="TAA", mode="PRESERVE", position=6, strand="forward")
    avoid = MotifConstraint(id="a", motif="GAA", mode="AVOID", strand="forward")
    reference = Construct(sequence="ATGGAATAA", features=(cds, lock), constraints=(preserve, avoid))
    assert validate_candidate(reference, "ATGGAGTAA").valid
    assert not validate_candidate(reference, "GTGGAGTAG").valid


def test_structured_report_and_export_gate(coding_construct: Construct) -> None:
    reference = coding_construct.model_copy(update={"import_warnings": ("Synthetic warning",)})
    good = validate_candidate(reference, "ATGGAGTAA")
    assert good.valid
    assert {c.status for c in good.checks} == {CheckStatus.PASS, CheckStatus.WARNING}
    good.require_valid()
    bad = validate_candidate(reference, "ATGGACTAA")
    assert CheckStatus.ERROR in {c.status for c in bad.checks}
    with pytest.raises(InvalidDesignError):
        bad.require_valid()
    assert '"ERROR"' in bad.model_dump_json()


@pytest.mark.parametrize("candidate", ["", "acgt", "ACGN", "ATG", "ATGGAATAAAA"])
def test_candidate_alphabet_and_length(coding_construct: Construct, candidate: str) -> None:
    assert not validate_candidate(coding_construct, candidate).valid


@given(st.text(alphabet="ACGT", min_size=1, max_size=500))
def test_sequence_import_round_trips(dna: str) -> None:
    assert parse_raw(dna.lower()).sequence == dna
    assert parse_fasta(">synthetic\n" + dna + "\n").sequence == dna
