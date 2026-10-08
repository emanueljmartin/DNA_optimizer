from io import StringIO
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pydantic import ValidationError

from dna_design import (
    Construct,
    Feature,
    MotifConstraint,
    ParseError,
    check_feasibility,
    parse_genbank,
    translate_feature,
    validate_candidate,
)

from .conftest import genbank_text


@pytest.mark.parametrize("molecule_type", ["RNA", "mRNA", "protein"])
def test_genbank_non_dna_molecule_with_dna_alphabet_is_rejected(molecule_type: str) -> None:
    record = SeqRecord(Seq("ACGT"), id="synthetic", name="synthetic")
    record.annotations["molecule_type"] = molecule_type
    output = StringIO()
    SeqIO.write(record, output, "genbank")
    with pytest.raises(ParseError, match="Unsupported GenBank molecule type"):
        parse_genbank(output.getvalue())


def test_genbank_coordinates_in_literal_independent_fixture() -> None:
    text = (Path(__file__).parent / "fixtures/synthetic.gb").read_text()
    reference = parse_genbank(text)
    promoter, forward, reverse = reference.features
    assert (promoter.start, promoter.end) == (0, 3)
    assert (forward.start, forward.end, forward.strand) == (3, 12, 1)
    assert (reverse.start, reverse.end, reverse.strand) == (12, 21, -1)
    assert translate_feature(reference.sequence, forward) == "ME"
    assert translate_feature(reference.sequence, reverse) == "ME"
    assert validate_candidate(reference, "GGGATGGAGTAATTATTCCAT").valid
    assert validate_candidate(reference, "GGGATGGAATAATTACTCCAT").valid
    assert not validate_candidate(reference, "AGGATGGAATAATTATTCCAT").valid


def test_genbank_omitted_genetic_code_means_standard_code() -> None:
    reference = parse_genbank(genbank_text(qualifiers={"translation": ["ME"]}))
    assert reference.features[0].translation_table == 1
    # GTG is an alternative start in table 11, but not in standard table 1.
    with pytest.raises(ParseError, match="start codon"):
        parse_genbank(genbank_text(sequence="GTGGAATAA", qualifiers={"translation": ["ME"]}))


@pytest.mark.parametrize("extra", ["before", "after", "no-terminator"])
def test_genbank_no_silent_surrounding_data_loss(extra: str) -> None:
    text = genbank_text()
    if extra == "before":
        text = "unparsed prefix\n" + text
    elif extra == "after":
        text += "unparsed suffix\n"
    else:
        text = text.removesuffix("//\n")
    with pytest.raises(ParseError):
        parse_genbank(text)


def test_genbank_gene_overlap_is_conservatively_locked() -> None:
    gene = SeqFeature(SimpleLocation(0, 9, strand=1), type="gene", qualifiers={"gene": ["test"]})
    reference = parse_genbank(genbank_text(extra_features=(gene,)))
    assert not validate_candidate(reference, "ATGGAGTAA").valid


def test_serialization_bypass_unserializable_settings_are_reported(
    coding_construct: Construct,
) -> None:
    candidate = coding_construct.model_copy(update={"rbs_settings": {"invalid": object()}})
    # Depending on serializer mode the object survives until schema checking or
    # raises a serialization error; both must be reported rather than accepted.
    assert not validate_candidate(coding_construct, candidate).valid


def test_bypassed_enum_values_revalidate_without_warning(coding_construct: Construct) -> None:
    candidate = coding_construct.model_copy(update={"topology": "circular"})
    assert not validate_candidate(coding_construct, candidate).valid


def test_required_motif_in_locked_wrong_sequence_is_proven_impossible() -> None:
    required = MotifConstraint(id="r", motif="GAA", mode="REQUIRED", strand="forward")
    reference = Construct(sequence="CCC", constraints=(required,))
    assert check_feasibility(reference).infeasible


def test_required_scope_with_no_complete_window_is_proven_impossible() -> None:
    required = MotifConstraint(id="r", motif="GAA", mode="REQUIRED", region={"start": 0, "end": 2})
    free = Feature(
        id="f", name="f", feature_type="spacer", start=0, end=3, mutability="FREE_MUTATION"
    )
    reference = Construct(sequence="CCC", features=(free,), constraints=(required,))
    assert check_feasibility(reference).infeasible


def test_required_motif_can_be_introduced_in_designable_sequence() -> None:
    required = MotifConstraint(id="r", motif="GAA", mode="REQUIRED", strand="forward")
    free = Feature(
        id="f", name="f", feature_type="spacer", start=0, end=3, mutability="FREE_MUTATION"
    )
    reference = Construct(sequence="CCC", features=(free,), constraints=(required,))
    assert not check_feasibility(reference).infeasible
    assert validate_candidate(reference, "GAA").valid


@pytest.mark.parametrize("enzyme", ["EcoRI", "DoesNotExist"])
def test_direct_constraint_cannot_mislabel_restriction_recognition(enzyme: str) -> None:
    with pytest.raises(ValidationError):
        MotifConstraint(id="m", motif="AAA", mode="AVOID", enzyme=enzyme)
