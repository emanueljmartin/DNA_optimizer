"""Permanent safeguards for aggregate public-reference derivation exclusions."""

import hashlib
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from scripts.derive_codon_profile import derive


def test_derivation_excludes_invalid_cds_instead_of_crashing_or_flattening(tmp_path: Path) -> None:
    valid = "ATGGCTTGTGATGAATTTGGTCATATTAAACTGATGAATCCTCAACGTTCTACTGTTTGGTATTAA"
    record = SeqRecord(Seq(valid + "AAAGAATAA" + "ATGGAAGAA"), id="synthetic.1", name="synthetic")
    record.annotations["molecule_type"] = "DNA"
    good = SeqFeature(
        SimpleLocation(0, len(valid), strand=1), type="CDS", qualifiers={"transl_table": ["11"]}
    )
    record.features = [
        good,
        SeqFeature(
            SimpleLocation(len(valid), len(valid) + 9, strand=1),
            type="CDS",
            qualifiers={"transl_table": ["11"]},
        ),
        SeqFeature(
            SimpleLocation(len(valid) + 9, len(record), strand=1),
            type="CDS",
            qualifiers={"transl_table": ["11"]},
        ),
        SeqFeature(
            CompoundLocation(
                [SimpleLocation(0, 3, strand=1), SimpleLocation(3, len(valid), strand=1)]
            ),
            type="CDS",
            qualifiers={"transl_table": ["11"]},
        ),
        SeqFeature(good.location, type="CDS", qualifiers={"pseudo": [""], "transl_table": ["11"]}),
        SeqFeature(good.location, type="CDS", qualifiers={"transl_table": ["1"]}),
        SeqFeature(
            good.location,
            type="CDS",
            qualifiers={"transl_table": ["11"], "translation": ["INCORRECT"]},
        ),
    ]
    source = tmp_path / "public-test.gb"
    SeqIO.write(record, source, "genbank")
    profile = derive(source, "synthetic", "synthetic", "test", "2026-10-08")
    assert profile.included_cds == 1
    assert profile.excluded_cds == {
        "invalid_complete_cds": 2,
        "unsupported_location": 1,
        "exception_or_pseudogene": 1,
        "other_code_or_frame": 1,
        "translation_disagreement": 1,
    }
    assert sum(profile.counts.values()) == (len(valid) // 3) - 2
    assert profile.counts["TAA"] == 0
    assert profile.source_sha256 == hashlib.sha256(source.read_bytes()).hexdigest()


def test_derivation_rejects_non_dna_source(tmp_path: Path) -> None:
    record = SeqRecord(Seq("ACGT"), id="synthetic.1", name="synthetic")
    record.annotations["molecule_type"] = "RNA"
    source = tmp_path / "rna.gb"
    SeqIO.write(record, source, "genbank")
    with pytest.raises(ValueError, match="DNA"):
        derive(source, "synthetic", "synthetic", "test", "2026-10-08")
