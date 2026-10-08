"""Export gates, project reference authority and annotated GenBank fidelity."""

import csv
import io

import pytest
from Bio import SeqIO

from dna_design import Construct, Feature, InvalidDesignError, parse_genbank
from dna_design.gibson import GibsonSettings, plan_gibson
from dna_design.optimization import AdditionalLock
from dna_design.projects import (
    DesignProject,
    export_project,
    load_project,
    save_draft,
    validate_project,
)

from .conftest import genbank_text
from .test_gibson import synthetic_dna


def test_project_roundtrip_and_validated_exports(coding_construct: Construct) -> None:
    project = DesignProject(
        original=coding_construct,
        sequence="ATGGAGTAA",
        settings={"seed": 7},
        history=({"operation": "codon", "seed": 7},),
    )
    assert load_project(export_project(project, "json")) == project
    assert (
        str(SeqIO.read(io.StringIO(export_project(project, "fasta")), "fasta").seq)
        == project.sequence
    )
    genbank = export_project(project, "genbank")
    restored = parse_genbank(genbank)
    assert restored.sequence == project.sequence
    assert restored.features[0].start == coding_construct.features[0].start
    assert restored.features[0].translation_table == 11
    assert restored.features[0].qualifiers["translation"] == ("ME",)


@pytest.mark.parametrize("format", ["json", "fasta", "genbank", "csv"])
def test_invalid_final_design_cannot_export_but_draft_can_be_saved(
    coding_construct: Construct, format: str
) -> None:
    project = DesignProject(original=coding_construct, sequence="ATGGACTAA")
    assert not validate_project(project).valid
    with pytest.raises(InvalidDesignError):
        export_project(project, format)
    assert load_project(save_draft(project)) == project


def test_carried_initiation_locks_are_export_gates(coding_construct: Construct) -> None:
    project = DesignProject(
        original=coding_construct,
        sequence="ATGGAGTAA",
        initiation_locks=(
            AdditionalLock(start=3, end=6, sequence="GAA", reason="designed context"),
        ),
    )
    with pytest.raises(InvalidDesignError):
        export_project(project, "fasta")


def test_fragment_csv_and_forged_plan_cannot_bypass_reconstruction() -> None:
    dna = synthetic_dna(250)
    original = Construct(sequence=dna)
    plan = plan_gibson(original, dna, GibsonSettings(maximum_ordered_length=100))
    project = DesignProject(original=original, gibson=plan)
    rows = list(csv.DictReader(io.StringIO(export_project(project, "csv"))))
    assert len(rows) == len(plan.fragments)
    assert all(int(row["ordered_length"]) <= 100 for row in rows)
    forged = plan.model_copy(update={"source_sha256": "0" * 64})
    with pytest.raises(InvalidDesignError):
        export_project(project.model_copy(update={"gibson": forged}), "csv")


def test_fragment_genbank_roundtrip_and_reverse_feature() -> None:
    feature = Feature(
        id="f", name="fragment", feature_type="CDS", start=0, end=9, strand=-1, cds_mode="fragment"
    )
    original = Construct(sequence="TTATTCCAT", features=(feature,))
    imported = parse_genbank(export_project(DesignProject(original=original), "genbank"))
    assert imported.features[0].strand == -1
    assert imported.features[0].cds_mode == "fragment"
    assert imported.features[0].qualifiers["translation"] == ("ME",)


def test_source_and_literature_metadata_survive_genbank_export() -> None:
    text = genbank_text()
    original = parse_genbank(text)
    project = DesignProject(original=original)
    emitted = export_project(project, "genbank")
    record = SeqIO.read(io.StringIO(emitted), "genbank")
    original_record = SeqIO.read(io.StringIO(text), "genbank")
    assert [(str(f.location), f.qualifiers) for f in record.features if f.type == "source"] == [
        (str(f.location), f.qualifiers) for f in original_record.features if f.type == "source"
    ]


def test_genbank_roundtrip_preserves_explicit_ids_mutability_and_rbs_mask() -> None:
    from .test_rbs import rbs_construct

    original, _ = rbs_construct()
    imported = parse_genbank(export_project(DesignProject(original=original), "genbank"))
    for before, after in zip(original.features, imported.features, strict=True):
        assert (after.id, after.mutability, after.rbs_mask, after.strand) == (
            before.id,
            before.mutability,
            before.rbs_mask,
            before.strand,
        )
