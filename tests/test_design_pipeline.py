"""Cross-phase original-reference validation; synthetic inputs stay local."""

import csv
import io

import pytest
from Bio import SeqIO

from dna_design import validate_candidate
from dna_design.gibson import GibsonSettings, plan_gibson, reconstruct
from dna_design.optimization import SearchConfig, optimize
from dna_design.projects import (
    DesignProject,
    export_project,
    load_project,
    save_draft,
    validate_project,
)
from dna_design.rbs import RBSDesignSettings, design_rbs, predict_rbs
from dna_design.synthesis import SynthesisProfile, repair_synthesis

from .test_optimization import synthetic_profile
from .test_rbs import CountingTestBackend, rbs_construct


@pytest.mark.parametrize("strand", [1, -1])
def test_phases_two_through_four_keep_original_controls_and_prediction(strand: int) -> None:
    original, context = rbs_construct(strand)
    codon_profile = synthetic_profile()
    config = SearchConfig(seed=1, iterations=100)
    codon = optimize(original, config, codon_profile)
    synthesis = repair_synthesis(
        original,
        SynthesisProfile(homopolymer_limit=4),
        config,
        codon_profile=codon_profile,
        initial_sequence=codon.designs[0].sequence,
    )
    backend = CountingTestBackend()
    rbs = design_rbs(
        original,
        context,
        400,
        RBSDesignSettings(seed=3, iterations=80, joint_cds_nt=12),
        initial_sequence=synthesis.search.designs[0].sequence,
        backend=backend,
    )
    final = optimize(
        original,
        SearchConfig(
            seed=4,
            iterations=100,
            n_terminal_lock_nt=0,
            additional_locks=(rbs.initiation_locks[0],),
        ),
        codon_profile,
        initial_sequence=rbs.search.designs[0].sequence,
    )
    for sequence in (
        codon.designs[0].sequence,
        synthesis.search.designs[0].sequence,
        rbs.search.designs[0].sequence,
        final.designs[0].sequence,
    ):
        assert validate_candidate(original, sequence).valid
    assert predict_rbs(original, final.designs[0].sequence, context, backend) == rbs.predictions[0]
    assert final.reference_sha256 == codon.reference_sha256 == rbs.search.reference_sha256
    sequence = final.designs[0].sequence
    plan = plan_gibson(
        original,
        sequence,
        GibsonSettings(
            maximum_ordered_length=70,
            minimum_ordered_length=20,
            overlap_length=6,
        ),
    )
    assert plan.validation.valid
    assert len(plan.fragments) == 2
    assert reconstruct(plan) == sequence
    project = DesignProject(
        original=original,
        sequence=sequence,
        initiation_locks=(rbs.initiation_locks[0],),
        gibson=plan,
    )
    restored = load_project(save_draft(project))
    assert restored == project and validate_project(restored).valid
    assert predict_rbs(restored.original, restored.sequence, context, backend) == rbs.predictions[0]
    assert load_project(export_project(restored, "json")) == restored
    for format in ("fasta", "genbank"):
        assert (
            str(SeqIO.read(io.StringIO(export_project(restored, format)), format).seq) == sequence
        )
    rows = list(csv.DictReader(io.StringIO(export_project(restored, "csv"))))
    assert [row["sequence"] for row in rows] == [fragment.sequence for fragment in plan.fragments]
