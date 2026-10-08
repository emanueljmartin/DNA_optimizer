"""RBS control tests with a named test double and real local backend integration."""

import importlib.util
import math
from collections import Counter

import pytest
from Bio.Seq import Seq
from pydantic import ValidationError

from dna_design import Construct, Feature, InvalidDesignError, Region, validate_candidate
from dna_design.optimization import SearchConfig, optimize
from dna_design.rbs import (
    OSTIRBackend,
    RBSBackendUnavailable,
    RBSContext,
    RBSDesignSettings,
    RBSHost,
    RBSModelError,
    RBSPrediction,
    design_rbs,
    design_rbs_library,
    predict_rbs,
)

from .test_optimization import synthetic_profile

HAS_LOCAL_BACKEND = (
    importlib.util.find_spec("ostir") is not None and importlib.util.find_spec("RNA") is not None
)


def rbs_construct(strand: int = 1) -> tuple[Construct, RBSContext]:
    upstream = "AACCAACCAACCAACCAACC" + "AAAGGAGGAAA" + "AACC"
    dna = upstream + "ATG" + "GAA" * 15 + "TAA"
    cds_start, cds_end = len(upstream), len(dna)
    rbs_start, rbs_end = 20, 31
    region = Region(start=0, end=len(upstream))
    if strand == -1:
        dna = str(Seq(dna).reverse_complement())
        cds_start, cds_end = len(dna) - cds_end, len(dna) - cds_start
        rbs_start, rbs_end = len(dna) - rbs_end, len(dna) - rbs_start
        region = Region(start=cds_end, end=len(dna))
    cds = Feature(
        id="cds",
        name="synthetic CDS",
        feature_type="CDS",
        start=cds_start,
        end=cds_end,
        strand=strand,
    )
    rbs = Feature(
        id="rbs",
        name="synthetic RBS",
        feature_type="RBS",
        start=rbs_start,
        end=rbs_end,
        strand=strand,
        mutability="DESIGNABLE_RBS",
        rbs_mask="NNAGGAGGNNN",
    )
    return Construct(sequence=dna, features=(rbs, cds)), RBSContext(
        cds_id="cds", rbs_id="rbs", upstream_region=region
    )


class CountingTestBackend:
    """Deliberately fictional unit-test scorer; no biological prediction claims."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, str]] = []

    def predict(
        self, upstream: str, downstream: str, host: RBSHost, *, allow_uncalibrated: bool = False
    ) -> RBSPrediction:
        self.calls.append((upstream, downstream))
        return RBSPrediction(
            relative_tir=100 + 10 * upstream.count("G") + 5 * downstream.count("G"),
            start_codon=downstream[:3],
            spacing_bp=6,
            energy_components={"test_only": 0},
            unbound_mrna_structure="." * (len(upstream) + len(downstream)),
            unbound_mrna_mfe_kcal_mol=0,
            model_versions={"test_double": "v1"},
            host=host,
            warnings=("Fictional test double; not a scientific model",),
        )


@pytest.mark.parametrize("strand", [1, -1])
def test_masked_design_preserves_all_original_biology_and_improves_target(strand: int) -> None:
    reference, context = rbs_construct(strand)
    backend = CountingTestBackend()
    result = design_rbs(
        reference, context, 400, RBSDesignSettings(seed=4, iterations=100), backend=backend
    )
    design = result.search.designs[0]
    assert design.validation.valid
    assert validate_candidate(reference, design.sequence).valid
    initial_error = abs(math.log(result.initial_prediction.relative_tir / 400))
    assert result.log_fold_errors[0] < initial_error
    assert (
        max(Counter(backend.calls).values()) >= 2
    )  # Final prediction did not trust the score cache.
    assert all(
        i in result.search.config.allowed_positions
        for i, (a, b) in enumerate(zip(reference.sequence, design.sequence, strict=True))
        if a != b
    )
    assert (
        result.search.designs[0].sequence[reference.features[1].start : reference.features[1].end]
        == reference.sequence[reference.features[1].start : reference.features[1].end]
    )
    lock = result.initiation_locks[0]
    later = optimize(
        reference,
        SearchConfig(seed=2, iterations=50, n_terminal_lock_nt=0, additional_locks=(lock,)),
        synthetic_profile(),
        initial_sequence=design.sequence,
    )
    assert later.designs[0].sequence[lock.start : lock.end] == lock.sequence
    assert validate_candidate(reference, later.designs[0].sequence).valid


def test_joint_cds_design_stays_within_requested_region_and_preserves_protein() -> None:
    reference, context = rbs_construct()
    result = design_rbs(
        reference,
        context,
        400,
        RBSDesignSettings(seed=8, iterations=140, joint_cds_nt=12),
        backend=CountingTestBackend(),
    )
    cds = reference.features[1]
    for design in result.search.designs:
        assert validate_candidate(reference, design.sequence).valid
        assert design.sequence[cds.start + 12 :] == reference.sequence[cds.start + 12 :]
        assert design.sequence[cds.start : cds.start + 3] == "ATG"
    assert result.initiation_locks[0].end == cds.start + context.downstream_nt


@pytest.mark.parametrize("strand", [1, -1])
def test_carried_lock_protects_every_base_used_by_prediction(strand: int) -> None:
    reference, context = rbs_construct(strand)
    backend = CountingTestBackend()
    result = design_rbs(
        reference, context, 400, RBSDesignSettings(seed=5, iterations=40), backend=backend
    )
    original_prediction = result.predictions[0]
    later = optimize(
        reference,
        SearchConfig(
            seed=3,
            iterations=160,
            n_terminal_lock_nt=0,
            additional_locks=(result.initiation_locks[0],),
        ),
        synthetic_profile(),
        initial_sequence=result.search.designs[0].sequence,
    )
    assert (
        predict_rbs(reference, later.designs[0].sequence, context, backend) == original_prediction
    )
    assert any(
        a != b
        for a, b in zip(result.search.designs[0].sequence, later.designs[0].sequence, strict=True)
    )


def test_original_reference_gate_runs_before_prediction() -> None:
    reference, context = rbs_construct()
    backend = CountingTestBackend()
    cds = reference.features[1]
    invalid = reference.sequence[: cds.start + 3] + "GAC" + reference.sequence[cds.start + 6 :]
    with pytest.raises(InvalidDesignError):
        predict_rbs(reference, invalid, context, backend)
    assert not backend.calls


def test_rbs_library_returns_distinct_rbs_sequences_without_rewriting_search_history() -> None:
    reference, context = rbs_construct()
    results = design_rbs_library(
        reference,
        context,
        (300.0, 400.0),
        RBSDesignSettings(seed=7, iterations=90),
        backend=CountingTestBackend(),
    )
    rbs = reference.features[0]
    selected = [result.search.designs[result.selected_design_index] for result in results]
    assert len({design.sequence[rbs.start : rbs.end] for design in selected}) == 2
    assert all(validate_candidate(reference, design.sequence).valid for design in selected)
    assert results[1].search.config.alternatives == 2
    with pytest.raises(ValueError, match="distinct"):
        design_rbs_library(
            reference, context, (300, 300), RBSDesignSettings(seed=1), backend=CountingTestBackend()
        )


def test_pao1_parameters_are_explicit_and_uncalibrated_use_is_opt_in() -> None:
    with pytest.raises(ValidationError, match="explicit anti-SD"):
        RBSHost(id="pseudomonas-pao1", calibrated_ecoli=False)
    with pytest.raises(ValidationError, match="cannot label"):
        RBSHost(
            id="pseudomonas-pao1",
            anti_sd="ACCTCCTTA",
            source="synthetic-test",
            calibrated_ecoli=True,
        )
    host = RBSHost(
        id="pseudomonas-pao1",
        anti_sd="ACCTCCTTA",
        source="Synthetic parameter fixture; not a PAO1 calibration",
        calibrated_ecoli=False,
    )
    with pytest.raises(RBSModelError, match="allow_uncalibrated"):
        OSTIRBackend().predict("AGGAGGAAAAAAA", "ATGGAATAA", host)


def test_missing_backend_does_not_fall_back_to_a_heuristic(monkeypatch: pytest.MonkeyPatch) -> None:
    from importlib.metadata import PackageNotFoundError

    import dna_design.rbs as rbs_module

    def missing(_: str) -> str:
        raise PackageNotFoundError("optional-test")

    monkeypatch.setattr(rbs_module, "version", missing)
    with pytest.raises(RBSBackendUnavailable, match=r"\[rbs\]"):
        OSTIRBackend().predict("AGGAGGAAAAAAA", "ATGGAATAA", RBSHost())


@pytest.mark.skipif(
    not HAS_LOCAL_BACKEND, reason="Optional local rbs dependencies are not installed"
)
@pytest.mark.parametrize("strand", [1, -1])
def test_real_ostir_prediction_reports_energy_structure_and_versions(strand: int) -> None:
    reference, context = rbs_construct(strand)
    prediction = predict_rbs(reference, reference.sequence, context)
    assert prediction.relative_tir > 0
    assert prediction.model_versions == {"ostir": "1.1.3", "ViennaRNA": "2.7.2"}
    assert "dG_total" in prediction.energy_components
    assert len(prediction.unbound_mrna_structure) == 70
    assert any("not experimentally measured" in note for note in prediction.warnings)
    counterpart, opposite = rbs_construct(-strand)
    assert prediction.relative_tir == pytest.approx(
        predict_rbs(counterpart, counterpart.sequence, opposite).relative_tir
    )


@pytest.mark.skipif(
    not HAS_LOCAL_BACKEND, reason="Optional local rbs dependencies are not installed"
)
def test_real_local_rbs_design_preserves_masks_and_validates_final_result() -> None:
    reference, context = rbs_construct()
    baseline = predict_rbs(reference, reference.sequence, context)
    result = design_rbs(
        reference, context, baseline.relative_tir / 2, RBSDesignSettings(seed=12, iterations=30)
    )
    assert result.search.designs
    assert result.log_fold_errors[0] <= abs(math.log(2))
    assert validate_candidate(reference, result.search.designs[0].sequence).valid
    assert result.predictions[0] == predict_rbs(
        reference, result.search.designs[0].sequence, context
    )


def test_bad_context_and_non_designable_rbs_fail_explicitly() -> None:
    reference, context = rbs_construct()
    wrong = context.model_copy(update={"upstream_region": Region(start=0, end=34)})
    with pytest.raises(RBSModelError, match="end at"):
        predict_rbs(reference, reference.sequence, wrong, CountingTestBackend())
    rbs = reference.features[0].model_copy(update={"mutability": "LOCKED", "rbs_mask": None})
    locked = Construct(sequence=reference.sequence, features=(rbs, reference.features[1]))
    with pytest.raises(RBSModelError, match="DESIGNABLE_RBS"):
        design_rbs(locked, context, 100, RBSDesignSettings(seed=1), backend=CountingTestBackend())


def test_library_keeps_previously_designed_cds_as_the_starting_sequence() -> None:
    reference, context = rbs_construct()
    initial = reference.sequence[:82] + "G" + reference.sequence[83:]
    assert initial != reference.sequence
    assert validate_candidate(reference, initial).valid
    results = design_rbs_library(
        reference,
        context,
        (200,),
        RBSDesignSettings(seed=3, iterations=0),
        backend=CountingTestBackend(),
        initial_sequence=initial,
    )
    assert results[0].search.designs[0].sequence == initial
