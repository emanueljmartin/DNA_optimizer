"""Synthetic design invariants; these tests are not expression benchmarks."""

import math

import pytest
from Bio.Seq import Seq
from hypothesis import given, settings
from hypothesis import strategies as st
from pydantic import ValidationError

from dna_design import Construct, Feature, MotifConstraint, validate_candidate
from dna_design.metrics import codon_metrics, gc_fraction, local_gc
from dna_design.optimization import (
    AdditionalLock,
    ObjectiveWeights,
    SearchConfig,
    optimize,
    validate_search_candidate,
)
from dna_design.profiles import CodonProfile, load_codon_profile, relative_codon_weights


def synthetic_profile() -> CodonProfile:
    data = load_codon_profile("ecoli").model_dump()
    data.update(
        id="synthetic-test",
        version="test-v1",
        reference_kind="custom",
        source_accession="synthetic-test",
        source_url="https://example.org/synthetic-test",
        source_sha256="0" * 64,
        included_cds=1,
        excluded_cds={},
        counts={codon: 1 for codon in data["counts"]},
    )
    data["counts"]["GAA"] = 1
    data["counts"]["GAG"] = 10
    return CodonProfile.model_validate(data)


@pytest.mark.parametrize(
    "host,accession,included", [("ecoli", "NC_000913.3", 4294), ("pao1", "NC_002516.2", 5570)]
)
def test_packaged_profiles_have_versioned_provenance(
    host: str, accession: str, included: int
) -> None:
    profile = load_codon_profile(host)
    assert profile.source_accession == accession
    assert profile.included_cds == included
    assert profile.reference_kind == "genome_wide"
    assert len(profile.counts) == 64
    assert len(profile.fingerprint) == 64
    assert all(0 <= weight <= 1 for weight in relative_codon_weights(profile).values())
    assert CodonProfile.model_validate_json(profile.model_dump_json()) == profile


def test_profile_unknown_host_and_invalid_counts_are_rejected() -> None:
    with pytest.raises(ValueError, match="Unknown"):
        load_codon_profile("unlisted")
    data = synthetic_profile().model_dump()
    data["counts"].pop("AAA")
    with pytest.raises(ValidationError, match="64"):
        CodonProfile.model_validate(data)


def test_genome_adaptation_is_not_high_expression_cai(coding_construct: Construct) -> None:
    profile = synthetic_profile()
    metrics = codon_metrics(
        coding_construct, coding_construct.sequence, profile, rare_threshold=0.2
    )[0]
    assert metrics.adaptation_score == pytest.approx(0.1)
    assert metrics.cai is None
    assert metrics.rare_codon_indices == (1,)
    assert metrics.rare_fraction == 1
    assert metrics.codon_counts == {"GAA": 1}
    assert metrics.warnings
    assert codon_metrics(coding_construct, "ATGGAGTAA", profile)[0].adaptation_score == 1
    zero = profile.model_dump()
    zero["counts"]["GAA"] = 0
    assert (
        codon_metrics(
            coding_construct, coding_construct.sequence, CodonProfile.model_validate(zero)
        )[0].adaptation_score
        == 0
    )
    high = profile.model_dump()
    high["reference_kind"] = "high_expression"
    assert codon_metrics(
        coding_construct, coding_construct.sequence, CodonProfile.model_validate(high)
    )[0].cai == pytest.approx(0.1)


def test_codon_metrics_respect_genetic_code(coding_construct: Construct) -> None:
    feature = coding_construct.features[0].model_copy(update={"translation_table": 4})
    reference = Construct(sequence=coding_construct.sequence, features=(feature,))
    with pytest.raises(ValueError, match="genetic code"):
        codon_metrics(reference, reference.sequence, synthetic_profile())


def test_gc_statistics_and_explicit_parameters() -> None:
    assert gc_fraction("AGGC") == 0.75
    assert local_gc("AGGC", 2) == ((0, 0.5), (1, 1.0), (2, 1.0))
    assert local_gc("AC", 3) == ()
    with pytest.raises(ValidationError, match="explicit"):
        SearchConfig(seed=1, weights=ObjectiveWeights(global_gc=1))
    with pytest.raises(ValidationError, match="explicit"):
        SearchConfig(seed=1, weights=ObjectiveWeights(rare_codons=1))


def test_seeded_search_improves_and_returns_distinct_valid_alternatives(
    coding_construct: Construct,
) -> None:
    config = SearchConfig(seed=22, iterations=80, n_terminal_lock_nt=0, alternatives=3)
    first = optimize(coding_construct, config, synthetic_profile())
    assert first == optimize(coding_construct, config, synthetic_profile())
    assert first.status == "success"
    assert len({d.sequence for d in first.designs}) == 3
    assert first.designs[0].objective_loss == 0
    assert all(validate_candidate(coding_construct, d.sequence).valid for d in first.designs)
    assert first.designs[0].changes
    assert first.profile_fingerprint == synthetic_profile().fingerprint


def test_search_repairs_forbidden_site_before_scoring(coding_construct: Construct) -> None:
    avoid = MotifConstraint(id="avoid", motif="GAA", mode="AVOID", strand="forward")
    reference = coding_construct.model_copy(update={"constraints": (avoid,)})
    result = optimize(
        reference, SearchConfig(seed=2, iterations=100, n_terminal_lock_nt=0), synthetic_profile()
    )
    assert result.status == "success"
    assert "GAA" not in result.designs[0].sequence
    assert result.designs[0].validation.valid


def test_proven_infeasibility_is_distinct_from_search_exhaustion() -> None:
    avoid = MotifConstraint(id="avoid", motif="AAA", mode="AVOID", strand="forward")
    impossible = Construct(sequence="AAA", constraints=(avoid,))
    config = SearchConfig(seed=1, weights=ObjectiveWeights(codon_adaptation=0))
    assert optimize(impossible, config).status == "infeasible"
    free = Feature(
        id="f", name="f", feature_type="spacer", start=0, end=3, mutability="FREE_MUTATION"
    )
    uncertain = Construct(sequence="AAA", features=(free,), constraints=(avoid,))
    exhausted = optimize(uncertain, config.model_copy(update={"iterations": 0}))
    assert exhausted.status == "budget_exhausted"
    assert not exhausted.designs
    assert not exhausted.validation.infeasible


def test_search_does_not_weaken_overlap_or_locked_controls(cds: Feature) -> None:
    overlap = Feature(
        id="overlap",
        name="overlap",
        feature_type="CDS",
        start=1,
        end=7,
        strand=1,
        cds_mode="fragment",
    )
    reference = Construct(sequence="ATGGAATAA", features=(cds, overlap))
    result = optimize(
        reference, SearchConfig(seed=1, iterations=60, n_terminal_lock_nt=0), synthetic_profile()
    )
    assert all(d.sequence[3:6] == "GAA" for d in result.designs)
    locked = Construct(
        sequence="ATGGAATAA",
        features=(cds, Feature(id="l", name="l", feature_type="protected region", start=0, end=9)),
    )
    result = optimize(locked, SearchConfig(seed=1, iterations=20), synthetic_profile())
    assert result.designs[0].sequence == locked.sequence
    assert result.proposals_evaluated == 0


def test_default_initiation_lock_and_additional_design_lock(coding_construct: Construct) -> None:
    result = optimize(coding_construct, SearchConfig(seed=1, iterations=20), synthetic_profile())
    assert result.designs[0].sequence == coding_construct.sequence
    config = SearchConfig(
        seed=1,
        iterations=40,
        n_terminal_lock_nt=0,
        additional_locks=(
            AdditionalLock(start=3, end=6, sequence="GAG", reason="designed initiation"),
        ),
    )
    result = optimize(coding_construct, config, synthetic_profile(), initial_sequence="ATGGAGTAA")
    assert all(d.sequence[3:6] == "GAG" for d in result.designs)
    assert optimize(coding_construct, config, synthetic_profile()).status == "invalid_reference"


def test_multiobjective_gc_tradeoff_and_custom_scorer(coding_construct: Construct) -> None:
    config = SearchConfig(
        seed=1,
        iterations=50,
        n_terminal_lock_nt=0,
        gc_target=2 / 9,
        weights=ObjectiveWeights(codon_adaptation=1, global_gc=10),
    )
    result = optimize(coding_construct, config, synthetic_profile())
    assert set(result.designs[0].components) == {"codon_adaptation", "global_gc"}
    custom = config.model_copy(update={"weights": ObjectiveWeights(codon_adaptation=0, custom=1)})
    with pytest.raises(ValueError, match="scorer"):
        optimize(coding_construct, custom)
    with pytest.raises(ValueError, match="finite"):
        optimize(coding_construct, custom, custom_objective=lambda _: math.nan)


def test_invalid_reference_schema_is_a_structured_search_failure(
    coding_construct: Construct,
) -> None:
    malformed = coding_construct.model_copy(update={"rbs_settings": {"invalid": object()}})
    result = optimize(malformed, SearchConfig(seed=1), synthetic_profile())
    assert result.status == "invalid_reference"
    assert result.reference_sha256 is None
    assert not result.designs


def test_rare_objective_requires_a_reference_even_without_codon_weight() -> None:
    config = SearchConfig(
        seed=1, rare_threshold=0.2, weights=ObjectiveWeights(codon_adaptation=0, rare_codons=1)
    )
    with pytest.raises(ValueError, match="profile"):
        optimize(Construct(sequence="AAA"), config)


def test_partial_lock_uses_the_valid_starting_codon_instead_of_original() -> None:
    reference = Construct(
        sequence="ATGTTATAA",
        features=(Feature(id="c", name="c", feature_type="CDS", start=0, end=9, strand=1),),
    )
    data = synthetic_profile().model_dump()
    data["counts"]["CTC"] = 100
    config = SearchConfig(
        seed=2,
        iterations=100,
        n_terminal_lock_nt=0,
        additional_locks=(AdditionalLock(start=3, end=4, sequence="C", reason="partial lock"),),
    )
    result = optimize(
        reference, config, CodonProfile.model_validate(data), initial_sequence="ATGCTTTAA"
    )
    assert result.designs[0].sequence[3:6] == "CTC"
    assert validate_candidate(reference, result.designs[0].sequence).valid


def test_fresh_search_validation_enforces_region_and_initiation_settings(
    coding_construct: Construct,
) -> None:
    synonymous = "ATGGAGTAA"
    assert validate_candidate(coding_construct, synonymous).valid
    for config in (
        SearchConfig(seed=1),
        SearchConfig(seed=1, n_terminal_lock_nt=0, allowed_positions=(6, 7, 8)),
    ):
        report = validate_search_candidate(coding_construct, synonymous, config)
        assert not report.valid
        assert any(
            c.code == "SEARCH_REGION_PROTECTION" and c.positions == (5,) for c in report.checks
        )


@settings(max_examples=15, deadline=None)
@given(st.integers(0, 2**16), st.sampled_from((1, -1)))
def test_seeded_optimizer_preserves_both_strands(seed: int, strand: int) -> None:
    sequence = "ATG" + "GAA" * 4 + "TAA"
    if strand == -1:
        sequence = str(Seq(sequence).reverse_complement())
    feature = Feature(
        id="c", name="c", feature_type="CDS", start=0, end=len(sequence), strand=strand
    )
    reference = Construct(sequence=sequence, features=(feature,))
    result = optimize(
        reference, SearchConfig(seed=seed, iterations=30, n_terminal_lock_nt=0), synthetic_profile()
    )
    assert result.designs
    assert all(validate_candidate(reference, d.sequence).valid for d in result.designs)
