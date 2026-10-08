"""Inserted synthetic defects and repairs; no empirical vendor-performance claims."""

import pytest
from pydantic import ValidationError

from dna_design import Construct, Feature, InvalidDesignError, MotifConstraint, validate_candidate
from dna_design.optimization import ObjectiveWeights, SearchConfig
from dna_design.synthesis import (
    FractionRange,
    SynthesisProfile,
    load_synthesis_profile,
    repair_synthesis,
    scan_synthesis,
)


def designable(sequence: str, *, topology: str = "linear") -> Construct:
    return Construct(
        sequence=sequence,
        topology=topology,
        features=(
            Feature(
                id="free",
                name="synthetic spacer",
                feature_type="spacer",
                start=0,
                end=len(sequence),
                mutability="FREE_MUTATION",
            ),
        ),
    )


@pytest.mark.parametrize(
    "sequence,profile,metric",
    [
        ("AAAAAAAAAACGT", SynthesisProfile(homopolymer_limit=4), "homopolymer"),
        (
            "GGGGGGAAAAAA",
            SynthesisProfile(global_gc_range=FractionRange(minimum=0.2, maximum=0.4)),
            "global_gc",
        ),
        (
            "GGGGGGAAAAAA",
            SynthesisProfile(gc_window=4, local_gc_range=FractionRange(minimum=0.2, maximum=0.8)),
            "local_gc",
        ),
        ("ACGTACGTACGT", SynthesisProfile(repeat_min_length=4), "direct_repeat"),
        ("AAAACCCCGGGGTTTT", SynthesisProfile(inverted_repeat_min_length=4), "inverted_repeat"),
        ("AAAAAAAAAACGT", SynthesisProfile(complexity_window=6, minimum_entropy=1.2), "complexity"),
    ],
)
def test_scanner_identifies_inserted_defects(
    sequence: str, profile: SynthesisProfile, metric: str
) -> None:
    reference = designable(sequence)
    report = scan_synthesis(reference, sequence, profile)
    assert any(f.metric == metric and f.penalty > 0 for f in report.findings)
    assert report.biological_validation.valid
    with pytest.raises(InvalidDesignError):
        report.require_policy()


def test_scanner_counts_origin_spanning_homopolymer_once() -> None:
    reference = designable("AAACGAAA", topology="circular")
    report = scan_synthesis(reference, reference.sequence, SynthesisProfile(homopolymer_limit=4))
    bad = [f for f in report.findings if f.metric == "homopolymer" and f.penalty]
    assert len(bad) == 1
    assert bad[0].observed == 6
    assert set(bad[0].positions) == {0, 1, 2, 5, 6, 7}
    entire = designable("AAAA", topology="circular")
    assert (
        scan_synthesis(entire, entire.sequence, SynthesisProfile(homopolymer_limit=2)).risk_loss
        == 2
    )


def test_generic_does_not_invent_thresholds_or_certify_vendor_acceptance() -> None:
    profile = load_synthesis_profile("generic")
    assert profile.homopolymer_limit is None
    reference = designable("ACGTACGT")
    report = scan_synthesis(reference, reference.sequence, profile)
    assert report.global_gc == 0.5
    assert all(f.observed is None for f in report.findings)
    with pytest.raises(InvalidDesignError, match="Unconfigured"):
        report.require_policy()
    assert load_synthesis_profile("genscript-public").repeat_min_length == 20
    twist = load_synthesis_profile("twist-public")
    assert twist.maximum_local_gc_spread == 0.5
    assert twist.global_gc_range == FractionRange(minimum=0.25, maximum=0.65)
    assert twist.source.startswith("https://www.twistbioscience.com/")
    assert "proprietary" in twist.approximation


def test_local_gc_spread_uses_complete_windows() -> None:
    reference = designable("G" * 50 + "A" * 50)
    report = scan_synthesis(reference, reference.sequence, load_synthesis_profile("twist-public"))
    finding = next(f for f in report.findings if f.metric == "local_gc_spread")
    assert finding.observed == 1
    assert finding.penalty == 0.5


def test_repair_reduces_risk_in_explicit_free_spacer_and_keeps_motifs() -> None:
    avoid = MotifConstraint(id="avoid", mode="AVOID", motif="GGGG", strand="forward")
    reference = designable("AAAAAAAAAACGT").model_copy(update={"constraints": (avoid,)})
    config = SearchConfig(seed=4, iterations=150, weights=ObjectiveWeights(codon_adaptation=0))
    result = repair_synthesis(reference, SynthesisProfile(homopolymer_limit=4), config)
    assert result.after[0].risk_loss < result.before.risk_loss
    assert result.after[0].risk_loss == 0
    assert validate_candidate(reference, result.search.designs[0].sequence).valid
    assert result == repair_synthesis(reference, SynthesisProfile(homopolymer_limit=4), config)


def test_synonymous_repair_validates_against_original_and_preserves_lock() -> None:
    sequence = "CCC" + "ATG" + "AAA" * 5 + "TAA"
    cds = Feature(id="cds", name="c", feature_type="CDS", start=3, end=len(sequence), strand=1)
    lock = Feature(id="lock", name="l", feature_type="protected region", start=0, end=3)
    reference = Construct(sequence=sequence, features=(lock, cds))
    config = SearchConfig(
        seed=7, iterations=180, n_terminal_lock_nt=0, weights=ObjectiveWeights(codon_adaptation=0)
    )
    result = repair_synthesis(reference, SynthesisProfile(homopolymer_limit=4), config)
    assert result.after[0].risk_loss < result.before.risk_loss
    designed = result.search.designs[0].sequence
    assert designed[:3] == "CCC"
    assert validate_candidate(reference, designed).valid
    subsequent = repair_synthesis(
        reference, SynthesisProfile(homopolymer_limit=4), config, initial_sequence=designed
    )
    assert validate_candidate(reference, subsequent.search.designs[0].sequence).valid


def test_immutable_risk_remains_visible() -> None:
    reference = Construct(sequence="AAAAAAAAAA")
    config = SearchConfig(seed=1, weights=ObjectiveWeights(codon_adaptation=0))
    result = repair_synthesis(reference, SynthesisProfile(homopolymer_limit=4), config)
    assert result.after[0].risk_loss == result.before.risk_loss == 6
    assert result.search.designs[0].sequence == reference.sequence
    with pytest.raises(InvalidDesignError):
        result.after[0].require_policy()


def test_incomplete_or_reversed_policy_is_rejected() -> None:
    with pytest.raises(ValidationError):
        SynthesisProfile(gc_window=20)
    with pytest.raises(ValidationError):
        SynthesisProfile(complexity_window=20)
    with pytest.raises(ValidationError):
        FractionRange(minimum=0.8, maximum=0.2)
