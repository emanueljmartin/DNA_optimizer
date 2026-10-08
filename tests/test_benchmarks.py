"""Synthetic arithmetic/provenance tests; never experimental benchmark evidence."""

import hashlib

import pytest
from pydantic import ValidationError

from dna_design.benchmarks import BenchmarkManifest, evaluate_benchmark, load_benchmark


def manifest(text, **updates):
    return BenchmarkManifest.model_validate(
        {
            "dataset_id": "synthetic-arithmetic",
            "dataset_version": "v1",
            "source": "Synthetic test fixture",
            "license": "Project-generated synthetic observations",
            "license_source": "tests/test_benchmarks.py",
            "accessed_date": "2026-10-08",
            "csv_sha256": hashlib.sha256(text.encode()).hexdigest(),
            "evidence": "synthetic",
            "measurement": "Synthetic value",
            "measured_units": "test-unit",
            "prediction": "Fictional test values",
            "predicted_units": "test-unit",
            "model_version": "test-double-v1",
            "split_policy": "Explicit synthetic disjoint groups",
            "group_definition": "One group per synthetic sequence",
            "scale": "ratio",
            "comparable_units": True,
            **updates,
        }
    )


HEADER = "id,sequence,group,host,condition,split,measured,predicted\n"
TEXT = (
    HEADER + "a,AAA,a,test,c,heldout,1,2\nb,AAC,b,test,c,heldout,2,4\nc,AAG,c,test,c,heldout,4,8\n"
)


def test_exact_metrics_and_explicit_evidence_provenance():
    report = evaluate_benchmark(TEXT, manifest(TEXT))
    metrics = report["groups"][0]["metrics"]
    assert metrics["spearman"] == pytest.approx(1)
    assert metrics["pearson_log"] == pytest.approx(1)
    assert metrics["median_fold_error"] == pytest.approx(2)
    assert metrics["fractions_within_fold"] == {"2": 1, "5": 1, "10": 1}
    assert report["independence_status"] == "unverified"
    assert report["manifest"]["evidence"] == "synthetic"
    assert any("Synthetic" in warning for warning in report["warnings"])


def test_tied_ranks_zeros_and_constant_predictions_are_reported_without_pseudocounts():
    text = (
        HEADER
        + "a,AAA,a,test,c,heldout,0,1\nb,AAC,b,test,c,heldout,1,1\nc,AAG,c,test,c,heldout,2,3\n"
    )
    metrics = evaluate_benchmark(text, manifest(text))["groups"][0]["metrics"]
    assert metrics["spearman"] == pytest.approx(0.8660254037844387)
    assert metrics["nonpositive_pairs_excluded_from_logs"] == 1
    constant = text.replace(",3\n", ",1\n")
    assert (
        evaluate_benchmark(constant, manifest(constant))["groups"][0]["metrics"]["spearman"] is None
    )


def test_multiplicative_calibration_only_uses_separate_positive_training_rows():
    text = TEXT + "d,ACA,d,test,c,calibration,10,20\ne,ACC,e,test,c,calibration,20,40\n"
    result = evaluate_benchmark(text, manifest(text, calibration="multiplicative_log"))
    assert result["calibration_records"] == 2
    assert result["groups"][0]["metrics"]["median_fold_error"] == pytest.approx(1)
    assert result["heldout_records"] == 3


@pytest.mark.parametrize(
    "extra", ["d,AAA,d,test,c,calibration,10,20\n", "d,ACA,a,test,c,calibration,10,20\n"]
)
def test_replicate_groups_and_sequence_leakage_are_rejected(extra):
    text = TEXT + extra
    with pytest.raises(ValueError, match="overlap"):
        load_benchmark(text, manifest(text))


def test_checksums_nonfinite_values_training_overlap_and_duplicate_ids_are_rejected():
    with pytest.raises(ValueError, match="checksum"):
        load_benchmark(TEXT + "\n", manifest(TEXT))
    with pytest.raises(ValueError, match="training"):
        load_benchmark(TEXT, manifest(TEXT, training_groups=("a",)))
    broken = TEXT.replace(",1,2\n", ",nan,2\n")
    with pytest.raises(ValueError, match="row 2"):
        load_benchmark(broken, manifest(broken))
    duplicate = TEXT.replace("b,AAC", "a,AAC")
    with pytest.raises(ValueError, match="IDs"):
        load_benchmark(duplicate, manifest(duplicate))


def test_ordinal_scores_and_incomparable_units_never_report_fold_error():
    result = evaluate_benchmark(
        TEXT, manifest(TEXT, dataset_id="boel2016", scale="ordinal", comparable_units=False)
    )
    metrics = result["groups"][0]["metrics"]
    assert metrics["spearman"] == pytest.approx(1)
    assert metrics["pearson_log"] is None and metrics["median_fold_error"] is None
    with pytest.raises(ValidationError, match="ordinal"):
        manifest(TEXT, dataset_id="boel2016")
    with pytest.raises(ValidationError, match="units"):
        manifest(TEXT, predicted_units="different")


def test_host_condition_groups_are_not_pooled_and_unknown_upstream_overlap_is_visible():
    mixed = TEXT + "d,ACA,d,other,c,heldout,20,1\ne,ACC,e,other,c,heldout,10,2\n"
    report = evaluate_benchmark(mixed, manifest(mixed))
    assert len(report["groups"]) == 2
    assert report["groups"][0]["metrics"]["spearman"] == pytest.approx(-1)
    assert report["groups"][1]["metrics"]["spearman"] == pytest.approx(1)


def test_mapped_public_columns_and_zero_exclusions_remain_auditable():
    public_shape = TEXT.replace("measured", "published score").replace("predicted", "model value")
    mapped = manifest(
        public_shape, columns={"measured": "published score", "predicted": "model value"}
    )
    assert len(load_benchmark(public_shape, mapped)) == 3
    missing = TEXT.replace("id,sequence", "id,other")
    with pytest.raises(ValueError, match="columns"):
        load_benchmark(missing, manifest(missing))
