"""Local versioned data loaders and explicit held-out experimental evaluation.

No experimental observations are bundled. File splits and model-training overlap
are distinct; a supplied CSV cannot by itself prove independent model validation.
"""

import csv
import hashlib
import io
import math
import statistics
from collections import defaultdict
from typing import Annotated, Any, Literal, Self

from pydantic import Field, model_validator

from .models import FrozenModel, Identifier

PUBLIC_DATASETS = {
    "kudla2009": {
        "doi": "10.1126/science.1170160",
        "purpose": "Synonymous GFP expression ranking",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3902468/",
    },
    "boel2016": {
        "doi": "10.1038/nature16509",
        "purpose": "Heterologous protein expression ranking",
        "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5054687/",
        "data": "https://pmc.ncbi.nlm.nih.gov/articles/instance/5054687/bin/NIHMS806290-supplement-Protein_Data_Set.csv",
        "caution": "Published scores are ordinal; do not interpret them as fold changes.",
    },
    "salis2009": {
        "doi": "10.1038/nbt.1568",
        "purpose": "RBS model calibration or independently held-out subsets",
        "source": "https://doi.org/10.1038/nbt.1568",
        "caution": "OSTIR uses Salis2009 calibration; reuse is not independent validation.",
    },
    "farasat2014": {
        "doi": "10.15252/msb.20134955",
        "purpose": "Host-specific RBS validation",
        "source": "https://link.springer.com/article/10.15252/msb.20134955",
    },
    "johns2018": {
        "doi": "10.1038/nmeth.4633",
        "purpose": "Species-selective regulatory expression",
        "source": "https://doi.org/10.1038/nmeth.4633",
        "caution": "Regulatory activity includes transcription; it is not isolated TIR.",
    },
}


class ColumnMapping(FrozenModel):
    id: Identifier = "id"
    sequence: Identifier = "sequence"
    group: Identifier = "group"
    host: Identifier = "host"
    condition: Identifier = "condition"
    split: Identifier = "split"
    measured: Identifier = "measured"
    predicted: Identifier = "predicted"


class BenchmarkManifest(FrozenModel):
    schema_version: Literal[1] = 1
    dataset_id: Identifier
    dataset_version: Identifier
    source: Identifier
    license: Identifier
    license_source: Identifier
    accessed_date: Identifier
    csv_sha256: Annotated[str, Field(pattern=r"^[a-f0-9]{64}$")]
    evidence: Literal["experimental", "synthetic"]
    measurement: Identifier
    measured_units: Identifier
    prediction: Identifier
    predicted_units: Identifier
    model_version: Identifier
    split_policy: Identifier
    group_definition: Identifier
    transforms: tuple[str, ...] = ()
    exclusions: tuple[str, ...] = ()
    columns: ColumnMapping = Field(default_factory=ColumnMapping)
    scale: Literal["ordinal", "ratio"] = "ordinal"
    comparable_units: bool = False
    calibration: Literal["none", "multiplicative_log"] = "none"
    training_overlap: Literal["unknown", "verified_disjoint"] = "unknown"
    training_overlap_source: str | None = None
    training_sequence_sha256: tuple[str, ...] = ()
    training_groups: tuple[str, ...] = ()

    @model_validator(mode="after")
    def assumptions(self) -> Self:
        if self.comparable_units and self.measured_units != self.predicted_units:
            raise ValueError("Comparable units must agree explicitly")
        if self.calibration != "none" and self.scale != "ratio":
            raise ValueError("Multiplicative calibration requires ratio-scale measurements")
        if self.training_overlap == "verified_disjoint" and not self.training_overlap_source:
            raise ValueError(
                "Disjoint model-training claims require a recorded verification source"
            )
        if self.dataset_id == "boel2016" and self.scale != "ordinal":
            raise ValueError(
                "Boël expression scores are ordinal; fold-error interpretation is unsupported"
            )
        return self


class BenchmarkRow(FrozenModel):
    id: Identifier
    sequence: Annotated[str, Field(pattern=r"^[ACGT]+$")]
    group: Identifier
    host: Identifier
    condition: Identifier
    split: Literal["calibration", "heldout", "excluded"]
    measured: Annotated[float, Field(ge=0, allow_inf_nan=False)]
    predicted: Annotated[float, Field(ge=0, allow_inf_nan=False)]


def load_benchmark(text: str, manifest: BenchmarkManifest) -> tuple[BenchmarkRow, ...]:
    """Mapped public CSV v1; exact UTF-8 text hash and explicit column meanings."""
    manifest = BenchmarkManifest.model_validate(manifest.model_dump())
    if hashlib.sha256(text.encode("utf-8")).hexdigest() != manifest.csv_sha256:
        raise ValueError("Benchmark CSV checksum does not match its provenance manifest")
    reader = csv.DictReader(io.StringIO(text))
    mapping = manifest.columns.model_dump()
    if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ValueError("Benchmark CSV requires unique column headings")
    if not set(mapping.values()).issubset(reader.fieldnames):
        raise ValueError("CSV lacks explicitly mapped benchmark columns")
    rows = []
    for number, row in enumerate(reader, 2):
        if number > 100_001:
            raise ValueError("Benchmark loader supports up to 100,000 records")
        if None in row:
            raise ValueError(f"CSV row {number} contains extra unmapped fields")
        try:
            rows.append(
                BenchmarkRow.model_validate(
                    {field: row[column] for field, column in mapping.items()}
                )
            )
        except ValueError as exc:
            raise ValueError(
                f"Invalid benchmark row {number}; check mapped columns and values"
            ) from exc
    if not rows:
        raise ValueError("Benchmark CSV contains no records")
    if len({row.id for row in rows}) != len(rows):
        raise ValueError("Benchmark observation IDs must be unique; retain replicate groups")
    calibration = [row for row in rows if row.split == "calibration"]
    heldout = [row for row in rows if row.split == "heldout"]
    if not heldout:
        raise ValueError("Explicit held-out records are required")
    if {row.group for row in calibration} & {row.group for row in heldout} or (
        {row.sequence for row in calibration} & {row.sequence for row in heldout}
    ):
        raise ValueError("Calibration and held-out groups/sequences overlap; leakage is prohibited")
    for observation in heldout:
        digest = hashlib.sha256(observation.sequence.encode()).hexdigest()
        if (
            digest in manifest.training_sequence_sha256
            or observation.group in manifest.training_groups
        ):
            raise ValueError("Held-out records overlap declared model-training data")
    return tuple(rows)


def _correlation(x: list[float], y: list[float]) -> float | None:
    if len(x) < 2 or len(set(x)) < 2 or len(set(y)) < 2:
        return None
    return statistics.correlation(x, y)


def _ranks(values: list[float]) -> list[float]:
    positions: dict[float, list[int]] = defaultdict(list)
    for index, value in enumerate(sorted(values), 1):
        positions[value].append(index)
    return [statistics.mean(positions[value]) for value in values]


def _metrics(
    rows: list[BenchmarkRow], manifest: BenchmarkManifest, log_scale: float
) -> dict[str, Any]:
    measured = [row.measured for row in rows]
    predicted = [row.predicted for row in rows]
    positive = [row for row in rows if row.measured > 0 and row.predicted > 0]
    logs_m = [math.log(row.measured) for row in positive]
    logs_p = [math.log(row.predicted) + log_scale for row in positive]
    errors = [abs(a - b) for a, b in zip(logs_m, logs_p, strict=True)]
    folds_assessed = manifest.scale == "ratio" and (
        manifest.comparable_units or manifest.calibration != "none"
    )
    median_log = statistics.median(errors) if errors and folds_assessed else None
    return {
        "records": len(rows),
        "positive_pairs": len(positive),
        "nonpositive_pairs_excluded_from_logs": len(rows) - len(positive),
        "spearman": _correlation(_ranks(measured), _ranks(predicted)),
        "pearson_linear": _correlation(measured, predicted),
        "pearson_log": _correlation(logs_m, logs_p) if manifest.scale == "ratio" else None,
        "median_log_fold_error": median_log,
        "median_fold_error": math.exp(median_log)
        if median_log is not None and median_log < 709
        else None,
        "fractions_within_fold": {
            str(limit): sum(error <= math.log(limit) + 1e-12 for error in errors) / len(errors)
            for limit in (2, 5, 10)
        }
        if errors and folds_assessed
        else None,
        "fold_metrics_assessed": folds_assessed,
    }


def evaluate_benchmark(text: str, manifest: BenchmarkManifest) -> dict[str, Any]:
    rows = load_benchmark(text, manifest)
    heldout = [row for row in rows if row.split == "heldout"]
    calibration = [row for row in rows if row.split == "calibration"]
    calibration_groups: dict[tuple[str, str], list[BenchmarkRow]] = defaultdict(list)
    groups: dict[tuple[str, str], list[BenchmarkRow]] = defaultdict(list)
    for row in calibration:
        calibration_groups[(row.host, row.condition)].append(row)
    for row in heldout:
        groups[(row.host, row.condition)].append(row)
    reports = []
    for (host, condition), selected in sorted(groups.items()):
        calibrators = calibration_groups[(host, condition)]
        log_scale = 0.0
        if manifest.calibration != "none":
            if len(calibrators) < 2 or any(
                row.measured <= 0 or row.predicted <= 0 for row in calibrators
            ):
                raise ValueError(
                    "Each host/condition requires at least two positive calibration records"
                )
            log_scale = statistics.mean(
                math.log(row.measured) - math.log(row.predicted) for row in calibrators
            )
        reports.append(
            {
                "host": host,
                "condition": condition,
                "calibration_records": len(calibrators),
                "calibration_log_scale": log_scale,
                "metrics": _metrics(selected, manifest, log_scale),
            }
        )
    warnings = [
        "A held-out file split does not establish independence from upstream model calibration",
        "Correlations assess association; they do not prove causal expression improvement",
    ]
    if manifest.evidence == "synthetic":
        warnings.append("Synthetic data validate arithmetic only; no experimental accuracy claim")
    if manifest.training_overlap == "unknown":
        warnings.append(
            "Model-training overlap is unknown; independent validation cannot be claimed"
        )
    if manifest.dataset_id == "salis2009" and "ostir" in manifest.model_version.lower():
        warnings.append(
            "OSTIR uses Salis2009 upstream calibration; verify exact excluded training records"
        )
    if any(row.measured == 0 or row.predicted == 0 for row in heldout):
        warnings.append(
            "Zero pairs contribute to rankings but are excluded from log/fold metrics; "
            "no pseudocounts"
        )
    return {
        "report_version": "benchmark-evaluation-v1",
        "manifest": manifest.model_dump(mode="json"),
        "loader_version": "public-mapped-csv-v1",
        "heldout_records": len(heldout),
        "calibration_records": len(calibration),
        "excluded_records": sum(row.split == "excluded" for row in rows),
        "groups": reports,
        "warnings": warnings,
        "independence_status": "user_documented_disjoint"
        if manifest.training_overlap == "verified_disjoint"
        else "unverified",
    }


class BenchmarkRequest(FrozenModel):
    csv: Annotated[str, Field(min_length=1, max_length=2_000_000)]
    manifest: BenchmarkManifest
