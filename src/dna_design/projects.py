"""Reproducible project snapshots and mandatory fresh export validation."""

import csv
import io
import json
import re
from typing import Any, Literal, Self

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pydantic import Field, JsonValue, model_validator

from .gibson import GibsonPlan, validate_gibson
from .models import Construct, FrozenModel
from .optimization import AdditionalLock, SearchConfig, validate_search_candidate
from .synthesis import SynthesisProfile, scan_synthesis
from .translation import translate_feature
from .validation import CheckStatus, Finding, ValidationReport


class DesignProject(FrozenModel):
    schema_version: Literal[1] = 1
    original: Construct
    sequence: str | None = None
    initiation_locks: tuple[AdditionalLock, ...] = ()
    synthesis_policy: SynthesisProfile | None = None
    require_synthesis_policy: bool = False
    gibson: GibsonPlan | None = None
    settings: dict[str, JsonValue] = Field(default_factory=dict)
    history: tuple[dict[str, JsonValue], ...] = ()
    model_versions: dict[str, str] = Field(default_factory=dict)

    @model_validator(mode="after")
    def initial_sequence(self) -> Self:
        if self.sequence is None:
            object.__setattr__(self, "sequence", self.original.sequence)
        return self


def validate_project(project: DesignProject) -> ValidationReport:
    project = DesignProject.model_validate(project.model_dump())
    assert project.sequence is not None
    # N-terminal search policy is not a new original control. Carried RBS locks are.
    config = SearchConfig(seed=0, n_terminal_lock_nt=0, additional_locks=project.initiation_locks)
    report = validate_search_candidate(
        project.original, project.sequence, config, baseline_sequence=project.sequence
    )
    checks = list(report.checks)
    if project.require_synthesis_policy:
        if project.synthesis_policy is None:
            checks.append(
                Finding(
                    status=CheckStatus.ERROR,
                    code="SYNTHESIS_POLICY_MISSING",
                    message="Required synthesis policy is not configured",
                )
            )
        elif report.valid:
            synthesis = scan_synthesis(project.original, project.sequence, project.synthesis_policy)
            passed = synthesis.risk_loss == 0 and all(
                f.observed is not None for f in synthesis.findings
            )
            checks.append(
                Finding(
                    status=CheckStatus.PASS if passed else CheckStatus.ERROR,
                    code="SYNTHESIS_POLICY",
                    message="Selected synthesis criteria must be assessed and passing",
                )
            )
    if project.gibson is not None and report.valid:
        checks.extend(validate_gibson(project.original, project.sequence, project.gibson).checks)
    return ValidationReport(checks=tuple(checks))


def _record(project: DesignProject) -> SeqRecord:
    assert project.sequence is not None
    original = project.original
    identifier = re.sub(r"[^A-Za-z0-9_.-]", "_", original.name)[:16] or "design"
    record = SeqRecord(
        Seq(project.sequence),
        id=identifier,
        name=identifier,
        description=f"{original.name}; validated DNA Design Platform design",
    )
    annotations: dict[str, Any] = dict(original.record_annotations)
    record.annotations = annotations
    record.annotations.pop("source_features", None)
    record.annotations.pop("references", None)
    if original.original_record_text and original.input_format == "genbank":
        source_record = SeqIO.read(io.StringIO(original.original_record_text), "genbank")
        if "references" in source_record.annotations:
            record.annotations["references"] = source_record.annotations["references"]
        record.features.extend(f for f in source_record.features if f.type == "source")
    record.annotations.update(molecule_type="DNA", topology=original.topology.value)
    record.annotations["comment"] = (
        "Validated against the original construct; "
        "full settings and provenance are in the JSON project."
    )
    for feature in original.features:
        qualifiers = {key: list(values) for key, values in feature.qualifiers.items()}
        qualifiers.pop("dna_design_mutability", None)
        qualifiers["label"] = [feature.name]
        qualifiers["dna_design_id"] = [feature.id]
        qualifiers["dna_mutability"] = [str(feature.mutability)]
        if feature.description:
            qualifiers["note"] = [feature.description]
        if feature.rbs_mask:
            qualifiers["dna_design_rbs_mask"] = [feature.rbs_mask]
        if feature.is_cds:
            qualifiers["transl_table"] = [str(feature.translation_table)]
            qualifiers["codon_start"] = [str(feature.frame_offset + 1)]
            qualifiers["translation"] = [
                translate_feature(project.sequence, feature).removesuffix("*")
            ]
            qualifiers["dna_design_cds_mode"] = [feature.cds_mode.value]
        record.features.append(
            SeqFeature(
                SimpleLocation(feature.start, feature.end, strand=int(feature.strand) or None),
                type=feature.feature_type,
                qualifiers=qualifiers,
            )
        )
    return record


def export_project(
    project: DesignProject, format: Literal["fasta", "genbank", "json", "csv"]
) -> str:
    """Validate original controls, carried locks and any assembly again before export."""
    validate_project(project).require_valid()
    if format == "json":
        return project.model_dump_json(indent=2) + "\n"
    if format == "csv":
        if project.gibson is None:
            raise ValueError("Fragment ordering CSV requires a validated Gibson plan")
        output = io.StringIO(newline="")
        writer = csv.writer(output)
        writer.writerow(
            (
                "fragment",
                "start_1based",
                "end_1based_unwrapped",
                "sequence",
                "unique_length",
                "ordered_length",
                "gc_fraction",
                "left_overlap",
                "right_overlap",
            )
        )
        for fragment in project.gibson.fragments:
            writer.writerow(
                (
                    fragment.number,
                    fragment.start + 1,
                    fragment.end,
                    fragment.sequence,
                    fragment.unique_length,
                    fragment.ordered_length,
                    fragment.gc,
                    fragment.left_overlap,
                    fragment.right_overlap,
                )
            )
        return output.getvalue()
    if format not in {"fasta", "genbank"}:
        raise ValueError("Unsupported export format")
    output = io.StringIO()
    SeqIO.write(_record(project), output, format)
    return output.getvalue()


def load_project(text: str) -> DesignProject:
    """Allow unfinished projects to load; export always enforces fresh validation."""
    payload = json.loads(text)
    if isinstance(payload, dict) and payload.get("format") == "dna-design-draft-v1":
        payload = payload["project"]
    return DesignProject.model_validate(payload)


def save_draft(project: DesignProject) -> str:
    return (
        json.dumps(
            {
                "format": "dna-design-draft-v1",
                "project": project.model_dump(mode="json"),
                "validation": validate_project(project).model_dump(mode="json"),
                "notice": "Work in progress; this snapshot is not a validated final design",
            },
            indent=2,
        )
        + "\n"
    )
