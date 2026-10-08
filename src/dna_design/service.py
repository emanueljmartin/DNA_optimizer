"""Local application orchestration; the domain modules own scientific calculations."""

from typing import Annotated, Any, Literal

from pydantic import Field

from .coordinates import internal_to_ui, ui_to_internal
from .errors import DNAError
from .execution import ExecutionControl
from .gibson import GibsonSettings, plan_gibson
from .metrics import codon_metrics, gc_fraction, local_gc
from .models import (
    Construct,
    Feature,
    FrozenModel,
    MotifConstraint,
    MotifMode,
    Region,
    Strand,
    Topology,
)
from .motifs import find_motif, restriction_constraint
from .optimization import SearchConfig, SearchResult, optimize
from .parsing import parse_fasta, parse_genbank, parse_raw
from .permissions import effective_permissions
from .profiles import CodonProfile, load_codon_profile
from .projects import DesignProject, load_project, validate_project
from .rbs import RBSContext, RBSDesignSettings, RBSHost, design_rbs, design_rbs_library, predict_rbs
from .reverse_translation import reverse_translate
from .synthesis import SynthesisProfile, repair_synthesis, scan_synthesis
from .translation import extract_feature, translate_feature


class ImportRequest(FrozenModel):
    format: Literal["raw", "fasta", "genbank", "protein", "project"]
    text: Annotated[str, Field(min_length=1, max_length=2_000_000)]
    name: str = "untitled"
    topology: Topology = Topology.LINEAR
    annotate_cds: bool = False
    seed: int = 1
    host: str = "ecoli"


def import_input(request: ImportRequest) -> DesignProject:
    if request.format == "project":
        return load_project(request.text)
    if request.format == "genbank":
        original = parse_genbank(request.text)
    elif request.format == "fasta":
        original = parse_fasta(request.text, topology=request.topology)
    elif request.format == "protein":
        original = reverse_translate(
            request.text, load_codon_profile(request.host), seed=request.seed, name=request.name
        )
    else:
        original = parse_raw(request.text, name=request.name, topology=request.topology)
    if request.annotate_cds and request.format in {"raw", "fasta"}:
        feature = Feature(
            id="cds",
            name="CDS",
            feature_type="CDS",
            start=0,
            end=len(original.sequence),
            strand=Strand.FORWARD,
        )
        original = Construct.model_validate({**original.model_dump(), "features": (feature,)})
    return DesignProject(original=original)


class AnnotationRequest(FrozenModel):
    project: DesignProject
    id: str
    name: str
    feature_type: str
    start_1based: Annotated[int, Field(ge=1, strict=True)]
    end_1based: Annotated[int, Field(ge=1, strict=True)]
    strand: Strand = Strand.UNKNOWN
    mutability: str | None = None
    translation_table: int = 11
    cds_mode: str = "complete"
    frame_offset: int = 0
    rbs_mask: str | None = None
    delete: bool = False


def annotate(request: AnnotationRequest) -> DesignProject:
    original = request.project.original
    existing = next((feature for feature in original.features if feature.id == request.id), None)
    features = [feature for feature in original.features if feature.id != request.id]
    if not request.delete:
        start, end = ui_to_internal(
            request.start_1based, request.end_1based, len(original.sequence)
        )
        feature = Feature.model_validate(
            {
                **(existing.model_dump() if existing else {}),
                "id": request.id,
                "name": request.name,
                "feature_type": request.feature_type,
                "start": start,
                "end": end,
                "strand": request.strand,
                "mutability": request.mutability,
                "translation_table": request.translation_table,
                "cds_mode": request.cds_mode,
                "frame_offset": request.frame_offset,
                "rbs_mask": request.rbs_mask,
            }
        )
        features.append(feature)
    original = Construct.model_validate({**original.model_dump(), "features": features})
    # Explicit curation establishes a new reference; never silently retain prior designs.
    return DesignProject(original=original)


class MotifRequest(FrozenModel):
    project: DesignProject
    id: str
    mode: str = "AVOID"
    motif: str | None = None
    enzyme: str | None = None
    strand: str = "both"
    start_1based: int | None = None
    end_1based: int | None = None
    position_1based: int | None = None
    delete: bool = False


def set_motif(request: MotifRequest) -> DesignProject:
    original = request.project.original
    constraints = [constraint for constraint in original.constraints if constraint.id != request.id]
    if not request.delete:
        pattern = request.motif
        if request.enzyme:
            pattern = restriction_constraint(
                request.enzyme, mode=MotifMode(request.mode), id=request.id
            ).motif
        if not pattern:
            raise ValueError("Provide an enzyme name or IUPAC motif")
        if (request.start_1based is None) != (request.end_1based is None):
            raise ValueError("Motif scope requires both start and end")
        region = None
        if request.start_1based is not None and request.end_1based is not None:
            start, end = ui_to_internal(
                request.start_1based, request.end_1based, len(original.sequence)
            )
            region = Region(start=start, end=end)
        constraint = MotifConstraint.model_validate(
            {
                "id": request.id,
                "mode": request.mode,
                "motif": pattern,
                "enzyme": request.enzyme,
                "strand": request.strand,
                "region": region,
                "position": request.position_1based - 1
                if request.position_1based is not None
                else None,
            }
        )
        constraints.append(constraint)
    original = Construct.model_validate({**original.model_dump(), "constraints": constraints})
    return DesignProject(original=original)


class AnalysisRequest(FrozenModel):
    project: DesignProject
    host: str = "ecoli"
    codon_profile: CodonProfile | None = None
    rare_threshold: float | None = None
    gc_window: Annotated[int, Field(ge=1, strict=True)] = 50
    synthesis_profile: SynthesisProfile | None = None
    motif_search: str | None = None


def analyze(request: AnalysisRequest) -> dict[str, Any]:
    project = request.project
    sequence = project.sequence
    assert sequence is not None
    validation = validate_project(project)
    original = project.original
    features: list[dict[str, Any]] = []
    codon_changes: list[dict[str, Any]] = []
    for feature in original.features:
        start, end = internal_to_ui(feature.start, feature.end, len(original.sequence))
        row = {**feature.model_dump(mode="json"), "start_1based": start, "end_1based": end}
        if feature.is_cds:
            try:
                row.update(
                    original_protein=translate_feature(original.sequence, feature),
                    designed_protein=translate_feature(sequence, feature),
                )
                old = extract_feature(original.sequence, feature)[feature.frame_offset :]
                new = extract_feature(sequence, feature)[feature.frame_offset :]
                protein = row["original_protein"]
                for offset in range(0, min(len(old), len(new)), 3):
                    if old[offset : offset + 3] == new[offset : offset + 3]:
                        continue
                    coding_offset = feature.frame_offset + offset
                    positions = tuple(
                        feature.start + coding_offset + i + 1
                        if feature.strand == Strand.FORWARD
                        else feature.end - coding_offset - i
                        for i in range(3)
                    )
                    codon_changes.append(
                        {
                            "feature_id": feature.id,
                            "codon_index_1based": offset // 3 + 1,
                            "positions_1based": positions,
                            "original": old[offset : offset + 3],
                            "designed": new[offset : offset + 3],
                            "amino_acid": protein[offset // 3]
                            if offset // 3 < len(protein)
                            else "*",
                            "category": "synonymous"
                            if row["original_protein"] == row["designed_protein"]
                            else "protein mismatch",
                        }
                    )
            except DNAError as exc:
                row["translation_error"] = str(exc)
        features.append(row)
    permissions: list[dict[str, Any]] = []
    for permission in effective_permissions(original):
        state = permission.state.value
        if permissions and permissions[-1]["state"] == state:
            permissions[-1]["end"] = permission.position + 1
        else:
            permissions.append(
                {"start": permission.position, "end": permission.position + 1, "state": state}
            )
    codons: list[dict[str, Any]] = []
    metric_warning = None
    if set(sequence) <= set("ACGT") and len(sequence) == len(original.sequence):
        try:
            codons = [
                metric.model_dump(mode="json")
                for metric in codon_metrics(
                    original,
                    sequence,
                    request.codon_profile or load_codon_profile(request.host),
                    rare_threshold=request.rare_threshold,
                )
            ]
        except DNAError as exc:
            metric_warning = str(exc)
        except ValueError as exc:
            metric_warning = str(exc)
        synthesis = scan_synthesis(
            original,
            sequence,
            request.synthesis_profile or project.synthesis_policy or SynthesisProfile(),
        )
        gc = gc_fraction(sequence)
        windows = local_gc(sequence, request.gc_window)
    else:
        synthesis, gc, windows = None, None, ()
    matches: tuple[dict[str, Any], ...] = ()
    if request.motif_search and set(sequence) <= set("ACGT"):
        search = MotifConstraint(
            id="search", mode=MotifMode.REQUIRED, motif=request.motif_search.upper()
        )
        matches = tuple(
            {**match.model_dump(mode="json"), "positions": match.positions(len(sequence))}
            for match in find_motif(sequence, search, original.topology)
        )
    return {
        "validation": validation.model_dump(mode="json"),
        "valid": validation.valid,
        "features": features,
        "permissions": permissions,
        "gc": gc,
        "local_gc": windows,
        "codon_metrics": codons,
        "codon_changes": codon_changes,
        "metric_warning": metric_warning,
        "synthesis": synthesis.model_dump(mode="json") if synthesis else None,
        "changes": [
            {"position_1based": i + 1, "original": a, "designed": b}
            for i, (a, b) in enumerate(zip(original.sequence, sequence, strict=False))
            if a != b
        ],
        "motif_matches": matches,
    }


class OperationRequest(FrozenModel):
    operation: Literal["codon", "synthesis", "rbs", "rbs-library", "gibson", "rbs-predict"]
    project: DesignProject
    host: str = "ecoli"
    codon_profile: CodonProfile | None = None
    search: SearchConfig | None = None
    synthesis_profile: SynthesisProfile | None = None
    synthesis_weight: Annotated[float, Field(ge=0, allow_inf_nan=False)] = 0
    rbs_context: RBSContext | None = None
    rbs_settings: RBSDesignSettings | None = None
    target_relative_tir: float | None = None
    library_targets: tuple[float, ...] = ()
    gibson_settings: GibsonSettings | None = None


class RBSContextRequest(FrozenModel):
    project: DesignProject
    cds_id: str
    rbs_id: str
    upstream_nt: Annotated[int, Field(ge=1, strict=True)] = 35
    downstream_nt: Annotated[int, Field(ge=3, strict=True)] = 35
    host: str = "ecoli-mg1655"
    anti_sd: str | None = None
    source: str | None = None
    allow_uncalibrated: bool = False


def rbs_context(request: RBSContextRequest) -> RBSContext:
    original = request.project.original
    features = {feature.id: feature for feature in original.features}
    if request.cds_id not in features or request.rbs_id not in features:
        raise ValueError("Select original CDS and RBS features")
    cds, rbs = features[request.cds_id], features[request.rbs_id]
    if cds.strand == Strand.FORWARD:
        region = Region(
            start=max(0, min(cds.start - request.upstream_nt, rbs.start)), end=cds.start
        )
    else:
        region = Region(
            start=cds.end,
            end=min(len(original.sequence), max(cds.end + request.upstream_nt, rbs.end)),
        )
    host = RBSHost(
        id=request.host,
        anti_sd=request.anti_sd,
        source=request.source,
        calibrated_ecoli=request.host == "ecoli-mg1655" and not request.allow_uncalibrated,
    )
    return RBSContext(
        cds_id=cds.id,
        rbs_id=rbs.id,
        upstream_region=region,
        downstream_nt=request.downstream_nt,
        host=host,
        allow_uncalibrated=request.allow_uncalibrated,
    )


def execute_operation(request: OperationRequest, control: ExecutionControl) -> dict[str, Any]:
    project, operation = request.project, request.operation
    original, sequence = project.original, project.sequence
    assert sequence is not None
    payload: dict[str, Any] = {"operation": operation}
    updated = project
    search: SearchResult | None = None
    selected_profile: CodonProfile | None = None
    if operation in {"codon", "synthesis"}:
        if request.search is None:
            raise ValueError("A seeded search configuration is required")
        search_config = SearchConfig.model_validate(
            {
                **request.search.model_dump(),
                "additional_locks": request.search.additional_locks + project.initiation_locks,
            }
        )
        profile = request.codon_profile or load_codon_profile(request.host)
        selected_profile = profile
        if operation == "codon":
            policy = request.synthesis_profile or project.synthesis_policy or SynthesisProfile()
            if request.synthesis_weight:
                search_config = SearchConfig.model_validate(
                    {
                        **search_config.model_dump(),
                        "weights": {
                            **search_config.weights.model_dump(),
                            "custom": request.synthesis_weight,
                        },
                    }
                )
            search = optimize(
                original,
                search_config,
                profile,
                initial_sequence=sequence,
                execution=control,
                custom_objective=(lambda dna: scan_synthesis(original, dna, policy).risk_loss)
                if request.synthesis_weight
                else None,
            )
            payload["result"] = search.model_dump(mode="json")
            updated = project.model_copy(update={"synthesis_policy": policy})
        else:
            policy = request.synthesis_profile or project.synthesis_policy or SynthesisProfile()
            repair = repair_synthesis(
                original,
                policy,
                search_config,
                initial_sequence=sequence,
                codon_profile=profile,
                execution=control,
            )
            search = repair.search
            payload["result"] = repair.model_dump(mode="json")
            updated = project.model_copy(update={"synthesis_policy": policy})
        if search.designs:
            updated = updated.model_copy(
                update={"sequence": search.designs[0].sequence, "gibson": None}
            )
        payload["alternatives"] = [design.sequence for design in search.designs]
    elif operation.startswith("rbs"):
        if request.rbs_context is None:
            raise ValueError("Explicit oriented RBS context is required")
        if operation == "rbs-predict":
            payload["result"] = predict_rbs(original, sequence, request.rbs_context).model_dump(
                mode="json"
            )
        else:
            if request.rbs_settings is None:
                raise ValueError("Seeded RBS settings are required")
            settings = RBSDesignSettings.model_validate(
                {
                    **request.rbs_settings.model_dump(),
                    "additional_locks": request.rbs_settings.additional_locks
                    + project.initiation_locks,
                }
            )
            if operation == "rbs-library":
                results = design_rbs_library(
                    original,
                    request.rbs_context,
                    request.library_targets,
                    settings,
                    execution=control,
                    initial_sequence=sequence,
                )
                payload["result"] = [result.model_dump(mode="json") for result in results]
                designs = [
                    (
                        result.search.designs[result.selected_design_index].sequence,
                        result.initiation_locks[result.selected_design_index],
                    )
                    for result in results
                ]
            else:
                if request.target_relative_tir is None:
                    raise ValueError("Positive relative TIR target is required")
                result = design_rbs(
                    original,
                    request.rbs_context,
                    request.target_relative_tir,
                    settings,
                    initial_sequence=sequence,
                    execution=control,
                )
                payload["result"] = result.model_dump(mode="json")
                designs = [
                    (design.sequence, result.initiation_locks[i])
                    for i, design in enumerate(result.search.designs)
                ]
            payload["alternatives"] = [dna for dna, _ in designs]
            payload["alternative_locks"] = [
                [
                    *[existing.model_dump(mode="json") for existing in project.initiation_locks],
                    lock.model_dump(mode="json"),
                ]
                for _, lock in designs
            ]
            if designs:
                updated = project.model_copy(
                    update={
                        "sequence": designs[0][0],
                        "gibson": None,
                        "initiation_locks": project.initiation_locks + (designs[0][1],),
                    }
                )
    else:
        if request.gibson_settings is None:
            raise ValueError("Explicit ordered-fragment settings are required")
        plan = plan_gibson(original, sequence, request.gibson_settings, execution=control)
        payload["result"] = plan.model_dump(mode="json")
        updated = project.model_copy(update={"gibson": plan if plan.validation.valid else None})
    settings_data = request.model_dump(mode="json", exclude={"project"})
    if selected_profile is not None:
        settings_data["codon_profile"] = selected_profile.model_dump(mode="json")
    history = {"operation": operation, "settings": settings_data}
    from importlib.metadata import version

    versions = {
        **updated.model_versions,
        "dna-design-platform": version("dna-design-platform"),
        "biopython": version("biopython"),
    }
    if search is not None and search.profile_fingerprint:
        versions["codon_profile_sha256"] = search.profile_fingerprint
    if operation.startswith("rbs"):
        versions.update({"ostir": version("ostir"), "ViennaRNA": version("ViennaRNA")})
    updated = DesignProject.model_validate(
        {
            **updated.model_dump(),
            "settings": settings_data,
            "model_versions": versions,
            "history": updated.history + (history,),
        }
    )
    payload["project"] = updated.model_dump(mode="json")
    control.checkpoint(1, 1, "Completed")
    return payload
