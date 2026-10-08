"""Local thermodynamic RBS backend and original-reference masked design."""

import importlib
import math
import warnings
from collections import OrderedDict
from importlib.metadata import PackageNotFoundError, version
from typing import Annotated, Protocol, Self

from Bio.Seq import Seq
from pydantic import Field, field_validator, model_validator

from .errors import DNAError
from .execution import ExecutionControl
from .models import CDSMode, Construct, Feature, FrozenModel, Identifier, Mutability, Region, Strand
from .optimization import AdditionalLock, ObjectiveWeights, SearchConfig, SearchResult, optimize
from .translation import extract_feature
from .validation import validate_candidate


class RBSBackendUnavailable(DNAError):
    """Install the optional local backend; never substitute a heuristic silently."""


class RBSModelError(DNAError):
    """The selected model cannot interpret this host, input, or result."""


class RBSHost(FrozenModel):
    id: Identifier = "ecoli-mg1655"
    anti_sd: str | None = None
    source: Identifier | None = None
    calibrated_ecoli: bool = True

    @field_validator("anti_sd")
    @classmethod
    def valid_anti_sd(cls, value: str | None) -> str | None:
        if value is not None and (len(value) != 9 or set(value) - set("ACGT")):
            raise ValueError("OSTIR anti-SD must be nine explicitly sourced DNA bases")
        return value

    @model_validator(mode="after")
    def calibration_claim(self) -> Self:
        if self.id == "ecoli-mg1655" and self.calibrated_ecoli:
            if self.anti_sd is None:
                object.__setattr__(self, "anti_sd", "ACCTCCTTA")
            if self.source is None:
                object.__setattr__(
                    self, "source", "https://github.com/barricklab/ostir/blob/v1.1.3/ostir/ostir.py"
                )
        elif self.anti_sd is None or self.source is None:
            raise ValueError("Uncalibrated hosts require an explicit anti-SD and its source")
        if self.calibrated_ecoli and (self.id != "ecoli-mg1655" or self.anti_sd != "ACCTCCTTA"):
            raise ValueError(
                "OSTIR's published E. coli calibration cannot label another host/anti-SD"
            )
        return self


class RBSPrediction(FrozenModel):
    relative_tir: Annotated[float, Field(gt=0, allow_inf_nan=False)]
    start_codon: str
    spacing_bp: int
    energy_components: dict[str, Annotated[float, Field(allow_inf_nan=False)]]
    unbound_mrna_structure: str
    unbound_mrna_mfe_kcal_mol: Annotated[float, Field(allow_inf_nan=False)]
    model_versions: dict[str, str]
    host: RBSHost
    warnings: tuple[str, ...]


class RBSBackend(Protocol):
    def predict(
        self, upstream: str, downstream: str, host: RBSHost, *, allow_uncalibrated: bool = False
    ) -> RBSPrediction: ...


class OSTIRBackend:
    """OSTIR 1.1.3 GPLv3 plus ViennaRNA Python API; no online prediction service."""

    def predict(
        self, upstream: str, downstream: str, host: RBSHost, *, allow_uncalibrated: bool = False
    ) -> RBSPrediction:
        host = RBSHost.model_validate(host.model_dump())
        if not host.calibrated_ecoli and not allow_uncalibrated:
            raise RBSModelError("Uncalibrated host requires explicit allow_uncalibrated=True")
        if set(upstream) - set("ACGT") or not downstream or set(downstream) - set("ACGT"):
            raise RBSModelError("RBS context requires concrete uppercase DNA")
        if downstream[:3] not in {"ATG", "GTG", "TTG"}:
            raise RBSModelError("OSTIR supports ATG/GTG/TTG starts; this start is unsupported")
        try:
            installed = {"ostir": version("ostir"), "ViennaRNA": version("ViennaRNA")}
            if installed != {"ostir": "1.1.3", "ViennaRNA": "2.7.2"}:
                raise RBSBackendUnavailable(
                    "Install the declared rbs extra with supported pinned versions"
                )
            with warnings.catch_warnings():
                # 1.1.3 uses the Python RNA API, but retains an obsolete CLI-binary
                # presence probe at import. Ignore only that exact stale warning.
                warnings.filterwarnings(
                    "ignore", message=r"^RBS Calculator Vienna is missing dependency ViennaRNA!$"
                )
                ostir = importlib.import_module("ostir")
                factory_module = importlib.import_module("ostir.ostir_factory")
                vienna = importlib.import_module("ostir.ViennaRNA")
                rna = importlib.import_module("RNA")
        except (ModuleNotFoundError, PackageNotFoundError) as exc:
            raise RBSBackendUnavailable(
                "Install dna-design-platform[rbs] for local RBS prediction"
            ) from exc
        sequence = upstream + downstream
        start = len(upstream) + 1  # Explicit conversion for OSTIR's 1-based API.
        try:
            rows = ostir.run_ostir(
                sequence, start=start, end=start, aSD=host.anti_sd, threads=1, decimal_places=8
            )
            if len(rows) != 1 or rows[0]["start_position"] != start:
                raise RBSModelError("OSTIR returned no prediction for the specified start")
            row = rows[0]
            # Obtain folding settings from the upstream model, not a second set
            # of fabricated scientific parameters. This is the UNBOUND structure.
            factory = factory_module.OSTIRFactory(sequence, [start, start], host.anti_sd, None)
            params = vienna.get_paramater_object(
                factory.RNA_model, factory.temp, factory.dangles_default
            )
            structure, energy = rna.fold_compound(sequence.replace("T", "U"), params).mfe()
            notes = [
                "Relative model TIR, not experimentally measured protein expression",
                "Upstream OSTIR calibration is not independent validation in this project",
                "Structure shown is the unbound context, not the ribosome-bound conformation",
            ]
            if len(upstream) < factory.cutoff or len(downstream) < factory.cutoff:
                notes.append(
                    "Context truncates the 35-nt window; boundary effects can affect prediction"
                )
            if not host.calibrated_ecoli:
                notes.append(
                    "Host transfer is uncalibrated; no PAO1 experimental accuracy is claimed"
                )
            return RBSPrediction(
                relative_tir=float(row["expression"]),
                start_codon=str(row["start_codon"]),
                spacing_bp=int(row["RBS_distance_bp"]),
                energy_components={
                    key: float(value) for key, value in row.items() if key.startswith("dG_")
                },
                unbound_mrna_structure=str(structure),
                unbound_mrna_mfe_kcal_mol=float(energy),
                model_versions=installed,
                host=host,
                warnings=tuple(notes),
            )
        except (ValueError, TypeError, KeyError, IndexError, RuntimeError) as exc:
            raise RBSModelError(f"Local OSTIR prediction failed: {exc}") from exc


class RBSContext(FrozenModel):
    cds_id: Identifier
    rbs_id: Identifier
    upstream_region: Region
    downstream_nt: Annotated[int, Field(ge=3, strict=True)] = 35
    host: RBSHost = Field(default_factory=RBSHost)
    allow_uncalibrated: bool = False


def _context_features(reference: Construct, context: RBSContext) -> tuple[Feature, Feature]:
    features = {f.id: f for f in reference.features}
    if context.cds_id not in features or context.rbs_id not in features:
        raise RBSModelError("RBS context must refer to original feature IDs")
    cds, rbs = features[context.cds_id], features[context.rbs_id]
    region = context.upstream_region
    if not cds.is_cds or cds.cds_mode != CDSMode.COMPLETE or cds.translation_table != 11:
        raise RBSModelError("This RBS backend requires an explicit complete bacterial table-11 CDS")
    if rbs.feature_type != "RBS" or rbs.strand != cds.strand:
        raise RBSModelError("RBS feature must have the same explicit strand as its CDS")
    if region.end > len(reference.sequence) or not (
        region.start <= rbs.start < rbs.end <= region.end
    ):
        raise RBSModelError("Upstream region must contain the RBS feature and remain in bounds")
    adjacent = region.end == cds.start if cds.strand == Strand.FORWARD else region.start == cds.end
    if not adjacent:
        raise RBSModelError(
            "Upstream region must end at the oriented CDS start, including any spacer"
        )
    if reference.host_profile and reference.host_profile.id != context.host.id:
        raise RBSModelError("RBS host disagrees with the original construct's host profile")
    return cds, rbs


def _context_sequences(reference: Construct, sequence: str, context: RBSContext) -> tuple[str, str]:
    cds, _ = _context_features(reference, context)
    region = context.upstream_region
    upstream = sequence[region.start : region.end]
    if cds.strand == Strand.REVERSE:
        upstream = str(Seq(upstream).reverse_complement())
    downstream = extract_feature(sequence, cds)[: context.downstream_nt]
    return upstream, downstream


def predict_rbs(
    reference: Construct, sequence: str, context: RBSContext, backend: RBSBackend | None = None
) -> RBSPrediction:
    validate_candidate(reference, sequence).require_valid()
    upstream, downstream = _context_sequences(reference, sequence, context)
    return (backend or OSTIRBackend()).predict(
        upstream, downstream, context.host, allow_uncalibrated=context.allow_uncalibrated
    )


class RBSDesignSettings(FrozenModel):
    seed: Annotated[int, Field(strict=True)]
    iterations: Annotated[int, Field(ge=0, strict=True)] = 200
    alternatives: Annotated[int, Field(ge=1, strict=True)] = 1
    joint_cds_nt: Annotated[int, Field(ge=0, strict=True)] = 0
    protect_downstream_nt: Annotated[int, Field(ge=3, strict=True)] = 30
    additional_locks: tuple[AdditionalLock, ...] = ()


class RBSDesignResult(FrozenModel):
    target_relative_tir: Annotated[float, Field(gt=0, allow_inf_nan=False)]
    initial_prediction: RBSPrediction
    predictions: tuple[RBSPrediction, ...]
    log_fold_errors: tuple[float, ...]
    initiation_locks: tuple[AdditionalLock, ...]
    search: SearchResult
    selected_design_index: Annotated[int, Field(ge=0, strict=True)] = 0


def design_rbs(
    reference: Construct,
    context: RBSContext,
    target_relative_tir: float,
    settings: RBSDesignSettings,
    *,
    initial_sequence: str | None = None,
    backend: RBSBackend | None = None,
    execution: ExecutionControl | None = None,
) -> RBSDesignResult:
    if not math.isfinite(target_relative_tir) or target_relative_tir <= 0:
        raise ValueError("RBS target must be a positive finite relative model TIR")
    context = RBSContext.model_validate(context.model_dump())
    settings = RBSDesignSettings.model_validate(settings.model_dump())
    cds, rbs = _context_features(reference, context)
    if rbs.mutability != Mutability.DESIGNABLE_RBS or rbs.rbs_mask is None:
        raise RBSModelError("RBS design requires an original DESIGNABLE_RBS annotation and mask")
    if settings.joint_cds_nt > cds.end - cds.start:
        raise RBSModelError("Joint CDS design region exceeds the original CDS")
    initial = reference.sequence if initial_sequence is None else initial_sequence
    initial_prediction = predict_rbs(reference, initial, context, backend)
    allowed = set(range(rbs.start, rbs.end))
    oriented = tuple(range(cds.start, cds.end))
    if cds.strand == Strand.REVERSE:
        oriented = tuple(reversed(oriented))
    allowed.update(oriented[: settings.joint_cds_nt])
    config = SearchConfig(
        seed=settings.seed,
        iterations=settings.iterations,
        alternatives=settings.alternatives,
        n_terminal_lock_nt=0,
        weights=ObjectiveWeights(codon_adaptation=0, custom=1),
        allowed_positions=tuple(sorted(allowed)),
        additional_locks=settings.additional_locks,
    )
    cache: OrderedDict[tuple[str, str], RBSPrediction] = OrderedDict(
        [(_context_sequences(reference, initial, context), initial_prediction)]
    )

    def objective(sequence: str) -> float:
        key = _context_sequences(reference, sequence, context)
        if key not in cache:
            cache[key] = predict_rbs(reference, sequence, context, backend)
            if len(cache) > 512:
                cache.popitem(last=False)
        cache.move_to_end(key)
        return abs(math.log(cache[key].relative_tir / target_relative_tir))

    search = optimize(
        reference, config, initial_sequence=initial, custom_objective=objective, execution=execution
    )
    # Final predictions are fresh, just like the final biological validation.
    predictions = tuple(
        predict_rbs(reference, d.sequence, context, backend) for d in search.designs
    )
    protected = min(
        cds.end - cds.start,
        max(settings.protect_downstream_nt, settings.joint_cds_nt, context.downstream_nt),
    )
    start, end = (
        (context.upstream_region.start, cds.start + protected)
        if cds.strand == Strand.FORWARD
        else (cds.end - protected, context.upstream_region.end)
    )
    locks = tuple(
        AdditionalLock(
            start=start,
            end=end,
            sequence=d.sequence[start:end],
            reason="Preserve designed initiation context during subsequent codon optimization",
        )
        for d in search.designs
    )
    return RBSDesignResult(
        target_relative_tir=target_relative_tir,
        initial_prediction=initial_prediction,
        predictions=predictions,
        log_fold_errors=tuple(
            abs(math.log(p.relative_tir / target_relative_tir)) for p in predictions
        ),
        initiation_locks=locks,
        search=search,
    )


def design_rbs_library(
    reference: Construct,
    context: RBSContext,
    targets: tuple[float, ...],
    settings: RBSDesignSettings,
    *,
    backend: RBSBackend | None = None,
    execution: ExecutionControl | None = None,
    initial_sequence: str | None = None,
) -> tuple[RBSDesignResult, ...]:
    if not targets or len(set(targets)) != len(targets):
        raise ValueError("RBS library requires distinct explicit positive relative targets")
    results: list[RBSDesignResult] = []
    used: set[str] = set()
    _, rbs_feature = _context_features(reference, context)
    for index, target in enumerate(targets):
        # Per-target seeds derive reproducibly from the explicit library seed.
        candidate_settings = RBSDesignSettings.model_validate(
            {
                **settings.model_dump(),
                "seed": settings.seed + index,
                "alternatives": max(settings.alternatives, len(used) + 1),
            }
        )
        result = design_rbs(
            reference,
            context,
            target,
            candidate_settings,
            backend=backend,
            execution=execution,
            initial_sequence=initial_sequence,
        )
        available = [
            i
            for i, d in enumerate(result.search.designs)
            if d.sequence[rbs_feature.start : rbs_feature.end] not in used
        ]
        if not available:
            raise RBSModelError("Library search budget exhausted without a distinct valid design")
        selected = min(available, key=lambda i: result.log_fold_errors[i])
        used.add(result.search.designs[selected].sequence[rbs_feature.start : rbs_feature.end])
        results.append(result.model_copy(update={"selected_design_index": selected}))
    return tuple(results)
