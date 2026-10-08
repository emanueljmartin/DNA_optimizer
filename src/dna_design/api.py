"""Optional loopback FastAPI application; sequence bodies never reach remote services."""

import asyncio
import importlib.util
import threading
import time
import uuid
from collections.abc import AsyncIterator, Awaitable, Callable
from concurrent.futures import ThreadPoolExecutor
from contextlib import asynccontextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Annotated, Any, Literal
from urllib.parse import urlsplit

from Bio.Data import CodonTable
from fastapi import FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import Field, ValidationError
from starlette.middleware.base import BaseHTTPMiddleware

from .benchmarks import PUBLIC_DATASETS, BenchmarkRequest, evaluate_benchmark
from .errors import DNAError
from .execution import ExecutionControl, OperationCancelled
from .help import HELP
from .models import FrozenModel
from .profiles import load_codon_profile
from .projects import DesignProject, export_project, save_draft, validate_project
from .service import (
    AnalysisRequest,
    AnnotationRequest,
    ImportRequest,
    MotifRequest,
    OperationRequest,
    RBSContextRequest,
    analyze,
    annotate,
    execute_operation,
    import_input,
    rbs_context,
    set_motif,
)
from .synthesis import load_synthesis_profile


@dataclass
class Job:
    id: str
    status: str = "queued"
    progress: float = 0
    stage: str = "Queued"
    result: dict[str, Any] | None = None
    error: str | None = None
    cancel: threading.Event = field(default_factory=threading.Event)
    created: float = field(default_factory=time.monotonic)
    lock: threading.Lock = field(default_factory=threading.Lock)

    def snapshot(self) -> dict[str, Any]:
        with self.lock:
            return {
                "id": self.id,
                "status": self.status,
                "progress": self.progress,
                "stage": self.stage,
                "result": self.result,
                "error": self.error,
            }


class JobManager:
    def __init__(self) -> None:
        self.pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="dna-local-job")
        self.jobs: dict[str, Job] = {}
        self.lock = threading.Lock()

    def submit(self, request: OperationRequest) -> Job:
        with self.lock:
            for key in list(self.jobs):
                old = self.jobs[key]
                if old.status in {"completed", "failed", "cancelled"} and (
                    time.monotonic() - old.created > 1800 or len(self.jobs) >= 16
                ):
                    del self.jobs[key]
            if sum(job.status in {"queued", "running"} for job in self.jobs.values()) >= 4:
                raise HTTPException(
                    429, "Local job queue is full; cancel or wait for an existing job"
                )
            job = Job(id=uuid.uuid4().hex)
            self.jobs[job.id] = job
        self.pool.submit(self._run, job, request)
        return job

    def _run(self, job: Job, request: OperationRequest) -> None:
        def progress(completed: int, total: int, stage: str) -> None:
            with job.lock:
                job.progress = min(1, completed / max(1, total))
                job.stage = stage

        with job.lock:
            job.status = "running"
        control = ExecutionControl(progress=progress, cancelled=job.cancel.is_set)
        try:
            control.checkpoint(0, 1, "Starting local computation")
            result = execute_operation(request, control)
            control.checkpoint(1, 1, "Completed")
            with job.lock:
                job.result, job.status = result, "completed"
        except OperationCancelled:
            with job.lock:
                job.status, job.stage = "cancelled", "Cancelled; previous project preserved"
        except (DNAError, ValueError, ValidationError) as exc:
            with job.lock:
                job.status, job.error = "failed", str(exc)
        except Exception:
            # Do not log private payloads or expose stack-local sequence data.
            with job.lock:
                job.status, job.error = (
                    "failed",
                    "Unexpected local computation error; previous project preserved",
                )

    def get(self, identifier: str) -> Job:
        with self.lock:
            if identifier not in self.jobs:
                raise HTTPException(
                    404, "Job not found or expired; project files remain in your browser"
                )
            return self.jobs[identifier]

    def close(self) -> None:
        with self.lock:
            for job in self.jobs.values():
                job.cancel.set()
        self.pool.shutdown(wait=True, cancel_futures=True)


class LocalRequests(BaseHTTPMiddleware):
    async def dispatch(
        self, request: Request, call_next: Callable[[Request], Awaitable[Response]]
    ) -> Response:
        if request.method in {"POST", "DELETE", "PUT", "PATCH"}:
            if request.headers.get("x-dna-design-client") != "local-ui":
                return JSONResponse(
                    {"detail": "Local application client header is required"}, status_code=403
                )
            origin = request.headers.get("origin")
            if origin and urlsplit(origin).netloc != request.headers.get("host"):
                return JSONResponse(
                    {"detail": "Cross-origin sequence requests are disabled"}, status_code=403
                )
            body = bytearray()
            async for part in request.stream():
                body.extend(part)
                if len(body) > 2_000_000:
                    return JSONResponse(
                        {"detail": "Request exceeds the 2 MB local application limit"},
                        status_code=413,
                    )
            request._body = bytes(body)
        response = await call_next(request)
        response.headers["Cache-Control"] = "no-store"
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Referrer-Policy"] = "no-referrer"
        response.headers["Content-Security-Policy"] = (
            "default-src 'self'; connect-src 'self'; script-src 'self'; "
            "style-src 'self' 'unsafe-inline'; img-src 'self' data:; "
            "object-src 'none'; base-uri 'self'; frame-ancestors 'none'"
        )
        return response


class ExportRequest(FrozenModel):
    project: DesignProject
    format: Literal["fasta", "genbank", "csv", "json"] = "json"


class SelectionRequest(FrozenModel):
    index: Annotated[int, Field(ge=0, strict=True)]


def create_app() -> FastAPI:
    manager = JobManager()

    @asynccontextmanager
    async def lifespan(app: FastAPI) -> AsyncIterator[None]:
        yield
        await asyncio.to_thread(manager.close)

    app = FastAPI(
        title="DNA Design Platform",
        version="0.1.0",
        docs_url=None,
        redoc_url=None,
        openapi_url="/api/schema",
        lifespan=lifespan,
    )
    app.state.jobs = manager
    app.add_middleware(LocalRequests)

    @app.exception_handler(DNAError)
    async def biological_error(request: Request, exc: DNAError) -> JSONResponse:
        return JSONResponse({"detail": str(exc)}, status_code=422)

    @app.exception_handler(ValueError)
    async def value_error(request: Request, exc: ValueError) -> JSONResponse:
        return JSONResponse({"detail": str(exc)}, status_code=422)

    @app.exception_handler(RequestValidationError)
    async def schema_error(request: Request, exc: RequestValidationError) -> JSONResponse:
        errors = [
            {"location": list(error["loc"]), "message": error["msg"]} for error in exc.errors()
        ]
        return JSONResponse({"detail": errors}, status_code=422)

    @app.get("/api/health")
    def health() -> dict[str, Any]:
        return {"status": "ok", "local_only": True, "version": "0.1.0"}

    @app.get("/api/catalog")
    def catalog() -> dict[str, Any]:
        return {
            "hosts": [
                load_codon_profile(host).model_dump(mode="json") for host in ("ecoli", "pao1")
            ],
            "synthesis_profiles": [
                load_synthesis_profile(key).model_dump(mode="json")
                for key in ("generic", "twist-public", "genscript-public")
            ],
            "genetic_codes": [
                {"id": key, "name": table.names[0]}
                for key, table in CodonTable.unambiguous_dna_by_id.items()
                if not set(table.stop_codons) & set(table.forward_table)
            ],
            "rbs_available": importlib.util.find_spec("ostir") is not None
            and importlib.util.find_spec("RNA") is not None,
            "help": HELP,
        }

    @app.post("/api/import")
    def import_sequence(request: ImportRequest) -> dict[str, Any]:
        project = import_input(request)
        if len(project.original.sequence) > 100_000:
            raise ValueError("The interactive application supports constructs up to 100,000 bases")
        return project.model_dump(mode="json")

    @app.post("/api/annotations")
    def annotation(request: AnnotationRequest) -> dict[str, Any]:
        return annotate(request).model_dump(mode="json")

    @app.post("/api/motifs")
    def motif(request: MotifRequest) -> dict[str, Any]:
        return set_motif(request).model_dump(mode="json")

    @app.post("/api/analyze")
    def analysis(request: AnalysisRequest) -> dict[str, Any]:
        return analyze(request)

    @app.post("/api/rbs-context")
    def context(request: RBSContextRequest) -> dict[str, Any]:
        return rbs_context(request).model_dump(mode="json")

    @app.get("/api/demo")
    def demo() -> dict[str, Any]:
        from .models import Construct, Feature, Strand

        dna = "AACCAACCAACCAACCAACC" + "AAAGGAGGAAA" + "AACC" + "ATG" + "GAA" * 20 + "TAA"
        original = Construct(
            name="Synthetic demo",
            sequence=dna,
            features=(
                Feature(
                    id="upstream",
                    name="Protected upstream",
                    feature_type="promoter",
                    start=0,
                    end=20,
                ),
                Feature.model_validate(
                    {
                        "id": "rbs",
                        "name": "RBS",
                        "feature_type": "RBS",
                        "start": 20,
                        "end": 31,
                        "strand": 1,
                        "mutability": "DESIGNABLE_RBS",
                        "rbs_mask": "NNAGGAGGNNN",
                    }
                ),
                Feature(
                    id="cds",
                    name="Synthetic CDS",
                    feature_type="CDS",
                    start=35,
                    end=len(dna),
                    strand=Strand.FORWARD,
                ),
            ),
        )
        return DesignProject(original=original).model_dump(mode="json")

    @app.post("/api/jobs", status_code=202)
    def start_job(request: OperationRequest) -> dict[str, Any]:
        if len(request.project.original.sequence) > 100_000:
            raise ValueError("Interactive jobs support up to 100,000 bases")
        if request.search and (
            request.search.iterations > 20_000 or request.search.alternatives > 32
        ):
            raise ValueError("Interactive search allows up to 20,000 proposals and 32 alternatives")
        if request.rbs_settings and (
            request.rbs_settings.iterations > 5000
            or len(request.library_targets) > 32
            or request.rbs_settings.alternatives > 32
        ):
            raise ValueError(
                "Interactive RBS jobs allow up to 5,000 proposals and 32 library targets"
            )
        return manager.submit(request).snapshot()

    @app.get("/api/jobs/{identifier}")
    def job_status(identifier: str) -> dict[str, Any]:
        return manager.get(identifier).snapshot()

    @app.delete("/api/jobs/{identifier}")
    def cancel_job(identifier: str) -> dict[str, Any]:
        job = manager.get(identifier)
        job.cancel.set()
        return job.snapshot()

    @app.post("/api/jobs/{identifier}/select")
    def select_design(identifier: str, request: SelectionRequest) -> dict[str, Any]:
        result = manager.get(identifier).snapshot()["result"]
        if not result or request.index >= len(result.get("alternatives", [])):
            raise ValueError("Completed job has no such alternative")
        project = DesignProject.model_validate(result["project"])
        locks = result.get("alternative_locks")
        data = {
            **project.model_dump(),
            "sequence": result["alternatives"][request.index],
            "gibson": None,
            "history": project.history
            + ({"operation": "select_alternative", "index": request.index},),
        }
        if locks:
            data["initiation_locks"] = locks[request.index]
        project = DesignProject.model_validate(data)
        validate_project(project).require_valid()
        return project.model_dump(mode="json")

    @app.post("/api/export")
    def export(request: ExportRequest) -> Response:
        text = export_project(request.project, request.format)
        media = (
            "application/json"
            if request.format == "json"
            else "text/csv"
            if request.format == "csv"
            else "text/plain"
        )
        return Response(text, media_type=media)

    @app.post("/api/save-draft")
    def draft(project: DesignProject) -> Response:
        return Response(save_draft(project), media_type="application/json")

    @app.get("/api/benchmarks/catalog")
    def benchmark_catalog() -> dict[str, Any]:
        return {"loader_version": "public-mapped-csv-v1", "datasets": PUBLIC_DATASETS}

    @app.post("/api/benchmarks/evaluate")
    def benchmark(request: BenchmarkRequest) -> dict[str, Any]:
        return evaluate_benchmark(request.csv, request.manifest)

    web = Path(__file__).parent / "web"
    if (web / "index.html").exists():
        app.mount("/", StaticFiles(directory=web, html=True), name="browser")
    else:

        @app.get("/")
        def missing_browser() -> dict[str, str]:
            return {
                "message": "Build the browser with npm ci and npm run build in frontend, "
                "then run dna-design serve"
            }

    return app
