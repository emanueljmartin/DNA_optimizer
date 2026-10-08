"""Local launcher and reproducible project exports; never publish automatically."""

import argparse
import importlib
import json
from pathlib import Path

from .benchmarks import BenchmarkManifest, evaluate_benchmark
from .projects import export_project, load_project, validate_project


def main() -> None:
    parser = argparse.ArgumentParser(prog="dna-design")
    commands = parser.add_subparsers(dest="command", required=True)
    serve = commands.add_parser("serve", help="Open the local application backend/browser")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8000)
    validate = commands.add_parser("validate", help="Independently validate a saved JSON project")
    validate.add_argument("project", type=Path)
    export = commands.add_parser("export", help="Export only a freshly validated project")
    export.add_argument("project", type=Path)
    export.add_argument("--format", choices=("fasta", "genbank", "json", "csv"), required=True)
    export.add_argument("--output", type=Path, required=True)
    benchmark = commands.add_parser(
        "benchmark", help="Evaluate local CSV predictions with explicit provenance/splits"
    )
    benchmark.add_argument("csv", type=Path)
    benchmark.add_argument("manifest", type=Path)
    benchmark.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.command == "serve":
        try:
            uvicorn = importlib.import_module("uvicorn")
        except ImportError:
            parser.error("Install dna-design-platform[app] to run the local browser application")
        uvicorn.run(
            "dna_design.api:create_app",
            factory=True,
            host=arguments.host,
            port=arguments.port,
            access_log=False,
        )
    elif arguments.command == "validate":
        report = validate_project(load_project(arguments.project.read_text()))
        print(report.model_dump_json(indent=2))
        raise SystemExit(0 if report.valid else 1)
    elif arguments.command == "benchmark":
        manifest = BenchmarkManifest.model_validate_json(arguments.manifest.read_text())
        benchmark_report = evaluate_benchmark(arguments.csv.read_text(encoding="utf-8"), manifest)
        with arguments.output.open("x", encoding="utf-8") as output:
            json.dump(benchmark_report, output, indent=2, allow_nan=False)
            output.write("\n")
    else:
        text = export_project(load_project(arguments.project.read_text()), arguments.format)
        # Explicit output path, refusing replacement of any existing file.
        with arguments.output.open("x") as output:
            output.write(text)
