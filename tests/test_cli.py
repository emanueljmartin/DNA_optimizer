"""Command-line validation, local evaluation and refusal to replace existing files."""

import json
import sys

import pytest

from dna_design.cli import main
from dna_design.projects import DesignProject, save_draft

from .test_benchmarks import TEXT, manifest


def test_cli_validation_and_existing_export_are_safe(tmp_path, monkeypatch, coding_construct):
    path = tmp_path / "project.json"
    path.write_text(save_draft(DesignProject(original=coding_construct)))
    monkeypatch.setattr(sys, "argv", ["dna-design", "validate", str(path)])
    with pytest.raises(SystemExit) as exit:
        main()
    assert exit.value.code == 0
    output = tmp_path / "design.fasta"
    monkeypatch.setattr(
        sys,
        "argv",
        ["dna-design", "export", str(path), "--format", "fasta", "--output", str(output)],
    )
    main()
    before = output.read_text()
    with pytest.raises(FileExistsError):
        main()
    assert output.read_text() == before


def test_cli_benchmark_uses_only_local_files(tmp_path, monkeypatch):
    data = tmp_path / "data.csv"
    data.write_text(TEXT)
    provenance = tmp_path / "manifest.json"
    provenance.write_text(manifest(TEXT).model_dump_json())
    report = tmp_path / "report.json"
    monkeypatch.setattr(
        sys,
        "argv",
        ["dna-design", "benchmark", str(data), str(provenance), "--output", str(report)],
    )
    main()
    assert json.loads(report.read_text())["manifest"]["evidence"] == "synthetic"
