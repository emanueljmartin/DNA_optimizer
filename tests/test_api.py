"""Local API contracts, privacy, curation and fresh export gating."""

import json
import threading
import time

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx2")
from fastapi.testclient import TestClient

from dna_design.api import create_app
from dna_design.projects import DesignProject

HEADERS = {"X-DNA-Design-Client": "local-ui"}


@pytest.fixture
def client():
    with TestClient(create_app(), headers=HEADERS) as local:
        yield local


def wait_job(client, identifier):
    for _ in range(300):
        job = client.get("/api/jobs/" + identifier).json()
        if job["status"] not in {"queued", "running"}:
            return job
        time.sleep(0.01)
    pytest.fail("Synthetic local job did not complete")


def test_privacy_limits_and_redacted_schema_errors(client):
    assert client.get("/api/health").json()["local_only"]
    assert client.get("/api/schema").status_code == 200
    assert (
        client.post("/api/import", headers={"X-DNA-Design-Client": ""}, json={}).status_code == 403
    )
    assert (
        client.post("/api/import", headers={"Origin": "https://other.example"}, json={}).status_code
        == 403
    )
    assert client.post("/api/import", content="A" * 2_000_001).status_code == 413
    response = client.post("/api/import", json={"format": "raw", "text": "CONFIDENTIAL_NOT_DNA"})
    assert response.status_code == 422
    assert "CONFIDENTIAL_NOT_DNA" not in response.text
    assert client.get("/api/health").headers["cache-control"] == "no-store"
    assert "connect-src 'self'" in client.get("/").headers["content-security-policy"]


def test_import_default_locks_ui_coordinates_and_annotation_metadata_preservation(client):
    project = client.post("/api/import", json={"format": "raw", "text": "ATGGAATAA"}).json()
    assert not project["original"]["features"]
    report = client.post("/api/analyze", json={"project": project}).json()
    assert report["permissions"] == [{"start": 0, "end": 9, "state": "LOCKED"}]
    request = {
        "project": project,
        "id": "cds",
        "name": "Test",
        "feature_type": "CDS",
        "start_1based": 1,
        "end_1based": 9,
        "strand": 1,
    }
    project = client.post("/api/annotations", json=request).json()
    project["original"]["features"][0]["qualifiers"] = {"gene": ["synthetic"]}
    project["original"]["features"][0]["description"] = "Preserved imported note"
    project["sequence"] = "ATGGAGTAA"
    project["history"] = [{"operation": "previous"}]
    response = client.post(
        "/api/annotations", json={**request, "project": project, "name": "Edited"}
    )
    assert response.status_code == 200, response.text
    edited = response.json()
    assert edited["sequence"] == "ATGGAATAA"
    assert edited["history"] == []
    feature = edited["original"]["features"][0]
    assert feature["start"] == 0 and feature["end"] == 9
    assert feature["qualifiers"] == {"gene": ["synthetic"]}
    assert feature["description"] == "Preserved imported note"


@pytest.mark.parametrize("format", ["fasta", "genbank", "json", "csv"])
def test_invalid_final_rejected_but_draft_restores_original(client, format):
    project = client.get("/api/demo").json()
    project["sequence"] = "T" + project["sequence"][1:]
    response = client.post("/api/export", json={"project": project, "format": format})
    assert response.status_code == 422
    draft = client.post("/api/save-draft", json=project)
    assert draft.status_code == 200
    assert (
        not draft.json()["validation"]["valid"]
        if "valid" in draft.json()["validation"]
        else any(check["status"] == "ERROR" for check in draft.json()["validation"]["checks"])
    )
    restored = client.post("/api/import", json={"format": "project", "text": draft.text}).json()
    assert restored["original"] == project["original"]
    assert restored["sequence"] == project["sequence"]


def test_seeded_job_alternatives_history_and_context(client):
    project = client.get("/api/demo").json()
    context = client.post(
        "/api/rbs-context", json={"project": project, "cds_id": "cds", "rbs_id": "rbs"}
    )
    assert context.status_code == 200, context.text
    assert context.json()["upstream_region"] == {"start": 0, "end": 35}
    request = {
        "operation": "codon",
        "project": project,
        "search": {"seed": 12, "iterations": 30, "alternatives": 2},
    }
    first = client.post("/api/jobs", json=request)
    assert first.status_code == 202
    job = wait_job(client, first.json()["id"])
    assert job["status"] == "completed", job
    result = job["result"]
    assert len(result["alternatives"]) == 2
    assert result["project"]["original"] == project["original"]
    assert result["project"]["history"][0]["settings"]["search"]["seed"] == 12
    assert "codon_profile_sha256" in result["project"]["model_versions"]
    selected = client.post("/api/jobs/" + job["id"] + "/select", json={"index": 1})
    assert selected.status_code == 200, selected.text
    assert selected.json()["sequence"] == result["alternatives"][1]
    assert (
        client.post(
            "/api/export", json={"project": selected.json(), "format": "genbank"}
        ).status_code
        == 200
    )
    second = client.post("/api/jobs", json=request).json()
    assert wait_job(client, second["id"])["result"]["alternatives"] == result["alternatives"]
    assert (
        client.post(
            "/api/jobs", json={**request, "search": {"seed": 1, "iterations": 20001}}
        ).status_code
        == 422
    )


def test_cancelled_job_has_no_result_or_project_mutation(client, monkeypatch):
    started = threading.Event()

    def cancellable(request, control):
        started.set()
        while True:
            control.checkpoint(1, 10, "Synthetic cancellation test")
            time.sleep(0.001)

    monkeypatch.setattr("dna_design.api.execute_operation", cancellable)
    project = client.get("/api/demo").json()
    before = json.dumps(project, sort_keys=True)
    job = client.post(
        "/api/jobs", json={"operation": "codon", "project": project, "search": {"seed": 1}}
    ).json()
    assert started.wait(1)
    assert client.delete("/api/jobs/" + job["id"]).status_code == 200
    cancelled = wait_job(client, job["id"])
    assert cancelled["status"] == "cancelled" and cancelled["result"] is None
    assert json.dumps(project, sort_keys=True) == before
    assert client.get("/api/jobs/missing").status_code == 404


def test_catalog_and_backend_metric_help(client):
    catalog = client.get("/api/catalog").json()
    assert len(catalog["hosts"]) == 2
    assert len(catalog["help"]) >= 10
    assert all(
        set(entry) >= {"definition", "effect", "limitations", "reference"}
        for entry in catalog["help"].values()
    )
    assert (
        DesignProject.model_validate(client.get("/api/demo").json()).original.name
        == "Synthetic demo"
    )


@pytest.mark.parametrize(
    "strand,sequence,positions", [(1, "ATGGAGTAA", [4, 5, 6]), (-1, "TTACTCCAT", [6, 5, 4])]
)
def test_codon_report_has_oriented_positions_and_synonymous_category(
    client, strand, sequence, positions
):
    original = "ATGGAATAA" if strand == 1 else "TTATTCCAT"
    project = client.post("/api/import", json={"format": "raw", "text": original}).json()
    project = client.post(
        "/api/annotations",
        json={
            "project": project,
            "id": "cds",
            "name": "CDS",
            "feature_type": "CDS",
            "start_1based": 1,
            "end_1based": 9,
            "strand": strand,
        },
    ).json()
    project["sequence"] = sequence
    report = client.post("/api/analyze", json={"project": project}).json()
    assert report["valid"]
    assert report["codon_changes"][0] == {
        "feature_id": "cds",
        "codon_index_1based": 2,
        "positions_1based": positions,
        "original": "GAA",
        "designed": "GAG",
        "amino_acid": "E",
        "category": "synonymous",
    }


def test_combined_preferences_and_profile_are_recorded_and_backend_benchmark_is_local(client):
    from .test_benchmarks import TEXT, manifest

    project = client.get("/api/demo").json()
    job = client.post(
        "/api/jobs",
        json={
            "operation": "codon",
            "project": project,
            "search": {"seed": 8, "iterations": 10},
            "synthesis_weight": 2,
            "synthesis_profile": {"homopolymer_limit": 3},
        },
    ).json()
    completed = wait_job(client, job["id"])
    assert completed["status"] == "completed", completed
    result = completed["result"]
    assert "custom" in result["result"]["designs"][0]["components"]
    assert result["project"]["settings"]["codon_profile"]["counts"]
    response = client.post(
        "/api/benchmarks/evaluate",
        json={"csv": TEXT, "manifest": manifest(TEXT).model_dump(mode="json")},
    )
    assert response.status_code == 200, response.text
    assert response.json()["manifest"]["evidence"] == "synthetic"
    assert "boel2016" in client.get("/api/benchmarks/catalog").json()["datasets"]
