"""Endpoint tests for the FastAPI service (offline via conftest)."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("fastapi")

from fastapi.testclient import TestClient  # noqa: E402

from checkmsg import minerals  # noqa: E402
from checkmsg.service.app import create_app  # noqa: E402
from checkmsg.synthetic import PeakSpec, generate  # noqa: E402


@pytest.fixture(scope="module")
def client():
    return TestClient(create_app())


def _diamond_raman():
    axis = np.linspace(100.0, 1700.0, 1601)
    spec = generate([PeakSpec(1332.0, 1.0, 1.5, 0.5)], axis,
                    technique="raman", units="cm-1", noise=0.005, seed=42)
    return spec.axis.tolist(), spec.intensity.tolist()


def test_meta_endpoints(client):
    assert client.get("/healthz").json() == {"status": "ok"}
    assert client.get("/version").json()["name"] == "checkmsg"
    root = client.get("/").json()
    assert "/diagnose" in root["endpoints"]
    assert root["disclaimer"]


def test_techniques(client):
    data = client.get("/techniques").json()
    assert data["count"] == 11  # 7 core + PL/FTIR/Mössbauer/CL
    for item in data["items"]:
        assert item["accuracy_tier"] in ("info", "approximation", "formula-only", "pedagogical")
        assert "caveat" in item


def test_catalog(client):
    data = client.get("/catalog").json()
    assert data["count"] == len(minerals.names())
    assert client.get("/catalog/diamond").json()["species"] == "diamond"
    assert client.get("/catalog/zzznope").status_code == 404
    # species filter
    corundum = client.get("/catalog", params={"species": "corundum"}).json()
    assert corundum["count"] >= 1


def test_catalog_alias_resolution(client):
    # Find any catalog entry that has an alias, and confirm the alias resolves.
    alias = name = None
    for n in minerals.names():
        p = minerals.get(n)
        if p.aliases:
            alias, name = p.aliases[0], n
            break
    if alias:
        resolved = client.get(f"/catalog/{alias}")
        assert resolved.status_code == 200
        assert resolved.json()["name"] == name


def test_glossary(client):
    assert client.get("/glossary").json()["count"] >= 50
    epr = client.get("/glossary/EPR").json()
    assert epr["expanded"] == "EPR (electron paramagnetic resonance)"
    assert client.get("/glossary/zzznope").status_code == 404


def test_diagnose_novice_vs_expert(client):
    axis, intensity = _diamond_raman()
    base = {"technique": "raman", "axis": axis, "intensity": intensity}

    novice = client.post("/diagnose", json={"tier": "novice", "spectra": [base]})
    assert novice.status_code == 200
    nj = novice.json()
    assert nj["verdict"]["name"] == "diamond"
    assert nj["confidence"]["band"] in ("high", "medium", "low", "inconclusive")
    assert nj["disclaimer"]
    assert "candidate_scores" not in nj  # gated out at novice

    expert = client.post("/diagnose", json={"tier": "expert", "spectra": [base]})
    ej = expert.json()
    assert "candidate_scores" in ej
    assert "evidence" in ej and "reasoning_trace" in ej
    assert ej["evidence"][0].get("weight") is not None  # weight surfaced at expert


def test_analyze_single_technique(client):
    axis, intensity = _diamond_raman()
    resp = client.post("/analyze/raman",
                       json={"axis": axis, "intensity": intensity, "tier": "practitioner"})
    assert resp.status_code == 200
    assert resp.json()["verdict"]["name"] == "diamond"


def test_validation_and_hardening(client):
    axis, intensity = _diamond_raman()
    base = {"technique": "raman", "axis": axis, "intensity": intensity}

    # Too many spectra (limit is 12).
    too_many = client.post("/diagnose", json={"tier": "novice", "spectra": [base] * 13})
    assert too_many.status_code == 422

    # muon-xray is not an allowed technique.
    bad = client.post("/diagnose", json={"tier": "novice",
                                         "spectra": [{**base, "technique": "muon-xray"}]})
    assert bad.status_code == 422

    # EPR without frequency_GHz is rejected.
    epr = client.post("/diagnose", json={"tier": "novice",
                                         "spectra": [{**base, "technique": "epr"}]})
    assert epr.status_code == 422

    # Unknown technique in the path.
    assert client.post("/analyze/bogus",
                       json={"axis": axis, "intensity": intensity}).status_code == 422
