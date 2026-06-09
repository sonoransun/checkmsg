"""Upload-path test: 2-column CSV multipart -> diagnosis."""

from __future__ import annotations

import json

import numpy as np
import pytest

pytest.importorskip("fastapi")

from fastapi.testclient import TestClient  # noqa: E402

from checkmsg import io as io_mod  # noqa: E402
from checkmsg.service.app import create_app  # noqa: E402
from checkmsg.synthetic import PeakSpec, generate  # noqa: E402


def test_diagnose_upload_csv(tmp_path):
    client = TestClient(create_app())
    axis = np.linspace(100.0, 1700.0, 1601)
    spec = generate([PeakSpec(1332.0, 1.0, 1.5, 0.5)], axis,
                    technique="raman", units="cm-1", noise=0.005, seed=42)
    path = tmp_path / "diamond.csv"
    io_mod.write_csv(spec, path)

    with open(path, "rb") as fh:
        resp = client.post(
            "/diagnose/upload",
            files=[("files", ("diamond.csv", fh, "text/csv"))],
            data={"specs": json.dumps([{"technique": "raman"}]), "tier": "expert"},
        )
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["verdict"]["name"] == "diamond"
    assert "candidate_scores" in body
    assert body["disclaimer"]
