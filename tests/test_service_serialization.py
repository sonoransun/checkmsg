"""Round-trip tests for the dataclass -> dict -> schema serialization."""

from __future__ import annotations

import pytest

pytest.importorskip("fastapi")

from checkmsg import minerals  # noqa: E402
from checkmsg.diagnose import diagnose_profile  # noqa: E402
from checkmsg.service.schemas import DiagnoseResponse  # noqa: E402
from checkmsg.service.serialization import report_to_dict  # noqa: E402


@pytest.fixture(scope="module")
def ruby_report():
    return diagnose_profile(minerals.get("ruby"))


@pytest.mark.parametrize("tier", ["novice", "practitioner", "expert"])
def test_round_trip_validates_against_schema(ruby_report, tier):
    payload = report_to_dict(ruby_report, tier=tier)
    model = DiagnoseResponse.model_validate(payload)  # raises on schema mismatch
    assert model.disclaimer
    assert model.verdict.name == "ruby"
    assert model.confidence.band in ("high", "medium", "low", "inconclusive")


def test_tier_gating(ruby_report):
    novice = report_to_dict(ruby_report, tier="novice")
    practitioner = report_to_dict(ruby_report, tier="practitioner")
    expert = report_to_dict(ruby_report, tier="expert")

    assert "candidate_scores" not in novice
    assert "evidence" not in novice

    assert "evidence" in practitioner
    assert "candidate_scores" not in practitioner
    # weights are an expert-only detail
    assert all("weight" not in e for e in practitioner["evidence"])

    assert "candidate_scores" in expert
    assert any(e.get("weight") is not None for e in expert["evidence"])
    assert any(e.get("rationale") for e in expert["evidence"])


def test_verdict_enriched_and_disclaimer_present(ruby_report):
    payload = report_to_dict(ruby_report, tier="novice")
    assert payload["verdict"]["species"] == "corundum"
    assert payload["verdict"]["chemical_formula"]
    assert payload["disclaimer"]
