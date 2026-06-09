"""Keep docs/accuracy.md in lockstep with the scoring registries."""

from __future__ import annotations

from pathlib import Path

import pytest

from checkmsg import scoring

DOC = Path(__file__).resolve().parents[1] / "docs" / "accuracy.md"


@pytest.fixture(scope="module")
def doc_text() -> str:
    return DOC.read_text(encoding="utf-8")


def test_doc_exists(doc_text):
    assert doc_text.strip()


def test_disclaimer_is_verbatim(doc_text):
    assert scoring.SERVICE_DISCLAIMER in doc_text


def test_every_weight_key_is_documented(doc_text):
    for key in scoring.WEIGHTS:
        assert key in doc_text, f"weight key {key!r} missing from accuracy.md"


def test_every_caveat_technique_is_documented(doc_text):
    for tech in scoring.CAVEATS:
        assert tech in doc_text, f"caveat technique {tech!r} missing from accuracy.md"


def test_band_labels_documented(doc_text):
    for band in scoring.BANDS:
        assert band in doc_text
