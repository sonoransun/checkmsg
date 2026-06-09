"""Tests for the cross-cutting scoring registries, bands, and conflict logic."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

from checkmsg import scoring

SRC_ROOT = Path(__file__).resolve().parents[1] / "src" / "checkmsg"


@dataclass
class _Ev:
    """Lightweight stand-in for diagnose.Evidence (structural duck-type)."""

    technique: str
    weight: float
    favors: tuple[str, ...] = ()
    rules_out: tuple[str, ...] = ()


# ---------------------------------------------------------------------------
# Confidence bands
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value,expected",
    [
        (1.0, "high"),
        (0.85, "high"),
        (0.8499, "medium"),
        (0.65, "medium"),
        (0.6499, "low"),
        (0.40, "low"),
        (0.3999, "inconclusive"),
        (0.0, "inconclusive"),
    ],
)
def test_confidence_band_boundaries(value, expected):
    assert scoring.confidence_band(value) == expected


def test_downgrade_band_steps_and_clamp():
    assert scoring.downgrade_band("high") == "medium"
    assert scoring.downgrade_band("medium") == "low"
    assert scoring.downgrade_band("low") == "inconclusive"
    assert scoring.downgrade_band("inconclusive") == "inconclusive"
    assert scoring.downgrade_band("high", 2) == "low"


# ---------------------------------------------------------------------------
# Weight registry parity (guards the diagnose.py refactor)
# ---------------------------------------------------------------------------

# These values are the literals that previously lived inline in diagnose.py.
_EXPECTED_WEIGHTS = {
    "raman.dominant_peak": 1.0,
    "raman.amorphous": 1.5,
    "raman.match_unit": 1.0,
    "raman.no_peaks": 0.5,
    "uvvis.chromophore": 0.6,
    "uvvis.no_chromophore": 0.3,
    "xrf.detected": 0.4,
    "xrf.major_set_unit": 0.4,
    "xrf.per_element": 0.3,
    "libs.detected": 0.3,
    "libs.major_set_unit": 0.3,
    "libs.per_element": 0.25,
    "epr.center_match": 0.7,
    "laicpms.isotope_set": 0.5,
    "squid.ordering": 0.7,
    "squid.tc_tn": 0.5,
    "squid.saturation": 0.4,
    "pl.line_match": 0.7,
    "pl.synthetic_marker": 0.5,
    "ftir.diamond_type": 0.7,
    "ftir.band_match": 0.5,
    "ftir.polymer_flag": 0.5,
    "mossbauer.valence": 0.7,
    "mossbauer.site_match": 0.5,
    "cl.band_match": 0.5,
}


def test_weight_registry_values_match_legacy_literals():
    assert set(scoring.WEIGHTS) == set(_EXPECTED_WEIGHTS)
    for key, value in _EXPECTED_WEIGHTS.items():
        assert scoring.WEIGHTS[key].value == value, key
        assert scoring.WEIGHTS[key].rationale  # every weight is documented


# ---------------------------------------------------------------------------
# Caveat coverage audit
# ---------------------------------------------------------------------------

_HONESTY_KEYWORDS = (
    "pedagogical", "approximation", "approximate", "heuristic",
    "not measurement", "out of scope", "relative indicator", "toy",
    "physically-motivated", "formula-correct", "not absolute", "degrad",
)


def test_every_caveat_points_at_a_real_honesty_marker():
    for tech, caveat in scoring.CAVEATS.items():
        file_part = caveat.source_marker.split(":")[0]
        path = SRC_ROOT / file_part
        assert path.exists(), f"{tech}: missing source file {file_part}"
        text = path.read_text(encoding="utf-8").lower()
        assert any(k in text for k in _HONESTY_KEYWORDS), (
            f"{tech}: no honesty marker found in {file_part}"
        )
        assert caveat.severity in (
            "info", "approximation", "formula-only", "pedagogical")


def test_caveats_for_techniques_normalises_squid():
    cavs = scoring.caveats_for_techniques({"raman", "squid-mh"})
    techs = {c.technique for c in cavs}
    assert techs == {"raman", "squid"}


# ---------------------------------------------------------------------------
# Conflict detection
# ---------------------------------------------------------------------------


def test_rules_out_verdict_conflict():
    ev = [_Ev("squid", 0.7, favors=("magnetite",), rules_out=("diamond",))]
    conflicts = scoring.detect_conflicts(ev, "diamond", {"diamond": 1.0, "magnetite": 0.7})
    kinds = {c.kind for c in conflicts}
    assert "rules_out_verdict" in kinds


def test_favors_non_winner_conflict():
    ev = [_Ev("uvvis", 0.6, favors=("red_spinel",))]
    # red_spinel scores within 75% of the verdict's score.
    conflicts = scoring.detect_conflicts(ev, "ruby", {"ruby": 1.0, "red_spinel": 0.9})
    assert any(c.kind == "favors_non_winner" and c.against == "red_spinel" for c in conflicts)


def test_no_favors_non_winner_when_runner_up_is_far():
    ev = [_Ev("uvvis", 0.6, favors=("red_spinel",))]
    conflicts = scoring.detect_conflicts(ev, "ruby", {"ruby": 1.0, "red_spinel": 0.2})
    assert not any(c.kind == "favors_non_winner" for c in conflicts)


def test_technique_disagreement_conflict():
    ev = [
        _Ev("raman", 1.0, favors=("diamond",)),
        _Ev("squid", 0.7, rules_out=("diamond",)),
    ]
    conflicts = scoring.detect_conflicts(ev, "diamond", {"diamond": 1.0})
    assert any(c.kind == "technique_disagreement" for c in conflicts)


def test_clean_evidence_has_no_conflicts():
    ev = [_Ev("raman", 1.0, favors=("ruby",))]
    conflicts = scoring.detect_conflicts(ev, "ruby", {"ruby": 1.0, "red_spinel": 0.0})
    assert conflicts == []


def test_no_conflicts_when_verdict_none():
    ev = [_Ev("raman", 1.0, rules_out=("ruby",))]
    assert scoring.detect_conflicts(ev, None, {}) == []


def test_band_after_conflicts_downgrades_only_on_rules_out():
    ro = [scoring.Conflict("rules_out_verdict", "squid", "x", 0.7, "diamond")]
    nw = [scoring.Conflict("favors_non_winner", "uvvis", "x", 0.6, "red_spinel")]
    assert scoring.band_after_conflicts("high", ro) == "medium"
    assert scoring.band_after_conflicts("high", nw) == "high"
    assert scoring.band_after_conflicts("high", []) == "high"
