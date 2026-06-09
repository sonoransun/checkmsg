"""Tests for the optional calibrated-confidence layer."""

from __future__ import annotations

from dataclasses import dataclass, field

from checkmsg import calibrate, minerals
from checkmsg.diagnose import diagnose_profile


@dataclass
class _Ev:
    technique: str
    weight: float
    favors: tuple[str, ...] = ()


@dataclass
class _Report:
    """Minimal duck-typed stand-in for DiagnosticReport features."""
    verdict: str | None
    confidence: float
    candidate_scores: dict
    evidence: list = field(default_factory=list)
    evidence_agreement: int = 0
    conflicts: list = field(default_factory=list)


def _stub(conf: float) -> _Report:
    return _Report(verdict="ruby", confidence=conf,
                   candidate_scores={"ruby": 3.0, "red_spinel": 1.0},
                   evidence=[_Ev("raman", 1.0, ("ruby",))], evidence_agreement=2)


def test_feature_extraction():
    f = calibrate.features_from_report(_stub(0.8))
    assert set(f) == set(calibrate.FEATURE_NAMES)
    assert f["separation_ratio"] == 0.8
    assert f["n_techniques"] == 1.0


def test_platt_predict_proba_monotonic():
    cal = calibrate.Calibrator(method="platt", platt=calibrate.PlattParams(A=4.0, B=-2.0))
    probs = [cal.predict_proba(_stub(r)) for r in (0.2, 0.5, 0.8, 1.0)]
    assert probs == sorted(probs)  # non-decreasing in separation ratio
    assert all(0.0 <= p <= 1.0 for p in probs)


def test_calibrator_json_round_trip():
    cal = calibrate.Calibrator(method="platt", platt=calibrate.PlattParams(1.5, -0.5),
                               metadata={"ece": 0.1})
    back = calibrate.Calibrator.from_dict(cal.to_dict())
    assert back.method == "platt"
    assert abs(back.platt.A - 1.5) < 1e-9 and abs(back.platt.B + 0.5) < 1e-9
    assert back.predict_proba(_stub(0.7)) == cal.predict_proba(_stub(0.7))


def test_fit_platt_recovers_increasing_sigmoid():
    # Higher ratios are more often correct → A should be positive.
    ratios = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    labels = [0, 0, 0, 1, 1, 1, 1, 1]
    params = calibrate.fit_platt(ratios, labels)
    assert params.A > 0.0


def test_metrics_bounds():
    assert 0.0 <= calibrate.brier_score([0.2, 0.8], [0, 1]) <= 1.0
    assert 0.0 <= calibrate.expected_calibration_error([0.2, 0.8, 0.5], [0, 1, 1]) <= 1.0


def test_default_diagnose_path_is_uncalibrated():
    report = diagnose_profile(minerals.get("ruby"))
    assert report.calibrated_confidence is None
    assert report.calibration_method is None


def test_opt_in_calibration_when_artifact_present():
    cal = calibrate.load_calibrator("platt")
    if cal is None:
        return  # artifact not built in this environment — opt-in degrades to no-op
    from checkmsg.diagnose import diagnose
    spectra = [minerals.synthesize_raman(minerals.get("ruby"), seed=0)]
    report = diagnose(spectra, calibrated=True)
    assert report.calibration_method == "platt"
    assert 0.0 <= report.calibrated_confidence <= 1.0


def test_gnb_round_trip_and_monotonic():
    # Only separation_ratio (feature 0) discriminates; the rest match the _stub.
    gnb = calibrate.GNBParams(
        class_means={"correct": [0.9, 2, 3, 2, 1, 1, 1, 0],
                     "incorrect": [0.3, 2, 3, 2, 1, 1, 1, 0]},
        class_vars={"correct": [0.05] * 8, "incorrect": [0.05] * 8},
        log_priors={"correct": -0.3, "incorrect": -1.3})
    cal = calibrate.Calibrator(method="gnb", gnb=gnb)
    back = calibrate.Calibrator.from_dict(cal.to_dict())
    assert back.method == "gnb" and back.gnb is not None
    assert 0.0 <= back.predict_proba(_stub(0.95)) <= 1.0
    # The "correct" class has higher separation_ratio mean → higher P for high ratio.
    assert back.predict_proba(_stub(0.95)) > back.predict_proba(_stub(0.2))


def test_fit_gnb_separates_classes():
    X = [[0.9] * 8] * 5 + [[0.1] * 8] * 5
    y = [1] * 5 + [0] * 5
    gnb = calibrate.fit_gnb(X, y)
    assert set(gnb.class_means) == {"correct", "incorrect"}
    assert gnb.class_means["correct"][0] > gnb.class_means["incorrect"][0]


def test_opt_in_gnb_when_artifact_present():
    cal = calibrate.load_calibrator("gnb")
    if cal is None:
        return
    from checkmsg.diagnose import diagnose
    spectra = [minerals.synthesize_raman(minerals.get("ruby"), seed=0)]
    report = diagnose(spectra, calibrated="gnb")
    assert report.calibration_method == "gnb"
    assert 0.0 <= report.calibrated_confidence <= 1.0
