"""Optional calibrated confidence — turns the separation ratio into a probability.

The default ``diagnose`` confidence is a *separation ratio*, not a calibrated
P(correct). This module fits a Platt-scaling sigmoid on leave-mineral-out
cross-validated synthetic self-diagnoses so an opt-in caller can read an
estimated probability instead. It is scipy/numpy only (no scikit-learn, no
pickle), the fitted parameters live in a human-readable JSON artifact, and the
default diagnosis path is unchanged — calibration is a parallel readout.

HONEST LIMIT: the calibrator is a pedagogical approximation trained on
noise-clean synthetic self-diagnoses; it estimates P(the pipeline's verdict
equals the true catalog mineral on synthetic data), NOT a real-world guarantee
of correctness.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from functools import lru_cache
from importlib.resources import files
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from checkmsg.diagnose import DiagnosticReport

FEATURE_NAMES: tuple[str, ...] = (
    "separation_ratio", "margin", "top_score", "agreement",
    "n_techniques", "n_evidence", "verdict_favor_sum", "has_conflict",
)


def features_from_report(report: DiagnosticReport) -> dict[str, float]:
    """Extract the calibration feature vector from a report (no new computation)."""
    scores = sorted(report.candidate_scores.values(), reverse=True)
    top = scores[0] if scores else 0.0
    second = scores[1] if len(scores) > 1 else 0.0
    techs: set[str] = set()
    favor_sum = 0.0
    for ev in report.evidence:
        t = ev.technique
        techs.add("squid" if t in ("squid-mh", "squid-chi") else t)
        if report.verdict and report.verdict in ev.favors:
            favor_sum += ev.weight
    return {
        "separation_ratio": float(report.confidence),
        "margin": float(top - second),
        "top_score": float(top),
        "agreement": float(report.evidence_agreement),
        "n_techniques": float(len(techs)),
        "n_evidence": float(len(report.evidence)),
        "verdict_favor_sum": float(favor_sum),
        "has_conflict": 1.0 if report.conflicts else 0.0,
    }


def feature_vector(report: DiagnosticReport) -> np.ndarray:
    f = features_from_report(report)
    return np.array([f[k] for k in FEATURE_NAMES])


@dataclass(frozen=True)
class PlattParams:
    A: float
    B: float
    feature: str = "separation_ratio"


@dataclass(frozen=True)
class GNBParams:
    """Gaussian Naive Bayes over the full feature vector (correct vs incorrect)."""

    class_means: dict      # "correct"/"incorrect" -> per-feature mean list
    class_vars: dict       # per-feature variance list (floored)
    log_priors: dict       # log class prior


@dataclass(frozen=True)
class Calibrator:
    method: str
    platt: PlattParams | None = None
    gnb: GNBParams | None = None
    feature_names: tuple[str, ...] = FEATURE_NAMES
    metadata: dict = field(default_factory=dict)

    def predict_proba(self, report: DiagnosticReport) -> float:
        f = features_from_report(report)
        if self.method == "platt" and self.platt is not None:
            z = self.platt.A * f[self.platt.feature] + self.platt.B
            z = max(-60.0, min(60.0, z))
            return float(1.0 / (1.0 + math.exp(-z)))
        if self.method == "gnb" and self.gnb is not None:
            from scipy.special import logsumexp
            x = np.array([f[k] for k in self.feature_names])
            logp = {}
            for cls in ("correct", "incorrect"):
                m = np.array(self.gnb.class_means[cls])
                v = np.array(self.gnb.class_vars[cls])
                logp[cls] = (-0.5 * np.sum(np.log(2.0 * np.pi * v) + (x - m) ** 2 / v)
                             + self.gnb.log_priors[cls])
            return float(np.exp(logp["correct"]
                                - logsumexp([logp["correct"], logp["incorrect"]])))
        return float(f["separation_ratio"])

    def to_dict(self) -> dict:
        d: dict = {
            "schema_version": 1, "method": self.method,
            "feature_names": list(self.feature_names), "metadata": self.metadata,
        }
        if self.platt is not None:
            d["params"] = {"A": self.platt.A, "B": self.platt.B, "feature": self.platt.feature}
        if self.gnb is not None:
            d["gnb"] = {"class_means": self.gnb.class_means, "class_vars": self.gnb.class_vars,
                        "log_priors": self.gnb.log_priors}
        return d

    @classmethod
    def from_dict(cls, d: dict) -> Calibrator:
        platt = None
        if "params" in d:
            p = d["params"]
            platt = PlattParams(float(p["A"]), float(p["B"]), p.get("feature", "separation_ratio"))
        gnb = None
        if "gnb" in d:
            g = d["gnb"]
            gnb = GNBParams(g["class_means"], g["class_vars"], g["log_priors"])
        return cls(method=d["method"], platt=platt, gnb=gnb,
                   feature_names=tuple(d.get("feature_names", FEATURE_NAMES)),
                   metadata=d.get("metadata", {}))


# ---------------------------------------------------------------------------
# Fitting + metrics (used by tools/build_calibration.py and tests)
# ---------------------------------------------------------------------------


def fit_platt(ratios, labels) -> PlattParams:
    """Fit P(correct)=sigmoid(A*r+B) with Platt 1999 target smoothing."""
    from scipy.optimize import minimize
    r = np.asarray(ratios, dtype=float)
    y = np.asarray(labels, dtype=float)
    n_pos = float(y.sum())
    n_neg = float(len(y) - n_pos)
    t_plus = (n_pos + 1.0) / (n_pos + 2.0)
    t_minus = 1.0 / (n_neg + 2.0)
    t = np.where(y > 0.5, t_plus, t_minus)

    def nll(params):
        a, b = params
        z = np.clip(a * r + b, -60.0, 60.0)
        p = 1.0 / (1.0 + np.exp(-z))
        p = np.clip(p, 1e-12, 1.0 - 1e-12)
        return -np.sum(t * np.log(p) + (1.0 - t) * np.log(1.0 - p))

    res = minimize(nll, np.array([1.0, 0.0]), method="L-BFGS-B")
    return PlattParams(float(res.x[0]), float(res.x[1]))


def fit_gnb(X, labels, *, eps: float = 1e-6) -> GNBParams:
    """Fit a Gaussian Naive Bayes calibrator over the feature matrix (rows = samples)."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(labels, dtype=float)
    n_feat = X.shape[1]
    means: dict = {}
    variances: dict = {}
    log_priors: dict = {}
    for cls, val in (("correct", 1.0), ("incorrect", 0.0)):
        mask = y == val
        if mask.sum() == 0:
            means[cls] = [0.0] * n_feat
            variances[cls] = [1.0] * n_feat
            log_priors[cls] = -60.0
            continue
        means[cls] = X[mask].mean(axis=0).tolist()
        variances[cls] = (X[mask].var(axis=0) + eps).tolist()
        log_priors[cls] = float(np.log(mask.sum() / len(y)))
    return GNBParams(means, variances, log_priors)


def brier_score(probs, labels) -> float:
    p = np.asarray(probs, dtype=float)
    y = np.asarray(labels, dtype=float)
    return float(np.mean((p - y) ** 2))


def expected_calibration_error(probs, labels, n_bins: int = 10) -> float:
    p = np.asarray(probs, dtype=float)
    y = np.asarray(labels, dtype=float)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0
    n = len(p)
    for i in range(n_bins):
        hi_inclusive = i == n_bins - 1
        mask = (p >= edges[i]) & (p <= edges[i + 1] if hi_inclusive else p < edges[i + 1])
        if mask.sum() == 0:
            continue
        ece += mask.sum() / n * abs(float(y[mask].mean()) - float(p[mask].mean()))
    return float(ece)


# ---------------------------------------------------------------------------
# Bundled-artifact loader
# ---------------------------------------------------------------------------


@lru_cache(maxsize=4)
def load_calibrator(method: str = "platt") -> Calibrator | None:
    """Load the bundled calibration artifact, or None if it is not present."""
    path = files("checkmsg.refdata.data") / f"calibration_{method}.json"
    if not path.is_file():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    if tuple(data.get("feature_names", ())) != FEATURE_NAMES:
        return None  # schema drift → refuse a stale artifact
    return Calibrator.from_dict(data)
