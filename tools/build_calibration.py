#!/usr/bin/env python
"""Fit the Platt confidence calibrator on synthetic leave-mineral-out CV.

Generates a training set by self-diagnosing every catalog entry across a few
seeds, labels each by whether the verdict matched the true mineral, fits a Platt
sigmoid, reports ECE/Brier from leave-mineral-out cross-validation, and writes
the human-readable artifact `src/checkmsg/refdata/data/calibration_platt.json`.

Deterministic (fixed seeds; no wall-clock), so `--check` re-fits in memory and
compares to the committed artifact — the CI guard against calibration drift.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from checkmsg import calibrate, minerals
from checkmsg.diagnose import diagnose_profile

ARTIFACT = (Path(__file__).resolve().parents[1]
            / "src" / "checkmsg" / "refdata" / "data" / "calibration_platt.json")
SEEDS = (0, 1, 2)


def build_dataset() -> list[tuple[str, dict, float]]:
    """Return (mineral, feature_dict, label) rows from synthetic self-diagnoses."""
    rows: list[tuple[str, dict, float]] = []
    for name in minerals.names():
        profile = minerals.get(name)
        for seed in SEEDS:
            report = diagnose_profile(profile, seed=seed)
            if report.verdict is None:
                continue
            label = 1.0 if report.verdict == name else 0.0
            rows.append((name, calibrate.features_from_report(report), label))
    return rows


def lmo_cv(rows: list[tuple[str, float, float]]) -> tuple[float, float, int]:
    names = sorted({r[0] for r in rows})
    oof_p: list[float] = []
    oof_y: list[float] = []
    for held in names:
        train = [(f["separation_ratio"], y) for n, f, y in rows if n != held]
        test = [(f["separation_ratio"], y) for n, f, y in rows if n == held]
        if not test or len({y for _, y in train}) < 2:
            continue
        params = calibrate.fit_platt([r for r, _ in train], [y for _, y in train])
        for r, y in test:
            z = max(-60.0, min(60.0, params.A * r + params.B))
            oof_p.append(1.0 / (1.0 + np.exp(-z)))
            oof_y.append(y)
    if not oof_y:
        return 0.0, 0.0, 0
    return (calibrate.expected_calibration_error(oof_p, oof_y),
            calibrate.brier_score(oof_p, oof_y), len(oof_y))


def _meta(ece, brier, n_oof, n_rows):
    return {
        "ece": round(ece, 4), "brier": round(brier, 4),
        "n_samples": n_rows, "n_oof": n_oof,
        "catalog_size": len(minerals.names()), "seeds": list(SEEDS),
        "limits": ("Calibrated on noise-clean synthetic leave-mineral-out "
                   "self-diagnosis; not a real-world probability of correctness."),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fit/check the confidence calibrators.")
    ap.add_argument("--check", action="store_true",
                    help="re-fit Platt and compare to the committed artifact (no write).")
    args = ap.parse_args(argv)

    rows = build_dataset()
    ratios = [f["separation_ratio"] for _, f, _ in rows]
    labels = [y for _, _, y in rows]
    platt = calibrate.fit_platt(ratios, labels)

    if args.check:
        if not ARTIFACT.is_file():
            print("missing artifact", ARTIFACT)
            return 1
        old = json.loads(ARTIFACT.read_text(encoding="utf-8"))["params"]
        drift = abs(old["A"] - platt.A) > 1e-3 or abs(old["B"] - platt.B) > 1e-3
        print("DRIFT" if drift else "OK")
        return 1 if drift else 0

    ece, brier, n_oof = lmo_cv(rows)
    meta = _meta(ece, brier, n_oof, len(rows))

    # Platt artifact
    cal = calibrate.Calibrator(method="platt", platt=platt, metadata=meta)
    ARTIFACT.write_text(json.dumps(cal.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {ARTIFACT.name} | A={platt.A:.4f} B={platt.B:.4f} "
          f"ECE={ece:.4f} Brier={brier:.4f} n={len(rows)}")

    # GNB artifact (full feature vector)
    X = [[f[k] for k in calibrate.FEATURE_NAMES] for _, f, _ in rows]
    gnb = calibrate.fit_gnb(X, labels)
    gnb_cal = calibrate.Calibrator(method="gnb", gnb=gnb, metadata=meta)
    gnb_path = ARTIFACT.with_name("calibration_gnb.json")
    gnb_path.write_text(json.dumps(gnb_cal.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {gnb_path.name} (Gaussian Naive Bayes over {len(calibrate.FEATURE_NAMES)} features)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
