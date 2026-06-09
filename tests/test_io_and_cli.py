import subprocess
import sys

import numpy as np

from checkmsg.io import read_csv, write_csv
from checkmsg.spectrum import Spectrum


def test_read_write_round_trip(tmp_path):
    s = Spectrum(np.linspace(0, 100, 101), np.arange(101, dtype=float), "raman", "cm-1")
    p = tmp_path / "spec.csv"
    write_csv(s, p)
    s2 = read_csv(p, technique="raman", units="cm-1")
    assert np.allclose(s2.axis, s.axis)
    assert np.allclose(s2.intensity, s.intensity)


def test_cli_analyze_raman(tmp_path):
    # Build a tiny synthetic Raman file and ensure the CLI runs.
    s = Spectrum(np.linspace(100, 1700, 200), np.zeros(200), "raman", "cm-1")
    s.intensity[100] = 1.0
    p = tmp_path / "x.csv"
    write_csv(s, p)
    proc = subprocess.run(
        [sys.executable, "-m", "checkmsg.cli", "analyze", "raman", str(p)],
        capture_output=True, text=True, timeout=60,
    )
    assert proc.returncode == 0
    assert "Raman peaks" in proc.stdout


def test_cli_glossary_lookup():
    proc = subprocess.run(
        [sys.executable, "-m", "checkmsg.cli", "glossary", "EPR"],
        capture_output=True, text=True, timeout=60,
    )
    assert proc.returncode == 0
    assert "electron paramagnetic resonance" in proc.stdout

    bad = subprocess.run(
        [sys.executable, "-m", "checkmsg.cli", "glossary", "zzznope"],
        capture_output=True, text=True, timeout=60,
    )
    assert bad.returncode == 1


def test_cli_diagnose_novice_tier(tmp_path):
    # A diamond-like Raman line; diagnose at the novice tier prints a plain verdict.
    from checkmsg.synthetic import PeakSpec, generate
    axis = np.linspace(100, 1700, 1601)
    spec = generate([PeakSpec(1332.0, 1.0, 1.5, 0.5)], axis,
                    technique="raman", units="cm-1", noise=0.005, seed=42)
    p = tmp_path / "diamond.csv"
    write_csv(spec, p)
    proc = subprocess.run(
        [sys.executable, "-m", "checkmsg.cli", "diagnose", f"raman:{p}", "--tier", "novice"],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0
    assert "Most likely:" in proc.stdout
    assert "confidence" in proc.stdout.lower()
    assert "Reasoning trace" not in proc.stdout  # novice omits the expert scaffolding
