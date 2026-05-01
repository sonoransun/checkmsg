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
