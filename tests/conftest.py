from __future__ import annotations

import numpy as np
import pytest


@pytest.fixture(autouse=True)
def isolate_cache(tmp_path, monkeypatch):
    """Force every test into a private CHECKMSG_CACHE and forbid network by default."""
    monkeypatch.setenv("CHECKMSG_CACHE", str(tmp_path / "cache"))
    monkeypatch.setenv("CHECKMSG_OFFLINE", "1")
    yield


@pytest.fixture
def diamond_axis() -> np.ndarray:
    return np.linspace(100.0, 1700.0, 1601)


@pytest.fixture
def diamond_spectrum(diamond_axis):
    from checkmsg.synthetic import PeakSpec, generate
    return generate(
        [PeakSpec(1332.0, 1.0, 1.5, 0.5)],
        diamond_axis, technique="raman", units="cm-1", noise=0.005, seed=42,
    )


@pytest.fixture
def quartz_xrf_spectrum():
    from checkmsg.synthetic import PeakSpec, generate
    axis = np.linspace(0.5, 12.0, 4096)
    peaks = [
        PeakSpec(1.740, 1.0, 0.045, 0.025),  # Si Ka
        PeakSpec(8.048, 0.4, 0.045, 0.025),  # Cu Ka
    ]
    return generate(peaks, axis, technique="xrf", units="keV", noise=0.002, seed=7)
