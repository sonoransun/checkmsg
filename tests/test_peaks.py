import numpy as np

from checkmsg.peaks import detect, fit_voigt
from checkmsg.spectrum import Spectrum
from checkmsg.synthetic import PeakSpec, generate


def test_detect_finds_strong_peak():
    s = generate([PeakSpec(1332.0, 1.0, 1.5, 0.5)],
                 np.linspace(100, 1700, 1601), "raman", "cm-1", noise=0.005, seed=0)
    peaks = detect(s, min_snr=10.0)
    assert peaks
    biggest = max(peaks, key=lambda p: p.height)
    assert abs(biggest.position - 1332.0) < 1.0


def test_detect_returns_empty_for_flat_spectrum():
    s = Spectrum(np.linspace(0, 100, 200), np.zeros(200), "raman", "cm-1")
    assert detect(s) == []


def test_voigt_fit_subpixel_accuracy():
    s = generate([PeakSpec(685.4, 1.0, 2.0, 0.7)],
                 np.linspace(500, 900, 401), "raman", "cm-1", noise=0.003, seed=1)
    fit = fit_voigt(s, around=685.0, window=20.0)
    assert fit is not None
    assert abs(fit.position - 685.4) < 0.5


def test_voigt_fit_returns_none_for_short_window():
    s = Spectrum(np.array([1.0, 2.0, 3.0]), np.array([0, 1, 0], dtype=float), "raman", "cm-1")
    assert fit_voigt(s, around=2.0, window=0.1) is None
