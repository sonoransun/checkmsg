import numpy as np

from checkmsg.preprocess import als_baseline, area_normalize, min_max_normalize, savgol, snip_baseline


def test_als_baseline_recovers_linear_drift():
    x = np.linspace(0, 100, 1000)
    drift = 0.05 * x + 2.0
    peak = np.exp(-((x - 50) ** 2) / (2 * 1.5 ** 2))
    y = drift + peak
    bl = als_baseline(y, lam=1e5, p=0.001)
    # baseline should track the drift, far from the peak region
    bg_region = (np.abs(x - 50) > 20)
    assert np.allclose(bl[bg_region], drift[bg_region], atol=0.5)


def test_snip_baseline_under_peaks():
    x = np.linspace(0, 100, 2000)
    bg = 1.0 + 0.5 * np.exp(-((x - 60) ** 2) / (40 ** 2))
    peak = 5.0 * np.exp(-((x - 30) ** 2) / (2 * 0.6 ** 2))
    y = bg + peak
    bl = snip_baseline(y, iterations=40)
    # SNIP must not eat the peak: residual at peak is large
    assert (y - bl)[np.argmax(peak)] > 4.0


def test_savgol_preserves_size_and_smooths():
    rng = np.random.default_rng(0)
    y = rng.normal(0, 1, 200)
    smoothed = savgol(y, 11, 3)
    assert smoothed.shape == y.shape
    assert smoothed.std() < y.std()


def test_min_max_normalize_handles_constant():
    assert np.array_equal(min_max_normalize(np.ones(5)), np.zeros(5))


def test_area_normalize_unit_area():
    x = np.linspace(0, 1, 100)
    y = np.exp(-((x - 0.5) ** 2) / 0.01)
    n = area_normalize(y, x)
    assert abs(np.trapezoid(np.abs(n), x) - 1.0) < 1e-9
