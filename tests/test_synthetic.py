import numpy as np

from checkmsg.synthetic import PeakSpec, amorphous_halo, generate, linear_baseline, voigt_pseudo


def test_voigt_pseudo_centered_max():
    x = np.linspace(0, 100, 1001)
    y = voigt_pseudo(x, center=50.0, sigma=2.0, gamma=1.0, amplitude=1.0)
    assert abs(x[np.argmax(y)] - 50.0) < 0.1


def test_generate_reproducible_with_seed():
    args = ([PeakSpec(50, 1.0, 1.0, 0.5)], np.linspace(0, 100, 1001))
    a = generate(*args, technique="raman", units="cm-1", noise=0.01, seed=99)
    b = generate(*args, technique="raman", units="cm-1", noise=0.01, seed=99)
    assert np.array_equal(a.intensity, b.intensity)


def test_generate_with_baseline_offsets():
    axis = np.linspace(0, 10, 200)
    base = linear_baseline(axis, slope=0.5, intercept=2.0)
    s = generate([], axis, "raman", "cm-1", noise=0.0, baseline=base, seed=0)
    assert s.intensity[0] == 2.0
    assert s.intensity[-1] == 2.0 + 0.5 * (axis[-1] - axis[0])


def test_amorphous_halo_peak_position():
    axis = np.linspace(0, 1500, 1500)
    halo = amorphous_halo(axis, center=500, width=80, amplitude=1.0)
    assert abs(axis[np.argmax(halo)] - 500) < 1.5
