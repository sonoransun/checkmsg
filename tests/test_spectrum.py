import numpy as np
import pytest

from checkmsg.spectrum import Spectrum


def test_basic_construction():
    s = Spectrum(np.array([1.0, 2.0, 3.0]), np.array([10.0, 20.0, 30.0]), "raman", "cm-1")
    assert len(s) == 3
    assert s.technique == "raman"
    assert s.units == "cm-1"


def test_axis_intensity_shape_mismatch_raises():
    with pytest.raises(ValueError, match="shape mismatch"):
        Spectrum(np.array([1.0, 2.0]), np.array([1.0, 2.0, 3.0]), "raman", "cm-1")


def test_axis_resorted_when_unordered():
    s = Spectrum(np.array([3.0, 1.0, 2.0]), np.array([30.0, 10.0, 20.0]), "raman", "cm-1")
    assert list(s.axis) == [1.0, 2.0, 3.0]
    assert list(s.intensity) == [10.0, 20.0, 30.0]


def test_slice_returns_subset():
    x = np.linspace(0, 100, 101)
    y = np.arange(101, dtype=float)
    s = Spectrum(x, y, "raman", "cm-1").slice(20, 40)
    assert s.axis[0] == 20.0
    assert s.axis[-1] == 40.0
    assert len(s) == 21


def test_normalize_modes():
    x = np.linspace(0, 10, 11)
    y = np.arange(11, dtype=float)
    s = Spectrum(x, y, "raman", "cm-1")
    assert s.normalize("max").intensity.max() == pytest.approx(1.0)
    n = s.normalize("minmax").intensity
    assert n.min() == pytest.approx(0.0) and n.max() == pytest.approx(1.0)


def test_with_intensity_preserves_axis():
    x = np.linspace(0, 1, 5)
    s = Spectrum(x, np.zeros(5), "xrf", "keV", {"sample": "A"})
    s2 = s.with_intensity(np.ones(5))
    assert np.array_equal(s2.axis, x)
    assert s2.metadata == {"sample": "A"}
    assert np.array_equal(s2.intensity, np.ones(5))
