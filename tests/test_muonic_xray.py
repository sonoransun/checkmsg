"""Muonic-atom K_α tabulation + spectrum builder tests."""

import numpy as np
import pytest

from checkmsg.muon import VoxelGrid, simulate_muonic_xray
from checkmsg.muon.muonic_xray import muonic_kalpha_keV
from checkmsg.refdata.muon_data import MUONIC_KALPHA_keV


def test_au_kalpha_matches_table():
    assert muonic_kalpha_keV("Au") == pytest.approx(6019.0)


def test_kalpha_scales_with_Z_squared_approximately():
    """Z² scaling: K_α(Z=14) / K_α(Z=6) ≈ (14/6)² ≈ 5.4 (ignoring corrections)."""
    e_si = muonic_kalpha_keV("Si")
    e_c = muonic_kalpha_keV("C")
    ratio = e_si / e_c
    target = (14.0 / 6.0) ** 2
    assert 0.6 * target < ratio < 1.4 * target


def test_muonic_xray_unknown_element_raises():
    with pytest.raises(KeyError):
        muonic_kalpha_keV("Unobtainium")


def test_table_has_at_least_40_elements():
    assert len(MUONIC_KALPHA_keV) >= 40


def test_simulate_muonic_xray_peak_at_au_kalpha():
    """A bunch of muons stopping in gold should produce a spectrum with a peak near 6019 keV."""
    g = VoxelGrid.filled((4, 4, 4), "gold")
    stopping_voxels = [(2, 2, 2)] * 200  # 200 muons all stop in the gold voxel
    spec = simulate_muonic_xray(stopping_voxels, g, fwhm_keV=10.0, noise=0.0)
    assert spec.technique == "muon-xray"
    assert spec.units == "keV"
    peak_idx = int(np.argmax(spec.intensity))
    peak_energy = float(spec.axis[peak_idx])
    assert abs(peak_energy - 6019.0) < 30.0


def test_simulate_muonic_xray_empty_input_returns_flat_spectrum():
    g = VoxelGrid.filled((4, 4, 4), "gold")
    spec = simulate_muonic_xray([], g, fwhm_keV=10.0, noise=0.0)
    assert spec.intensity.shape[0] > 0
    assert float(np.max(np.abs(spec.intensity))) < 1e-6
