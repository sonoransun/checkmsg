"""Forward projection: voxel-path lengths, transmission, and scattering sanity checks."""

import numpy as np

from checkmsg.muon import MuonSource, VoxelGrid, simulate_scattering, simulate_transmission, trace_one
from checkmsg.muon.forward import _voxel_path


def _uniform_grid():
    return VoxelGrid.filled((10, 10, 10), "corundum", spacing_mm=(2.0, 2.0, 2.0))


def test_voxel_path_total_length_matches_geometric_distance():
    """Through a uniform grid, sum(segment_cm) along the entry direction must match the
    straight-line distance through the grid."""
    g = _uniform_grid()
    # Grid spans -1 to 19 mm on each axis; total length on axial ray is 20 mm = 2.0 cm
    entry = np.array([30.0, 9.0, 9.0])
    direction = np.array([-1.0, 0.0, 0.0])
    path = _voxel_path(entry, direction, g)
    assert path
    total_cm = sum(seg for _, _, _, seg in path)
    assert abs(total_cm - 2.0) < 0.01


def test_voxel_path_outside_returns_empty():
    """A ray that misses the grid entirely returns an empty path."""
    g = _uniform_grid()
    entry = np.array([100.0, 100.0, 100.0])
    direction = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)  # away from the grid
    path = _voxel_path(entry, direction, g)
    assert path == []


def test_trace_one_through_corundum_loses_expected_energy():
    """A 300 MeV/c muon through 2 cm of corundum (8 g/cm² mass thickness, MIP-ish dE/dx) should
    lose ~16 MeV (within 50 % — physics is approximate)."""
    g = _uniform_grid()
    track = trace_one(np.array([30.0, 9.0, 9.0]),
                      np.array([-1.0, 0.0, 0.0]),
                      300.0, g)
    assert track.transmitted is True
    assert 5.0 < track.energy_loss_MeV < 25.0


def test_trace_one_low_momentum_stops():
    """A 30 MeV/c muon (very short range) should stop in the first cm of corundum."""
    g = _uniform_grid()
    track = trace_one(np.array([30.0, 9.0, 9.0]),
                      np.array([-1.0, 0.0, 0.0]),
                      30.0, g)
    assert track.transmitted is False
    assert track.stopping_voxel is not None


def test_high_Z_insert_increases_scattering():
    """Replacing a corundum block with platinum should significantly raise the RMS scatter."""
    g_clean = _uniform_grid()
    g_pt = _uniform_grid()
    g_pt.set_box((4, 4, 4), (6, 6, 6), "platinum")
    src = MuonSource(mean_momentum_MeV=300.0, flux_per_s=1e9)
    s_clean = simulate_scattering(g_clean, src, n_projections=4, pixels_per_side=8,
                                  muons_per_ray=4, rng_seed=0)
    s_pt = simulate_scattering(g_pt, src, n_projections=4, pixels_per_side=8,
                               muons_per_ray=4, rng_seed=0)
    assert s_pt.rms_mrad.max() > s_clean.rms_mrad.max() * 1.5


def test_transmission_drops_with_low_momentum():
    """At 30 MeV/c, far more muons stop in the volume than at 500 MeV/c."""
    g = _uniform_grid()
    src_low = MuonSource(mean_momentum_MeV=30.0, flux_per_s=1e9)
    src_high = MuonSource(mean_momentum_MeV=500.0, flux_per_s=1e9)
    t_low = simulate_transmission(g, src_low, n_projections=4, pixels_per_side=8,
                                  muons_per_ray=4, rng_seed=0)
    t_high = simulate_transmission(g, src_high, n_projections=4, pixels_per_side=8,
                                   muons_per_ray=4, rng_seed=0)
    assert t_low.transmission().mean() < t_high.transmission().mean()
