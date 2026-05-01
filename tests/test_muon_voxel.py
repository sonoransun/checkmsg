import numpy as np
import pytest

from checkmsg.muon.voxel import VoxelGrid


def test_filled_grid_uniform_material():
    g = VoxelGrid.filled((4, 4, 4), "corundum")
    assert g.shape == (4, 4, 4)
    assert g.materials[0, 0, 0].name == "corundum"
    assert g.densities[0, 0, 0] == pytest.approx(4.0)


def test_set_box_overrides_material_and_density():
    g = VoxelGrid.filled((6, 6, 6), "corundum")
    g.set_box((1, 1, 1), (3, 3, 3), "platinum")
    assert g.materials[1, 1, 1].name == "platinum"
    assert g.materials[3, 3, 3].name == "corundum"   # exclusive upper bound


def test_set_sphere_localised():
    g = VoxelGrid.filled((10, 10, 10), "polymer")
    g.set_sphere((5, 5, 5), 2.0, "diamond")
    assert g.materials[5, 5, 5].name == "diamond"
    assert "polymer" in g.materials[0, 0, 0].name.lower()


def test_density_array_round_trip():
    g = VoxelGrid.filled((4, 4, 4), "corundum")
    arr = g.density_array()
    assert arr.shape == (4, 4, 4)
    assert np.allclose(arr, 4.0)


def test_z_eff_array_per_voxel():
    g = VoxelGrid.filled((4, 4, 4), "corundum")
    g.set_box((1, 1, 1), (2, 2, 2), "platinum")
    z = g.z_eff_array()
    assert z[1, 1, 1] == pytest.approx(78.0)
    assert z[0, 0, 0] == pytest.approx(11.31)


def test_void_density_is_zero():
    g = VoxelGrid.filled((4, 4, 4), "corundum")
    g.set_box((0, 0, 0), (2, 2, 2), "vacuum", density_g_cc=0.0)
    assert g.densities[0, 0, 0] == 0.0


def test_bounds_match_origin_and_spacing():
    g = VoxelGrid.filled((10, 10, 10), "water",
                         spacing_mm=(2.0, 2.0, 2.0),
                         origin_mm=(5.0, 5.0, 5.0))
    lo, hi = g.bounds_mm()
    assert lo[0] == pytest.approx(4.0)
    assert hi[0] == pytest.approx(24.0)


def test_shape_mismatch_raises():
    with pytest.raises(ValueError):
        VoxelGrid(materials=np.empty((3, 3, 3), dtype=object),
                  densities=np.zeros((4, 4, 4)))
