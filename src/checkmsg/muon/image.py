"""`MuonImage` — reconstructed 3-D density / scattering maps + stopping-muon spectrum."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from checkmsg.muon.voxel import VoxelGrid
from checkmsg.spectrum import Spectrum


@dataclass
class MuonImage:
    """Reconstructed muon image bundle.

    Any subset of (`density_map`, `scattering_density_map`, `muonic_xray_spectrum`)
    may be `None`. The grid records the reconstruction geometry — a slice through
    the grid using the inferred density should match the ground-truth density to
    within the simulator's error budget.
    """

    density_map: np.ndarray | None
    scattering_density_map: np.ndarray | None
    muonic_xray_spectrum: Spectrum | None
    grid: VoxelGrid
    n_muons: int
    exposure_s: float
    metadata: dict = field(default_factory=dict)

    def central_slice(self, axis: int = 2,
                       which: str = "density") -> np.ndarray | None:
        """Return the central 2-D slice of the chosen 3-D map.

        `which` ∈ {"density", "scattering"}.
        """
        if which == "density":
            arr = self.density_map
        elif which == "scattering":
            arr = self.scattering_density_map
        else:
            raise ValueError(f"unknown map {which!r}; choose 'density' or 'scattering'")
        if arr is None:
            return None
        if axis < 0 or axis > 2:
            raise ValueError("axis must be 0, 1, or 2")
        mid = arr.shape[axis] // 2
        return np.take(arr, mid, axis=axis)
