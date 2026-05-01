"""High-level pipeline: simulate transmission/scattering, reconstruct, build MuonImage."""

from __future__ import annotations

import numpy as np

from checkmsg.muon.forward import (
    ScatteringSinogram,
    TransmissionSinogram,
    simulate_scattering,
    simulate_transmission,
    trace_one,
)
from checkmsg.muon.image import MuonImage
from checkmsg.muon.muonic_xray import simulate_muonic_xray
from checkmsg.muon.reconstruct import reconstruct_sart
from checkmsg.muon.source import MuonSource
from checkmsg.muon.voxel import VoxelGrid


def _reconstruct_to_3d(sino: np.ndarray, angles_rad: np.ndarray,
                        grid_shape: tuple[int, int, int],
                        iterations: int = 5) -> np.ndarray:
    """Reconstruct a 3-D volume by SART-reconstructing each z-layer independently.

    The sinogram is shape (n_proj, n_pixels) and represents an integrated quantity
    along z; we use it as the slice-z=mid reconstruction and broadcast across z.
    For more rigour, multi-layer sinograms could be separated per detector row.
    """
    image_size = max(grid_shape[0], grid_shape[1])
    slice_2d = reconstruct_sart(sino, angles_rad, image_size, iterations=iterations)
    # Broadcast across the z dimension; this is sufficient for the symmetric
    # parallel-beam geometry used by the example scripts.
    nz = grid_shape[2]
    return np.broadcast_to(slice_2d[..., None], (image_size, image_size, nz)).copy()


def analyze(
    grid: VoxelGrid,
    source: MuonSource,
    *,
    transmission: bool = True,
    scattering: bool = False,
    muonic_xray: bool = False,
    n_projections: int = 30,
    pixels_per_side: int = 32,
    muons_per_ray: int = 16,
    sart_iterations: int = 5,
    rng_seed: int = 0,
) -> MuonImage:
    """Run the full muon-imaging pipeline against a known voxel grid.

    Returns a `MuonImage` with whichever subset of (density_map, scattering_density_map,
    muonic_xray_spectrum) the caller requested. Intended for benchmarking and for the
    curriculum example — real instrument data would feed sinograms directly into the
    reconstruction helpers, skipping the simulator.
    """
    density_map: np.ndarray | None = None
    scatter_map: np.ndarray | None = None
    xray: object = None

    if transmission:
        sino_t: TransmissionSinogram = simulate_transmission(
            grid, source,
            n_projections=n_projections,
            pixels_per_side=pixels_per_side,
            muons_per_ray=muons_per_ray,
            rng_seed=rng_seed,
        )
        # Convert survival fraction to opacity μx = -ln(survival).
        survival = sino_t.transmission()
        opacity = -np.log(np.clip(survival, 1e-6, 1.0))
        density_map = _reconstruct_to_3d(opacity, sino_t.angles_rad,
                                          grid.shape, iterations=sart_iterations)

    if scattering:
        sino_s: ScatteringSinogram = simulate_scattering(
            grid, source,
            n_projections=n_projections,
            pixels_per_side=pixels_per_side,
            muons_per_ray=muons_per_ray,
            rng_seed=rng_seed + 1,
        )
        scatter_map = _reconstruct_to_3d(sino_s.rms_mrad ** 2,  # squared-angle is additive
                                          sino_s.angles_rad,
                                          grid.shape, iterations=sart_iterations)

    if muonic_xray:
        # Stopping muons need to (a) actually enter the grid and (b) lose all
        # their energy inside it. Beam aimed at -x with random y/z offsets so
        # the muons sample the cross-section.
        rng = np.random.default_rng(rng_seed + 2)
        stopping_source = MuonSource(
            mean_momentum_MeV=min(source.mean_momentum_MeV, 50.0),
            momentum_FWHM_MeV=2.0,
            flux_per_s=source.flux_per_s,
            direction=(-1.0, 0.0, 0.0),
            angular_spread_mrad=source.angular_spread_mrad,
            polarity="negative",
        )
        n_to_simulate = max(200, source.expected_count(0.001))
        n_to_simulate = min(n_to_simulate, 2000)
        bounds_lo, bounds_hi = grid.bounds_mm()
        centre = 0.5 * (bounds_lo + bounds_hi)
        half = 0.5 * (bounds_hi - bounds_lo)
        diag = float(np.linalg.norm(half))
        momenta, dirs = stopping_source.sample(n_to_simulate, rng=rng)
        stopping_voxels: list[tuple[int, int, int]] = []
        for k in range(n_to_simulate):
            offset = (rng.random(2) - 0.5) * 2.0 * np.array([half[1], half[2]])
            entry = (centre + np.array([1.5 * diag, offset[0], offset[1]]))
            track = trace_one(entry, dirs[k], float(momenta[k]), grid)
            if track.stopping_voxel is not None:
                stopping_voxels.append(track.stopping_voxel)
        xray = simulate_muonic_xray(stopping_voxels, grid, fwhm_keV=10.0,
                                     noise=0.005, seed=rng_seed + 5)

    return MuonImage(
        density_map=density_map,
        scattering_density_map=scatter_map,
        muonic_xray_spectrum=xray,
        grid=grid,
        n_muons=source.expected_count(1.0),
        exposure_s=1.0,
        metadata={
            "n_projections": n_projections,
            "pixels_per_side": pixels_per_side,
            "muons_per_ray": muons_per_ray,
            "sart_iterations": sart_iterations,
        },
    )
