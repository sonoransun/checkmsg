"""Forward projection: simulate muons traversing a voxel grid.

Each muon is treated as a straight line through the grid (small-angle
approximation; valid when the total scattering angle is far below 100 mrad).
For each ray we accumulate:

  - **Mass thickness** Σ ρ_i × Δℓ_i   (g/cm²) — drives Bethe-Bloch energy loss.
  - **Radiation-length thickness** Σ Δℓ_i / X₀_i — drives Highland multi-scattering.

Stopping is detected by comparing total kinetic energy lost to the muon's incoming
kinetic energy. Stopped muons report their stopping voxel for muonic-X-ray mode.

Two top-level helpers:
  - `simulate_transmission` — counts surviving muons per pixel of an exit detector.
  - `simulate_scattering` — computes the RMS scattering angle along each track.
Both use the same per-ray Siddon-style traversal under the hood.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from checkmsg.muon.physics import (
    bethe_bloch_dE_dx,
    highland_scattering_rms_mrad,
)
from checkmsg.muon.source import MuonSource
from checkmsg.muon.voxel import VoxelGrid
from checkmsg.refdata.muon_data import M_MU_MeV


@dataclass
class MuonTrack:
    """A single muon trajectory through the grid."""
    entry_pos_mm: np.ndarray
    entry_momentum_MeV: float
    entry_dir: np.ndarray
    exit_pos_mm: np.ndarray | None
    exit_dir: np.ndarray | None
    energy_loss_MeV: float
    scattering_angle_mrad: float
    stopping_voxel: tuple[int, int, int] | None = None
    transmitted: bool = True


@dataclass
class TransmissionSinogram:
    """Per-projection-angle, per-pixel counts of transmitted muons.

    Used by `reconstruct_*` to build the density map.
    """
    angles_rad: np.ndarray  # shape (n_proj,)
    counts: np.ndarray  # shape (n_proj, n_pixels), int — surviving muons
    incident_per_ray: int  # number of muons per (angle, pixel) bin in the simulation
    pixel_size_mm: float
    grid_shape: tuple[int, int, int] = (0, 0, 0)
    grid_spacing_mm: tuple[float, float, float] = (1.0, 1.0, 1.0)

    def transmission(self) -> np.ndarray:
        """Surviving fraction per (angle, pixel)."""
        return self.counts.astype(float) / max(self.incident_per_ray, 1)


@dataclass
class ScatteringSinogram:
    """Per-pixel RMS scattering angle, used to build the scattering-density image."""
    angles_rad: np.ndarray
    rms_mrad: np.ndarray  # shape (n_proj, n_pixels)
    pixel_size_mm: float
    grid_shape: tuple[int, int, int] = (0, 0, 0)
    grid_spacing_mm: tuple[float, float, float] = (1.0, 1.0, 1.0)
    metadata: dict = field(default_factory=dict)


# --- Ray traversal --------------------------------------------------------------


def _voxel_path(entry: np.ndarray, direction: np.ndarray, grid: VoxelGrid
                ) -> list[tuple[int, int, int, float]]:
    """3D-DDA traversal returning (ix, iy, iz, segment_length_cm) per voxel hit.

    Lengths are returned in centimetres so they can be combined with material
    densities (g/cc) and X_0 (g/cm²) directly. Handles entry from outside the
    grid by first clipping the ray to the grid bounding box (slab method).
    """
    spacing = np.asarray(grid.spacing_mm, dtype=float)
    origin = np.asarray(grid.origin_mm, dtype=float)
    shape = np.asarray(grid.shape)
    bounds_lo, bounds_hi = grid.bounds_mm()

    # Slab-method ray-box intersection.
    with np.errstate(divide="ignore", invalid="ignore"):
        dir_safe = np.where(np.abs(direction) > 1e-12, direction, 1e-12)
        t_lo = (bounds_lo - entry) / dir_safe
        t_hi = (bounds_hi - entry) / dir_safe
    t_enter = float(np.maximum(np.minimum(t_lo, t_hi), 0.0).max())
    t_exit = float(np.minimum(np.maximum(t_lo, t_hi), 1e18).min())
    if t_exit <= t_enter:
        return []

    pos = entry + direction * (t_enter + 1e-6)
    vox = (pos - origin) / spacing + 0.5
    cur = np.floor(vox).astype(int)
    if not (np.all(cur >= 0) and np.all(cur < shape)):
        return []

    step = np.where(direction > 0, 1, -1)
    next_vox_boundary = cur + (step > 0).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_max = np.where(np.abs(direction) > 1e-12,
                         (next_vox_boundary - vox) * spacing / direction,
                         np.inf)
        t_delta = np.where(np.abs(direction) > 1e-12,
                           spacing / np.abs(direction),
                           np.inf)
    # Reject any negative t_max (numerical clamp at the entry edge).
    t_max = np.where(t_max < 0, np.inf, t_max)

    out: list[tuple[int, int, int, float]] = []
    prev_t = 0.0
    while np.all(cur >= 0) and np.all(cur < shape):
        axis = int(np.argmin(t_max))
        t_now = float(t_max[axis])
        if t_now == np.inf:
            break
        seg_mm = t_now - prev_t   # distance through THIS voxel
        out.append((int(cur[0]), int(cur[1]), int(cur[2]), seg_mm / 10.0))
        prev_t = t_now
        cur[axis] += step[axis]
        t_max[axis] += t_delta[axis]
    return out


# --- High-level simulators -------------------------------------------------------


def trace_one(entry: np.ndarray, direction: np.ndarray, momentum_MeV: float,
              grid: VoxelGrid) -> MuonTrack:
    """Trace a single muon through `grid`, accumulating energy loss + scattering.

    Returns a `MuonTrack` recording where (and whether) the muon exited.
    """
    direction = direction / np.linalg.norm(direction)
    path = _voxel_path(entry, direction, grid)
    p = momentum_MeV
    E_kin_in = math.hypot(p, M_MU_MeV) - M_MU_MeV
    E_kin_remaining = E_kin_in
    rms_total_squared = 0.0
    last_pos = entry.copy()
    stopping_vox: tuple[int, int, int] | None = None

    for ix, iy, iz, seg_cm in path:
        material = grid.materials[ix, iy, iz]
        density = grid.densities[ix, iy, iz]
        if density <= 0 or material.X_0_g_cm2 == float("inf"):
            last_pos = last_pos + direction * (seg_cm * 10.0)
            continue
        x_g_cm2 = density * seg_cm
        # Energy loss
        dEdx = bethe_bloch_dE_dx(p, material)
        loss_MeV = dEdx * x_g_cm2
        if loss_MeV >= E_kin_remaining:
            stopping_vox = (ix, iy, iz)
            E_kin_remaining = 0.0
            break
        E_kin_remaining -= loss_MeV
        p = math.sqrt((E_kin_remaining + M_MU_MeV) ** 2 - M_MU_MeV ** 2)
        # Scattering (root-sum-square accumulation)
        rms_local = highland_scattering_rms_mrad(p, x_g_cm2, material.X_0_g_cm2)
        rms_total_squared += rms_local ** 2
        last_pos = last_pos + direction * (seg_cm * 10.0)

    transmitted = stopping_vox is None
    return MuonTrack(
        entry_pos_mm=entry.copy(),
        entry_momentum_MeV=momentum_MeV,
        entry_dir=direction.copy(),
        exit_pos_mm=last_pos.copy() if transmitted else None,
        exit_dir=direction.copy() if transmitted else None,
        energy_loss_MeV=E_kin_in - E_kin_remaining,
        scattering_angle_mrad=math.sqrt(rms_total_squared),
        stopping_voxel=stopping_vox,
        transmitted=transmitted,
    )


def simulate_transmission(grid: VoxelGrid, source: MuonSource, *,
                          n_projections: int = 30,
                          pixels_per_side: int = 32,
                          muons_per_ray: int = 32,
                          rng_seed: int = 0) -> TransmissionSinogram:
    """Sweep the source around the grid (z-axis tomography) and count survivors per
    detector pixel for each projection angle.

    The grid is assumed to be roughly centred at `grid.origin_mm`; rays are launched
    from the +x boundary at each projection angle, a flat detector sits on the −x side.
    The `rng_seed` parameter is reserved for future stochastic noise additions.
    """
    angles = np.linspace(0.0, math.pi, n_projections, endpoint=False)
    bounds_lo, bounds_hi = grid.bounds_mm()
    centre = 0.5 * (bounds_lo + bounds_hi)
    half_size = 0.5 * (bounds_hi - bounds_lo)
    diag = float(np.linalg.norm(half_size))

    pixel_size = (2 * diag) / pixels_per_side
    counts = np.zeros((n_projections, pixels_per_side), dtype=np.int32)

    # The "v" axis runs across the projection — pixels along it.
    for ai, angle in enumerate(angles):
        # Beam direction in the xy-plane at this angle.
        n_dir = np.array([math.cos(angle), math.sin(angle), 0.0])
        # Pixel positions along the perpendicular axis.
        perp = np.array([-math.sin(angle), math.cos(angle), 0.0])
        for pi in range(pixels_per_side):
            # Pixel offset from the centre.
            offset = (pi - pixels_per_side / 2 + 0.5) * pixel_size
            entry = centre + n_dir * (1.5 * diag) + perp * offset
            # Direction is opposite of n_dir (going through the grid toward −n_dir).
            direction = -n_dir
            survivors = 0
            for _ in range(muons_per_ray):
                p_MeV = source.mean_momentum_MeV  # mean only; spread irrelevant for transmission
                track = trace_one(entry.copy(), direction, p_MeV, grid)
                if track.transmitted:
                    survivors += 1
            counts[ai, pi] = survivors
    return TransmissionSinogram(
        angles_rad=angles,
        counts=counts,
        incident_per_ray=muons_per_ray,
        pixel_size_mm=pixel_size,
        grid_shape=grid.shape,
        grid_spacing_mm=tuple(grid.spacing_mm),
    )


def simulate_scattering(grid: VoxelGrid, source: MuonSource, *,
                        n_projections: int = 30,
                        pixels_per_side: int = 32,
                        muons_per_ray: int = 32,
                        rng_seed: int = 0) -> ScatteringSinogram:
    """Per-pixel RMS scattering angle (mrad) for the same beam geometry as
    `simulate_transmission`. Sensitive to high-Z material via 1/X_0 scaling.
    `rng_seed` is reserved for future stochastic noise additions."""
    angles = np.linspace(0.0, math.pi, n_projections, endpoint=False)
    bounds_lo, bounds_hi = grid.bounds_mm()
    centre = 0.5 * (bounds_lo + bounds_hi)
    half_size = 0.5 * (bounds_hi - bounds_lo)
    diag = float(np.linalg.norm(half_size))
    pixel_size = (2 * diag) / pixels_per_side
    rms = np.zeros((n_projections, pixels_per_side), dtype=float)

    for ai, angle in enumerate(angles):
        n_dir = np.array([math.cos(angle), math.sin(angle), 0.0])
        perp = np.array([-math.sin(angle), math.cos(angle), 0.0])
        for pi in range(pixels_per_side):
            offset = (pi - pixels_per_side / 2 + 0.5) * pixel_size
            entry = centre + n_dir * (1.5 * diag) + perp * offset
            direction = -n_dir
            sumsq = 0.0
            n_used = 0
            for _ in range(muons_per_ray):
                p_MeV = source.mean_momentum_MeV
                track = trace_one(entry.copy(), direction, p_MeV, grid)
                if track.transmitted:
                    sumsq += track.scattering_angle_mrad ** 2
                    n_used += 1
            if n_used > 0:
                rms[ai, pi] = math.sqrt(sumsq / n_used)
    return ScatteringSinogram(
        angles_rad=angles,
        rms_mrad=rms,
        pixel_size_mm=pixel_size,
        grid_shape=grid.shape,
        grid_spacing_mm=tuple(grid.spacing_mm),
    )
