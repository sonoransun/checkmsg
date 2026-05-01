"""Idealised on-demand muon source.

The toolkit's `MuonSource` is a deliberately simple model: monoenergetic-with-
Gaussian-spread momentum, narrow angular cone, fixed direction, configurable
flux. Real surface-muon beams (PSI, J-PARC) reach ~10⁸ µ/s; cosmic-ray muography
operates near 1 µ cm⁻² min⁻¹. The bundled default flux of 10⁹ µ/s is an
aspirational design target representing a "theoretical on-demand high muon
source" — useful for exploring what high-flux muography could enable, without
implying any extant hardware can deliver it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np


@dataclass(frozen=True)
class MuonSource:
    """An idealised muon beam.

    Parameters:
        mean_momentum_MeV: central momentum of the beam (e.g. 30 MeV/c surface,
            300 MeV/c MIP, 1 GeV/c cosmic-equivalent).
        momentum_FWHM_MeV: Gaussian momentum spread (full width at half maximum).
        flux_per_s: total muons delivered per second through the collimator.
        direction: unit-vector beam direction (default: −z, downward).
        angular_spread_mrad: half-angle of the divergence cone (1 mrad = collimated).
        polarity: sign of the muon charge. Negative muons are required for muonic-
            atom (stopping-muon) imaging; positive muons typically used for µSR.
    """

    mean_momentum_MeV: float = 300.0
    momentum_FWHM_MeV: float = 5.0
    flux_per_s: float = 1.0e9
    direction: tuple[float, float, float] = (0.0, 0.0, -1.0)
    angular_spread_mrad: float = 1.0
    polarity: Literal["positive", "negative"] = "negative"
    seed: int | None = field(default=None, repr=False)

    def sample(self, n: int, *, rng: np.random.Generator | None = None
               ) -> tuple[np.ndarray, np.ndarray]:
        """Draw `n` muon (momentum, direction) pairs from the source.

        Returns:
            momenta: shape (n,) — sampled momenta in MeV/c.
            directions: shape (n, 3) — unit-vector directions.
        """
        if rng is None:
            rng = np.random.default_rng(self.seed)
        sigma_p = self.momentum_FWHM_MeV / 2.355
        momenta = np.maximum(
            rng.normal(self.mean_momentum_MeV, sigma_p, size=n),
            1.0,
        )
        # Sample directions inside a small cone around `self.direction`.
        nominal = np.asarray(self.direction, dtype=float)
        nominal = nominal / np.linalg.norm(nominal)
        # Pick two basis vectors orthogonal to `nominal`.
        if abs(nominal[2]) < 0.999:
            ortho = np.cross(nominal, np.array([0.0, 0.0, 1.0]))
        else:
            ortho = np.cross(nominal, np.array([1.0, 0.0, 0.0]))
        ortho = ortho / np.linalg.norm(ortho)
        ortho2 = np.cross(nominal, ortho)
        sigma = self.angular_spread_mrad * 1e-3
        theta_x = rng.normal(0.0, sigma, size=n)
        theta_y = rng.normal(0.0, sigma, size=n)
        # First-order small-angle expansion of the rotated direction.
        directions = (
            nominal[None, :]
            + theta_x[:, None] * ortho[None, :]
            + theta_y[:, None] * ortho2[None, :]
        )
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        return momenta, directions

    def expected_count(self, exposure_s: float) -> int:
        """Total muons delivered in the given exposure window (rounded)."""
        return int(round(self.flux_per_s * exposure_s))
