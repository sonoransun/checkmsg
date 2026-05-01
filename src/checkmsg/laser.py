"""Excitation laser catalogue and per-wavelength corrections.

Eight laser lines are supported, spanning deep-UV through NIR:

  275, 325, 405, 457, 488, 514, 633, 830 nm.

Raman shift in cm-1 is invariant with excitation wavelength — peaks of a given
material land at the same `omega` regardless of laser. What *does* change with
the laser is:

  - absolute scattered wavelength (a 1332 cm-1 line is at ~575 nm with 532 nm
    excitation, ~707 nm with 633 nm)
  - scattering intensity (Raman cross-section scales as 1/lambda^4)
  - resonance enhancement (10x-1000x when laser tunes to an electronic absorption)
  - fluorescence background (visible lasers excite Cr3+, Mn2+, etc.; UV and NIR
    largely escape that)

This module captures each effect through small, transparent functions used by
`checkmsg.synthetic` and the example scripts. They are not full first-principles
models — they are physically-motivated approximations sufficient to demonstrate
the wavelength dependence on synthetic gem-spectroscopy data.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

# Sorted UV -> NIR. These are the eight lasers explicitly supported.
SUPPORTED_LASERS: tuple[float, ...] = (275.0, 325.0, 405.0, 457.0, 488.0, 514.0, 633.0, 830.0)

LaserRegime = Literal["deep-UV", "UV", "violet", "blue", "green", "red", "NIR"]

# h*c / (1 nm * 1 eV) -> photon energy in eV when given lambda in nm.
_HC_eV_NM = 1239.84198


def _regime_for(wavelength_nm: float) -> LaserRegime:
    if wavelength_nm < 300:
        return "deep-UV"
    if wavelength_nm < 400:
        return "UV"
    if wavelength_nm < 430:
        return "violet"
    if wavelength_nm < 480:
        return "blue"
    if wavelength_nm < 580:
        return "green"
    if wavelength_nm < 700:
        return "red"
    return "NIR"


# Heuristic per-laser fluorescence interference factor for typical gem materials.
# 0 = clean Raman, 1 = severe fluorescence drowning the Raman signal. Empirical
# trends from gemological practice: Cr3+ in beryl/corundum emits in the red,
# excited efficiently by 488/514 nm; UV (275/325) escapes most visible PL bands;
# 830 nm is below most absorption edges.
_GEM_FLUO_BY_LASER: dict[float, float] = {
    275.0: 0.10,
    325.0: 0.15,
    405.0: 0.45,
    457.0: 0.65,
    488.0: 0.85,
    514.0: 1.00,
    633.0: 0.40,
    830.0: 0.05,
}


@dataclass(frozen=True)
class LaserConfig:
    """A single excitation wavelength with derived properties."""

    wavelength_nm: float

    def __post_init__(self) -> None:
        if self.wavelength_nm <= 0:
            raise ValueError(f"laser wavelength must be positive, got {self.wavelength_nm}")

    @classmethod
    def supported(cls) -> tuple[LaserConfig, ...]:
        return tuple(cls(wl) for wl in SUPPORTED_LASERS)

    @property
    def regime(self) -> LaserRegime:
        return _regime_for(self.wavelength_nm)

    @property
    def photon_eV(self) -> float:
        return _HC_eV_NM / self.wavelength_nm

    @property
    def lambda4_scale(self) -> float:
        """Relative Raman scattering efficiency vs. a 532 nm reference (proportional to 1/lambda^4)."""
        return (532.0 / self.wavelength_nm) ** 4

    @property
    def fluorescence_factor(self) -> float:
        """Fraction in [0, 1] estimating fluorescence interference for typical gems.

        Lasers not in the bundled table fall back by linear interpolation between
        neighbouring entries.
        """
        if self.wavelength_nm in _GEM_FLUO_BY_LASER:
            return _GEM_FLUO_BY_LASER[self.wavelength_nm]
        keys = sorted(_GEM_FLUO_BY_LASER)
        if self.wavelength_nm <= keys[0]:
            return _GEM_FLUO_BY_LASER[keys[0]]
        if self.wavelength_nm >= keys[-1]:
            return _GEM_FLUO_BY_LASER[keys[-1]]
        for i in range(len(keys) - 1):
            a, b = keys[i], keys[i + 1]
            if a <= self.wavelength_nm <= b:
                t = (self.wavelength_nm - a) / (b - a)
                return (1 - t) * _GEM_FLUO_BY_LASER[a] + t * _GEM_FLUO_BY_LASER[b]
        return 0.5  # unreachable

    def shift_to_wavelength(self, raman_shift_cm: float) -> float:
        """Convert a Stokes Raman shift (cm-1) into the absolute scattered wavelength (nm)."""
        # 1/lambda_s = 1/lambda_0 - shift_cm * 1e-7
        inv_l0_per_nm = 1.0 / self.wavelength_nm
        inv_ls_per_nm = inv_l0_per_nm - raman_shift_cm * 1e-7
        if inv_ls_per_nm <= 0:
            raise ValueError(
                f"Stokes shift {raman_shift_cm} cm-1 exceeds 1/lambda_0 for {self.wavelength_nm} nm"
            )
        return 1.0 / inv_ls_per_nm

    def wavelength_to_shift(self, scattered_nm: float) -> float:
        """Inverse of `shift_to_wavelength` — a wavelength in nm to Stokes shift (cm-1)."""
        return (1.0 / self.wavelength_nm - 1.0 / scattered_nm) * 1e7

    def resonance_enhancement(
        self,
        absorption_band_nm: float,
        sigma_nm: float = 30.0,
        max_enhancement: float = 50.0,
    ) -> float:
        """Multiplicative Raman enhancement when the laser tunes near an electronic absorption.

        Modeled as a Gaussian centered on the absorption with the given sigma, scaled to
        `max_enhancement` at exact resonance and 1.0 far away.
        """
        delta = self.wavelength_nm - absorption_band_nm
        gaussian = math.exp(-(delta ** 2) / (2.0 * sigma_nm ** 2))
        return 1.0 + (max_enhancement - 1.0) * gaussian

    def above_bandgap(self, bandgap_eV: float) -> bool:
        """True if photon energy exceeds the host material's bandgap (Raman is then weak)."""
        return self.photon_eV >= bandgap_eV

    def __str__(self) -> str:
        return f"{self.wavelength_nm:g} nm ({self.regime})"


def by_wavelength(wavelength_nm: float) -> LaserConfig:
    """Construct a `LaserConfig`, validating against the supported catalogue.

    Raises ValueError for wavelengths not in `SUPPORTED_LASERS`. Use
    `LaserConfig(wavelength_nm)` directly to bypass the check (e.g. for 532 nm references).
    """
    if wavelength_nm not in SUPPORTED_LASERS:
        raise ValueError(
            f"{wavelength_nm} nm is not a supported laser; choose one of {SUPPORTED_LASERS}"
        )
    return LaserConfig(wavelength_nm)
