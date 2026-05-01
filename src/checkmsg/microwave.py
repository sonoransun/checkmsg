"""Microwave band catalogue for CW EPR / ESR.

Seven canonical bands span 1 GHz through 100 GHz; each defines a typical
resonance-field range for g ~ 2 paramagnetic centers:

    L (1.1 GHz, ~39 mT)   S (3.0 GHz, ~107 mT)   C (5.0 GHz, ~178 mT)
    X (9.5 GHz, ~339 mT)  K (24 GHz, ~856 mT)    Q (35 GHz, ~1248 mT)
    W (94 GHz, ~3355 mT)

X-band is the workhorse in gemological practice. Q and W bands resolve smaller
g-anisotropies; L/S/C see use for biological samples and large-volume cavities.

The conversion between resonance field B (mT) and frequency nu (GHz) is

    h * nu = g * mu_B * B    =>    B [mT] = (h / mu_B) * nu / g
                                          ~= 71.4477 * nu_GHz / g

This module exposes a small `MicrowaveBand` dataclass plus a registry of the
seven canonical bands. The frequency on a real spectrometer can be anything
inside its tunable range, so the constants are the *typical* center frequency
and EPR analyses still take the actual operating frequency as input.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

# Bohr magneton / Planck constant in mT/GHz.
# h * 1 GHz = mu_B * (h/mu_B) * 1 GHz / g; (h/mu_B) ~= 71.4477 mT/GHz at g=1.
_MT_PER_GHZ_AT_G1 = 71.44773  # mT * (g=1) per GHz; B = this * freq / g

BandRegime = Literal["LF", "X-class", "high-field"]


@dataclass(frozen=True)
class MicrowaveBand:
    """A named CW EPR microwave band."""

    name: str
    frequency_GHz: float

    def __post_init__(self) -> None:
        if self.frequency_GHz <= 0:
            raise ValueError(f"frequency must be > 0 GHz, got {self.frequency_GHz}")

    @property
    def regime(self) -> BandRegime:
        if self.frequency_GHz < 8:
            return "LF"
        if self.frequency_GHz < 50:
            return "X-class"
        return "high-field"

    def resonance_field_mT(self, g: float = 2.0) -> float:
        """Resonance field for a given g-factor at this band's frequency."""
        if g <= 0:
            raise ValueError(f"g must be > 0, got {g}")
        return _MT_PER_GHZ_AT_G1 * self.frequency_GHz / g

    def g_at(self, field_mT: float) -> float:
        """g-factor implied by resonance at the given field (assuming this band's frequency)."""
        if field_mT <= 0:
            raise ValueError(f"field must be > 0 mT, got {field_mT}")
        return _MT_PER_GHZ_AT_G1 * self.frequency_GHz / field_mT

    def __str__(self) -> str:
        return f"{self.name}-band ({self.frequency_GHz:g} GHz, {self.regime})"


# Canonical band frequencies — these are the values quoted by every commercial
# EPR vendor. Real instruments vary by ~10% within each band's tunable range.
SUPPORTED_BANDS: tuple[MicrowaveBand, ...] = (
    MicrowaveBand("L", 1.1),
    MicrowaveBand("S", 3.0),
    MicrowaveBand("C", 5.0),
    MicrowaveBand("X", 9.5),
    MicrowaveBand("K", 24.0),
    MicrowaveBand("Q", 35.0),
    MicrowaveBand("W", 94.0),
)

_BY_NAME: dict[str, MicrowaveBand] = {b.name: b for b in SUPPORTED_BANDS}


def by_name(name: str) -> MicrowaveBand:
    """Return the canonical band with the given letter name (e.g. 'X', 'Q')."""
    key = name.strip().upper()
    if key not in _BY_NAME:
        raise KeyError(f"unknown band {name!r}; known: {sorted(_BY_NAME)}")
    return _BY_NAME[key]


def resonance_field_mT(frequency_GHz: float, g: float = 2.0) -> float:
    """Standalone helper: resonance field at any frequency (not restricted to the catalogue)."""
    if frequency_GHz <= 0 or g <= 0:
        raise ValueError("frequency and g must be > 0")
    return _MT_PER_GHZ_AT_G1 * frequency_GHz / g


def g_factor(field_mT: float, frequency_GHz: float) -> float:
    """Standalone helper: g-factor implied by a resonance at (field, frequency)."""
    if field_mT <= 0 or frequency_GHz <= 0:
        raise ValueError("field and frequency must be > 0")
    return _MT_PER_GHZ_AT_G1 * frequency_GHz / field_mT
