from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

Technique = Literal["raman", "xrf", "libs", "uvvis", "epr", "laicpms", "muon-xray"]


@dataclass
class Spectrum:
    axis: np.ndarray
    intensity: np.ndarray
    technique: Technique
    units: str
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.axis = np.asarray(self.axis, dtype=float)
        self.intensity = np.asarray(self.intensity, dtype=float)
        if self.axis.shape != self.intensity.shape:
            raise ValueError(f"axis/intensity shape mismatch: {self.axis.shape} vs {self.intensity.shape}")
        if self.axis.ndim != 1:
            raise ValueError("Spectrum requires 1-D axis/intensity")
        order = np.argsort(self.axis)
        if not np.all(order == np.arange(self.axis.size)):
            self.axis = self.axis[order]
            self.intensity = self.intensity[order]

    def slice(self, low: float, high: float) -> Spectrum:
        mask = (self.axis >= low) & (self.axis <= high)
        return Spectrum(self.axis[mask], self.intensity[mask], self.technique, self.units, dict(self.metadata))

    def normalize(self, mode: str = "max") -> Spectrum:
        y = self.intensity
        if mode == "max":
            scale = float(np.max(np.abs(y))) or 1.0
            y = y / scale
        elif mode == "area":
            scale = float(np.trapezoid(np.abs(y), self.axis)) or 1.0
            y = y / scale
        elif mode == "minmax":
            lo, hi = float(np.min(y)), float(np.max(y))
            y = (y - lo) / ((hi - lo) or 1.0)
        else:
            raise ValueError(f"unknown normalize mode: {mode}")
        return Spectrum(self.axis.copy(), y, self.technique, self.units, dict(self.metadata))

    def with_intensity(self, new_intensity: np.ndarray) -> Spectrum:
        return Spectrum(self.axis.copy(), np.asarray(new_intensity), self.technique, self.units, dict(self.metadata))

    def __len__(self) -> int:
        return int(self.axis.size)

    def __repr__(self) -> str:
        return (
            f"Spectrum(technique={self.technique!r}, units={self.units!r}, "
            f"n={len(self)}, range=[{self.axis[0]:.3g}, {self.axis[-1]:.3g}])"
        )
