from __future__ import annotations

from pathlib import Path

import numpy as np

from checkmsg.spectrum import Spectrum, Technique


def read_csv(path: str | Path, technique: Technique, units: str, metadata: dict | None = None) -> Spectrum:
    """Read a two-column ASCII/CSV file (axis, intensity). Comment lines starting with '#' are ignored."""
    arr = np.loadtxt(str(path), delimiter=_sniff_delim(path), comments="#")
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError(f"{path}: expected 2-column data, got shape {arr.shape}")
    return Spectrum(arr[:, 0], arr[:, 1], technique, units, dict(metadata or {}))


def write_csv(spectrum: Spectrum, path: str | Path) -> None:
    arr = np.column_stack([spectrum.axis, spectrum.intensity])
    np.savetxt(str(path), arr, delimiter=",", header=f"{spectrum.technique},{spectrum.units}")


def _sniff_delim(path: str | Path) -> str:
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            if "," in line:
                return ","
            if "\t" in line:
                return "\t"
            return None  # whitespace
    return None
