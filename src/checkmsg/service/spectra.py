"""Build ``Spectrum`` objects from API inputs (inline arrays or uploaded CSV)."""

from __future__ import annotations

import os
import tempfile

import numpy as np

from checkmsg import io as io_mod
from checkmsg.spectrum import Spectrum

# Technique -> axis units (kept local so the service package is self-contained).
UNITS: dict[str, str] = {
    "raman": "cm-1", "xrf": "keV", "libs": "nm", "uvvis": "nm",
    "epr": "mT", "laicpms": "m/z", "squid-mh": "mT", "squid-chi": "K",
    "pl": "nm", "ftir": "cm-1", "mossbauer": "mm/s", "cl": "nm",
}


def _metadata_for(technique: str, *, frequency_GHz=None, temperature_K=None,
                  bias_mT=None, isotope_keys=None) -> dict:
    meta: dict = {}
    if technique == "epr" and frequency_GHz is not None:
        meta["frequency_GHz"] = float(frequency_GHz)
    elif technique == "squid-mh":
        meta["squid_mode"] = "dc-mh"
        meta["temperature_K"] = float(temperature_K if temperature_K is not None else 295.0)
    elif technique == "squid-chi":
        meta["squid_mode"] = "rf-chi-T"
        if bias_mT is not None:
            meta["applied_field_mT"] = float(bias_mT)
    elif technique == "laicpms" and isotope_keys:
        meta["isotope_keys"] = list(isotope_keys)
    return meta


def spectrum_from_arrays(technique: str, axis, intensity, **params) -> Spectrum:
    """Build a Spectrum from inline numeric arrays + per-technique params."""
    meta = _metadata_for(technique, **params)
    return Spectrum(np.asarray(axis, dtype=float), np.asarray(intensity, dtype=float),
                    technique, UNITS[technique], meta)


def spectrum_from_upload(technique: str, content: bytes, **params) -> Spectrum:
    """Build a Spectrum from an uploaded 2-column CSV.

    ``io.read_csv`` is path-based, so the bytes are written to a temp file that
    is removed afterwards.
    """
    meta = _metadata_for(technique, **params)
    tmp = tempfile.NamedTemporaryFile(suffix=".csv", delete=False)
    try:
        tmp.write(content)
        tmp.flush()
        tmp.close()
        return io_mod.read_csv(tmp.name, technique=technique, units=UNITS[technique], metadata=meta)
    finally:
        os.unlink(tmp.name)
