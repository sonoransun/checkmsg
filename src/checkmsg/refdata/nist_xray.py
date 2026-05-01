from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files

import pandas as pd


@dataclass(frozen=True)
class XrayLine:
    element: str
    Z: int
    line: str
    energy_keV: float
    relative_intensity: float


@lru_cache(maxsize=1)
def line_table() -> pd.DataFrame:
    """Bundled NIST X-ray K/L characteristic line table (keV)."""
    path = files("checkmsg.refdata.data") / "nist_xray_lines.csv"
    return pd.read_csv(str(path), comment="#")


def lines_for(element: str) -> list[XrayLine]:
    df = line_table()
    sub = df[df["element"] == element]
    return [XrayLine(**row) for row in sub.to_dict(orient="records")]


def candidates_at(energy_keV: float, tolerance_keV: float = 0.05) -> list[XrayLine]:
    """Return all NIST lines within `tolerance_keV` of the given energy."""
    df = line_table()
    mask = (df["energy_keV"] - energy_keV).abs() <= tolerance_keV
    return [XrayLine(**row) for row in df[mask].to_dict(orient="records")]
