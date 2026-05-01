from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files

import pandas as pd


@dataclass(frozen=True)
class AtomicLine:
    element: str
    wavelength_nm: float
    relative_intensity: float
    species: str


@lru_cache(maxsize=1)
def line_table() -> pd.DataFrame:
    path = files("checkmsg.refdata.data") / "nist_asd_subset.csv"
    return pd.read_csv(str(path), comment="#")


def lines_for(element: str) -> list[AtomicLine]:
    df = line_table()
    sub = df[df["element"] == element]
    return [AtomicLine(**row) for row in sub.to_dict(orient="records")]


def candidates_at(wavelength_nm: float, tolerance_nm: float = 0.3) -> list[AtomicLine]:
    df = line_table()
    mask = (df["wavelength_nm"] - wavelength_nm).abs() <= tolerance_nm
    return [AtomicLine(**row) for row in df[mask].to_dict(orient="records")]


def all_elements() -> list[str]:
    return sorted(line_table()["element"].unique().tolist())
