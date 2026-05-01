"""RRUFF Project — Raman reference spectra for minerals.

The RRUFF Project (https://rruff.info) provides freely-licensed Raman spectra of
mineral type specimens. We fetch raw 532 nm processed spectra on first use and
cache them locally; subsequent calls (and tests) hit the cache only.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import requests

from checkmsg.refdata.cache import cache_root, is_offline
from checkmsg.spectrum import Spectrum

RRUFF_BASE = "https://rruff.info/repository/sample/by_minerals"


@dataclass(frozen=True)
class RRUFFEntry:
    mineral: str
    rruff_id: str
    laser_nm: int
    orientation: str | None = None


# Curated handful of "ideal" reference spectra — enough for the example scenarios.
# Each is the RRUFF processed file URL (532 nm laser, 'oriented' or 'unoriented').
KNOWN: dict[str, str] = {
    "diamond": "https://rruff.info/repository/sample_child_record_raman/by_minerals/Diamond__R050204__Raman__532__unoriented__Raman_Data_Processed__19129.txt",
    "moissanite": "https://rruff.info/repository/sample_child_record_raman/by_minerals/Moissanite__R061083__Raman__532__unoriented__Raman_Data_Processed__9837.txt",
    "corundum": "https://rruff.info/repository/sample_child_record_raman/by_minerals/Corundum__R040112__Raman__532__unoriented__Raman_Data_Processed__24876.txt",
    "ruby": "https://rruff.info/repository/sample_child_record_raman/by_minerals/Corundum__R060020__Raman__532__unoriented__Raman_Data_Processed__9395.txt",
    "beryl": "https://rruff.info/repository/sample_child_record_raman/by_minerals/Beryl__R040146__Raman__532__unoriented__Raman_Data_Processed__9056.txt",
    "quartz": "https://rruff.info/repository/sample_child_record_raman/by_minerals/Quartz__R040031__Raman__532__unoriented__Raman_Data_Processed__9442.txt",
    "zircon": "https://rruff.info/repository/sample_child_record_raman/by_minerals/Zircon__R050050__Raman__532__unoriented__Raman_Data_Processed__9576.txt",
}


def _cache_path(name: str) -> Path:
    d = cache_root() / "rruff"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{name.lower()}.txt"


def _meta_path() -> Path:
    return cache_root() / "rruff" / "_meta.json"


def fetch(name: str, timeout: float = 30.0) -> Spectrum:
    """Fetch (or load from cache) a RRUFF reference Raman spectrum by mineral name."""
    name_l = name.lower()
    if name_l not in KNOWN:
        raise KeyError(f"no known RRUFF entry for {name!r}; known: {sorted(KNOWN)}")
    p = _cache_path(name_l)
    if not p.exists():
        if is_offline():
            raise RuntimeError(
                f"RRUFF cache miss for {name!r} and CHECKMSG_OFFLINE=1; "
                f"warm the cache with one online run first.")
        url = KNOWN[name_l]
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        p.write_text(r.text)
        meta = {}
        if _meta_path().exists():
            meta = json.loads(_meta_path().read_text())
        meta[name_l] = {"url": url}
        _meta_path().write_text(json.dumps(meta, indent=2))
    return _parse(p, name_l)


def cached_minerals() -> list[str]:
    d = cache_root() / "rruff"
    if not d.exists():
        return []
    return sorted(p.stem for p in d.glob("*.txt"))


def load_cached(name: str) -> Spectrum:
    """Load a cached RRUFF spectrum without attempting network."""
    p = _cache_path(name.lower())
    if not p.exists():
        raise FileNotFoundError(f"no cached RRUFF spectrum for {name!r}")
    return _parse(p, name.lower())


def _parse(path: Path, name: str) -> Spectrum:
    """RRUFF processed files are CSV with header lines starting with ##; data is x, y."""
    text = path.read_text()
    rows: list[tuple[float, float]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("##") or line.startswith("#"):
            continue
        parts = [p for p in line.replace(";", ",").split(",") if p.strip()]
        if len(parts) < 2:
            parts = line.split()
        try:
            x = float(parts[0])
            y = float(parts[1])
        except (ValueError, IndexError):
            continue
        rows.append((x, y))
    if not rows:
        raise ValueError(f"no parseable rows in {path}")
    arr = np.asarray(rows, dtype=float)
    return Spectrum(arr[:, 0], arr[:, 1], "raman", "cm-1", {"source": "RRUFF", "mineral": name})


def install_synthetic_fallback(mineral: str, spectrum: Spectrum) -> None:
    """Write a fabricated reference into the cache (used by tests / offline demos).

    The library treats this exactly like a real RRUFF entry on subsequent calls.
    """
    p = _cache_path(mineral.lower())
    arr = np.column_stack([spectrum.axis, spectrum.intensity])
    header = f"##NAME={mineral}\n##SOURCE=synthetic\n"
    body = "\n".join(f"{x:.4f},{y:.6f}" for x, y in arr)
    p.write_text(header + body)
