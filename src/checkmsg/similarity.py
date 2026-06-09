"""Spectral-embedding similarity search — nearest catalog references.

A standalone retrieval capability: it embeds a spectrum into a fixed-length unit
vector (resample to a per-technique canonical grid + baseline-subtract +
L2-normalise — the same transform ``match.cosine`` applies) and ranks the
catalog by cosine similarity. It is deliberately NOT wired into ``diagnose`` —
the additive verdict stays transparent and auditable; this is a parallel "what
does it look most like?" tool. Offline and pedagogical: the index is built
lazily from the catalog synthesizers and cached.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from checkmsg import minerals
from checkmsg.spectrum import Spectrum

# Canonical resampling grid per technique (mirrors the synthesizer axes so a
# query and the catalog references land on the same grid).
_GRIDS: dict[str, np.ndarray] = {
    "raman": np.linspace(100.0, 1500.0, 1401),
    "uvvis": np.linspace(380.0, 800.0, 421),
    "xrf": np.linspace(0.5, 15.0, 2048),
    "libs": np.linspace(200.0, 800.0, 3000),
    "pl": np.linspace(400.0, 800.0, 601),
    "cl": np.linspace(350.0, 750.0, 401),
    "ftir": np.linspace(400.0, 4000.0, 1801),
}


@dataclass(frozen=True)
class SimilarityHit:
    name: str
    score: float
    technique: str


def embed(spectrum: Spectrum) -> np.ndarray:
    """Resample to the technique's canonical grid, baseline-subtract, L2-normalise."""
    grid = _GRIDS.get(spectrum.technique)
    if grid is None:
        raise ValueError(f"similarity not supported for technique {spectrum.technique!r}")
    y = np.interp(grid, spectrum.axis, spectrum.intensity, left=0.0, right=0.0)
    y = y - float(np.min(y))
    norm = float(np.linalg.norm(y))
    return y / norm if norm > 0 else y


def _synthesize(profile: minerals.MineralProfile, technique: str) -> Spectrum | None:
    if technique == "raman" and profile.raman_peaks_cm:
        return minerals.synthesize_raman(profile, noise=0.005, seed=0)
    if technique == "uvvis" and profile.uvvis_bands_nm:
        return minerals.synthesize_uvvis(profile, seed=0)
    if technique == "xrf" and profile.xrf_signature:
        return minerals.synthesize_xrf(profile, seed=0)
    if technique == "libs" and profile.libs_signature:
        return minerals.synthesize_libs(profile, seed=0)
    if technique == "pl" and profile.pl_centers:
        return minerals.synthesize_pl(profile, seed=0)
    if technique == "cl" and profile.cl_bands:
        return minerals.synthesize_cl(profile, seed=0)
    if technique == "ftir" and (profile.ftir_bands or profile.diamond_type):
        return minerals.synthesize_ftir(profile, seed=0)
    return None


@lru_cache(maxsize=16)
def _index(technique: str) -> tuple[tuple[str, ...], np.ndarray]:
    """Build (and cache) the catalog embedding matrix for one technique."""
    names: list[str] = []
    rows: list[np.ndarray] = []
    for name, profile in minerals.CATALOG.items():
        spec = _synthesize(profile, technique)
        if spec is None:
            continue
        names.append(name)
        rows.append(embed(spec))
    matrix = np.array(rows) if rows else np.zeros((0, len(_GRIDS.get(technique, [0]))))
    return tuple(names), matrix


def nearest(spectrum: Spectrum, k: int = 5) -> list[SimilarityHit]:
    """Return the ``k`` catalog entries whose reference spectrum is most similar."""
    names, matrix = _index(spectrum.technique)
    if matrix.shape[0] == 0:
        return []
    q = embed(spectrum)
    scores = matrix @ q  # both unit vectors → cosine similarity
    order = np.argsort(scores)[::-1][:k]
    return [SimilarityHit(names[i], float(scores[i]), spectrum.technique) for i in order]
