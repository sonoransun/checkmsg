from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from checkmsg.peaks import Peak
from checkmsg.spectrum import Spectrum


@dataclass(frozen=True)
class MatchResult:
    name: str
    score: float
    detail: dict


def cosine(a: Spectrum, b: Spectrum, lo: float | None = None, hi: float | None = None) -> float:
    """Cosine similarity after interpolating both spectra onto a common grid.

    Both spectra must share a technique/units. Score in [-1, 1]; 1 = identical shape.
    """
    if a.technique != b.technique:
        raise ValueError(f"technique mismatch: {a.technique} vs {b.technique}")
    x_lo = max(a.axis[0], b.axis[0]) if lo is None else lo
    x_hi = min(a.axis[-1], b.axis[-1]) if hi is None else hi
    if x_hi <= x_lo:
        return 0.0
    n = max(256, min(len(a), len(b)))
    grid = np.linspace(x_lo, x_hi, n)
    ya = np.interp(grid, a.axis, a.intensity)
    yb = np.interp(grid, b.axis, b.intensity)
    ya = ya - ya.min()
    yb = yb - yb.min()
    na = np.linalg.norm(ya)
    nb = np.linalg.norm(yb)
    if na == 0 or nb == 0:
        return 0.0
    return float(np.dot(ya, yb) / (na * nb))


def peak_list_match(peaks_a: list[Peak], peaks_b: list[Peak], tolerance: float) -> dict:
    """Symmetric peak-list match: fraction of peaks within `tolerance` of a counterpart.

    Returns a dict with `score` (Jaccard-like, weighted by relative height) and the matched pairs.
    """
    if not peaks_a or not peaks_b:
        return {"score": 0.0, "matches": []}
    pa = sorted(peaks_a, key=lambda p: p.position)
    pb = sorted(peaks_b, key=lambda p: p.position)
    used_b: set[int] = set()
    matches: list[tuple[Peak, Peak]] = []
    for p in pa:
        best_j, best_d = None, tolerance + 1.0
        for j, q in enumerate(pb):
            if j in used_b:
                continue
            d = abs(p.position - q.position)
            if d < best_d:
                best_d = d
                best_j = j
        if best_j is not None and best_d <= tolerance:
            used_b.add(best_j)
            matches.append((p, pb[best_j]))
    union = len(pa) + len(pb) - len(matches)
    score = len(matches) / union if union else 0.0
    return {"score": float(score), "matches": matches}


def rank(unknown: Spectrum, candidates: dict[str, Spectrum], top: int = 5) -> list[MatchResult]:
    """Rank candidates by cosine similarity to the unknown."""
    out = []
    for name, c in candidates.items():
        s = cosine(unknown, c)
        out.append(MatchResult(name=name, score=s, detail={"method": "cosine"}))
    out.sort(key=lambda r: r.score, reverse=True)
    return out[:top]
