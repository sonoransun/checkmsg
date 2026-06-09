"""Generic spectral line/band reference primitive, shared by PL, FTIR, and CL.

These three techniques all reduce to "detect peaks/bands on a 1-D spectrum, then
match the observed positions against a reference table of named centres/bands
keyed to minerals." This module factors that out (a generalisation of
``refdata.chromophores.assign``) so each technique module stays thin.

Mössbauer is deliberately NOT built on this — its observable is a doublet/sextet
parameterised by isomer shift and quadrupole splitting, scored by residuals (see
``refdata.mossbauer_sites``).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BandSet:
    """A named set of reference line/band positions diagnostic of a centre.

    ``require_all`` semantics are decided by the caller of :func:`match_bands`,
    not stored here, so the same table can be matched strictly or loosely.
    """

    name: str
    centers: tuple[float, ...]
    tolerance: float
    typical_in: tuple[str, ...]
    notes: str = ""
    references: tuple[str, ...] = ()


def match_bands(
    positions: list[float],
    bandsets: tuple[BandSet, ...],
    *,
    require_all: bool = True,
    tolerance_scale: float = 1.0,
) -> list[tuple[float, BandSet]]:
    """Match observed ``positions`` against ``bandsets``.

    With ``require_all=True`` a bandset is reported only when *every* one of its
    centres has a matching observed position within tolerance (diagnostic
    multiplets — FTIR diamond types, CL fingerprints). With ``require_all=False``
    a bandset is reported when *any* centre matches (a single zero-phonon line is
    diagnostic on its own — PL). Each returned tuple pairs a matched observed
    position with its bandset; a bandset may appear multiple times (once per
    matched centre), mirroring ``chromophores.assign``.
    """
    out: list[tuple[float, BandSet]] = []
    for bs in bandsets:
        tol = bs.tolerance * tolerance_scale
        matched: list[float] = []
        complete = True
        for c in bs.centers:
            hits = [p for p in positions if abs(p - c) <= tol]
            if hits:
                matched.append(min(hits, key=lambda p: abs(p - c)))
            else:
                complete = False
        if require_all and not complete:
            continue
        for pos in matched:
            out.append((pos, bs))
    return out
