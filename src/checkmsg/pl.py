"""Photoluminescence (PL) analysis — sharp defect zero-phonon-line matching.

PL is the modern frontier for separating natural, HPHT-synthetic, and CVD-synthetic
diamond: NV⁻/NV⁰ and especially the Si-V centre (~737 nm) are growth markers.
This module is intentionally pedagogical — it matches detected emission lines
against tabulated room/low-T ZPL positions; it does not model the phonon
sidebands or temperature/strain shifts of a real spectrometer.
"""

from __future__ import annotations

from dataclasses import dataclass

from checkmsg import preprocess
from checkmsg.peaks import Peak, detect, prominent
from checkmsg.refdata.line_bands import BandSet, match_bands
from checkmsg.refdata.pl_lines import PL_CENTERS, SYNTHETIC_MARKER_NM
from checkmsg.spectrum import Spectrum


@dataclass
class PlResult:
    lines: list[Peak]
    assignments: list[tuple[float, BandSet]]
    cleaned: Spectrum

    def defects(self) -> list[BandSet]:
        seen: list[BandSet] = []
        for _pos, bs in self.assignments:
            if bs not in seen:
                seen.append(bs)
        return seen

    @property
    def best(self) -> BandSet | None:
        return self.defects()[0] if self.assignments else None

    def has_synthetic_marker(self) -> bool:
        return any(any(abs(pos - m) <= 4.0 for m in SYNTHETIC_MARKER_NM)
                   for pos, _bs in self.assignments)

    def headline(self) -> str:
        names = [b.name for b in self.defects()]
        tag = " [synthetic marker]" if self.has_synthetic_marker() else ""
        return f"PL centres: {', '.join(names) if names else 'none identified'}{tag}"


def preprocess_pl(spectrum: Spectrum) -> Spectrum:
    if spectrum.technique != "pl":
        raise ValueError(f"expected pl spectrum, got {spectrum.technique}")
    y = preprocess.savgol(spectrum.intensity, window=11, order=3)
    baseline = preprocess.als_baseline(y, lam=1e5, p=0.01)
    return spectrum.with_intensity(y - baseline)


def analyze(spectrum: Spectrum, min_snr: float = 6.0) -> PlResult:
    cleaned = preprocess_pl(spectrum)
    lines = prominent(detect(cleaned, min_snr=min_snr))
    assignments = match_bands([p.position for p in lines], PL_CENTERS, require_all=False)
    return PlResult(lines=lines, assignments=assignments, cleaned=cleaned)
