"""Cathodoluminescence (CL) analysis — activator-band emission matching.

CL reveals luminescent activator centres (diamond band-A, Cr³⁺ in corundum,
Mn²⁺ in carbonates, REE³⁺ in zircon). Natural-vs-synthetic growth-zoning
discrimination genuinely needs *imaging*; this 1-D module is a pedagogical
approximation that matches emission bands only and surfaces zoning context as
descriptive notes.
"""

from __future__ import annotations

from dataclasses import dataclass

from checkmsg import preprocess
from checkmsg.peaks import Peak, detect, prominent
from checkmsg.refdata.cl_bands import CL_BANDS
from checkmsg.refdata.line_bands import BandSet, match_bands
from checkmsg.spectrum import Spectrum


@dataclass
class ClResult:
    bands: list[Peak]
    assignments: list[tuple[float, BandSet]]
    cleaned: Spectrum

    def emitters(self) -> list[BandSet]:
        seen: list[BandSet] = []
        for _pos, bs in self.assignments:
            if bs not in seen:
                seen.append(bs)
        return seen

    @property
    def best(self) -> BandSet | None:
        return self.emitters()[0] if self.assignments else None

    def headline(self) -> str:
        names = [b.name for b in self.emitters()]
        return f"CL bands: {', '.join(names) if names else 'none identified'}"


def preprocess_cl(spectrum: Spectrum) -> Spectrum:
    if spectrum.technique != "cl":
        raise ValueError(f"expected cl spectrum, got {spectrum.technique}")
    y = preprocess.savgol(spectrum.intensity, window=15, order=3)
    baseline = preprocess.als_baseline(y, lam=1e5, p=0.01)
    return spectrum.with_intensity(y - baseline)


def analyze(spectrum: Spectrum, min_snr: float = 5.0) -> ClResult:
    cleaned = preprocess_cl(spectrum)
    bands = prominent(detect(cleaned, min_snr=min_snr))
    assignments = match_bands([b.position for b in bands], CL_BANDS, require_all=False)
    return ClResult(bands=bands, assignments=assignments, cleaned=cleaned)
