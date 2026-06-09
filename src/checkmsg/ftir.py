"""FTIR (infrared) analysis — diamond type + functional-group band matching.

FTIR classifies diamond by nitrogen aggregation (Type Ia/Ib/IIa/IIb) and detects
structural water in beryl and polymer impregnation in jade. Bands are broad
absorption features; the band centres here are textbook values — a pedagogical
approximation, not the quantitative one-phonon deconvolution a real lab uses for
nitrogen concentration.
"""

from __future__ import annotations

from dataclasses import dataclass

from checkmsg import preprocess
from checkmsg.peaks import Peak, detect, prominent
from checkmsg.refdata.ftir_bands import DIAMOND_TYPE_BANDS, FTIR_BANDS
from checkmsg.refdata.line_bands import BandSet, match_bands
from checkmsg.spectrum import Spectrum


@dataclass
class FtirResult:
    bands: list[Peak]
    assignments: list[tuple[float, BandSet]]
    cleaned: Spectrum
    polarity: str
    diamond_type: str  # "Ia" | "Ib" | "IIb" | "" (IIa shows no N/B bands)

    def band_names(self) -> list[str]:
        seen: list[str] = []
        for _pos, bs in self.assignments:
            if bs.name not in seen:
                seen.append(bs.name)
        return seen

    def has_polymer(self) -> bool:
        return any("polymer" in bs.name for _pos, bs in self.assignments)

    def headline(self) -> str:
        parts = []
        if self.diamond_type:
            parts.append(f"diamond Type {self.diamond_type}")
        parts.extend(n for n in self.band_names() if not n.startswith("diamond Type"))
        return f"FTIR: {', '.join(parts) if parts else 'no diagnostic bands'}"


def preprocess_ftir(spectrum: Spectrum, polarity: str = "absorbance") -> Spectrum:
    if spectrum.technique != "ftir":
        raise ValueError(f"expected ftir spectrum, got {spectrum.technique}")
    y = preprocess.savgol(spectrum.intensity, window=15, order=3)
    if polarity == "transmittance":
        y = -y
    baseline = preprocess.als_baseline(y, lam=1e6, p=0.001)
    return spectrum.with_intensity(y - baseline)


def _detect_diamond_type(positions: list[float], tol: float = 12.0) -> str:
    """Return the diamond IR type whose full band multiplet is present, else ''."""
    for tcode, centers in DIAMOND_TYPE_BANDS.items():
        if all(any(abs(p - c) <= tol for p in positions) for c in centers):
            return tcode
    return ""


def analyze(spectrum: Spectrum, polarity: str = "absorbance", min_snr: float = 4.0) -> FtirResult:
    cleaned = preprocess_ftir(spectrum, polarity=polarity)
    bands = prominent(detect(cleaned, min_snr=min_snr))
    positions = [b.position for b in bands]
    assignments = match_bands(positions, FTIR_BANDS, require_all=True)
    return FtirResult(bands=bands, assignments=assignments, cleaned=cleaned,
                      polarity=polarity, diamond_type=_detect_diamond_type(positions))
