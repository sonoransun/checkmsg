from __future__ import annotations

from dataclasses import dataclass

from checkmsg import preprocess
from checkmsg.peaks import Peak, detect
from checkmsg.refdata.chromophores import Chromophore, assign
from checkmsg.spectrum import Spectrum


@dataclass
class UvVisResult:
    bands: list[Peak]
    assignments: list[tuple[float, Chromophore]]
    cleaned: Spectrum
    polarity: str  # "absorbance" or "transmittance"

    def chromophores(self) -> list[Chromophore]:
        seen: list[Chromophore] = []
        for _pos, c in self.assignments:
            if c not in seen:
                seen.append(c)
        return seen


def preprocess_uvvis(spectrum: Spectrum, polarity: str = "absorbance") -> Spectrum:
    if spectrum.technique != "uvvis":
        raise ValueError(f"expected uvvis spectrum, got {spectrum.technique}")
    y = preprocess.savgol(spectrum.intensity, window=15, order=3)
    if polarity == "transmittance":
        y = -y
    baseline = preprocess.als_baseline(y, lam=1e6, p=0.001)
    return spectrum.with_intensity(y - baseline)


def assign_bands(
    spectrum: Spectrum,
    polarity: str = "absorbance",
    min_snr: float = 4.0,
) -> UvVisResult:
    cleaned = preprocess_uvvis(spectrum, polarity=polarity)
    bands = detect(cleaned, min_snr=min_snr)
    assignments = assign([b.position for b in bands])
    return UvVisResult(bands=bands, assignments=assignments, cleaned=cleaned, polarity=polarity)


def is_transparent(spectrum: Spectrum, min_snr: float = 4.0) -> bool:
    """Heuristic: featureless spectrum across visible range = transparent (e.g. diamond)."""
    cleaned = preprocess_uvvis(spectrum)
    peaks = detect(cleaned, min_snr=min_snr)
    return not peaks
