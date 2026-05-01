from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from checkmsg import preprocess
from checkmsg.peaks import Peak, detect
from checkmsg.refdata.nist_xray import XrayLine, candidates_at, lines_for
from checkmsg.spectrum import Spectrum


@dataclass
class ElementAssignment:
    element: str
    Z: int
    matched_lines: list[XrayLine] = field(default_factory=list)
    summed_height: float = 0.0
    confidence: float = 0.0  # heuristic


@dataclass
class XrfResult:
    peaks: list[Peak]
    elements: list[ElementAssignment]
    cleaned: Spectrum


def preprocess_xrf(spectrum: Spectrum) -> Spectrum:
    if spectrum.technique != "xrf":
        raise ValueError(f"expected xrf spectrum, got {spectrum.technique}")
    baseline = preprocess.snip_baseline(spectrum.intensity, iterations=40)
    y_clean = np.clip(spectrum.intensity - baseline, 0, None)
    return spectrum.with_intensity(y_clean)


def identify_elements(
    spectrum: Spectrum,
    tolerance_keV: float = 0.06,
    min_snr: float = 6.0,
) -> XrfResult:
    cleaned = preprocess_xrf(spectrum)
    peaks = detect(cleaned, min_snr=min_snr)
    by_element: dict[str, ElementAssignment] = {}
    for peak in peaks:
        for line in candidates_at(peak.position, tolerance_keV):
            ea = by_element.setdefault(line.element, ElementAssignment(element=line.element, Z=line.Z))
            ea.matched_lines.append(line)
            ea.summed_height += peak.height * line.relative_intensity
    # Confidence: 1.0 if Ka+Kb both found at right ratio, else partial
    for element, ea in by_element.items():
        expected = lines_for(element)
        n_expected = len([line for line in expected if line.line in {"Ka", "Kb", "La", "Lb"}])
        n_found = len({line.line for line in ea.matched_lines})
        ea.confidence = min(1.0, n_found / max(1, n_expected))
    elements = sorted(by_element.values(), key=lambda e: e.summed_height, reverse=True)
    return XrfResult(peaks=peaks, elements=elements, cleaned=cleaned)


def relative_quant(result: XrfResult) -> dict[str, float]:
    """Crude relative quantification: peak-area fractions normalized to 1.0.

    Real XRF needs fundamental-parameters correction (matrix effects, Z-dependent yields).
    This is a relative indicator, not concentration.
    """
    total = sum(e.summed_height for e in result.elements) or 1.0
    return {e.element: e.summed_height / total for e in result.elements}
