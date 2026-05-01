from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from checkmsg import preprocess
from checkmsg import temperature as temp_mod
from checkmsg.match import cosine, peak_list_match
from checkmsg.peaks import Peak, detect, fit_voigt
from checkmsg.refdata import rruff
from checkmsg.spectrum import Spectrum


@dataclass(frozen=True)
class RamanCandidate:
    mineral: str
    cosine: float
    peak_score: float
    matched_peaks: int

    @property
    def combined(self) -> float:
        return 0.6 * self.cosine + 0.4 * self.peak_score


@dataclass
class RamanResult:
    peaks: list[Peak]
    candidates: list[RamanCandidate]
    cleaned: Spectrum

    @property
    def best(self) -> RamanCandidate | None:
        return self.candidates[0] if self.candidates else None


def preprocess_raman(spectrum: Spectrum) -> Spectrum:
    if spectrum.technique != "raman":
        raise ValueError(f"expected raman spectrum, got {spectrum.technique}")
    y = preprocess.savgol(spectrum.intensity, window=11, order=3)
    baseline = preprocess.als_baseline(y, lam=1e5, p=0.001)
    y_clean = y - baseline
    y_clean = np.clip(y_clean, 0, None)
    y_norm = preprocess.min_max_normalize(y_clean)
    return spectrum.with_intensity(y_norm)


def analyze(
    spectrum: Spectrum,
    candidates: list[str] | None = None,
    peak_tolerance_cm: float = 8.0,
    top: int = 5,
) -> RamanResult:
    """Baseline + peak-pick a Raman spectrum and rank against RRUFF references."""
    cleaned = preprocess_raman(spectrum)
    peaks = detect(cleaned, min_snr=8.0)
    if candidates is None:
        # default: only what's already cached — no surprise network fetches.
        candidates = rruff.cached_minerals() or list(rruff.KNOWN.keys())
    out: list[RamanCandidate] = []
    for name in candidates:
        try:
            ref = rruff.load_cached(name)
        except FileNotFoundError:
            try:
                ref = rruff.fetch(name)
            except Exception:
                continue
        ref_clean = preprocess_raman(ref)
        cos = cosine(cleaned, ref_clean)
        ref_peaks = detect(ref_clean, min_snr=8.0)
        pm = peak_list_match(peaks, ref_peaks, tolerance=peak_tolerance_cm)
        out.append(RamanCandidate(
            mineral=name,
            cosine=cos,
            peak_score=pm["score"],
            matched_peaks=len(pm["matches"]),
        ))
    out.sort(key=lambda c: c.combined, reverse=True)
    return RamanResult(peaks=peaks, candidates=out[:top], cleaned=cleaned)


def is_amorphous(spectrum: Spectrum, broad_fwhm_threshold: float = 60.0) -> bool:
    """Heuristic: amorphous materials are dominated by a single broad envelope.

    True if the strongest peak has FWHM > broad_fwhm_threshold cm-1 AND no narrow peak
    rises within an order of magnitude of it (rules out crystalline + broad fluorescence).
    """
    cleaned = preprocess_raman(spectrum)
    peaks = detect(cleaned, min_snr=8.0)
    if not peaks:
        return False
    biggest = max(peaks, key=lambda p: p.height)
    if biggest.width < broad_fwhm_threshold:
        return False
    threshold = biggest.height * 0.50
    narrow_strong = [p for p in peaks if p.width < 15.0 and p.height > threshold]
    return len(narrow_strong) == 0


def infer_temperature(
    spectrum: Spectrum,
    mode_cm: float,
    window_cm: float = 30.0,
) -> float:
    """Estimate sample temperature from a Stokes/anti-Stokes pair.

    The spectrum's axis must span both `+mode_cm` (Stokes) and `-mode_cm`
    (anti-Stokes). Heights are taken from a Voigt fit around each side; the
    ratio is inverted via `temperature.infer_temperature`.

    Useful for verifying that an LN2 / room-temperature measurement actually
    landed at the expected temperature.
    """
    if spectrum.technique != "raman":
        raise ValueError(f"expected raman spectrum, got {spectrum.technique}")
    if spectrum.axis.min() > -mode_cm or spectrum.axis.max() < mode_cm:
        raise ValueError(
            f"spectrum must span [-{mode_cm}, +{mode_cm}] cm-1; "
            f"got [{spectrum.axis.min():.0f}, {spectrum.axis.max():.0f}]"
        )
    stokes = fit_voigt(spectrum, around=mode_cm, window=window_cm)
    anti = fit_voigt(spectrum, around=-mode_cm, window=window_cm)
    if stokes is None or anti is None or stokes.height <= 0 or anti.height <= 0:
        raise RuntimeError("could not fit both Stokes and anti-Stokes peaks")
    ratio = anti.height / stokes.height
    return temp_mod.infer_temperature(mode_cm, ratio)
