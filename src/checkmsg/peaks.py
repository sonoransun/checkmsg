from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from lmfit.models import VoigtModel
from scipy.signal import find_peaks, peak_widths

from checkmsg.spectrum import Spectrum


@dataclass(frozen=True)
class Peak:
    position: float
    height: float
    width: float
    snr: float


def detect(spectrum: Spectrum, prominence: float | None = None, distance: float | None = None,
           min_snr: float = 5.0) -> list[Peak]:
    y = np.asarray(spectrum.intensity, dtype=float)
    x = np.asarray(spectrum.axis, dtype=float)
    if y.size < 5:
        return []
    noise = _estimate_noise(y)
    if prominence is None:
        prominence = max(min_snr * noise, (np.max(y) - np.min(y)) * 0.02)
    dist_samples = None
    if distance is not None:
        step = float(np.median(np.diff(x))) or 1.0
        dist_samples = max(1, int(distance / step))
    idx, props = find_peaks(y, prominence=prominence, distance=dist_samples)
    if idx.size == 0:
        return []
    widths_samples, _, _, _ = peak_widths(y, idx, rel_height=0.5)
    step = float(np.median(np.diff(x)))
    out: list[Peak] = []
    for k, i in enumerate(idx):
        snr = float(props["prominences"][k] / max(noise, 1e-12))
        if snr < min_snr:
            continue
        out.append(Peak(
            position=float(x[i]),
            height=float(y[i]),
            width=float(widths_samples[k] * step),
            snr=snr,
        ))
    return out


def fit_voigt(spectrum: Spectrum, around: float, window: float) -> Peak | None:
    """Fit a single Voigt profile in [around-window, around+window] for sub-pixel position."""
    s = spectrum.slice(around - window, around + window)
    if len(s) < 5:
        return None
    model = VoigtModel()
    pars = model.guess(s.intensity, x=s.axis)
    res = model.fit(s.intensity, pars, x=s.axis)
    center = float(res.params["center"].value)
    height = float(np.max(model.eval(res.params, x=s.axis)))
    sigma = float(res.params["sigma"].value)
    gamma_param = res.params["gamma"] if "gamma" in res.params else res.params["sigma"]
    gamma = float(gamma_param.value)
    fwhm_g = 2.0 * sigma * np.sqrt(2.0 * np.log(2.0))
    fwhm_l = 2.0 * gamma
    fwhm = 0.5346 * fwhm_l + np.sqrt(0.2166 * fwhm_l ** 2 + fwhm_g ** 2)
    noise = _estimate_noise(spectrum.intensity)
    return Peak(position=center, height=height, width=float(fwhm), snr=float(height / max(noise, 1e-12)))


def _estimate_noise(y: np.ndarray) -> float:
    """Robust noise estimate via MAD of first differences."""
    d = np.diff(np.asarray(y, dtype=float))
    if d.size == 0:
        return 1.0
    mad = float(np.median(np.abs(d - np.median(d))))
    return mad * 1.4826 / np.sqrt(2.0) or 1e-6
