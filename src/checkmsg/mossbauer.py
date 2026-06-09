"""⁵⁷Fe Mössbauer analysis — Fe²⁺/Fe³⁺ valence + site assignment.

The observable is a velocity spectrum (mm/s) of resonant-absorption dips: a
quadrupole *doublet* (two lines split by ΔEQ, centred on the isomer shift δ) for
paramagnetic Fe, or a six-line *sextet* for a magnetically ordered oxide. This
module simulates and extracts those parameters and scores them against the
``refdata.mossbauer_sites`` table (squid-style residual matching), staying
``Spectrum``-native (the velocity axis is single-valued, so no custom adapter is
needed). It is pedagogical: no thickness/texture/recoil-free-fraction correction.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from checkmsg import preprocess
from checkmsg.peaks import Peak, detect, prominent
from checkmsg.refdata.mossbauer_sites import SITES, FeSite
from checkmsg.spectrum import Spectrum

# Canonical ⁵⁷Fe magnetic sextet: line positions (mm/s) for α-Fe (B_hf ≈ 33 T),
# scaled linearly with B_hf, and the 3:2:1:1:2:3 intensity pattern.
_ALPHA_FE_LINES = np.array([-5.3, -3.1, -0.9, 0.9, 3.1, 5.3])
_SEXTET_INTENS = np.array([3.0, 2.0, 1.0, 1.0, 2.0, 3.0])


def _lorentz_dip(v: np.ndarray, center: float, hwhm: float, depth: float) -> np.ndarray:
    return depth / (1.0 + ((v - center) / hwhm) ** 2)


def simulate_doublet(v: np.ndarray, delta: float, d_eq: float, *,
                     hwhm: float = 0.15, depth: float = 0.12) -> np.ndarray:
    """Absorption (positive) of a quadrupole doublet centred on ``delta``."""
    return (_lorentz_dip(v, delta - d_eq / 2.0, hwhm, depth)
            + _lorentz_dip(v, delta + d_eq / 2.0, hwhm, depth))


def simulate_sextet(v: np.ndarray, delta: float, b_hf_T: float, *,
                    hwhm: float = 0.2, depth: float = 0.1) -> np.ndarray:
    """Absorption (positive) of a magnetic sextet with hyperfine field ``b_hf_T``."""
    positions = delta + _ALPHA_FE_LINES * (b_hf_T / 33.0)
    y = np.zeros_like(v)
    for pos, w in zip(positions, _SEXTET_INTENS, strict=True):
        y = y + _lorentz_dip(v, pos, hwhm, depth * w / 3.0)
    return y


def simulate_mossbauer(site: FeSite, *, velocities_mm_s: np.ndarray | None = None,
                       noise: float = 0.004, seed: int | None = 0) -> Spectrum:
    """Synthesise a transmission Mössbauer spectrum for one Fe site (dips below 1)."""
    v = np.linspace(-11.0, 11.0, 881) if site.magnetic_sextet else np.linspace(-4.0, 4.0, 401)
    if site.magnetic_sextet:
        absorption = simulate_sextet(v, site.isomer_shift_mm_s, site.hyperfine_field_T)
    else:
        absorption = simulate_doublet(v, site.isomer_shift_mm_s, site.quadrupole_splitting_mm_s)
    y = 1.0 - absorption  # transmission
    if noise > 0:
        rng = np.random.default_rng(seed)
        y = y + rng.normal(0.0, noise, size=y.size)
    return Spectrum(v, y, "mossbauer", "mm/s", {"site": site.name})


@dataclass(frozen=True)
class MossbauerCandidate:
    name: str          # site key
    residual: float    # normalised (δ, ΔEQ) residual; lower is better
    combined: float    # 1/(1+residual)


@dataclass
class MossbauerResult:
    cleaned: Spectrum
    peaks: list[Peak]
    extracted: dict
    candidates: list[MossbauerCandidate]

    @property
    def best(self) -> MossbauerCandidate | None:
        return self.candidates[0] if self.candidates else None

    def headline(self) -> str:
        e = self.extracted
        if not e:
            return "Mössbauer: no doublet/sextet resolved"
        kind = "sextet" if e.get("is_sextet") else "doublet"
        top = self.best.name if self.best else "?"
        return (f"Mössbauer {kind}: δ={e['delta']:.2f} ΔEQ={e['d_eq']:.2f} mm/s "
                f"({e['valence']}); top={top}")


def preprocess_mossbauer(spectrum: Spectrum) -> Spectrum:
    if spectrum.technique != "mossbauer":
        raise ValueError(f"expected mossbauer spectrum, got {spectrum.technique}")
    # Transmission dips -> absorption peaks.
    y = preprocess.savgol(spectrum.intensity, window=9, order=3)
    y = -(y - float(np.median(y)))
    return spectrum.with_intensity(y)


def _valence_from_delta(delta: float) -> str:
    if delta >= 0.9:
        return "Fe2+"
    if delta <= 0.6:
        return "Fe3+"
    return "Fe(mixed)"


def analyze(spectrum: Spectrum, min_snr: float = 5.0, top: int = 5) -> MossbauerResult:
    cleaned = preprocess_mossbauer(spectrum)
    peaks = prominent(detect(cleaned, min_snr=min_snr), frac=0.25)
    extracted: dict = {}
    if len(peaks) >= 2:
        ordered = sorted(peaks, key=lambda p: p.height, reverse=True)
        is_sextet = len(peaks) >= 6 and (max(p.position for p in peaks)
                                         - min(p.position for p in peaks)) > 6.0
        if is_sextet:
            positions = [p.position for p in peaks]
            delta = float(np.mean(positions))
            span = max(positions) - min(positions)
            b_hf = span / (2.0 * 5.3) * 33.0  # invert the α-Fe outer-line scaling
            extracted = {"delta": delta, "d_eq": 0.0, "valence": "Fe3+",
                         "is_sextet": True, "b_hf_T": b_hf}
        else:
            two = sorted(ordered[:2], key=lambda p: p.position)
            delta = float((two[0].position + two[1].position) / 2.0)
            d_eq = float(abs(two[1].position - two[0].position))
            extracted = {"delta": delta, "d_eq": d_eq,
                         "valence": _valence_from_delta(delta), "is_sextet": False}

    candidates: list[MossbauerCandidate] = []
    if extracted:
        for key, site in SITES.items():
            if extracted["is_sextet"] != site.magnetic_sextet:
                continue
            if site.magnetic_sextet:
                resid = abs(extracted["b_hf_T"] - site.hyperfine_field_T) / 5.0
            else:
                resid = (abs(extracted["delta"] - site.isomer_shift_mm_s) / site.delta_tol
                         + abs(extracted["d_eq"] - site.quadrupole_splitting_mm_s) / site.qs_tol)
            candidates.append(MossbauerCandidate(key, resid, 1.0 / (1.0 + resid)))
        candidates.sort(key=lambda c: c.residual)
        candidates = candidates[:top]
    return MossbauerResult(cleaned=cleaned, peaks=peaks, extracted=extracted, candidates=candidates)
