"""Sample-temperature physics for Raman analysis.

Two operating points are first-class — liquid-nitrogen (77 K) and room
temperature (~295 K) — but the helpers accept any positive Kelvin temperature.

What temperature changes:

  - phonon population (Bose-Einstein) drives Stokes/anti-Stokes intensity ratio
  - peak FWHM scales roughly as sqrt(T) (Klemens decay model, simplified)
  - peak position blue-shifts at low T from lattice contraction (Grüneisen)
  - Cr3+/N-V/H3 photoluminescence lines sharpen by ~10x at LN2

These helpers are heuristic but physically grounded — sufficient to make the
synthetic spectra of `checkmsg.synthetic` look right and to demonstrate
thermometric inversion of a measured Stokes/AS ratio.
"""

from __future__ import annotations

import math

# h*c/k in cm * K — multiply by Raman shift in cm-1 to get hbar*omega/k in K.
_HCK_CM_K = 1.4387768775039337  # exact: 100 * h * c / k_B with c in m/s

ROOM_K: float = 295.0
LN2_K: float = 77.0
SUPPORTED_TEMPERATURES_K: tuple[float, ...] = (LN2_K, ROOM_K)


def bose_einstein(omega_cm: float, T_K: float) -> float:
    """Phonon occupation number n(omega, T) = 1 / (exp(hbar*omega/kT) - 1).

    `omega_cm` is the mode frequency in cm-1; `T_K` the absolute temperature.
    """
    if T_K <= 0:
        raise ValueError("temperature must be > 0 K")
    x = _HCK_CM_K * omega_cm / T_K
    if x > 50:  # avoid overflow for high freq / low T
        return 0.0
    return 1.0 / (math.expm1(x))


def stokes_antistokes_ratio(omega_cm: float, T_K: float) -> float:
    """Anti-Stokes/Stokes intensity ratio for a non-resonant Raman mode.

    Formal expression (within the harmonic, non-resonant approximation):
      I_AS / I_S = ((nu_0 + omega) / (nu_0 - omega))^4 * exp(-hc*omega / kT)
    For most Raman experiments the (nu_0 +/- omega) factors are within a few
    percent, so we drop them and return just the Boltzmann factor — that's
    what's used for thermometry calibration.
    """
    if T_K <= 0:
        raise ValueError("temperature must be > 0 K")
    return math.exp(-_HCK_CM_K * omega_cm / T_K)


def infer_temperature(omega_cm: float, antistokes_over_stokes: float) -> float:
    """Invert `stokes_antistokes_ratio` to estimate sample temperature."""
    if antistokes_over_stokes <= 0:
        raise ValueError("ratio must be > 0")
    if antistokes_over_stokes >= 1.0:
        raise ValueError("AS/S >= 1 implies T = infinity; check input")
    return -_HCK_CM_K * omega_cm / math.log(antistokes_over_stokes)


def fwhm_factor(T_K: float, T_ref: float = ROOM_K, floor: float = 0.30) -> float:
    """Multiplicative scaling of peak FWHM at temperature T relative to T_ref.

    Approximates Klemens-Hart-Loudon four-phonon decay: width grows as ~sqrt(T)
    over a wide range. Below ~50 K the width is dominated by spectrometer
    resolution and inhomogeneous broadening, so we clamp with a `floor`.
    """
    if T_K <= 0:
        raise ValueError("temperature must be > 0 K")
    return max(floor, math.sqrt(T_K / T_ref))


def phonon_shift(
    omega_cm: float,
    T_K: float,
    T_ref: float = ROOM_K,
    alpha_thermal: float = 1.5e-5,
    grueneisen: float = 1.0,
) -> float:
    """Approximate temperature-induced peak shift (cm-1) relative to `T_ref`.

    Cooling tightens the lattice, hardening the phonon frequency. Sign convention
    follows the experimental observation: at low T, peaks blue-shift (positive
    return value).
    """
    return -2.0 * grueneisen * alpha_thermal * omega_cm * (T_K - T_ref)
