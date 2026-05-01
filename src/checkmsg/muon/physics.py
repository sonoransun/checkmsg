"""Muon physics primitives: Bethe-Bloch, Highland, CSDA range.

Three short, transparent helpers — no external particle-physics libraries. All
use textbook formulas valid to ~5–10 % in the typical surface-muon to MIP regime
(30 MeV/c to 1 GeV/c) for compact compound materials. Hadronic interactions,
muon decay in flight, and pair-production losses are out of scope.
"""

from __future__ import annotations

import math

from checkmsg.refdata.muon_data import K_BB, M_E_MeV, M_MU_MeV, Material


def _beta_gamma(p_MeV: float) -> tuple[float, float]:
    """Return (beta, gamma) for a muon of momentum p (MeV/c)."""
    if p_MeV <= 0:
        raise ValueError("momentum must be > 0 MeV/c")
    E = math.hypot(p_MeV, M_MU_MeV)
    gamma = E / M_MU_MeV
    beta = p_MeV / E
    return beta, gamma


def bethe_bloch_dE_dx(p_MeV: float, material: Material) -> float:
    """Mean ionisation loss for a muon, in MeV per (g cm⁻²) of material.

    The relativistic Bethe-Bloch formula:

        −dE/(ρ dx) = K (Z/A) (1/β²) [ ln(2 m_e β² γ² T_max / I²) − β² − δ/2 ]

    Density-effect correction δ is set to zero (small for the regime of interest);
    K = 0.30707 MeV cm² mol⁻¹.
    """
    if material.X_0_g_cm2 <= 0 or material.density_g_cc <= 0:
        return 0.0
    beta, gamma = _beta_gamma(p_MeV)
    if beta <= 0:
        return 0.0
    # Maximum kinetic energy transfer in a single collision (PDG eq.):
    #   T_max = 2 m_e β² γ² / (1 + 2 γ m_e/m_μ + (m_e/m_μ)²)
    me_over_mmu = M_E_MeV / M_MU_MeV
    T_max = (2.0 * M_E_MeV * (beta * gamma) ** 2
             / (1.0 + 2.0 * gamma * me_over_mmu + me_over_mmu ** 2))
    I_MeV = material.I_eV * 1e-6
    log_term = math.log(2.0 * M_E_MeV * (beta * gamma) ** 2 * T_max / I_MeV ** 2)
    # PDG 28.31.5: -dE/(rho dx) = K (Z/A) (1/β²) [ (1/2) ln(...) - β² - δ/2 ]
    # Density-effect correction δ is set to zero for the regime of interest.
    return K_BB * (material.Z_eff / material.A_eff) * (1.0 / beta ** 2) * (0.5 * log_term - beta ** 2)


def highland_scattering_rms_mrad(p_MeV: float, x_g_cm2: float, X_0_g_cm2: float) -> float:
    """RMS plane-projected scattering angle in milliradians (Highland formula).

        θ₀ = (13.6 / β c p) × √(x / X₀) × (1 + 0.038 ln(x / X₀))

    The 13.6 has units MeV; the input p is in MeV/c, x in g/cm², X₀ in g/cm². The
    result is the standard deviation of the projected angle distribution after
    traversing thickness x.
    """
    if x_g_cm2 <= 0 or X_0_g_cm2 <= 0:
        return 0.0
    beta, _ = _beta_gamma(p_MeV)
    if beta <= 0:
        return 0.0
    ratio = x_g_cm2 / X_0_g_cm2
    return (13.6 / (beta * p_MeV)) * math.sqrt(ratio) * (1.0 + 0.038 * math.log(max(ratio, 1e-12))) * 1000.0


def csda_range_g_cm2(p_MeV: float, material: Material, *, n_steps: int = 200) -> float:
    """Continuous-slowing-down approximation range in g cm⁻² for a muon entering with
    momentum p_MeV.

    Numerical integration of 1 / (dE/dx) over the kinetic energy from T_in down to
    a cutoff (1 MeV by default). The result is the *mass thickness* needed to bring
    the muon to rest; divide by the material density for the linear range in cm.
    """
    beta, gamma = _beta_gamma(p_MeV)
    T_in = (gamma - 1.0) * M_MU_MeV
    if T_in <= 1.0:
        return 0.0
    energies = [T_in - i * (T_in - 1.0) / n_steps for i in range(n_steps + 1)]
    integral = 0.0
    for j in range(n_steps):
        T_lo, T_hi = energies[j + 1], energies[j]
        E_lo = T_lo + M_MU_MeV
        E_hi = T_hi + M_MU_MeV
        p_lo = math.sqrt(max(E_lo ** 2 - M_MU_MeV ** 2, 0.0))
        p_hi = math.sqrt(max(E_hi ** 2 - M_MU_MeV ** 2, 0.0))
        dEdx_lo = bethe_bloch_dE_dx(p_lo, material) if p_lo > 0 else float("inf")
        dEdx_hi = bethe_bloch_dE_dx(p_hi, material) if p_hi > 0 else float("inf")
        avg_inv = 0.5 * (1.0 / max(dEdx_lo, 1e-12) + 1.0 / max(dEdx_hi, 1e-12))
        integral += avg_inv * (T_hi - T_lo)
    return integral
