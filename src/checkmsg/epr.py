"""Electron paramagnetic resonance simulator and analysis.

Implements a bounded but real spin-Hamiltonian simulator suitable for the
gemologically relevant centers (S in {1/2, 1, 3/2, 2, 5/2}, up to two
hyperfine-coupled nuclei, anisotropic g and zero-field splitting). The
simulator builds the spin Hamiltonian in the |M_S, M_I, ...> product basis,
diagonalises it on a (theta, phi, B) grid, applies the CW transition selection
rule via the perpendicular component of the spin operator, and outputs either
the absorption spectrum or its first field derivative — the latter being what
a CW EPR spectrometer with field modulation actually records.

The unit convention is: magnetic field in mT, energy in MHz, frequency in GHz.
The conversion mu_B / h ~= 13.99624 MHz/mT is the load-bearing constant.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Literal

import numpy as np

from checkmsg.spectrum import Spectrum

# Bohr magneton over Planck constant: hbar * 2 pi * f = g * mu_B * B  =>  f[MHz] = g * 13.99624 * B[mT]
MU_B_OVER_H_MHZ_PER_MT = 13.99624491
# Nuclear magneton (about 1/1836 of Bohr).
MU_N_OVER_H_MHZ_PER_MT = 0.007622593285

LineShape = Literal["lorentzian", "gaussian"]


@dataclass(frozen=True)
class Hyperfine:
    """A nucleus coupled to the electronic spin."""

    nucleus: str
    I: float  # noqa: E741 — standard physics notation for nuclear spin
    A_iso_MHz: float
    A_aniso_MHz: tuple[float, float, float] = (0.0, 0.0, 0.0)
    gn: float = 0.0  # nuclear Zeeman g-factor; ~0 unless explicitly modeled


@dataclass(frozen=True)
class SpinSystem:
    """A paramagnetic center described by its spin Hamiltonian parameters."""

    name: str
    S: float
    g: float | tuple[float, float, float] = 2.0
    D_MHz: float = 0.0
    E_MHz: float = 0.0
    hyperfine: tuple[Hyperfine, ...] = ()
    linewidth_mT: float = 0.3
    lineshape: LineShape = "lorentzian"
    host: str = ""

    def g_tensor(self) -> tuple[float, float, float]:
        if isinstance(self.g, tuple):
            return self.g
        return (float(self.g), float(self.g), float(self.g))

    @property
    def is_isotropic(self) -> bool:
        gx, gy, gz = self.g_tensor()
        if abs(gx - gy) > 1e-9 or abs(gy - gz) > 1e-9:
            return False
        if self.D_MHz or self.E_MHz:
            return False
        for h in self.hyperfine:
            if any(abs(a) > 1e-9 for a in h.A_aniso_MHz):
                return False
        return True

    @property
    def is_axial(self) -> bool:
        gx, gy, gz = self.g_tensor()
        if self.E_MHz:
            return False
        if abs(gx - gy) > 1e-9:
            return False
        for h in self.hyperfine:
            if abs(h.A_aniso_MHz[0] - h.A_aniso_MHz[1]) > 1e-9:
                return False
        return True


@dataclass(frozen=True)
class EprCandidate:
    name: str
    cosine: float
    g_score: float

    @property
    def combined(self) -> float:
        return 0.6 * self.cosine + 0.4 * self.g_score


@dataclass
class EprResult:
    spectrum: Spectrum
    g_factors: list[float]
    candidates: list[EprCandidate]

    @property
    def best(self) -> EprCandidate | None:
        return self.candidates[0] if self.candidates else None


@dataclass(frozen=True)
class SpinCount:
    total: float
    per_gram: float | None
    ratio_to_reference: float


# ---------- Spin algebra ----------

@lru_cache(maxsize=64)
def spin_matrices(S: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (Sx, Sy, Sz) for spin S as complex (2S+1)x(2S+1) ndarrays.

    Basis is descending in m: |S>, |S-1>, ..., |-S>.
    """
    dim = int(round(2.0 * S + 1.0))
    if dim < 1 or abs(2.0 * S + 1.0 - dim) > 1e-9:
        raise ValueError(f"S={S} must be a non-negative half-integer")
    m = np.array([S - k for k in range(dim)], dtype=float)
    Sz = np.diag(m).astype(complex)
    Splus = np.zeros((dim, dim), dtype=complex)
    for k in range(dim - 1):
        m_low = m[k + 1]
        Splus[k, k + 1] = np.sqrt(S * (S + 1) - m_low * (m_low + 1))
    Sminus = Splus.conj().T
    Sx = 0.5 * (Splus + Sminus)
    Sy = -0.5j * (Splus - Sminus)
    return Sx, Sy, Sz


def _kron_chain(matrices: list[np.ndarray]) -> np.ndarray:
    out = matrices[0]
    for m in matrices[1:]:
        out = np.kron(out, m)
    return out


def _build_operators(spin_system: SpinSystem):
    """Return Sx, Sy, Sz (electron) and per-nucleus (Ix, Iy, Iz), all in the product basis."""
    Sx_e, Sy_e, Sz_e = spin_matrices(spin_system.S)
    nuc_dims = [int(round(2.0 * h.I + 1.0)) for h in spin_system.hyperfine]
    nuc_eyes = [np.eye(d, dtype=complex) for d in nuc_dims]

    Sx = _kron_chain([Sx_e, *nuc_eyes])
    Sy = _kron_chain([Sy_e, *nuc_eyes])
    Sz = _kron_chain([Sz_e, *nuc_eyes])

    e_eye = np.eye(int(round(2 * spin_system.S + 1)), dtype=complex)
    nuc_ops = []
    for i, hf in enumerate(spin_system.hyperfine):
        Ix_i, Iy_i, Iz_i = spin_matrices(hf.I)
        before = nuc_eyes[:i]
        after = nuc_eyes[i + 1:]
        Ix = _kron_chain([e_eye, *before, Ix_i, *after])
        Iy = _kron_chain([e_eye, *before, Iy_i, *after])
        Iz = _kron_chain([e_eye, *before, Iz_i, *after])
        nuc_ops.append((Ix, Iy, Iz))

    return Sx, Sy, Sz, nuc_ops


def _static_hamiltonian(spin_system, Sx, Sy, Sz, nuc_ops) -> np.ndarray:
    """Field-independent part: D/E zero-field splitting + S.A.I hyperfine terms."""
    S = spin_system.S
    H = np.zeros_like(Sx)
    if spin_system.D_MHz:
        H = H + spin_system.D_MHz * (Sz @ Sz - (S * (S + 1) / 3.0) * np.eye(Sx.shape[0]))
    if spin_system.E_MHz:
        H = H + spin_system.E_MHz * (Sx @ Sx - Sy @ Sy)
    for hf, (Ix, Iy, Iz) in zip(spin_system.hyperfine, nuc_ops, strict=True):
        Axx = hf.A_iso_MHz + hf.A_aniso_MHz[0]
        Ayy = hf.A_iso_MHz + hf.A_aniso_MHz[1]
        Azz = hf.A_iso_MHz + hf.A_aniso_MHz[2]
        H = H + Axx * (Sx @ Ix) + Ayy * (Sy @ Iy) + Azz * (Sz @ Iz)
    return H


def _zeeman_per_mT(spin_system, n_unit, Sx, Sy, Sz, nuc_ops) -> np.ndarray:
    """Coefficient of B (mT) in the Hamiltonian (units MHz) for the given orientation."""
    gx, gy, gz = spin_system.g_tensor()
    bohr = MU_B_OVER_H_MHZ_PER_MT
    Hz = bohr * (gx * n_unit[0] * Sx + gy * n_unit[1] * Sy + gz * n_unit[2] * Sz)
    for hf, (Ix, Iy, Iz) in zip(spin_system.hyperfine, nuc_ops, strict=True):
        if hf.gn:
            mu_n = MU_N_OVER_H_MHZ_PER_MT
            Hz = Hz - hf.gn * mu_n * (n_unit[0] * Ix + n_unit[1] * Iy + n_unit[2] * Iz)
    return Hz


# ---------- Powder simulation ----------

def _orientation_grid(spin_system: SpinSystem, n_theta: int, n_phi: int):
    """Return (thetas, phis, weights) given the system's symmetry.

    For axial systems (E=0, g_xx==g_yy) only theta is needed; phi is one point with weight 2 pi.
    For isotropic systems both thetas and phis collapse to a single (0, 0) point.
    """
    if spin_system.is_isotropic:
        return np.array([0.0]), np.array([0.0]), np.array([4 * np.pi])
    if spin_system.is_axial:
        thetas = np.linspace(0.0, np.pi / 2, max(n_theta, 5))
        dt = np.pi / 2 / max(len(thetas) - 1, 1)
        weights = np.sin(thetas) * dt * 2 * (2 * np.pi)  # factor 2 for half-sphere, 2pi for phi
        weights[0] *= 0.5
        weights[-1] *= 0.5
        phis = np.zeros_like(thetas)
        return thetas, phis, weights
    # Rhombic / lower symmetry: integrate over half sphere (sufficient for E real).
    thetas = np.linspace(0.0, np.pi / 2, max(n_theta, 5))
    phis = np.linspace(0.0, 2 * np.pi, max(n_phi, 5), endpoint=False)
    sin_t = np.sin(thetas)
    dt = np.pi / 2 / max(len(thetas) - 1, 1)
    dp = 2 * np.pi / len(phis)
    w_t = sin_t * dt * 2  # factor 2 because we cover half sphere
    w_t[0] *= 0.5
    w_t[-1] *= 0.5
    weights = (w_t[:, None] * np.full((len(thetas), len(phis)), dp)).ravel()
    thetas_grid, phis_grid = np.meshgrid(thetas, phis, indexing="ij")
    return thetas_grid.ravel(), phis_grid.ravel(), weights


def simulate_field_sweep(
    spin_system: SpinSystem,
    frequency_GHz: float,
    fields_mT: np.ndarray,
    *,
    orientations: tuple[int, int] = (15, 15),
    derivative: bool = True,
    temperature_K: float = 295.0,
) -> Spectrum:
    """Simulate a CW EPR field-sweep at one frequency. Returns a derivative Spectrum by default."""
    fields_mT = np.asarray(fields_mT, dtype=float)
    Sx, Sy, Sz, nuc_ops = _build_operators(spin_system)
    H0 = _static_hamiltonian(spin_system, Sx, Sy, Sz, nuc_ops)

    target_MHz = frequency_GHz * 1000.0
    g_avg = sum(spin_system.g_tensor()) / 3.0
    sigma_MHz = spin_system.linewidth_mT * g_avg * MU_B_OVER_H_MHZ_PER_MT
    if sigma_MHz <= 0:
        sigma_MHz = 1.0  # safety floor

    thetas, phis, weights = _orientation_grid(spin_system, *orientations)

    dim = Sx.shape[0]
    mask = np.triu(np.ones((dim, dim), dtype=bool), k=1)
    absorption = np.zeros_like(fields_mT)

    for theta, phi, w in zip(thetas, phis, weights, strict=True):
        sin_t, cos_t = float(np.sin(theta)), float(np.cos(theta))
        sin_p, cos_p = float(np.sin(phi)), float(np.cos(phi))
        n = np.array([sin_t * cos_p, sin_t * sin_p, cos_t])
        if sin_t < 1e-10:
            e1 = np.array([1.0, 0.0, 0.0])
            e2 = np.array([0.0, 1.0, 0.0])
        else:
            e1 = np.array([cos_t * cos_p, cos_t * sin_p, -sin_t])
            e2 = np.array([-sin_p, cos_p, 0.0])

        Hz_unit = _zeeman_per_mT(spin_system, n, Sx, Sy, Sz, nuc_ops)
        S_e1 = e1[0] * Sx + e1[1] * Sy + e1[2] * Sz
        S_e2 = e2[0] * Sx + e2[1] * Sy + e2[2] * Sz

        H_all = H0[None, :, :] + fields_mT[:, None, None] * Hz_unit[None, :, :]
        eigvals, eigvecs = np.linalg.eigh(H_all)

        psi_dag = eigvecs.conj().swapaxes(-1, -2)
        Me1 = psi_dag @ S_e1 @ eigvecs
        Me2 = psi_dag @ S_e2 @ eigvecs
        mom2 = (np.abs(Me1) ** 2 + np.abs(Me2) ** 2)

        # eigvals[..., j] - eigvals[..., a]; arrange so dE[..., a, b] = E_b - E_a
        dE = eigvals[:, None, :] - eigvals[:, :, None]
        # Lorentzian on ΔE - hν
        if spin_system.lineshape == "gaussian":
            shape = (1.0 / (sigma_MHz * np.sqrt(2 * np.pi))) * np.exp(
                -0.5 * ((dE - target_MHz) / sigma_MHz) ** 2
            )
        else:
            shape = (sigma_MHz / np.pi) / ((dE - target_MHz) ** 2 + sigma_MHz ** 2)

        contrib = np.sum(mom2 * shape * mask[None, :, :], axis=(1, 2))
        absorption = absorption + w * contrib

    # Normalise (convenience for downstream cosine matching)
    peak = float(np.max(np.abs(absorption)))
    if peak > 0:
        absorption = absorption / peak

    if derivative:
        y = np.gradient(absorption, fields_mT)
    else:
        y = absorption

    return Spectrum(
        fields_mT, y, "epr", "mT",
        metadata={
            "frequency_GHz": frequency_GHz,
            "spin_system": spin_system.name,
            "derivative": derivative,
            "temperature_K": temperature_K,
        },
    )


# ---------- Analysis ----------

def infer_g_factors(spectrum: Spectrum, frequency_GHz: float | None = None,
                    prominence: float | None = None) -> list[float]:
    """Estimate g-values from zero-crossings of a derivative EPR spectrum."""
    if spectrum.technique != "epr":
        raise ValueError(f"expected epr spectrum, got {spectrum.technique}")
    if frequency_GHz is None:
        frequency_GHz = spectrum.metadata.get("frequency_GHz")
        if frequency_GHz is None:
            raise ValueError("frequency_GHz is not in metadata; pass it explicitly")
    y = np.asarray(spectrum.intensity, dtype=float)
    x = np.asarray(spectrum.axis, dtype=float)
    # zero crossings from + to - (centre of derivative line)
    sign_changes = np.where((y[:-1] > 0) & (y[1:] <= 0))[0]
    if prominence is None:
        prominence = 0.05 * float(np.max(np.abs(y)) or 1.0)
    out: list[float] = []
    for i in sign_changes:
        # Magnitude at this crossing — drop tiny derivative wiggles.
        local_amp = max(abs(y[max(i - 5, 0):i + 5]).max(), 1e-12)
        if local_amp < prominence:
            continue
        # Linear interpolation to find precise zero crossing
        if y[i] - y[i + 1] != 0:
            frac = y[i] / (y[i] - y[i + 1])
            x0 = x[i] + frac * (x[i + 1] - x[i])
        else:
            x0 = x[i]
        if x0 > 0:
            # g = nu_MHz / (mu_B/h * B_mT)
            out.append(float(frequency_GHz * 1000.0 / (MU_B_OVER_H_MHZ_PER_MT * x0)))
    return out


def _common_grid_cosine(a: Spectrum, b: Spectrum) -> float:
    """Cosine similarity between two derivative EPR spectra after restricting to overlap."""
    lo = max(a.axis[0], b.axis[0])
    hi = min(a.axis[-1], b.axis[-1])
    if hi <= lo:
        return 0.0
    grid = np.linspace(lo, hi, max(256, min(len(a), len(b))))
    ya = np.interp(grid, a.axis, a.intensity)
    yb = np.interp(grid, b.axis, b.intensity)
    na = float(np.linalg.norm(ya))
    nb = float(np.linalg.norm(yb))
    if na == 0 or nb == 0:
        return 0.0
    return float(np.dot(ya, yb) / (na * nb))


def analyze(
    spectrum: Spectrum,
    frequency_GHz: float | None = None,
    candidates: dict[str, SpinSystem] | None = None,
    fields_mT: np.ndarray | None = None,
    *,
    top: int = 5,
    orientations: tuple[int, int] = (11, 11),
) -> EprResult:
    """Identify the paramagnetic center(s) by simulating each candidate at the same conditions."""
    if spectrum.technique != "epr":
        raise ValueError(f"expected epr spectrum, got {spectrum.technique}")
    if frequency_GHz is None:
        frequency_GHz = spectrum.metadata.get("frequency_GHz")
        if frequency_GHz is None:
            raise ValueError("frequency_GHz must be supplied (or be in spectrum.metadata)")
    if candidates is None:
        from checkmsg.refdata import epr_centers
        candidates = epr_centers.CENTERS

    if fields_mT is None:
        fields_mT = spectrum.axis

    g_obs = infer_g_factors(spectrum, frequency_GHz=frequency_GHz)
    out: list[EprCandidate] = []
    for name, sys_ in candidates.items():
        try:
            sim = simulate_field_sweep(
                sys_, frequency_GHz=frequency_GHz, fields_mT=fields_mT,
                orientations=orientations,
            )
        except Exception:
            continue
        g_sim = infer_g_factors(sim, frequency_GHz=frequency_GHz)
        cos = _common_grid_cosine(spectrum, sim)
        g_score = _g_match_score(g_obs, g_sim, tolerance=0.005)
        out.append(EprCandidate(name=name, cosine=cos, g_score=g_score))
    out.sort(key=lambda c: c.combined, reverse=True)
    return EprResult(spectrum=spectrum, g_factors=g_obs, candidates=out[:top])


def _g_match_score(obs: list[float], sim: list[float], tolerance: float) -> float:
    """Symmetric Jaccard-like overlap of g-factor lists within `tolerance`."""
    if not obs or not sim:
        return 0.0
    used = set()
    matched = 0
    for g in obs:
        best_j, best_d = None, tolerance + 1
        for j, gs in enumerate(sim):
            if j in used:
                continue
            d = abs(g - gs)
            if d < best_d:
                best_d = d
                best_j = j
        if best_j is not None and best_d <= tolerance:
            used.add(best_j)
            matched += 1
    union = len(obs) + len(sim) - matched
    return float(matched / union) if union else 0.0


# ---------- Quantification ----------

def double_integral(spectrum: Spectrum) -> float:
    """Double-integrate a derivative EPR spectrum to get a number proportional to total spins."""
    if spectrum.technique != "epr":
        raise ValueError(f"expected epr spectrum, got {spectrum.technique}")
    x = np.asarray(spectrum.axis, dtype=float)
    y = np.asarray(spectrum.intensity, dtype=float)
    # First integral: cumulative -> absorption
    absorbance = np.cumulative_trapezoid(y, x, initial=0.0) if hasattr(np, "cumulative_trapezoid") \
        else _cumulative_trapezoid(y, x)
    # Subtract any DC offset/baseline drift in the first integral (linear detrend on tails)
    n = len(x)
    base = np.linspace(absorbance[0], absorbance[-1], n)
    absorbance = absorbance - base
    # Second integral: total area under absorption.
    return float(np.trapezoid(absorbance, x))


def _cumulative_trapezoid(y, x):
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


def count_spins(
    sample: Spectrum,
    reference: Spectrum,
    reference_spins: float,
    *,
    sample_mass_g: float | None = None,
    reference_mass_g: float | None = None,
) -> SpinCount:
    """Reference-standard spin quantification.

    Compares the double-integral of `sample` and `reference` (both must share frequency,
    modulation amplitude, and microwave power for the comparison to be meaningful).
    Returns absolute spin count and (when `sample_mass_g` given) spins per gram.
    """
    if sample.technique != "epr" or reference.technique != "epr":
        raise ValueError("both spectra must be EPR")
    fs = sample.metadata.get("frequency_GHz")
    fr = reference.metadata.get("frequency_GHz")
    if fs is not None and fr is not None and abs(fs - fr) / max(fr, 1e-9) > 0.05:
        raise ValueError(f"frequency mismatch: sample {fs} vs reference {fr} GHz")
    di_s = double_integral(sample)
    di_r = double_integral(reference)
    if di_r == 0:
        raise ValueError("reference double integral is zero")
    ratio = di_s / di_r
    # Adjust for reference mass scaling — the user supplies absolute spins in `reference_spins`,
    # so no additional scaling is needed unless reference_mass_g is given (then we report total).
    total = ratio * reference_spins
    per_gram = total / sample_mass_g if sample_mass_g else None
    return SpinCount(total=total, per_gram=per_gram, ratio_to_reference=ratio)


def absolute_spin_count(
    spectrum: Spectrum,
    *,
    cavity_Q: float,
    modulation_amplitude_mT: float,
    microwave_power_mW: float,
    calibration_factor: float = 1.0,
    sample_mass_g: float | None = None,
) -> SpinCount:
    """Absolute-mode quantification using instrument calibration.

    The expected signal scales as Q * modulation_amplitude * sqrt(power) * N_spins;
    so N_spins ~ double_integral / (Q * h_mod * sqrt(P) * calibration_factor).

    NOTE: real-instrument absolute counts require careful calibration of B1, sample
    geometry, dielectric losses, etc. This function exposes the formula-correct path
    for users who have all the calibration constants; tests cover formula correctness
    only, not measurement-level fidelity.
    """
    if spectrum.technique != "epr":
        raise ValueError("expected epr spectrum")
    if cavity_Q <= 0 or modulation_amplitude_mT <= 0 or microwave_power_mW <= 0:
        raise ValueError("Q, modulation amplitude and power must be > 0")
    di = double_integral(spectrum)
    n = di / (cavity_Q * modulation_amplitude_mT * np.sqrt(microwave_power_mW) * calibration_factor)
    per_gram = float(n) / sample_mass_g if sample_mass_g else None
    return SpinCount(total=float(n), per_gram=per_gram, ratio_to_reference=float("nan"))
