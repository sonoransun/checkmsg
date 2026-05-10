"""SQUID magnetometry — bulk magnetic moment and susceptibility for mineral ID.

Two complementary acquisition modes are bundled:

* **dc-SQUID** — quasi-static field sweep, two-junction interferometer.
  Best moment sensitivity; produces M(H) hysteresis loops at fixed temperature.
* **rf-SQUID** — single-junction tank-circuit lock-in readout.
  Better suited to χ(T) thermal sweeps and AC susceptibility χ'(ω) + iχ''(ω).

The module follows the LA-ICP-MS precedent for non-1-D-spectrum data: a custom
`SquidMeasurement` dataclass holds the axis/signal pair together with mode,
temperature, applied bias field, frequency, and a free-form metadata dict. A
`to_spectrum()` adapter wraps the measurement into a generic `Spectrum` so the
existing CLI/IO infrastructure works unchanged.

Forward simulators (`simulate_mh`, `simulate_chi_T`, `simulate_chi_ac`) are
physics-rich enough to capture the diagnostic features that distinguish minerals:

  - tanh-saturation hysteresis with explicit coercivity for ferro/ferri-magnets;
  - linear paramagnetic / diamagnetic loops with Curie-law slope tracking;
  - Curie-Weiss χ(T) with sharp cusp at Tc/TN for ordered phases;
  - Casimir-du Pré (Debye-relaxation) χ'(ω) + iχ''(ω) with a loss peak at ωτ=1.

The high-level entry point `analyze(measurement, candidates=...)` extracts the
ordering type plus relevant scalar parameters (coercivity, saturation moment,
Curie temperature, Weiss intercept, AC loss peak), scores each candidate
`MagneticMineral` from `refdata.squid_signatures`, and returns a ranked
`SquidResult`. Both `analyze_dc(...)` and `analyze_rf(...)` are thin wrappers
that route to `analyze` based on the measurement mode.

Out of scope: anything that needs an actual magnetometer (drift correction,
gradient-pickup geometry, He-3 cryogenics) — the module is pedagogical and the
forward models are deterministic given a `seed`.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from checkmsg.refdata.squid_signatures import (
    ORDERINGS,
    SIGNATURES,
    MagneticMineral,
    Ordering,
)
from checkmsg.spectrum import Spectrum

# ---------------------------------------------------------------------------
# Data primitives
# ---------------------------------------------------------------------------

Mode = Literal["dc-mh", "rf-chi-T", "rf-chi-ac"]


@dataclass(frozen=True)
class SquidMeasurement:
    """A SQUID magnetometer acquisition.

    `axis` and `signal` semantics depend on `mode`:

    * **dc-mh** — axis = applied field (mT), signal = magnetisation (emu/g).
    * **rf-chi-T** — axis = sample temperature (K), signal = SI volume susceptibility.
    * **rf-chi-ac** — axis = drive frequency (Hz), signal = complex susceptibility
      packed as a 1-D ``[chi_real_0, chi_imag_0, chi_real_1, ...]`` array (length 2N).

    Static fields (`temperature_K`, `applied_field_mT`, `frequency_Hz`) hold the
    parameters that were *fixed* during the sweep, and are 0.0 if not relevant
    to the mode.
    """

    mode: Mode
    axis: np.ndarray
    signal: np.ndarray
    temperature_K: float = 0.0
    applied_field_mT: float = 0.0
    frequency_Hz: float = 0.0
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Cast to ndarray. Use object.__setattr__ because the dataclass is frozen.
        ax = np.asarray(self.axis, dtype=float)
        sig = np.asarray(self.signal, dtype=float)
        if ax.ndim != 1 or sig.ndim != 1:
            raise ValueError("axis and signal must be 1-D arrays")
        if self.mode == "rf-chi-ac":
            if sig.size != 2 * ax.size:
                raise ValueError(
                    "rf-chi-ac signal must hold real+imag pairs (size = 2 * axis.size)"
                )
        else:
            if ax.size != sig.size:
                raise ValueError(
                    f"axis ({ax.size}) and signal ({sig.size}) must match for mode={self.mode}"
                )
        object.__setattr__(self, "axis", ax)
        object.__setattr__(self, "signal", sig)

    def chi_complex(self) -> tuple[np.ndarray, np.ndarray]:
        """Unpack rf-chi-ac signal into (chi_real, chi_imag). Other modes raise."""
        if self.mode != "rf-chi-ac":
            raise ValueError("chi_complex() is only defined for rf-chi-ac measurements")
        s = self.signal
        return s[0::2].copy(), s[1::2].copy()


@dataclass(frozen=True)
class SquidCandidate:
    """One ranked mineral / signature against an observed measurement."""

    name: str
    ordering_match: bool
    tc_residual_K: float
    moment_residual_emu_g: float
    susceptibility_residual: float
    combined: float


@dataclass
class SquidResult:
    """Output of `analyze()` — the extracted parameters and ranked candidates."""

    measurement: SquidMeasurement
    candidates: list[SquidCandidate]
    extracted: dict
    notes: list[str] = field(default_factory=list)

    @property
    def best(self) -> SquidCandidate | None:
        return self.candidates[0] if self.candidates else None

    def headline(self) -> str:
        ord_ = self.extracted.get("ordering", "?")
        parts = [f"ordering={ord_}"]
        if self.extracted.get("curie_K"):
            parts.append(f"Tc={self.extracted['curie_K']:.0f} K")
        if self.extracted.get("neel_K"):
            parts.append(f"TN={self.extracted['neel_K']:.0f} K")
        if self.extracted.get("saturation_emu_g"):
            parts.append(f"Ms={self.extracted['saturation_emu_g']:.2f} emu/g")
        if self.extracted.get("coercivity_mT"):
            parts.append(f"Hc={self.extracted['coercivity_mT']:.0f} mT")
        if self.best:
            parts.append(f"top={self.best.name}")
        return ", ".join(parts)


# ---------------------------------------------------------------------------
# Forward simulators
# ---------------------------------------------------------------------------


def simulate_mh(
    ordering: Ordering,
    *,
    fields_mT: np.ndarray | None = None,
    saturation_emu_g: float = 0.0,
    coercivity_mT: float = 0.0,
    susceptibility_si: float = 0.0,
    temperature_K: float = 295.0,
    noise: float = 0.005,
    seed: int | None = 0,
) -> SquidMeasurement:
    """Synthesise a dc-SQUID M(H) hysteresis loop for the given ordering type.

    Models, by ordering:

      - **ferromagnetic / ferrimagnetic** — tanh saturation with explicit
        coercivity: M(H) = Ms · tanh((H ∓ Hc) / Hk), with Hk ≈ 0.7·max(Hc, 1).
        The two branches are stitched at H=0 to give the canonical hysteretic
        loop; remanence Mr ≈ Ms · tanh(Hc/Hk).
      - **canted-afm** — small remanent moment (saturation_emu_g·0.005) + linear
        susceptibility slope. No appreciable coercivity above the Morin
        transition unless the literature value is supplied.
      - **antiferromagnetic** — small linear susceptibility, no remanence.
      - **paramagnetic** — linear M(H), slope = susceptibility · H, scaled by
        (1/T) Curie law relative to 295 K.
      - **diamagnetic** — linear with negative susceptibility (signal is
        intrinsically tiny: |χ| ~ 10⁻⁵ · ρ · H).

    The axis is in mT and forms a closed sweep H_max → −H_max → H_max so the
    returned `signal` traces the full loop.
    """
    rng = np.random.default_rng(seed)

    if fields_mT is None:
        # Two-direction sweep: descending then ascending, so a hysteresis loop is closed.
        # 4× coercivity gives clean saturation on both branches.
        h_max = max(4.0 * max(coercivity_mT, 50.0), 200.0)
        forward = np.linspace(h_max, -h_max, 201)
        backward = np.linspace(-h_max, h_max, 201)
        fields_mT = np.concatenate([forward, backward])
    fields = np.asarray(fields_mT, dtype=float)

    if ordering in ("ferromagnetic", "ferrimagnetic"):
        Hk = max(0.7 * max(coercivity_mT, 1.0), 1.0)
        # Identify which half of the loop we're on by detecting direction of travel.
        dh = np.gradient(fields)
        # Use a continuous switch: descending branch shifts the loop in +H, ascending in -H.
        offset = np.where(dh < 0.0, -coercivity_mT, +coercivity_mT)
        m = saturation_emu_g * np.tanh((fields + offset) / Hk)
    elif ordering == "canted-afm":
        # Small canted ferromagnet sitting on a linear AFM slope. Use a narrow
        # switching window (Hk = 5 % of Hc) so the loop crosses near ±Hc — this
        # is the right model for high-coercivity weak ferromagnets like hematite.
        Hk = max(0.05 * max(coercivity_mT, 1.0), 0.5)
        dh = np.gradient(fields)
        offset = np.where(dh < 0.0, -coercivity_mT, +coercivity_mT)
        canted_sat = saturation_emu_g * np.tanh((fields + offset) / Hk)
        slope_emu_g_per_mT = max(susceptibility_si, 1.0e-4) * 0.05
        m = canted_sat + slope_emu_g_per_mT * fields
    elif ordering == "antiferromagnetic":
        slope = max(susceptibility_si, 1.0e-5) * 0.05
        m = slope * fields
    elif ordering == "paramagnetic":
        chi = max(susceptibility_si, 1.0e-7)
        # Curie scaling: chi(T) ∝ 1/T.
        chi_at_T = chi * 295.0 / max(temperature_K, 1.0)
        m = chi_at_T * 0.05 * fields  # 0.05 = mT-to-emu/g conversion factor (toy)
    elif ordering == "diamagnetic":
        chi = -abs(susceptibility_si if susceptibility_si != 0 else 1.5e-5)
        m = chi * 0.05 * fields
    else:
        raise ValueError(f"unknown ordering {ordering!r}; expected one of {ORDERINGS}")

    if noise > 0:
        peak = float(np.max(np.abs(m)) or 1.0)
        m = m + rng.normal(0.0, noise * peak, size=m.size)

    return SquidMeasurement(
        mode="dc-mh",
        axis=fields,
        signal=m,
        temperature_K=temperature_K,
        applied_field_mT=0.0,
        frequency_Hz=0.0,
        metadata={
            "ordering": ordering,
            "saturation_emu_g": saturation_emu_g,
            "coercivity_mT_input": coercivity_mT,
            "noise": noise,
        },
    )


def simulate_chi_T(
    ordering: Ordering,
    *,
    temperatures_K: np.ndarray | None = None,
    curie_K: float = 0.0,
    neel_K: float = 0.0,
    weiss_K: float = 0.0,
    susceptibility_si: float = 0.0,
    saturation_emu_g: float = 0.0,
    applied_field_mT: float = 10.0,
    morin_K: float = 0.0,
    noise: float = 0.005,
    seed: int | None = 0,
) -> SquidMeasurement:
    """Synthesise an rf-SQUID χ(T) thermal sweep.

    Behaviour by ordering:

      - **paramagnetic** — Curie-Weiss χ(T) = C / (T − θ) with C calibrated so
        χ(295 K) ≈ susceptibility_si.
      - **ferromagnetic / ferrimagnetic** — Curie-Weiss above Tc, sharp cusp at
        Tc, decline below as 1/(T₀−T) (low-T saturation); peak at T=Tc.
      - **antiferromagnetic** — Curie-Weiss above TN, broad cusp at TN, decline
        below.
      - **canted-afm** — like AFM but with optional Morin transition: a step
        change at morin_K when supplied (signature of α-Fe2O3 single crystals).
      - **diamagnetic** — flat negative susceptibility (no T-dependence).
    """
    rng = np.random.default_rng(seed)
    if temperatures_K is None:
        T_max = max(2.0 * max(curie_K, neel_K, 295.0), 600.0)
        temperatures_K = np.linspace(2.0, T_max, 601)
    T = np.asarray(temperatures_K, dtype=float)

    if ordering == "paramagnetic":
        chi_295 = max(susceptibility_si, 1.0e-7)
        # Pick C so that chi(295 K) = chi_295 with given Weiss intercept.
        C = chi_295 * (295.0 - weiss_K)
        chi = C / np.maximum(T - weiss_K, 0.5)
    elif ordering == "ferromagnetic" or ordering == "ferrimagnetic":
        Tc = max(curie_K, 10.0)
        chi_295 = max(susceptibility_si, 1.0e-4)
        # Above Tc: Curie-Weiss with positive Weiss = Tc.
        C = chi_295 * (295.0 - Tc) if 295.0 > Tc else chi_295 * Tc
        above = T > Tc + 1.0
        below = T < Tc - 1.0
        chi = np.zeros_like(T)
        chi[above] = C / np.maximum(T[above] - Tc, 0.5)
        # Below Tc: increase from low-T saturation toward divergence at Tc.
        peak = chi[above].max() if np.any(above) else chi_295 * 10.0
        if np.any(below):
            chi[below] = peak * (T[below] / Tc) ** 2
        # Cusp band near Tc: average of the two formulas.
        cusp = ~(above | below)
        if np.any(cusp):
            chi[cusp] = peak
    elif ordering == "antiferromagnetic":
        TN = max(neel_K, 5.0)
        chi_295 = max(susceptibility_si, 1.0e-6)
        weiss = weiss_K if weiss_K != 0.0 else -TN
        C = chi_295 * (295.0 - weiss)
        above = T > TN + 1.0
        below = T < TN - 1.0
        chi = np.zeros_like(T)
        chi[above] = C / np.maximum(T[above] - weiss, 0.5)
        # Below TN: linear decline to 0.5 of cusp value at low T.
        cusp_value = C / max(TN - weiss, 0.5)
        if np.any(below):
            chi[below] = cusp_value * (0.5 + 0.5 * T[below] / TN)
        cusp = ~(above | below)
        if np.any(cusp):
            chi[cusp] = cusp_value
    elif ordering == "canted-afm":
        # Like AFM, but layered with the Morin transition step.
        TN = max(neel_K, 5.0)
        chi_295 = max(susceptibility_si, 1.0e-6)
        weiss = weiss_K if weiss_K != 0.0 else -TN
        C = chi_295 * (295.0 - weiss)
        above = T > TN + 1.0
        below = T < TN - 1.0
        chi = np.zeros_like(T)
        chi[above] = C / np.maximum(T[above] - weiss, 0.5)
        cusp_value = C / max(TN - weiss, 0.5)
        if np.any(below):
            chi[below] = cusp_value * (0.5 + 0.5 * T[below] / TN)
        cusp = ~(above | below)
        if np.any(cusp):
            chi[cusp] = cusp_value
        # Morin transition: spin-flop drops the in-plane component.
        if morin_K > 0.0:
            mask = T < morin_K
            chi[mask] = chi[mask] * 0.2
    elif ordering == "diamagnetic":
        chi_val = -abs(susceptibility_si if susceptibility_si != 0 else 1.5e-5)
        chi = np.full_like(T, chi_val)
    else:
        raise ValueError(f"unknown ordering {ordering!r}; expected one of {ORDERINGS}")

    if noise > 0:
        peak = float(np.max(np.abs(chi)) or 1.0)
        chi = chi + rng.normal(0.0, noise * peak, size=chi.size)

    return SquidMeasurement(
        mode="rf-chi-T",
        axis=T,
        signal=chi,
        temperature_K=0.0,
        applied_field_mT=applied_field_mT,
        frequency_Hz=0.0,
        metadata={
            "ordering": ordering,
            "curie_K_input": curie_K,
            "neel_K_input": neel_K,
            "weiss_K_input": weiss_K,
            "morin_K_input": morin_K,
            "noise": noise,
        },
    )


def simulate_chi_ac(
    ordering: Ordering,
    *,
    frequencies_Hz: np.ndarray | None = None,
    blocking_temp_K: float = 100.0,
    temperature_K: float = 295.0,
    chi_T: float = 1.0e-3,
    chi_S: float = 1.0e-5,
    tau_s: float | None = None,
    noise: float = 0.005,
    seed: int | None = 0,
) -> SquidMeasurement:
    """Synthesise an rf-SQUID AC susceptibility χ'(ω) + iχ''(ω) frequency sweep.

    Implements the Casimir-du Pré (Debye) relaxation model:

        χ(ω) = χ_S + (χ_T − χ_S) / (1 + iωτ)

    so that χ'(ω) is the in-phase response and χ''(ω) peaks at ω·τ = 1. The
    relaxation time τ follows an Arrhenius law τ(T) = τ₀·exp(blocking/T) when
    `tau_s` is left unset, mimicking superparamagnetic blocking. Setting
    `tau_s` directly gives a fixed-τ run.
    """
    rng = np.random.default_rng(seed)
    if frequencies_Hz is None:
        frequencies_Hz = np.logspace(-1.0, 5.0, 121)  # 0.1 Hz to 100 kHz
    f = np.asarray(frequencies_Hz, dtype=float)
    omega = 2.0 * np.pi * f

    if tau_s is None:
        tau_0 = 1.0e-9  # attempt frequency 1 GHz
        tau_s = tau_0 * np.exp(blocking_temp_K / max(temperature_K, 1.0))

    if ordering == "diamagnetic":
        chi_re = np.full_like(f, -abs(chi_S))
        chi_im = np.zeros_like(f)
    else:
        denom = 1.0 + (omega * tau_s) ** 2
        chi_re = chi_S + (chi_T - chi_S) / denom
        chi_im = (chi_T - chi_S) * (omega * tau_s) / denom

    if noise > 0:
        peak = float(max(np.max(np.abs(chi_re)), np.max(np.abs(chi_im))))
        # Use the actual signal scale — never a hard floor — so noise stays
        # proportional to the susceptibility (paramagnets can be 1e-6 SI or smaller).
        if peak > 0.0:
            chi_re = chi_re + rng.normal(0.0, noise * peak, size=chi_re.size)
            chi_im = chi_im + rng.normal(0.0, noise * peak, size=chi_im.size)

    # Pack as interleaved real+imag pairs for SquidMeasurement.
    packed = np.empty(chi_re.size + chi_im.size, dtype=float)
    packed[0::2] = chi_re
    packed[1::2] = chi_im

    return SquidMeasurement(
        mode="rf-chi-ac",
        axis=f,
        signal=packed,
        temperature_K=temperature_K,
        applied_field_mT=0.0,
        frequency_Hz=0.0,
        metadata={
            "ordering": ordering,
            "tau_s": float(tau_s),
            "chi_T": chi_T,
            "chi_S": chi_S,
            "blocking_temp_K_input": blocking_temp_K,
            "noise": noise,
        },
    )


# ---------------------------------------------------------------------------
# Parameter extraction
# ---------------------------------------------------------------------------


def extract_coercivity(measurement: SquidMeasurement) -> float:
    """Estimate Hc as the magnitude of the field at which M crosses zero (mT)."""
    if measurement.mode != "dc-mh":
        return 0.0
    fields = measurement.axis
    m = measurement.signal
    # Interpolate sign changes; report the median |H| at the crossings.
    sign = np.sign(m)
    crossings = np.where(np.diff(sign) != 0)[0]
    if crossings.size == 0:
        return 0.0
    h_at_cross = []
    for i in crossings:
        # Linear interpolate between i and i+1.
        f0, f1 = fields[i], fields[i + 1]
        m0, m1 = m[i], m[i + 1]
        if m1 == m0:
            continue
        frac = -m0 / (m1 - m0)
        h_at_cross.append(f0 + frac * (f1 - f0))
    if not h_at_cross:
        return 0.0
    return float(np.median(np.abs(h_at_cross)))


def extract_saturation(measurement: SquidMeasurement) -> float:
    """Estimate Ms as the average |M| across the upper-decile of |H| (emu/g)."""
    if measurement.mode != "dc-mh":
        return 0.0
    fields = measurement.axis
    m = measurement.signal
    h_max = float(np.max(np.abs(fields)))
    threshold = 0.9 * h_max
    mask = np.abs(fields) > threshold
    if not np.any(mask):
        return float(np.max(np.abs(m)))
    return float(np.mean(np.abs(m[mask])))


def extract_remanence(measurement: SquidMeasurement) -> float:
    """Estimate Mr as the average |M| at the H=0 branches of the loop (emu/g)."""
    if measurement.mode != "dc-mh":
        return 0.0
    fields = measurement.axis
    m = measurement.signal
    # Find the two sample indices closest to H=0 (one on each branch).
    abs_h = np.abs(fields)
    if abs_h.size < 4:
        return 0.0
    order = np.argsort(abs_h)
    nearest = order[:4]
    return float(np.mean(np.abs(m[nearest])))


def extract_curie_temperature(measurement: SquidMeasurement) -> float:
    """Locate Tc / TN as the temperature at which dχ/dT is most negative.

    For ferro/ferrimagnets χ(T) drops sharply *just above* Tc; for AFM the cusp
    is at TN. In both cases the steepest negative slope is co-located with the
    transition, so this works as a unified estimator.
    """
    if measurement.mode != "rf-chi-T":
        return 0.0
    T = measurement.axis
    chi = measurement.signal
    if T.size < 3:
        return 0.0
    # Smooth derivative via central differences. Susceptibility falls sharply
    # *above* an ordering transition, so the minimum derivative locates Tc.
    dchi_dT = np.gradient(chi, T)
    idx = int(np.argmin(dchi_dT))
    return float(T[idx])


def fit_curie_weiss(measurement: SquidMeasurement,
                    T_min_K: float | None = None) -> tuple[float, float]:
    """Fit χ(T) = C / (T − θ) above a chosen low-T cutoff.

    Returns (C, theta_K). When `T_min_K` is None, uses the upper half of the
    temperature range (which always lives above the ordering transition).
    """
    if measurement.mode != "rf-chi-T":
        return (0.0, 0.0)
    T = measurement.axis
    chi = measurement.signal
    if T_min_K is None:
        T_min_K = float(np.percentile(T, 60.0))
    mask = (T > T_min_K) & (chi > 0.0)
    if mask.sum() < 4:
        return (0.0, 0.0)
    Tm = T[mask]
    chim = chi[mask]
    inv_chi = 1.0 / chim
    # Linear fit: 1/chi = (1/C) * (T - theta) -> slope=1/C, intercept=-theta/C
    A = np.vstack([Tm, np.ones_like(Tm)]).T
    slope, intercept = np.linalg.lstsq(A, inv_chi, rcond=None)[0]
    if slope == 0.0:
        return (0.0, 0.0)
    C = 1.0 / slope
    theta = -intercept * C
    return (float(C), float(theta))


def extract_loss_peak(measurement: SquidMeasurement) -> float:
    """Locate the χ''(ω) peak frequency (Hz) — proxy for 1/τ relaxation rate."""
    if measurement.mode != "rf-chi-ac":
        return 0.0
    freqs = measurement.axis
    _re, im = measurement.chi_complex()
    if im.size == 0:
        return 0.0
    idx = int(np.argmax(im))
    return float(freqs[idx])


def infer_ordering(measurement: SquidMeasurement) -> Ordering:
    """Heuristic: classify the ordering type from a single measurement.

    Decision rules:

      - dc-mh: nonzero coercivity + large |Ms| → ferro/ferri (favour ferri unless
        Ms is metallic-large); negative slope only → diamagnetic; tiny linear
        positive slope → paramagnetic; small remanence + linear slope →
        canted-afm.
      - rf-chi-T: dχ/dT minimum above a meaningful threshold → ordered
        (ferro/ferri/AFM); flat negative → diamagnetic; smooth Curie-Weiss →
        paramagnetic.
      - rf-chi-ac: presence of a χ'' peak → relaxing (treat as superparamagnet
        which we map to paramagnetic for catalog matching); flat → use χ' sign.
    """
    if measurement.mode == "dc-mh":
        Hc = extract_coercivity(measurement)
        Ms = extract_saturation(measurement)
        Mr = extract_remanence(measurement)
        # Slope at high field (linear susceptibility component)
        fields = measurement.axis
        m = measurement.signal
        h_max = float(np.max(np.abs(fields)))
        outer = np.abs(fields) > 0.9 * h_max
        slope = float(np.polyfit(fields[outer], m[outer], 1)[0]) if outer.sum() > 1 else 0.0

        if Hc > 1.0 and Ms > 1.0 and Mr > 0.05 * Ms:
            return "ferromagnetic" if Ms > 100.0 else "ferrimagnetic"
        if Hc > 5.0 and Ms < 5.0 and Mr > 0.0:
            return "canted-afm"
        if slope < -1.0e-7:
            return "diamagnetic"
        if Ms < 0.05 and abs(slope) > 0.0:
            return "paramagnetic"
        # Fallback: small linear slope, classify by sign.
        return "paramagnetic" if slope > 0.0 else "diamagnetic"

    if measurement.mode == "rf-chi-T":
        chi = measurement.signal
        if np.all(chi < 0):
            return "diamagnetic"
        T = measurement.axis
        dchi_dT = np.gradient(chi, T)
        peak_neg = float(-np.min(dchi_dT))
        baseline_neg = float(np.median(np.abs(dchi_dT)))
        # Locate the temperature of steepest negative slope.
        idx = int(np.argmin(dchi_dT))
        T_min = float(T.min())
        T_max = float(T.max())
        T_at_drop = float(T[idx])

        sharp_cusp = peak_neg > 5.0 * max(baseline_neg, 1.0e-9)
        # If the drop sits at the low edge of the sweep, it's just the Curie-Weiss
        # divergence near θ — not a real ordering transition.
        on_low_edge = T_at_drop < T_min + 0.05 * (T_max - T_min)
        if not sharp_cusp or on_low_edge:
            return "paramagnetic"

        # Use the Curie-Weiss intercept sign to separate ferro/ferri (θ > 0)
        # from canted-AFM / AFM (θ < 0). Fit on the upper half of the sweep
        # where Curie-Weiss is reliably valid.
        _C, theta = fit_curie_weiss(measurement)
        if theta > 0.0:
            chi_at_max = float(np.max(chi))
            return "ferrimagnetic" if chi_at_max > 1.0e-3 else "ferromagnetic"
        # Negative θ → AFM family. A clean Morin step (chi drops to <40% of cusp
        # well below TN) flags canted-AFM (hematite-style); otherwise plain AFM.
        cusp_value = float(chi[idx])
        low_T_chi = float(np.mean(chi[T < 0.3 * T_at_drop])) if np.any(T < 0.3 * T_at_drop) else cusp_value
        if low_T_chi < 0.4 * cusp_value:
            return "canted-afm"
        return "antiferromagnetic"

    if measurement.mode == "rf-chi-ac":
        re, im = measurement.chi_complex()
        if np.mean(re) < 0:
            return "diamagnetic"
        if np.max(im) > 0.1 * np.max(np.abs(re)):
            return "paramagnetic"  # superparamagnetic relaxation
        return "paramagnetic"

    raise ValueError(f"unknown mode {measurement.mode!r}")


# ---------------------------------------------------------------------------
# Identification
# ---------------------------------------------------------------------------


def _score_candidate(extracted: dict, sig: MagneticMineral) -> SquidCandidate:
    obs_ord = extracted.get("ordering", "")
    ord_ok = (obs_ord == sig.ordering)
    Tc_obs = extracted.get("curie_K", 0.0) or extracted.get("neel_K", 0.0)
    Tc_sig = sig.curie_K or sig.neel_K
    tc_residual = abs(Tc_obs - Tc_sig) if Tc_obs and Tc_sig else 999.0
    Ms_obs = extracted.get("saturation_emu_g", 0.0)
    Ms_sig = sig.saturation_emu_g
    moment_residual = abs(Ms_obs - Ms_sig) if Ms_obs and Ms_sig else 999.0
    chi_obs = extracted.get("susceptibility_295_si", 0.0)
    chi_lo, chi_hi = sig.susceptibility_si
    if chi_lo == 0.0 and chi_hi == 0.0:
        chi_residual = 999.0
    elif chi_lo <= chi_obs <= chi_hi:
        chi_residual = 0.0
    else:
        target = 0.5 * (chi_lo + chi_hi)
        chi_residual = abs(chi_obs - target) / max(abs(target), 1.0e-9)

    # Combined score: ordering match dominates, then Tc, then moment.
    combined = 0.0
    if ord_ok:
        combined += 1.0
    if Tc_sig and tc_residual <= 0.05 * max(Tc_sig, 1.0):
        combined += 0.7
    elif Tc_sig and tc_residual <= 0.20 * max(Tc_sig, 1.0):
        combined += 0.3
    if Ms_sig and moment_residual <= 0.20 * max(Ms_sig, 1.0):
        combined += 0.5
    if chi_residual <= 1.0:
        combined += 0.2
    return SquidCandidate(
        name=sig.name,
        ordering_match=ord_ok,
        tc_residual_K=float(tc_residual),
        moment_residual_emu_g=float(moment_residual),
        susceptibility_residual=float(chi_residual),
        combined=float(combined),
    )


def analyze(
    measurement: SquidMeasurement,
    *,
    candidates: dict[str, MagneticMineral] | None = None,
    top: int = 5,
) -> SquidResult:
    """Identify the magnetic mineral signature behind a SQUID measurement.

    Extracts ordering type, Curie/Néel temperature, saturation moment, coercivity
    (where applicable), and 295-K susceptibility. Each `MagneticMineral` is then
    scored by ordering match + parameter residuals. Returns up to `top` ranked
    candidates and the extracted parameter dict.
    """
    if candidates is None:
        candidates = SIGNATURES

    ordering = infer_ordering(measurement)
    extracted: dict[str, float | str] = {"ordering": ordering}

    if measurement.mode == "dc-mh":
        extracted["coercivity_mT"] = extract_coercivity(measurement)
        extracted["saturation_emu_g"] = extract_saturation(measurement)
        extracted["remanence_emu_g"] = extract_remanence(measurement)
        # Linear-slope susceptibility from outer-decile fit.
        fields = measurement.axis
        m = measurement.signal
        h_max = float(np.max(np.abs(fields)))
        outer = np.abs(fields) > 0.9 * h_max
        if outer.sum() > 1:
            slope = float(np.polyfit(fields[outer], m[outer], 1)[0])
            # Convert mT-emu_per_g to dimensionless SI volume susceptibility (toy factor).
            extracted["susceptibility_295_si"] = slope / 0.05
        else:
            extracted["susceptibility_295_si"] = 0.0
    elif measurement.mode == "rf-chi-T":
        T = measurement.axis
        chi = measurement.signal
        tc = extract_curie_temperature(measurement)
        # Distinguish Curie from Néel by ordering classification.
        if ordering in ("ferromagnetic", "ferrimagnetic"):
            extracted["curie_K"] = tc
        elif ordering in ("antiferromagnetic", "canted-afm"):
            extracted["neel_K"] = tc
        C, theta = fit_curie_weiss(measurement)
        extracted["curie_constant"] = C
        extracted["weiss_K"] = theta
        # χ at 295 K via interpolation.
        if T.size > 0 and 295.0 >= float(T.min()) and 295.0 <= float(T.max()):
            extracted["susceptibility_295_si"] = float(np.interp(295.0, T, chi))
        else:
            extracted["susceptibility_295_si"] = float(chi[-1])
    elif measurement.mode == "rf-chi-ac":
        re, _im = measurement.chi_complex()
        extracted["susceptibility_295_si"] = float(np.mean(re))
        extracted["loss_peak_Hz"] = extract_loss_peak(measurement)

    scored = [_score_candidate(extracted, sig) for sig in candidates.values()]
    scored.sort(key=lambda c: c.combined, reverse=True)
    notes: list[str] = []
    if not scored:
        notes.append("no candidates supplied")
    elif scored[0].combined <= 0:
        notes.append("no candidate matches the observed ordering or parameters")
    return SquidResult(
        measurement=measurement,
        candidates=scored[: max(top, 1)],
        extracted=extracted,
        notes=notes,
    )


def analyze_dc(measurement: SquidMeasurement, **kwargs) -> SquidResult:
    """dc-SQUID-only entry point. Errors if the measurement isn't M(H)."""
    if measurement.mode != "dc-mh":
        raise ValueError(f"analyze_dc requires mode='dc-mh', got {measurement.mode!r}")
    return analyze(measurement, **kwargs)


def analyze_rf(measurement: SquidMeasurement, **kwargs) -> SquidResult:
    """rf-SQUID-only entry point. Errors if the measurement isn't χ(T) or χ_ac."""
    if measurement.mode not in ("rf-chi-T", "rf-chi-ac"):
        raise ValueError(
            f"analyze_rf requires χ(T) or χ_ac mode, got {measurement.mode!r}"
        )
    return analyze(measurement, **kwargs)


# ---------------------------------------------------------------------------
# Spectrum adapter (so existing CLI/IO/diagnose works unchanged)
# ---------------------------------------------------------------------------


def to_spectrum(measurement: SquidMeasurement) -> Spectrum:
    """Wrap a SquidMeasurement in the generic Spectrum container.

    Notes:
      * dc-mh hysteresis loops are not single-valued in H (forward and backward
        branches share field values), and `Spectrum.__post_init__` sorts on
        axis, which would scramble the loop. To preserve branch order we use
        a monotone sample-index axis (np.arange) and stash the actual H array
        in `metadata['squid_field_mT']`.
      * rf-chi-ac stores the imaginary part of the susceptibility in
        `metadata['imag_signal']`.
    """
    meta = dict(measurement.metadata)
    meta["squid_mode"] = measurement.mode
    meta["temperature_K"] = float(measurement.temperature_K)
    meta["applied_field_mT"] = float(measurement.applied_field_mT)
    meta["frequency_Hz"] = float(measurement.frequency_Hz)
    if measurement.mode == "dc-mh":
        meta["squid_field_mT"] = measurement.axis.copy()
        index = np.arange(measurement.axis.size, dtype=float)
        return Spectrum(index, measurement.signal.copy(), "squid-mh", "mT", meta)
    if measurement.mode == "rf-chi-T":
        return Spectrum(
            measurement.axis.copy(), measurement.signal.copy(), "squid-chi", "K", meta,
        )
    # rf-chi-ac: store real part as Spectrum.intensity, imag in metadata.
    re, im = measurement.chi_complex()
    meta["imag_signal"] = im
    return Spectrum(measurement.axis.copy(), re, "squid-chi", "Hz", meta)


def from_spectrum(spectrum: Spectrum) -> SquidMeasurement:
    """Inverse of `to_spectrum` — recover the SquidMeasurement from a Spectrum."""
    meta = dict(spectrum.metadata)
    mode = meta.pop("squid_mode", None)
    if mode is None:
        # Fall back to inferring from technique + units.
        if spectrum.technique == "squid-mh":
            mode = "dc-mh"
        elif spectrum.technique == "squid-chi" and spectrum.units == "K":
            mode = "rf-chi-T"
        elif spectrum.technique == "squid-chi" and spectrum.units == "Hz":
            mode = "rf-chi-ac"
        else:
            raise ValueError(
                f"cannot infer SQUID mode from technique={spectrum.technique!r} units={spectrum.units!r}"
            )
    T = float(meta.pop("temperature_K", 0.0))
    H = float(meta.pop("applied_field_mT", 0.0))
    f = float(meta.pop("frequency_Hz", 0.0))
    if mode == "dc-mh":
        # Recover the original branch-ordered field axis from metadata.
        if "squid_field_mT" in meta:
            axis = np.asarray(meta.pop("squid_field_mT"), dtype=float)
        else:
            axis = spectrum.axis.copy()
        return SquidMeasurement(
            mode=mode, axis=axis, signal=spectrum.intensity.copy(),
            temperature_K=T, applied_field_mT=H, frequency_Hz=f, metadata=meta,
        )
    if mode == "rf-chi-ac" and "imag_signal" in meta:
        im = np.asarray(meta.pop("imag_signal"), dtype=float)
        re = np.asarray(spectrum.intensity, dtype=float)
        packed = np.empty(re.size + im.size, dtype=float)
        packed[0::2] = re
        packed[1::2] = im
        return SquidMeasurement(
            mode=mode, axis=spectrum.axis.copy(), signal=packed,
            temperature_K=T, applied_field_mT=H, frequency_Hz=f, metadata=meta,
        )
    return SquidMeasurement(
        mode=mode, axis=spectrum.axis.copy(), signal=spectrum.intensity.copy(),
        temperature_K=T, applied_field_mT=H, frequency_Hz=f, metadata=meta,
    )


# ---------------------------------------------------------------------------
# Public re-exports
# ---------------------------------------------------------------------------


__all__ = [
    "Mode",
    "SquidMeasurement",
    "SquidCandidate",
    "SquidResult",
    "simulate_mh",
    "simulate_chi_T",
    "simulate_chi_ac",
    "extract_coercivity",
    "extract_saturation",
    "extract_remanence",
    "extract_curie_temperature",
    "fit_curie_weiss",
    "extract_loss_peak",
    "infer_ordering",
    "analyze",
    "analyze_dc",
    "analyze_rf",
    "to_spectrum",
    "from_spectrum",
]


def _orderings_check() -> Iterable[str]:
    return ORDERINGS
