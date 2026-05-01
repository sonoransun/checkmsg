"""LA-ICP-MS analysis: calibrated quantitation, isotope ratios, U-Pb dating, REE patterns.

This module turns a raw LA-ICP-MS acquisition into actionable quantities:

  - **Concentrations** in ppm via the Longerich et al. 1996 internal-standard equation,
    calibrated against bundled NIST SRM 612/610 reference glasses.
  - **Isotope ratios** for Pb (²⁰⁶/²⁰⁷/²⁰⁸/²⁰⁴) and Sr (⁸⁷/⁸⁶) with mass-bias correction.
  - **U-Pb concordant ages** from ²⁰⁶Pb/²³⁸U + ²⁰⁷Pb/²³⁵U using the Steiger & Jäger 1977
    decay constants.
  - **REE patterns** chondrite-normalised (McDonough & Sun 1995).
  - **Time-resolved transients** with gas-blank subtraction, sample-window integration,
    and simple change-point segmentation for depth-profile work.

The data primitive is `IcpmsRun` — a dict of `IcpmsTransient` channels (one per isotope)
plus blank/sample window definitions. The same module also accepts a bulk integrated
`Spectrum` (technique='laicpms', axis = isotope masses) for callers that already
processed their transient data externally.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass, field

import numpy as np

from checkmsg.refdata.icpms_data import (
    CHONDRITE_REE_PPM,
    ISOTOPES,
    LAMBDA_235,
    LAMBDA_238,
    NIST_REFERENCE_VALUES,
    REE_ELEMENTS,
    STACEY_KRAMERS_PB,
    IsotopeRecord,
)
from checkmsg.refdata.icpms_data import (
    isotope as _lookup_isotope,
)
from checkmsg.spectrum import Spectrum

# ---------- Data primitives ----------

@dataclass(frozen=True)
class Isotope:
    """A monitored mass channel (element + mass number + natural abundance)."""

    element: str
    mass: int
    natural_abundance: float = 1.0

    @property
    def key(self) -> str:
        return f"{self.element}{self.mass}"

    @classmethod
    def lookup(cls, key: str) -> Isotope:
        rec: IsotopeRecord = _lookup_isotope(key)
        return cls(element=rec.element, mass=rec.mass, natural_abundance=rec.natural_abundance)


@dataclass
class IcpmsTransient:
    """Counts-per-second time series for one isotope channel."""

    isotope: Isotope
    time_s: np.ndarray
    intensity_cps: np.ndarray

    def integrate(self, window: tuple[float, float]) -> float:
        """Mean cps over the given (start, end) seconds — robust against varying dwell rate."""
        t = np.asarray(self.time_s, dtype=float)
        y = np.asarray(self.intensity_cps, dtype=float)
        mask = (t >= window[0]) & (t <= window[1])
        if not np.any(mask):
            return 0.0
        return float(np.mean(y[mask]))


@dataclass
class IcpmsRun:
    """A complete LA-ICP-MS acquisition: per-channel transients plus integration windows."""

    transients: dict[str, IcpmsTransient]
    blank_window_s: tuple[float, float]
    sample_window_s: tuple[float, float]
    metadata: dict = field(default_factory=dict)

    def integrate(self, isotope_key: str, window: tuple[float, float] | None = None) -> float:
        if isotope_key not in self.transients:
            raise KeyError(f"no transient for {isotope_key}")
        return self.transients[isotope_key].integrate(window or self.sample_window_s)

    def blank_subtracted_cps(self) -> dict[str, float]:
        out: dict[str, float] = {}
        for key, tr in self.transients.items():
            sample = tr.integrate(self.sample_window_s)
            blank = tr.integrate(self.blank_window_s)
            out[key] = max(sample - blank, 0.0)
        return out

    def to_spectrum(self) -> Spectrum:
        """Bulk-integrated Spectrum: axis = isotope mass, intensity = blank-corrected cps."""
        cleaned = self.blank_subtracted_cps()
        # Sort by mass for a sensible axis ordering.
        items = sorted(cleaned.items(), key=lambda kv: self.transients[kv[0]].isotope.mass)
        masses = np.array([self.transients[k].isotope.mass for k, _ in items], dtype=float)
        intensities = np.array([v for _, v in items], dtype=float)
        meta = dict(self.metadata)
        meta.setdefault("integrated", True)
        meta["isotope_keys"] = [k for k, _ in items]
        return Spectrum(masses, intensities, "laicpms", "m/z", meta)

    def detect_segments(self, key: str, threshold_factor: float = 6.0,
                        min_gap_s: float = 1.0,
                        merge_fraction: float = 0.30) -> list[tuple[float, float]]:
        """Split the sample window into plateaus based on step changes in the chosen channel.

        Returns a list of (start, end) tuples in seconds. Designed for depth profiles
        where a coating and bulk material produce distinct plateaus.

        Algorithm:
          1. Find candidate change points where |dy/dt| > threshold_factor * MAD(dy).
          2. Enforce a minimum spacing of `min_gap_s` between change points.
          3. Merge adjacent segments whose means agree within `merge_fraction` of the
             larger value (kills noise-driven spurious splits within a single plateau).
        """
        if key not in self.transients:
            raise KeyError(f"no transient for {key}")
        tr = self.transients[key]
        t = np.asarray(tr.time_s, dtype=float)
        y = np.asarray(tr.intensity_cps, dtype=float)
        mask = (t >= self.sample_window_s[0]) & (t <= self.sample_window_s[1])
        if not np.any(mask):
            return []
        ts = t[mask]
        ys = y[mask]
        if len(ts) < 5:
            return [(ts[0], ts[-1])]
        dy = np.gradient(ys, ts)
        mad = float(np.median(np.abs(dy - np.median(dy)))) or 1e-9
        threshold = threshold_factor * mad
        change_idx = np.where(np.abs(dy) > threshold)[0]

        boundaries = [ts[0]]
        last_t = -math.inf
        for i in change_idx:
            ti = float(ts[i])
            if ti - last_t > min_gap_s:
                boundaries.append(ti)
                last_t = ti
        boundaries.append(ts[-1])
        boundaries = sorted(set(boundaries))

        segments = [(boundaries[i], boundaries[i + 1]) for i in range(len(boundaries) - 1)
                    if boundaries[i + 1] - boundaries[i] > min_gap_s / 2]
        if not segments:
            return [(ts[0], ts[-1])]

        # Merge adjacent segments with similar means.
        means = []
        for seg in segments:
            seg_mask = (ts >= seg[0]) & (ts <= seg[1])
            means.append(float(np.mean(ys[seg_mask])) if np.any(seg_mask) else 0.0)
        merged: list[tuple[float, float]] = [segments[0]]
        merged_means: list[float] = [means[0]]
        for seg, m in zip(segments[1:], means[1:], strict=True):
            prev_mean = merged_means[-1]
            scale = max(abs(prev_mean), abs(m), 1.0)
            if abs(m - prev_mean) / scale < merge_fraction:
                # merge into previous
                merged[-1] = (merged[-1][0], seg[1])
                merged_means[-1] = 0.5 * (prev_mean + m)
            else:
                merged.append(seg)
                merged_means.append(m)
        return merged


@dataclass(frozen=True)
class Concentration:
    element: str
    ppm: float
    detection_limit_ppm: float = 0.0


@dataclass
class IcpmsResult:
    concentrations: dict[str, Concentration]
    isotope_ratios: dict[str, float] = field(default_factory=dict)
    u_pb_age_Ma: float | None = None
    ree_pattern: dict[str, float] | None = None
    notes: list[str] = field(default_factory=list)

    def headline(self) -> str:
        parts = []
        if self.u_pb_age_Ma is not None:
            parts.append(f"U-Pb age = {self.u_pb_age_Ma:.0f} Ma")
        if "207/206" in self.isotope_ratios:
            parts.append(f"²⁰⁷Pb/²⁰⁶Pb = {self.isotope_ratios['207/206']:.4f}")
        if "87/86" in self.isotope_ratios:
            parts.append(f"⁸⁷Sr/⁸⁶Sr = {self.isotope_ratios['87/86']:.5f}")
        top = sorted(self.concentrations.values(), key=lambda c: c.ppm, reverse=True)[:5]
        if top:
            parts.append("top: " + ", ".join(f"{c.element}={c.ppm:.2f}ppm" for c in top))
        return "; ".join(parts) if parts else "no quantitation"


# ---------- Quantitation (Longerich 1996) ----------

def quantify(
    sample: IcpmsRun,
    calibration: IcpmsRun,
    *,
    internal_standard_element: str | None = None,
    internal_standard_ppm: float | None = None,
    reference: str = "NIST612",
) -> dict[str, Concentration]:
    """Convert a sample LA-ICP-MS run into concentrations (ppm).

    Implements the Longerich et al. 1996 quantitation equation:

        C_sample = (I_sample / I_cal) * (C_cal,ref) * (C_IS,sample / C_IS_apparent_sample)

    where the last bracket is the internal-standard correction. If `internal_standard_element`
    is None, the IS correction is skipped (sensitivity-only mode); accuracy degrades for
    matrices significantly different from the calibration glass.
    """
    if reference not in NIST_REFERENCE_VALUES:
        raise KeyError(f"unknown reference {reference!r}; choose from {list(NIST_REFERENCE_VALUES)}")
    ref_values = NIST_REFERENCE_VALUES[reference]

    sample_cps = sample.blank_subtracted_cps()
    cal_cps = calibration.blank_subtracted_cps()

    # First pass: apparent concentrations using sensitivity from the calibration alone.
    apparent: dict[str, float] = {}
    detection_limits: dict[str, float] = {}
    for key, sample_signal in sample_cps.items():
        tr = sample.transients[key]
        elem = tr.isotope.element
        ref_ppm = ref_values.get(elem)
        if ref_ppm is None or ref_ppm <= 0:
            continue
        cal_signal = cal_cps.get(key, 0.0)
        if cal_signal <= 0:
            continue
        sensitivity = cal_signal / ref_ppm   # cps per ppm
        apparent[elem] = sample_signal / sensitivity
        # 3 * sigma of blank divided by sensitivity = LOD (minimum detectable conc)
        blank_tr = sample.transients[key]
        blank_mask = ((blank_tr.time_s >= sample.blank_window_s[0])
                      & (blank_tr.time_s <= sample.blank_window_s[1]))
        if np.any(blank_mask):
            blank_sigma = float(np.std(blank_tr.intensity_cps[blank_mask]))
        else:
            blank_sigma = 0.0
        detection_limits[elem] = 3.0 * blank_sigma / sensitivity if sensitivity > 0 else 0.0

    # Optional internal-standard correction.
    if internal_standard_element is not None and internal_standard_ppm is not None:
        is_apparent = apparent.get(internal_standard_element)
        if is_apparent and is_apparent > 0:
            scale = internal_standard_ppm / is_apparent
            apparent = {el: c * scale for el, c in apparent.items()}
            detection_limits = {el: lod * scale for el, lod in detection_limits.items()}

    return {
        el: Concentration(element=el, ppm=ppm, detection_limit_ppm=detection_limits.get(el, 0.0))
        for el, ppm in apparent.items()
    }


# ---------- Isotope ratios ----------

def pb_ratios(run: IcpmsRun) -> dict[str, float]:
    """Standard Pb isotope ratios from a sample run (no mass-bias correction).

    Returns whichever ratios are computable given the channels present:
    206/204, 207/204, 208/204, 207/206, 208/206.
    """
    cps = run.blank_subtracted_cps()
    out: dict[str, float] = {}
    pb204 = cps.get("Pb204", 0.0)
    pb206 = cps.get("Pb206", 0.0)
    pb207 = cps.get("Pb207", 0.0)
    pb208 = cps.get("Pb208", 0.0)
    if pb204 > 0:
        if pb206:
            out["206/204"] = pb206 / pb204
        if pb207:
            out["207/204"] = pb207 / pb204
        if pb208:
            out["208/204"] = pb208 / pb204
    if pb206 > 0:
        if pb207:
            out["207/206"] = pb207 / pb206
        if pb208:
            out["208/206"] = pb208 / pb206
    return out


def sr_ratio(run: IcpmsRun, *, mass_bias_correct: bool = True) -> float:
    """⁸⁷Sr/⁸⁶Sr with internal mass-bias correction using ⁸⁶Sr/⁸⁸Sr = 0.1194.

    Returns NaN if Sr88 isn't measured (mass-bias correction disabled in that case
    when `mass_bias_correct=False`).
    """
    cps = run.blank_subtracted_cps()
    sr86 = cps.get("Sr86", 0.0)
    sr87 = cps.get("Sr87", 0.0)
    if sr86 <= 0 or sr87 <= 0:
        return float("nan")
    raw = sr87 / sr86
    if not mass_bias_correct:
        return raw
    sr88 = cps.get("Sr88", 0.0)
    if sr88 <= 0:
        return raw
    # Linear mass-bias law: m_corrected = m_raw * (1 + epsilon * (mass - reference_mass))
    # Solve epsilon from the known 86/88 = 0.1194 (Steiger & Jäger 1977).
    raw_86_over_88 = sr86 / sr88
    target_86_over_88 = 0.1194
    if raw_86_over_88 <= 0:
        return raw
    # 1 mass unit between 86 and 88: target/raw = (1 + 2*eps)
    epsilon = ((target_86_over_88 / raw_86_over_88) - 1.0) / 2.0
    # Apply same epsilon to 87/86 (1 mass unit difference)
    return raw * (1.0 + epsilon)


# ---------- U-Pb geochronology ----------

def _u_pb_age_from_ratio(ratio: float, lambda_decay: float) -> float:
    """Solve t in years from a daughter/parent ratio:  ratio = exp(lambda * t) - 1."""
    if ratio <= 0:
        return 0.0
    return math.log1p(ratio) / lambda_decay


def u_pb_age(
    run: IcpmsRun,
    *,
    common_pb_correction: bool = False,
    common_pb_206_204: float = STACEY_KRAMERS_PB["206/204"],
    common_pb_207_204: float = STACEY_KRAMERS_PB["207/204"],
    discordance_tolerance: float = 0.10,
) -> float:
    """Concordant U-Pb age in Ma.

    Computes both ²⁰⁶Pb/²³⁸U and ²⁰⁷Pb/²³⁵U ages and returns their mean if they agree
    within `discordance_tolerance`. Raises if the two diverge by more than that.
    """
    cps = run.blank_subtracted_cps()
    pb206 = cps.get("Pb206", 0.0)
    pb207 = cps.get("Pb207", 0.0)
    u238 = cps.get("U238", 0.0)
    u235 = cps.get("U235", 0.0)
    if u238 <= 0 or pb206 <= 0:
        raise ValueError("U-Pb age requires Pb206 and U238 channels")
    if u235 <= 0:
        # If U235 not directly measured, derive from U238 via natural ratio.
        from checkmsg.refdata.icpms_data import U238_OVER_U235
        u235 = u238 / U238_OVER_U235

    if common_pb_correction:
        pb204 = cps.get("Pb204", 0.0)
        if pb204 > 0:
            pb206 = max(pb206 - pb204 * common_pb_206_204, 0.0)
            if pb207 > 0:
                pb207 = max(pb207 - pb204 * common_pb_207_204, 0.0)

    # Sensitivity ratios: counts per atom for Pb vs U are roughly equal in static mode
    # (small mass discrimination). For our purposes we treat the cps ratio as the atom ratio.
    age_206_yr = _u_pb_age_from_ratio(pb206 / u238, LAMBDA_238)
    age_207_yr = _u_pb_age_from_ratio(pb207 / u235, LAMBDA_235) if pb207 > 0 else age_206_yr

    if age_206_yr > 0:
        rel_diff = abs(age_207_yr - age_206_yr) / age_206_yr
        if rel_diff > discordance_tolerance:
            raise ValueError(
                f"discordant U-Pb ages: 206/238 = {age_206_yr/1e6:.1f} Ma vs "
                f"207/235 = {age_207_yr/1e6:.1f} Ma (delta {rel_diff*100:.1f}%)"
            )

    return float(0.5 * (age_206_yr + age_207_yr) / 1e6)


# ---------- REE patterns ----------

def ree_pattern(concentrations: dict[str, Concentration]) -> dict[str, float]:
    """Chondrite-normalised REE pattern (CI chondrite, McDonough & Sun 1995).

    Returns the 14 REEs that have both a measured concentration and a chondrite reference.
    """
    out: dict[str, float] = {}
    for el in REE_ELEMENTS:
        c = concentrations.get(el)
        ref = CHONDRITE_REE_PPM.get(el)
        if c is None or ref is None or ref <= 0:
            continue
        out[el] = c.ppm / ref
    return out


def eu_anomaly(pattern: dict[str, float]) -> float:
    """Eu/Eu* anomaly: Eu_norm / sqrt(Sm_norm * Gd_norm). <1 = negative, >1 = positive."""
    sm, eu, gd = pattern.get("Sm"), pattern.get("Eu"), pattern.get("Gd")
    if sm is None or eu is None or gd is None or sm <= 0 or gd <= 0:
        return float("nan")
    return float(eu / math.sqrt(sm * gd))


def ce_anomaly(pattern: dict[str, float]) -> float:
    """Ce/Ce* anomaly: Ce_norm / sqrt(La_norm * Pr_norm)."""
    la, ce, pr = pattern.get("La"), pattern.get("Ce"), pattern.get("Pr")
    if la is None or ce is None or pr is None or la <= 0 or pr <= 0:
        return float("nan")
    return float(ce / math.sqrt(la * pr))


# ---------- High-level pipeline ----------

def analyze(
    sample: IcpmsRun,
    *,
    calibration: IcpmsRun | None = None,
    internal_standard: tuple[str, float] | None = None,
    reference: str = "NIST612",
) -> IcpmsResult:
    """Run the full LA-ICP-MS analysis pipeline on a sample run."""
    notes: list[str] = []

    if calibration is not None:
        is_el, is_ppm = (None, None) if internal_standard is None else internal_standard
        concentrations = quantify(
            sample, calibration,
            internal_standard_element=is_el, internal_standard_ppm=is_ppm,
            reference=reference,
        )
    else:
        concentrations = {}
        notes.append("No calibration spectrum provided — concentrations not computed.")

    isotope_ratios: dict[str, float] = {}
    isotope_ratios.update(pb_ratios(sample))
    sr = sr_ratio(sample)
    if not math.isnan(sr):
        isotope_ratios["87/86"] = sr

    age = None
    try:
        age = u_pb_age(sample, common_pb_correction=False)
    except ValueError:
        pass

    pattern = ree_pattern(concentrations) if concentrations else None
    if pattern:
        eu = eu_anomaly(pattern)
        ce = ce_anomaly(pattern)
        if not math.isnan(eu):
            notes.append(f"Eu/Eu* = {eu:.2f} ({'negative' if eu < 0.85 else 'no'} anomaly)")
        if not math.isnan(ce):
            notes.append(f"Ce/Ce* = {ce:.2f} ({'positive' if ce > 1.15 else 'no'} anomaly)")

    if age is not None:
        notes.append(f"U-Pb concordant age = {age:.0f} Ma")

    return IcpmsResult(
        concentrations=concentrations,
        isotope_ratios=isotope_ratios,
        u_pb_age_Ma=age,
        ree_pattern=pattern,
        notes=notes,
    )


# ---------- Convenience helpers ----------

def all_isotopes_for(elements: Iterable[str]) -> list[Isotope]:
    """Return all bundled isotopes for the given element list."""
    out: list[Isotope] = []
    for el in elements:
        for rec in ISOTOPES.values():
            if rec.element == el:
                out.append(Isotope(element=rec.element, mass=rec.mass,
                                   natural_abundance=rec.natural_abundance))
    return out


def run_from_spectrum(spec: Spectrum) -> IcpmsRun:
    """Wrap a bulk-integrated `Spectrum` (axis=mass, intensity=cps) as an `IcpmsRun`.

    The metadata key `isotope_keys` (a list aligned with the axis) lets the loader
    resolve every channel back to its element. If absent, the loader falls back to
    matching axis values against bundled `ISOTOPES` masses, picking the isotope of
    highest natural abundance for that mass when several elements share it.
    """
    if spec.technique != "laicpms":
        raise ValueError(f"expected laicpms spectrum, got {spec.technique}")
    keys = spec.metadata.get("isotope_keys")
    if keys is None:
        # Reverse-lookup masses against bundled ISOTOPES; pick highest abundance per mass.
        keys = []
        for m in spec.axis:
            mass_int = int(round(float(m)))
            best: tuple[str, float] | None = None
            for k, rec in ISOTOPES.items():
                if rec.mass == mass_int and (best is None or rec.natural_abundance > best[1]):
                    best = (k, rec.natural_abundance)
            keys.append(best[0] if best else f"M{mass_int}")
    transients: dict[str, IcpmsTransient] = {}
    t = np.array([0.0, 1.0])
    for key, cps in zip(keys, spec.intensity, strict=True):
        if key not in ISOTOPES:
            continue
        iso = Isotope.lookup(key)
        # Two identical samples (bulk-integrated values) at t=0 and t=1 — sample window covers both.
        transients[key] = IcpmsTransient(
            isotope=iso, time_s=t, intensity_cps=np.array([float(cps), float(cps)]),
        )
    return IcpmsRun(
        transients=transients,
        blank_window_s=(-1.0, -0.5),  # empty window -> blank subtraction returns sample as-is
        sample_window_s=(0.0, 1.0),
        metadata=dict(spec.metadata),
    )
