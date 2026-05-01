from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace

import numpy as np

from checkmsg.laser import LaserConfig
from checkmsg.spectrum import Spectrum, Technique
from checkmsg.temperature import ROOM_K, fwhm_factor, phonon_shift, stokes_antistokes_ratio


@dataclass(frozen=True)
class PeakSpec:
    position: float
    intensity: float = 1.0
    sigma: float = 1.0  # gaussian width
    gamma: float = 0.5  # lorentzian width
    # Optional electronic-absorption band (nm). When the excitation laser tunes
    # near this band, this peak is resonance-enhanced. None disables resonance.
    resonance_band_nm: float | None = None
    resonance_sigma_nm: float = 30.0
    resonance_max: float = 30.0  # multiplicative cap at exact resonance


def voigt_pseudo(x: np.ndarray, center: float, sigma: float, gamma: float, amplitude: float) -> np.ndarray:
    """Pseudo-Voigt approximation (Thompson, Cox, Hastings 1987) — fast and stable."""
    fg = 2.355 * sigma
    fl = 2.0 * gamma
    f = (fg ** 5 + 2.69269 * fg ** 4 * fl + 2.42843 * fg ** 3 * fl ** 2
         + 4.47163 * fg ** 2 * fl ** 3 + 0.07842 * fg * fl ** 4 + fl ** 5) ** 0.2
    eta = 1.36603 * (fl / f) - 0.47719 * (fl / f) ** 2 + 0.11116 * (fl / f) ** 3
    g = np.exp(-4.0 * np.log(2.0) * ((x - center) / f) ** 2)
    lo = 1.0 / (1.0 + 4.0 * ((x - center) / f) ** 2)
    return amplitude * (eta * lo + (1.0 - eta) * g)


def _apply_temperature(peak: PeakSpec, T_K: float) -> PeakSpec:
    """Narrow widths and blue-shift positions to reflect sample temperature."""
    width_scale = fwhm_factor(T_K)
    shift = phonon_shift(peak.position, T_K)
    return replace(
        peak,
        position=peak.position + shift,
        sigma=peak.sigma * width_scale,
        gamma=peak.gamma * width_scale,
    )


def _apply_laser(peak: PeakSpec, laser: LaserConfig) -> PeakSpec:
    """Apply 1/lambda^4 scaling and resonance enhancement to a peak."""
    intensity = peak.intensity * laser.lambda4_scale
    if peak.resonance_band_nm is not None:
        intensity *= laser.resonance_enhancement(
            peak.resonance_band_nm,
            sigma_nm=peak.resonance_sigma_nm,
            max_enhancement=peak.resonance_max,
        )
    return replace(peak, intensity=intensity)


def _fluorescence_baseline(axis: np.ndarray, technique: Technique, laser: LaserConfig,
                           amplitude_at_max_factor: float) -> np.ndarray:
    """Smooth fluorescence/PL background that grows linearly across the Raman axis.

    Real instrument backgrounds have complex shape; for the synthetic case we use
    a soft monotonic ramp scaled by the laser's gem-fluorescence factor. Only
    applied to Raman spectra (cm-1 axis).
    """
    x = np.asarray(axis, dtype=float)
    if technique != "raman" or x.size == 0:
        return np.zeros_like(x)
    # Normalise axis to [0, 1] over its observed extent.
    span = x.max() - x.min() or 1.0
    t = (x - x.min()) / span
    # Visible lasers tend to ramp toward higher Raman shifts (red side of the
    # absolute spectrum). UV/NIR show flatter, weaker backgrounds.
    profile = 0.4 + 0.6 * t
    return amplitude_at_max_factor * laser.fluorescence_factor * profile


def generate(
    peaks: Sequence[PeakSpec],
    axis: np.ndarray,
    technique: Technique,
    units: str,
    noise: float = 0.01,
    baseline: np.ndarray | None = None,
    seed: int | None = 0,
    metadata: dict | None = None,
    laser_nm: float | None = None,
    temperature_K: float = ROOM_K,
    fluorescence_amplitude: float = 0.5,
    include_antistokes: bool = False,
) -> Spectrum:
    """Generate a noisy synthetic spectrum from a set of peaks on the given axis.

    When `laser_nm` is set (a wavelength in nm), the spectrum gets:
      - 1/lambda^4 scaling on every peak (relative to 532 nm),
      - per-peak resonance enhancement when `PeakSpec.resonance_band_nm` is set,
      - a smooth fluorescence baseline whose amplitude follows the laser's
        gem-fluorescence factor.

    `temperature_K` (default 295) narrows widths and blue-shifts positions, mimicking
    cooling. The applied laser/temperature are recorded in `metadata`.
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(axis, dtype=float)
    y = np.zeros_like(x)

    laser = LaserConfig(laser_nm) if laser_nm is not None else None

    for p in peaks:
        adjusted = _apply_temperature(p, temperature_K)
        if laser is not None and technique == "raman":
            adjusted = _apply_laser(adjusted, laser)
        y = y + voigt_pseudo(x, adjusted.position, adjusted.sigma, adjusted.gamma, adjusted.intensity)
        if include_antistokes and technique == "raman" and p.position > 0:
            ratio = stokes_antistokes_ratio(p.position, temperature_K)
            y = y + voigt_pseudo(
                x, -adjusted.position, adjusted.sigma, adjusted.gamma, adjusted.intensity * ratio,
            )

    if baseline is not None:
        y = y + np.asarray(baseline, dtype=float)
    if laser is not None:
        y = y + _fluorescence_baseline(x, technique, laser, fluorescence_amplitude)

    if noise > 0:
        y = y + rng.normal(0.0, noise * (np.max(y) if np.max(y) > 0 else 1.0), size=y.size)

    meta = dict(metadata or {})
    if laser is not None:
        meta.setdefault("laser_nm", laser.wavelength_nm)
    meta.setdefault("temperature_K", float(temperature_K))
    return Spectrum(x, y, technique, units, meta)


def linear_baseline(axis: np.ndarray, slope: float = 0.0, intercept: float = 0.0) -> np.ndarray:
    x = np.asarray(axis, dtype=float)
    return slope * (x - x[0]) + intercept


def generate_laicpms_run(
    sample_signals_cps: dict[str, float],
    *,
    blank_cps: dict[str, float] | None = None,
    sample_duration_s: float = 30.0,
    blank_duration_s: float = 10.0,
    transition_s: float = 2.0,
    dwell_ms: float = 10.0,
    noise_factor: float = 0.05,
    seed: int | None = 0,
    metadata: dict | None = None,
):
    """Generate a synthetic LA-ICP-MS time-resolved run.

    `sample_signals_cps` maps isotope keys (e.g. "Pb206") to mean cps during the
    sample window. `blank_cps` maps the same keys to the gas-blank baseline; if
    omitted, defaults to a small uniform baseline (5 cps) per channel.

    The simulated transient is: blank window (length blank_duration_s) -> short
    transition with a sigmoidal ramp -> sample window (length sample_duration_s)
    where the cps fluctuates around the prescribed mean with Gaussian noise of
    std `noise_factor * mean`.
    """
    from checkmsg.laicpms import IcpmsRun, IcpmsTransient, Isotope

    rng = np.random.default_rng(seed)
    blank_cps = dict(blank_cps or {})
    total_s = blank_duration_s + transition_s + sample_duration_s
    n = int(total_s * 1000.0 / dwell_ms)
    t = np.linspace(0.0, total_s, n)

    transients: dict[str, IcpmsTransient] = {}
    for key, sample_signal in sample_signals_cps.items():
        iso = Isotope.lookup(key)
        baseline = blank_cps.get(key, 5.0)
        # sigmoidal step from baseline to sample_signal centered on the start of the sample window
        ramp_center = blank_duration_s + 0.5 * transition_s
        ramp = 1.0 / (1.0 + np.exp(-(t - ramp_center) * 6.0 / max(transition_s, 0.1)))
        mean = baseline + (sample_signal - baseline) * ramp
        noise_std = np.maximum(noise_factor * np.maximum(mean, 1.0), 1.0)
        intensity = mean + rng.normal(0.0, noise_std)
        intensity = np.clip(intensity, 0.0, None)
        transients[key] = IcpmsTransient(isotope=iso, time_s=t, intensity_cps=intensity)

    meta = dict(metadata or {})
    meta.setdefault("dwell_time_ms", dwell_ms)
    meta.setdefault("laser_wavelength_nm", 193)
    return IcpmsRun(
        transients=transients,
        blank_window_s=(0.5, blank_duration_s - 0.5),
        sample_window_s=(blank_duration_s + transition_s + 0.5, total_s - 0.5),
        metadata=meta,
    )


def generate_laicpms_spectrum(
    sample_signals_cps: dict[str, float],
    *,
    blank_cps: dict[str, float] | None = None,
    sample_duration_s: float = 30.0,
    blank_duration_s: float = 10.0,
    seed: int | None = 0,
    metadata: dict | None = None,
):
    """Bulk-integrated `Spectrum` form. Equivalent to generate_laicpms_run().to_spectrum()."""
    run = generate_laicpms_run(
        sample_signals_cps,
        blank_cps=blank_cps,
        sample_duration_s=sample_duration_s,
        blank_duration_s=blank_duration_s,
        seed=seed,
        metadata=metadata,
    )
    return run.to_spectrum()


def amorphous_halo(axis: np.ndarray, center: float, width: float, amplitude: float) -> np.ndarray:
    """Broad gaussian envelope used to fake amorphous (glass) Raman background."""
    x = np.asarray(axis, dtype=float)
    return amplitude * np.exp(-((x - center) ** 2) / (2.0 * width ** 2))


def generate_epr(
    spin_system,
    fields_mT: np.ndarray,
    frequency_GHz: float,
    *,
    noise: float = 0.005,
    modulation_amplitude_mT: float = 0.2,
    microwave_power_mW: float = 1.0,
    cavity_Q: float = 5000.0,
    sample_mass_g: float | None = None,
    derivative: bool = True,
    orientations: tuple[int, int] = (15, 15),
    seed: int | None = 0,
    temperature_K: float = ROOM_K,
    spin_count: float = 1.0,
):
    """Synthesize a CW EPR derivative spectrum for a `SpinSystem`.

    Wraps `epr.simulate_field_sweep`, applies modulation broadening (linewidth augmented by
    H_pp = sqrt(H_pp^2 + H_mod^2 / 2) form), Gaussian noise, and saturation roll-off
    `intensity / sqrt(1 + (P/Psat)^2)`. Records all relevant metadata on the returned Spectrum.
    """
    from dataclasses import replace as _replace

    from checkmsg import epr

    if modulation_amplitude_mT > 0 and spin_system.linewidth_mT > 0:
        broadened_lw = (spin_system.linewidth_mT ** 2 + 0.5 * modulation_amplitude_mT ** 2) ** 0.5
        sys_for_sim = _replace(spin_system, linewidth_mT=broadened_lw)
    else:
        sys_for_sim = spin_system

    spec = epr.simulate_field_sweep(
        sys_for_sim, frequency_GHz=frequency_GHz, fields_mT=np.asarray(fields_mT, dtype=float),
        orientations=orientations, derivative=derivative, temperature_K=temperature_K,
    )
    y = np.asarray(spec.intensity, dtype=float)

    # Saturation roll-off: simple model relative to a 1 mW reference saturation.
    Psat = 1.0
    y = y / np.sqrt(1.0 + (microwave_power_mW / Psat) ** 2)

    # Spin-count scaling (linear in N_spins)
    y = y * spin_count * cavity_Q * np.sqrt(microwave_power_mW)

    if noise > 0:
        rng = np.random.default_rng(seed)
        peak = float(np.max(np.abs(y)) or 1.0)
        y = y + rng.normal(0.0, noise * peak, size=y.size)

    metadata = dict(spec.metadata)
    metadata.update({
        "modulation_amplitude_mT": modulation_amplitude_mT,
        "microwave_power_mW": microwave_power_mW,
        "cavity_Q": cavity_Q,
        "spin_count_input": spin_count,
    })
    if sample_mass_g is not None:
        metadata["sample_mass_g"] = sample_mass_g
    return Spectrum(spec.axis.copy(), y, "epr", spec.units, metadata)
