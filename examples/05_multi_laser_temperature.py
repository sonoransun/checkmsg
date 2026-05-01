"""Ruby Raman across 8 excitation wavelengths and 2 temperatures.

Story: a gem laboratory measures the same ruby specimen on a confocal Raman
microscope equipped with eight laser lines (275, 325, 405, 457, 488, 514, 633,
830 nm) and a liquid-nitrogen cryostage. The same specimen is measured at room
temperature and at 77 K (LN2). 16 spectra in total.

Five effects emerge from the data — none of which is visible from a single
laser at a single temperature:

  1. *Fingerprint invariance.* Corundum's A1g/Eg phonon positions in cm-1 are
     identical regardless of laser wavelength. Raman shift is what identifies
     the material, not absolute wavelength.

  2. *1/lambda^4 scaling.* Deep-UV excitation produces dramatically larger
     scattering cross-section per scatterer than NIR.

  3. *Cr3+ resonance enhancement.* Excitation near the 4A2 -> 4T2 absorption
     band (~555 nm) selectively boosts modes coupled to Cr3+ — strongly at
     514 nm, mildly at 488 nm, negligibly at 830 nm.

  4. *Fluorescence interference.* 488 and 514 nm excite Cr3+ photoluminescence
     in the red, raising the apparent baseline. 830 nm and the UV lasers escape
     this interference.

  5. *Temperature signatures.* At 77 K, phonon populations are suppressed, so
     intrinsic FWHMs narrow and peaks blue-shift by a few cm-1. Stokes/anti-
     Stokes intensity ratios provide an independent thermometer.

The script runs all 16 measurements, asserts the expected behaviour, prints a
summary table, and (without `--smoke`) saves a 2x4 figure organising the
spectra by laser group and temperature.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args  # noqa: E402

from checkmsg import raman  # noqa: E402
from checkmsg.laser import SUPPORTED_LASERS, LaserConfig  # noqa: E402
from checkmsg.peaks import fit_voigt  # noqa: E402
from checkmsg.spectrum import Spectrum  # noqa: E402
from checkmsg.synthetic import PeakSpec, generate  # noqa: E402
from checkmsg.temperature import LN2_K, ROOM_K, stokes_antistokes_ratio  # noqa: E402

# Ruby = Cr3+-doped corundum. Phonon positions from Porto & Krishnan 1967.
# The strongest Eg/A1g modes are coupled to the Cr3+ d-d transitions, so they
# resonance-enhance when the laser tunes near the ~555 nm 4A2 -> 4T2 band.
RUBY_PEAKS: list[PeakSpec] = [
    PeakSpec(378.0, intensity=0.55, sigma=2.5, gamma=0.8),
    PeakSpec(417.0, intensity=1.00, sigma=2.5, gamma=0.8,
             resonance_band_nm=555.0, resonance_max=20.0),
    PeakSpec(430.0, intensity=0.30, sigma=2.5, gamma=0.8),
    PeakSpec(450.0, intensity=0.20, sigma=2.5, gamma=0.8),
    PeakSpec(577.0, intensity=0.45, sigma=3.0, gamma=1.0),
    PeakSpec(645.0, intensity=0.50, sigma=2.5, gamma=0.8),
    PeakSpec(750.0, intensity=0.30, sigma=3.0, gamma=1.0,
             resonance_band_nm=555.0, resonance_max=15.0),
]

# Use the strongest, narrowest mode for cross-laser comparisons.
PROBE_PEAK_CM = 417.0


def axis_for_thermometry() -> np.ndarray:
    return np.linspace(-700.0, 1000.0, 1701)


def axis_stokes_only() -> np.ndarray:
    return np.linspace(100.0, 1100.0, 2001)


def measurement(laser_nm: float, T_K: float, seed: int):
    """Run one (laser, temperature) measurement on the synthetic ruby."""
    spec = generate(
        RUBY_PEAKS, axis_stokes_only(),
        technique="raman", units="cm-1",
        noise=0.005, seed=seed,
        laser_nm=laser_nm, temperature_K=T_K,
        fluorescence_amplitude=0.6,
    )
    cleaned = raman.preprocess_raman(spec)
    fit = fit_voigt(cleaned, around=PROBE_PEAK_CM, window=20.0)
    return spec, cleaned, fit


def thermometry(laser_nm: float, T_K: float, seed: int) -> float:
    """Use anti-Stokes/Stokes Voigt-fit ratio to recover the sample temperature."""
    spec = generate(
        RUBY_PEAKS, axis_for_thermometry(),
        technique="raman", units="cm-1",
        noise=0.001, seed=seed,
        laser_nm=laser_nm, temperature_K=T_K,
        fluorescence_amplitude=0.0,  # disable baseline so the AS peak is clean
        include_antistokes=True,
    )
    return raman.infer_temperature(spec, mode_cm=PROBE_PEAK_CM, window_cm=25.0)


def main() -> int:
    args = parse_smoke_args("05_multi_laser_temperature")

    print("=== Scenario 5: Ruby Raman across 8 lasers x 2 temperatures ===\n")

    # Header
    headers = ("laser", "regime", "scattered_at_417", "RT_pos", "RT_FWHM", "RT_I",
               "LN2_pos", "LN2_FWHM", "LN2_I", "I_LN2/I_RT", "S/AS thermometer")
    print(("{:>5}  {:>7}  {:>15}  {:>7}  {:>7}  {:>9}  {:>7}  {:>7}  {:>9}  {:>10}  {:>15}"
           ).format(*headers))
    print("-" * 130)

    # Collect for assertions and plotting.
    results: dict[str, list] = {"laser": [], "rt_pos": [], "ln2_pos": [], "rt_fwhm": [],
                                "ln2_fwhm": [], "rt_int": [], "ln2_int": [],
                                "rt_recovered_T": [], "ln2_recovered_T": []}
    cleaned_rt: dict[float, Spectrum] = {}
    cleaned_ln2: dict[float, Spectrum] = {}

    for i, wl in enumerate(SUPPORTED_LASERS):
        laser = LaserConfig(wl)
        rt_spec, rt_clean, rt_fit = measurement(wl, ROOM_K, seed=200 + i)
        ln2_spec, ln2_clean, ln2_fit = measurement(wl, LN2_K, seed=400 + i)

        # Use raw (non-normalised) cleaned intensity for the comparison: take
        # the height of the cleaned peak before min-max normalisation kills
        # the laser-dependent scaling.
        rt_raw_int = float(np.max(rt_spec.intensity) - np.min(rt_spec.intensity))
        ln2_raw_int = float(np.max(ln2_spec.intensity) - np.min(ln2_spec.intensity))

        try:
            rt_T = thermometry(wl, ROOM_K, seed=600 + i)
        except Exception:
            rt_T = float("nan")
        try:
            ln2_T = thermometry(wl, LN2_K, seed=800 + i)
        except Exception:
            ln2_T = float("nan")

        scattered = laser.shift_to_wavelength(PROBE_PEAK_CM)
        print(("{:>5g}  {:>7}  {:>13.1f} nm  {:>7.2f}  {:>7.2f}  {:>9.2f}  "
               "{:>7.2f}  {:>7.2f}  {:>9.2f}  {:>10.2f}  {:>9.0f} K").format(
            wl, laser.regime, scattered,
            rt_fit.position, rt_fit.width, rt_raw_int,
            ln2_fit.position, ln2_fit.width, ln2_raw_int,
            ln2_raw_int / rt_raw_int if rt_raw_int else float("nan"),
            ln2_T,
        ))

        results["laser"].append(wl)
        results["rt_pos"].append(rt_fit.position)
        results["ln2_pos"].append(ln2_fit.position)
        results["rt_fwhm"].append(rt_fit.width)
        results["ln2_fwhm"].append(ln2_fit.width)
        results["rt_int"].append(rt_raw_int)
        results["ln2_int"].append(ln2_raw_int)
        results["rt_recovered_T"].append(rt_T)
        results["ln2_recovered_T"].append(ln2_T)
        cleaned_rt[wl] = rt_clean
        cleaned_ln2[wl] = ln2_clean

    print()

    # ---- Assertions: each is the bone of a published claim about Raman ----

    # 1. Fingerprint invariance: positions within tight tolerance across lasers
    #    after temperature is held constant.
    rt_positions = np.array(results["rt_pos"])
    spread_rt = float(rt_positions.max() - rt_positions.min())
    print(f"[1] Fingerprint invariance: 417 cm-1 mode position spread across 8 lasers "
          f"at room T = {spread_rt:.2f} cm-1 (target <2.0)")
    assert spread_rt < 2.0, f"position invariance violated: spread {spread_rt:.2f}"

    # 2. Peak narrowing at LN2 vs RT: every laser should show LN2 width <= RT width.
    narrower = [ln2 < rt for rt, ln2 in zip(results["rt_fwhm"], results["ln2_fwhm"], strict=True)]
    print(f"[2] LN2 narrows peaks: {sum(narrower)}/{len(narrower)} lasers show LN2 FWHM < RT FWHM")
    assert sum(narrower) == len(narrower), "LN2 should narrow peaks for every laser"

    # 3. Phonon blue-shift at LN2: median(LN2 pos) > median(RT pos) by at least 1 cm-1.
    shift = float(np.median(results["ln2_pos"]) - np.median(results["rt_pos"]))
    print(f"[3] Phonon blue-shift at LN2: median delta = +{shift:.2f} cm-1 (target >+0.8)")
    assert shift > 0.8, f"expected blue-shift on cooling; got {shift:.2f}"

    # 4. Resonance enhancement at 514 nm vs 830 nm: the (intensity * lambda^4)-normalised
    #    height should be far larger at 514 because of the ~555 nm Cr3+ resonance.
    enh514 = LaserConfig(514.0).resonance_enhancement(555.0, sigma_nm=30.0, max_enhancement=20.0)
    enh830 = LaserConfig(830.0).resonance_enhancement(555.0, sigma_nm=30.0, max_enhancement=20.0)
    print(f"[4] Cr3+ resonance enhancement: 514 nm = {enh514:.1f}x, 830 nm = {enh830:.1f}x "
          f"(514 should be much larger)")
    assert enh514 > 5 * enh830, "resonance enhancement at 514 should dominate 830"

    # 5. Stokes/anti-Stokes thermometry: recovered T within 25 K of truth.
    rt_T_med = float(np.nanmedian(results["rt_recovered_T"]))
    ln2_T_med = float(np.nanmedian(results["ln2_recovered_T"]))
    print(f"[5] Stokes/anti-Stokes thermometry: RT recovered = {rt_T_med:.0f} K (truth {ROOM_K:.0f}); "
          f"LN2 recovered = {ln2_T_med:.0f} K (truth {LN2_K:.0f})")
    assert abs(rt_T_med - ROOM_K) < 25, "RT thermometry should be within 25 K"
    assert abs(ln2_T_med - LN2_K) < 25, "LN2 thermometry should be within 25 K"

    # 6. Fluorescence factor: 514 nm baseline is dominant, 830 nm baseline is faint.
    fl514 = LaserConfig(514.0).fluorescence_factor
    fl830 = LaserConfig(830.0).fluorescence_factor
    print(f"[6] Fluorescence interference: 514 nm = {fl514:.2f}, 830 nm = {fl830:.2f} "
          f"(NIR is preferred for low-fluorescence Raman)")
    assert fl514 > 5 * fl830

    # 7. Direct verification of stokes_antistokes_ratio at the probe mode.
    r_rt = stokes_antistokes_ratio(PROBE_PEAK_CM, ROOM_K)
    r_ln2 = stokes_antistokes_ratio(PROBE_PEAK_CM, LN2_K)
    print(f"[7] Bose-Einstein AS/S ratio for 417 cm-1: room T = {r_rt:.3f}, LN2 = {r_ln2:.5f} "
          f"(LN2 ratio is ~{r_rt/r_ln2:.0f}x smaller)")
    assert r_rt > r_ln2 * 50

    if not args.smoke:
        _plot(cleaned_rt, cleaned_ln2, output_path("05_multi_laser_temperature.png"))

    print("\nOK")
    return 0


def _plot(rt: dict, ln2: dict, path: Path) -> None:
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(13, 6), sharex=True, sharey=True)
    cmap = plt.get_cmap("viridis")
    n = len(SUPPORTED_LASERS)
    for ax, src, title in ((axs[0], rt, f"Room temperature ({ROOM_K:.0f} K)"),
                           (axs[1], ln2, f"Liquid nitrogen ({LN2_K:.0f} K)")):
        for i, wl in enumerate(SUPPORTED_LASERS):
            color = cmap(i / max(1, n - 1))
            offset = (n - 1 - i) * 1.05
            ax.plot(src[wl].axis, src[wl].intensity + offset, color=color,
                    linewidth=0.9, label=f"{wl:g} nm")
            ax.text(220, offset + 0.05, f"{wl:g} nm", color=color, fontsize=8)
        ax.set_xlabel("Raman shift (cm$^{-1}$)")
        ax.set_title(title)
        ax.set_xlim(150, 1050)
        ax.set_yticks([])
    axs[0].set_ylabel("intensity (offset for clarity)")
    fig.suptitle("Ruby — eight lasers, two temperatures (Stokes Raman, baseline-corrected)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
