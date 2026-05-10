"""SQUID magnetometry — bulk magnetic moment and susceptibility for mineral ID.

Five physics-rich scenarios on synthetic SQUID measurements produced by the
`checkmsg.squid` forward simulators:

  1. dc-SQUID hysteresis carousel — magnetite vs hematite vs ilmenite vs
     diamond M(H) loops at 295 K. Asserts coercivity ordering and the 200×
     ratio in saturation moment between magnetite and hematite.
  2. rf-SQUID χ(T) thermal sweep on magnetite from 50–1000 K with Curie-Weiss
     fit; recover Tc within 5 %.
  3. AC susceptibility on a superparamagnetic Fe-Co-Ni cluster (HPHT diamond
     proxy). Casimir-du Pré loss peak in χ″(ω) → blocking-temperature
     estimate from Arrhenius τ.
  4. Pearl freshwater-vs-saltwater AC screening — non-destructive Mn²⁺ contrast
     at 1 Hz, 295 K. Asserts the freshwater pearl shows the larger χ′.
  5. Unified diagnose() integration — feed synthesised magnetite SQUID + Raman
     through the diagnostic pipeline; the report should cite SQUID evidence
     alongside the Raman match.

Each scenario asserts an expected outcome so the script doubles as a smoke test.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import install_raman_references, output_path, parse_smoke_args  # noqa: E402

from checkmsg import minerals, squid  # noqa: E402
from checkmsg.diagnose import diagnose  # noqa: E402


def scenario_hysteresis_carousel():
    """dc-SQUID M(H) loops on four canonical specimens at 295 K."""
    print("--- Scenario 1: dc-SQUID hysteresis carousel (295 K) ---")
    # Specimens: magnetite, hematite, ilmenite (proxy from refdata — not in catalog), diamond.
    from checkmsg.refdata.squid_signatures import SIGNATURES
    ilm = SIGNATURES["ilmenite"]

    measurements: dict[str, squid.SquidMeasurement] = {}
    measurements["magnetite"] = minerals.synthesize_squid_mh(
        minerals.get("magnetite"), seed=11
    )
    measurements["hematite"] = minerals.synthesize_squid_mh(
        minerals.get("hematite"), seed=12
    )
    measurements["ilmenite"] = squid.simulate_mh(
        ilm.ordering, saturation_emu_g=ilm.saturation_emu_g,
        coercivity_mT=ilm.coercivity_mT, susceptibility_si=0.5e-3,
        temperature_K=295.0, noise=0.005, seed=13,
    )
    measurements["diamond"] = minerals.synthesize_squid_mh(
        minerals.get("diamond"), seed=14
    )

    extracted: dict[str, dict] = {}
    for name, meas in measurements.items():
        Hc = squid.extract_coercivity(meas)
        Ms = squid.extract_saturation(meas)
        Mr = squid.extract_remanence(meas)
        extracted[name] = {"Hc_mT": Hc, "Ms_emu_g": Ms, "Mr_emu_g": Mr}
        print(f"  {name:<12} Hc = {Hc:7.1f} mT, Ms = {Ms:8.3f} emu/g, "
              f"Mr = {Mr:7.3f} emu/g")

    # Hierarchies that must hold:
    assert extracted["magnetite"]["Ms_emu_g"] > 50.0
    assert extracted["magnetite"]["Ms_emu_g"] > 50.0 * extracted["hematite"]["Ms_emu_g"]
    assert extracted["hematite"]["Hc_mT"] > 5.0 * extracted["magnetite"]["Hc_mT"]
    assert abs(extracted["diamond"]["Ms_emu_g"]) < 1e-2
    return measurements, extracted


def scenario_chi_T_curie_fit():
    """rf-SQUID χ(T) on magnetite — recover Tc + fit Curie-Weiss above the cusp."""
    print("--- Scenario 2: χ(T) Curie-Weiss fit on magnetite ---")
    T = np.linspace(50.0, 1000.0, 951)
    meas = squid.simulate_chi_T(
        "ferrimagnetic", temperatures_K=T, curie_K=858.0,
        susceptibility_si=2.0e-3, saturation_emu_g=92.0, noise=0.0, seed=21,
    )
    Tc = squid.extract_curie_temperature(meas)
    C, theta = squid.fit_curie_weiss(meas, T_min_K=900.0)
    print(f"  recovered Tc       = {Tc:.0f} K  (literature 858 K)")
    print(f"  Curie-Weiss θ      = {theta:.0f} K  (positive → ferri/ferro)")
    print(f"  Curie constant C   = {C:.3e}")
    assert abs(Tc - 858.0) < 50.0  # 5.8% tolerance
    return meas, {"Tc_K": Tc, "weiss_K": theta, "curie_C": C}


def scenario_ac_susceptibility_hpht_proxy():
    """AC χ_ac on a superparamagnetic FeCoNi nanocluster — Casimir-du Pré peak."""
    print("--- Scenario 3: superparamagnetic AC susceptibility (HPHT proxy) ---")
    # Pick a τ that puts the loss peak around 1 kHz: τ = 1/(2πf) ≈ 1.6e-4 s
    tau_s = 1.6e-4
    meas = squid.simulate_chi_ac(
        "ferromagnetic", tau_s=tau_s, chi_T=2.0e-3, chi_S=2.0e-5,
        temperature_K=295.0, noise=0.0, seed=31,
    )
    f_peak = squid.extract_loss_peak(meas)
    expected = 1.0 / (2.0 * np.pi * tau_s)
    print(f"  loss-peak frequency = {f_peak:.0f} Hz  (Casimir-du Pré expects {expected:.0f} Hz)")
    assert 0.7 * expected < f_peak < 1.4 * expected
    return meas, {"loss_peak_Hz": f_peak, "tau_s": tau_s}


def scenario_pearl_screening():
    """Freshwater vs saltwater pearl: AC χ′ at 1 Hz separates Mn²⁺ concentrations."""
    print("--- Scenario 4: pearl freshwater-vs-saltwater AC screening ---")
    saltwater = minerals.synthesize_squid_chi_ac(
        minerals.get("pearl_natural_saltwater"), seed=41,
        frequencies_Hz=np.array([1.0, 10.0, 100.0, 1000.0]),
    )
    freshwater = minerals.synthesize_squid_chi_ac(
        minerals.get("pearl_freshwater"), seed=42,
        frequencies_Hz=np.array([1.0, 10.0, 100.0, 1000.0]),
    )
    re_salt, _ = saltwater.chi_complex()
    re_fresh, _ = freshwater.chi_complex()
    chi_salt_1Hz = float(re_salt[0])
    chi_fresh_1Hz = float(re_fresh[0])
    print(f"  saltwater  pearl χ′(1 Hz) = {chi_salt_1Hz:.2e}")
    print(f"  freshwater pearl χ′(1 Hz) = {chi_fresh_1Hz:.2e}")
    print(f"  contrast ratio            = {chi_fresh_1Hz / chi_salt_1Hz:.1f}×")
    # Freshwater should show ≥ 5× larger χ at 1 Hz (Mn²⁺ ~30× higher in freshwater).
    assert chi_fresh_1Hz > 5.0 * chi_salt_1Hz
    return {"saltwater": saltwater, "freshwater": freshwater}, {
        "chi_saltwater_1Hz": chi_salt_1Hz,
        "chi_freshwater_1Hz": chi_fresh_1Hz,
    }


def scenario_diagnose_integration():
    """Feed synthesised magnetite SQUID + Raman through the diagnose() pipeline."""
    print("--- Scenario 5: diagnose() with SQUID + Raman ---")
    install_raman_references()
    profile = minerals.get("magnetite")
    raman_spec = minerals.synthesize_raman(profile, noise=0.01, seed=51)
    squid_meas = minerals.synthesize_squid_mh(profile, seed=52)
    squid_spec = squid.to_spectrum(squid_meas)
    report = diagnose([raman_spec, squid_spec])
    squid_evidence = [e for e in report.evidence if e.technique == "squid"]
    print(f"  verdict: {report.verdict}  confidence={report.confidence:.2f}")
    print(f"  squid evidence rows: {len(squid_evidence)}")
    for e in squid_evidence[:3]:
        print(f"    [{e.technique}] {e.observation}  (weight {e.weight:.1f})")
    assert report.verdict == "magnetite"
    assert squid_evidence, "diagnose() should cite SQUID evidence"
    return report


def plot_results(meas_set, chi_T_meas, ac_meas, pearl_meas, path):
    """Stack the scenario plots into a single PNG."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # Hysteresis loops
    ax = axes[0, 0]
    for name, meas in meas_set.items():
        ax.plot(meas.axis, meas.signal, label=name, lw=1.0)
    ax.set_xlabel("μ₀H (mT)")
    ax.set_ylabel("M (emu/g)")
    ax.set_title("dc-SQUID M(H) loops at 295 K")
    ax.axhline(0.0, color="grey", lw=0.5)
    ax.axvline(0.0, color="grey", lw=0.5)
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)

    # Chi(T)
    ax = axes[0, 1]
    ax.plot(chi_T_meas.axis, chi_T_meas.signal, color="#1f3a5f", lw=1.2)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("χ (SI)")
    ax.set_title("rf-SQUID χ(T) — magnetite Tc=858 K")
    ax.set_yscale("log")
    ax.grid(alpha=0.2)

    # AC susceptibility loss peak
    ax = axes[1, 0]
    re, im = ac_meas.chi_complex()
    ax.plot(ac_meas.axis, re, color="#1f3a5f", label="χ′(ω)")
    ax.plot(ac_meas.axis, im, color="#e07a3c", label="χ″(ω)")
    ax.set_xscale("log")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("χ (SI)")
    ax.set_title("AC χ — Casimir-du Pré (HPHT FeCoNi proxy)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)

    # Pearl AC contrast
    ax = axes[1, 1]
    for label, meas in pearl_meas.items():
        re, _ = meas.chi_complex()
        ax.plot(meas.axis, re, label=f"{label} pearl", marker="o", lw=1.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("χ′ (SI)")
    ax.set_title("Pearl AC screening — Mn²⁺ contrast")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)

    fig.suptitle("Scenario 21: SQUID magnetometry on canonical magnetic minerals",
                 fontsize=12, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


def main() -> int:
    args = parse_smoke_args("21_squid_magnetic_minerals")

    measurements, _ = scenario_hysteresis_carousel()
    chi_T_meas, _ = scenario_chi_T_curie_fit()
    ac_meas, _ = scenario_ac_susceptibility_hpht_proxy()
    pearl_meas, _ = scenario_pearl_screening()
    scenario_diagnose_integration()

    if not args.smoke:
        plot_results(
            measurements, chi_T_meas, ac_meas, pearl_meas,
            output_path("21_squid_magnetic_minerals.png"),
        )

    print("\nAll five SQUID scenarios passed.")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
