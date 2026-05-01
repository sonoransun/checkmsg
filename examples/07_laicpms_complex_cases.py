"""LA-ICP-MS — resolving the cases other techniques cannot.

Four scenarios where Raman / XRF / UV-VIS / EPR alone leave room for ambiguity,
solved by a single LA-ICP-MS measurement that delivers ppm-level concentrations,
isotope ratios, U-Pb age, and (where relevant) a depth profile:

  1. Pearl: natural Persian Gulf vs Akoya cultured vs Chinese freshwater.
     Mn/Sr quantitation discriminates seawater from freshwater origin; Pb isotope
     ratios distinguish the bead-nucleus signature in cultured stones.
  2. HPHT-treated vs natural type IIa diamond. Type IIa is too clean for normal
     trace techniques; HPHT growth leaves Fe/Co/Ni catalyst residues at low ppm.
  3. Surface-coating depth profile. A treated stone shows a transient that drops
     from coating composition into the bulk after the laser drills through.
  4. Zircon U-Pb age + REE pattern. Returns a geological age plus the chondrite-
     normalised REE pattern that confirms natural igneous zircon.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args  # noqa: E402

from checkmsg import laicpms  # noqa: E402
from checkmsg.refdata.icpms_data import LAMBDA_235, LAMBDA_238, U238_OVER_U235  # noqa: E402
from checkmsg.synthetic import generate_laicpms_run  # noqa: E402


def _calibration_run(seed: int = 100):
    """Synthetic NIST 612 calibration: ~38 ppm trace plus matrix."""
    from checkmsg.refdata.icpms_data import NIST_SRM_612
    # Use a uniform sensitivity of 1000 cps/ppm.
    sens = 1000.0
    sample = {}
    for el, ppm in NIST_SRM_612.items():
        for key, rec in laicpms.ISOTOPES.items():
            if rec.element == el:
                sample[key] = ppm * rec.natural_abundance * sens
    return generate_laicpms_run(sample, sample_duration_s=20, blank_duration_s=5,
                                noise_factor=0.03, seed=seed)


def scenario_pearl():
    print("--- Scenario 1: pearl natural / akoya / freshwater ---")
    sens = 1000.0  # cps/ppm

    def pearl(true_ppm: dict[str, float], pb_ratios: dict[str, float], seed: int):
        # Compose channel cps from concentrations + isotopic-fraction Pb signal.
        sample = {}
        # Use Ca44 as internal standard (calcite matrix; total Ca ~400,000 ppm)
        sample["Ca44"] = 400000.0 * laicpms.ISOTOPES["Ca44"].natural_abundance * sens
        for el, ppm in true_ppm.items():
            for key, rec in laicpms.ISOTOPES.items():
                if rec.element == el and el != "Pb":
                    sample[key] = ppm * rec.natural_abundance * sens
        if "Pb" in true_ppm:
            pb_ppm = true_ppm["Pb"]
            r206_204 = pb_ratios["206/204"]
            r207_204 = pb_ratios["207/204"]
            r208_204 = pb_ratios["208/204"]
            atoms_total = 1.0
            atoms_204 = atoms_total / (1 + r206_204 + r207_204 + r208_204)
            sample["Pb204"] = pb_ppm * atoms_204 * sens
            sample["Pb206"] = pb_ppm * atoms_204 * r206_204 * sens
            sample["Pb207"] = pb_ppm * atoms_204 * r207_204 * sens
            sample["Pb208"] = pb_ppm * atoms_204 * r208_204 * sens
        # Use very low blank for Pb channels so Pb204 (~tens of cps) isn't clipped.
        blank_cps = {k: 0.5 for k in sample if k.startswith("Pb")}
        return generate_laicpms_run(sample, blank_cps=blank_cps,
                                    sample_duration_s=20, blank_duration_s=5,
                                    noise_factor=0.04, seed=seed)

    pearls = {
        # Natural saltwater: marine Pb ~18.7 (Stacey-Kramers terrestrial average)
        "Persian Gulf natural": pearl(
            {"Mn": 25, "Sr": 1500, "Pb": 5.0, "Mg": 100},
            pb_ratios={"206/204": 18.70, "207/204": 15.63, "208/204": 38.65},
            seed=11,
        ),
        # Akoya cultured: bead from Mississippi freshwater mussel shell — radically different
        # Pb composition from old (~1.5 Ga) North-American craton material.
        "Akoya cultured": pearl(
            {"Mn": 28, "Sr": 1450, "Pb": 5.0, "Mg": 110},
            pb_ratios={"206/204": 16.40, "207/204": 15.40, "208/204": 36.10},
            seed=12,
        ),
        "Chinese freshwater cultured": pearl(
            {"Mn": 750, "Sr": 350, "Pb": 5.0, "Mg": 200},
            pb_ratios={"206/204": 19.05, "207/204": 15.71, "208/204": 39.20},
            seed=13,
        ),
    }
    cal = _calibration_run(seed=110)

    verdicts: dict[str, str] = {}
    for label, run in pearls.items():
        result = laicpms.analyze(
            run, calibration=cal,
            internal_standard=("Ca", 400000.0),
            reference="NIST612",
        )
        mn = result.concentrations.get("Mn")
        sr = result.concentrations.get("Sr")
        r206_204 = result.isotope_ratios.get("206/204", 0.0)
        mn_ppm = mn.ppm if mn else 0
        sr_ppm = sr.ppm if sr else 0
        print(f"  {label}:")
        print(f"     Mn = {mn_ppm:.1f} ppm,  Sr = {sr_ppm:.1f} ppm,  ²⁰⁶Pb/²⁰⁴Pb = {r206_204:.3f}")
        if mn_ppm > 200:
            verdicts[label] = "freshwater cultured"
        elif r206_204 < 17.5:
            verdicts[label] = "akoya cultured (bead Pb signature)"
        else:
            verdicts[label] = "natural saltwater"

    print(f"  verdicts: {verdicts}")
    assert "freshwater" in verdicts["Chinese freshwater cultured"]
    assert "akoya" in verdicts["Akoya cultured"]
    assert "natural" in verdicts["Persian Gulf natural"]
    return pearls


def scenario_diamond():
    print("\n--- Scenario 2: type IIa diamond — natural vs HPHT-treated ---")
    sens = 1000.0
    cal = _calibration_run(seed=200)

    def diamond_run(trace_ppm: dict[str, float], seed: int):
        sample = {"Ca44": 1.0 * laicpms.ISOTOPES["Ca44"].natural_abundance * sens}  # tiny matrix
        for el, ppm in trace_ppm.items():
            for key, rec in laicpms.ISOTOPES.items():
                if rec.element == el:
                    sample[key] = ppm * rec.natural_abundance * sens
        return generate_laicpms_run(sample, sample_duration_s=15, blank_duration_s=5,
                                    noise_factor=0.05, seed=seed)

    diamonds = {
        "Untreated type IIa": diamond_run({"B": 0.001, "Si": 0.05, "Fe": 0.05}, seed=21),
        "HPHT-treated type IIa": diamond_run(
            {"B": 0.005, "Si": 0.10, "Fe": 8.0, "Co": 1.5, "Ni": 3.5}, seed=22),
        # Natural diamond with iron-bearing inclusion: Fe is elevated but no
        # Co/Ni cluster pattern — characteristic of magmatic, not synthetic origin.
        "Natural with Fe-inclusion": diamond_run(
            {"B": 0.002, "Si": 0.10, "Fe": 12.0, "Co": 0.05, "Ni": 0.10}, seed=23),
    }

    verdicts: dict[str, str] = {}
    for label, run in diamonds.items():
        result = laicpms.analyze(run, calibration=cal,
                                 internal_standard=None, reference="NIST612")
        fe = result.concentrations.get("Fe", laicpms.Concentration("Fe", 0.0)).ppm
        co = result.concentrations.get("Co", laicpms.Concentration("Co", 0.0)).ppm
        ni = result.concentrations.get("Ni", laicpms.Concentration("Ni", 0.0)).ppm
        print(f"  {label}:  Fe={fe:.3f} Co={co:.3f} Ni={ni:.3f}  (ppm)")
        # Discrimination: HPHT growth uses Fe-Ni-Co metal flux, leaving all three
        # at correlated ppm levels. Natural Fe inclusions (sulphide, magnetite)
        # bring Fe but never the Co+Ni catalyst cluster.
        if fe > 1.0 and co > 0.5 and ni > 0.5:
            verdicts[label] = "HPHT-treated"
        elif fe > 1.0:
            verdicts[label] = "natural with Fe-inclusion"
        else:
            verdicts[label] = "natural type IIa"

    print(f"  verdicts: {verdicts}")
    assert verdicts["HPHT-treated type IIa"] == "HPHT-treated"
    assert verdicts["Untreated type IIa"] == "natural type IIa"
    assert verdicts["Natural with Fe-inclusion"] == "natural with Fe-inclusion"
    return diamonds


def scenario_depth_profile():
    print("\n--- Scenario 3: depth profile of a coated/treated stone ---")
    sens = 1000.0

    # Build a transient with two segments inside the sample window:
    #   surface coating (high trace) for first 10 s, bulk corundum after.
    t = np.linspace(0.0, 50.0, 5000)

    def channel(key, blank, surface, bulk):
        iso = laicpms.Isotope.lookup(key)
        y = np.full_like(t, blank)
        in_surface = (t >= 12.0) & (t < 22.0)
        in_bulk = (t >= 22.0) & (t < 45.0)
        y[in_surface] = surface
        y[in_bulk] = bulk
        rng = np.random.default_rng(hash(key) & 0xFFFF)
        y = y + rng.normal(0.0, 0.05 * np.maximum(y, 1.0))
        y = np.clip(y, 0.0, None)
        return key, laicpms.IcpmsTransient(isotope=iso, time_s=t, intensity_cps=y)

    transients = dict([
        channel("Al27", 5, 1000.0 * sens, 530000.0 * laicpms.ISOTOPES["Al27"].natural_abundance * sens / 10),
        channel("Cr52", 5, 5.0 * laicpms.ISOTOPES["Cr52"].natural_abundance * sens,
                0.5 * laicpms.ISOTOPES["Cr52"].natural_abundance * sens),
        channel("Pb208", 5, 50.0 * laicpms.ISOTOPES["Pb208"].natural_abundance * sens,
                0.05 * laicpms.ISOTOPES["Pb208"].natural_abundance * sens),
    ])
    run = laicpms.IcpmsRun(
        transients=transients, blank_window_s=(0.5, 10.0),
        sample_window_s=(12.5, 44.5),
        metadata={"laser_wavelength_nm": 193},
    )
    segments = run.detect_segments("Pb208", threshold_factor=4.0, min_gap_s=1.0)
    print(f"  detected segments: {len(segments)} -> {[(round(a,1), round(b,1)) for a,b in segments]}")
    # Compute Pb cps in each segment to confirm coating vs bulk
    pb_per_segment = [run.integrate("Pb208", seg) for seg in segments]
    print(f"  Pb208 cps per segment: {[round(v) for v in pb_per_segment]}")
    assert len(segments) >= 2, "should detect at least surface + bulk segments"
    # The first segment should be the surface (highest Pb)
    assert max(pb_per_segment) == pb_per_segment[0], (
        "surface segment should have highest Pb208 signal")
    return run


def scenario_zircon():
    print("\n--- Scenario 4: zircon U-Pb age + REE pattern ---")
    sens = 1000.0
    cal = _calibration_run(seed=300)

    # Pick a Cretaceous age: 100 Myr.  Compute expected daughter/parent counts.
    age_yr = 100e6
    pb206_per_u238 = np.expm1(LAMBDA_238 * age_yr)   # ≈ 0.01563
    pb207_per_u235 = np.expm1(LAMBDA_235 * age_yr)   # ≈ 0.10350
    u238_atoms = 5e6  # arbitrary: defines signal magnitude
    u235_atoms = u238_atoms / U238_OVER_U235
    pb206_atoms = u238_atoms * pb206_per_u238
    pb207_atoms = u235_atoms * pb207_per_u235
    pb204_atoms = pb206_atoms * 0.001  # very low common Pb in zircon (clean)

    # Zircon REE pattern: HREE-enriched, negative Eu anomaly, positive Ce anomaly.
    # Typical chondrite-normalised values for igneous zircon
    ree_chondrite_norm = {
        "La": 0.005, "Ce": 50.0, "Pr": 0.10, "Nd": 0.50, "Sm": 1.0, "Eu": 0.3,
        "Gd": 4.0, "Tb": 8.0, "Dy": 30.0, "Ho": 60.0, "Er": 200.0, "Tm": 350.0,
        "Yb": 700.0, "Lu": 900.0,
    }
    from checkmsg.refdata.icpms_data import CHONDRITE_REE_PPM
    ree_ppm = {el: norm * CHONDRITE_REE_PPM[el] for el, norm in ree_chondrite_norm.items()}

    sample = {
        "U238": u238_atoms * sens,
        "U235": u235_atoms * sens,
        "Pb204": pb204_atoms * sens,
        "Pb206": pb206_atoms * sens,
        "Pb207": pb207_atoms * sens,
        "Pb208": 1.0 * sens,
        "Zr90": 480000.0 * laicpms.ISOTOPES["Zr90"].natural_abundance * sens,
    }
    for el, ppm in ree_ppm.items():
        for key, rec in laicpms.ISOTOPES.items():
            if rec.element == el:
                sample[key] = ppm * rec.natural_abundance * sens
    sample.setdefault("Ca44", 100.0 * laicpms.ISOTOPES["Ca44"].natural_abundance * sens)

    run = generate_laicpms_run(sample, sample_duration_s=30, blank_duration_s=10,
                               noise_factor=0.03, seed=41)
    result = laicpms.analyze(run, calibration=cal,
                             internal_standard=("Zr", 480000.0),
                             reference="NIST612")
    print(f"  recovered U-Pb age: {result.u_pb_age_Ma:.1f} Ma  (truth: 100 Ma)")
    print(f"  Eu/Eu* anomaly: {laicpms.eu_anomaly(result.ree_pattern):.3f}  (zircon expects <0.5)")
    print(f"  Ce/Ce* anomaly: {laicpms.ce_anomaly(result.ree_pattern):.3f}  (zircon expects >5)")
    pat = result.ree_pattern
    print(f"  HREE/LREE: La={pat['La']:.2f}  Yb={pat['Yb']:.1f}  Lu={pat['Lu']:.1f}")
    assert abs(result.u_pb_age_Ma - 100.0) < 10.0
    eu = laicpms.eu_anomaly(pat)
    assert eu < 0.6, f"expected negative Eu anomaly (Eu/Eu* < 0.6), got {eu:.2f}"
    assert pat["Yb"] > pat["La"]
    return run


def main() -> int:
    args = parse_smoke_args("07_laicpms_complex_cases")
    print("=== Scenario 7: LA-ICP-MS for ambiguous & complex test subjects ===\n")
    pearls = scenario_pearl()
    diamonds = scenario_diamond()
    depth = scenario_depth_profile()
    zircon = scenario_zircon()

    if not args.smoke:
        _plot(pearls, diamonds, depth, zircon, output_path("07_laicpms_complex_cases.png"))
    print("\nOK")
    return 0


def _plot(pearls, diamonds, depth, zircon, path):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(13, 9))

    # Pearl: bar chart of Mn ppm
    ax = axs[0, 0]
    cal = _calibration_run(seed=999)
    for label, run in pearls.items():
        cps = run.blank_subtracted_cps()
        ax.bar(label, cps.get("Mn55", 0), label=label)
    ax.set_title("Pearl Mn signal — freshwater stands out")
    ax.set_ylabel("Mn55 cps")
    ax.tick_params(axis="x", rotation=15)

    # Diamond: trace metals
    ax = axs[0, 1]
    elements = ["Fe56", "Co59", "Ni60", "B11"]
    width = 0.25
    x = np.arange(len(elements))
    for i, (label, run) in enumerate(diamonds.items()):
        cps = run.blank_subtracted_cps()
        ax.bar(x + (i - 1) * width, [cps.get(e, 0) for e in elements], width, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(elements)
    ax.set_title("Diamond trace metals — HPHT catalyst residues")
    ax.set_yscale("log")
    ax.legend(fontsize=8)

    # Depth profile transient
    ax = axs[1, 0]
    for key in ("Al27", "Cr52", "Pb208"):
        tr = depth.transients[key]
        ax.plot(tr.time_s, tr.intensity_cps / max(tr.intensity_cps.max(), 1), label=key, linewidth=0.8)
    ax.axvspan(*depth.blank_window_s, color="lightgray", alpha=0.3)
    for seg in depth.detect_segments("Pb208"):
        ax.axvspan(*seg, color="orange", alpha=0.1)
    ax.set_xlabel("time (s)")
    ax.set_title("Depth profile: surface coating (orange) vs bulk")
    ax.set_yscale("log")
    ax.legend(fontsize=8)

    # Zircon REE pattern
    ax = axs[1, 1]
    cal = _calibration_run(seed=998)
    res = laicpms.analyze(zircon, calibration=cal,
                          internal_standard=("Zr", 480000.0))
    pat = res.ree_pattern
    elements_ordered = [e for e in ("La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb",
                                    "Dy", "Ho", "Er", "Tm", "Yb", "Lu") if e in pat]
    values = [pat[e] for e in elements_ordered]
    ax.semilogy(elements_ordered, values, "o-", color="navy")
    ax.set_title(f"Zircon REE pattern (chondrite-normalised, U-Pb age = {res.u_pb_age_Ma:.0f} Ma)")
    ax.set_ylabel("sample / chondrite")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
