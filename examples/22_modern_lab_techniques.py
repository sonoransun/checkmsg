"""Modern-lab techniques — PL, FTIR, Mössbauer, and CL for hard cases.

Five scenarios on synthetic spectra produced by the catalog synthesizers,
exercising the four techniques that real cutting-edge gem labs use:

  1. PL synthetic-diamond screen — a CVD diamond's Si-V (~737 nm) zero-phonon
     line is the near-definitive growth marker. Asserts the synthetic flag.
  2. FTIR diamond type — a natural diamond's aggregated-nitrogen bands classify
     it Type Ia. Asserts the type assignment.
  3. Mössbauer Fe valence — almandine garnet (Fe²⁺) vs blue sapphire (Fe³⁺);
     the isomer shift separates the oxidation states. Asserts both valences.
  4. CL activator bands — ruby's Cr³⁺ R-line red luminescence. Asserts the band.
  5. Unified diagnose() — natural vs CVD vs HPHT diamond. Raman is identical
     (1332.5), so PL/FTIR/EPR/SQUID carry the identification. Asserts that all
     three resolve to distinct verdicts.

Each scenario asserts an expected outcome so the script doubles as a smoke test.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args  # noqa: E402

from checkmsg import cl, ftir, minerals, mossbauer, pl  # noqa: E402
from checkmsg.diagnose import diagnose_profile  # noqa: E402


def scenario_pl():
    print("--- Scenario 1: PL synthetic-diamond screen ---")
    spec = minerals.synthesize_pl(minerals.get("diamond_cvd"), seed=1)
    res = pl.analyze(spec)
    print(f"  {res.headline()}")
    assert res.has_synthetic_marker(), "CVD diamond should show the Si-V synthetic marker"
    print("  ✓ Si-V centre flags CVD-synthetic growth")


def scenario_ftir():
    print("--- Scenario 2: FTIR diamond type ---")
    spec = minerals.synthesize_ftir(minerals.get("diamond"), seed=1)
    res = ftir.analyze(spec)
    print(f"  {res.headline()}")
    assert res.diamond_type == "Ia", f"expected Type Ia, got {res.diamond_type!r}"
    print("  ✓ aggregated-nitrogen bands → Type Ia")


def scenario_mossbauer():
    print("--- Scenario 3: Mössbauer Fe valence ---")
    alm = mossbauer.analyze(minerals.synthesize_mossbauer(minerals.get("almandine"), seed=1))
    sap = mossbauer.analyze(minerals.synthesize_mossbauer(minerals.get("sapphire_blue"), seed=1))
    print(f"  almandine: {alm.headline()}")
    print(f"  sapphire:  {sap.headline()}")
    assert alm.extracted["valence"] == "Fe2+", "almandine garnet is Fe2+"
    assert sap.extracted["valence"] == "Fe3+", "blue sapphire Fe is Fe3+"
    print("  ✓ isomer shift separates Fe²⁺ (almandine) from Fe³⁺ (sapphire)")


def scenario_cl():
    print("--- Scenario 4: CL activator bands ---")
    res = cl.analyze(minerals.synthesize_cl(minerals.get("ruby"), seed=1))
    print(f"  {res.headline()}")
    assert any("Cr3+" in b.name for b in res.emitters()), "ruby CL should show Cr3+ red"
    print("  ✓ Cr³⁺ R-line red luminescence")


def scenario_diagnose():
    print("--- Scenario 5: natural vs CVD vs HPHT diamond ---")
    verdicts = {}
    for name in ("diamond", "diamond_cvd", "diamond_hpht"):
        report = diagnose_profile(minerals.get(name))
        verdicts[name] = report.verdict
        techs = sorted({e.technique for e in report.evidence})
        print(f"  {name:14} -> {report.verdict:14} via {techs}")
    assert len(set(verdicts.values())) == 3, "natural / CVD / HPHT must be distinguished"
    print("  ✓ all three growth types resolve to distinct verdicts")


def main() -> None:
    args = parse_smoke_args("22_modern_lab_techniques")
    scenario_pl()
    scenario_ftir()
    scenario_mossbauer()
    scenario_cl()
    scenario_diagnose()

    if not args.smoke:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 2, figsize=(11, 7))
        specs = [
            ("PL — CVD diamond", minerals.synthesize_pl(minerals.get("diamond_cvd"), seed=1)),
            ("FTIR — natural diamond (Ia)", minerals.synthesize_ftir(minerals.get("diamond"), seed=1)),
            ("Mössbauer — almandine (Fe²⁺)",
             minerals.synthesize_mossbauer(minerals.get("almandine"), seed=1)),
            ("CL — ruby (Cr³⁺)", minerals.synthesize_cl(minerals.get("ruby"), seed=1)),
        ]
        for axis, (title, spec) in zip(axes.ravel(), specs, strict=True):
            axis.plot(spec.axis, spec.intensity, lw=0.8)
            axis.set_title(title, fontsize=10)
            axis.set_xlabel(spec.units)
        fig.suptitle("Modern lab techniques — PL / FTIR / Mössbauer / CL", weight="bold")
        fig.tight_layout()
        path = output_path("22_modern_lab_techniques.png")
        fig.savefig(path, dpi=130)
        print(f"saved plot to {path}")

    print("OK")


if __name__ == "__main__":
    main()
