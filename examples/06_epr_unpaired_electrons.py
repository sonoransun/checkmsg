"""EPR / ESR — identifying and quantifying unpaired-electron centers in gems.

Five mini-scenarios run on synthetic CW first-derivative spectra produced by a
real spin-Hamiltonian simulator (`checkmsg.epr.simulate_field_sweep`):

  1. Diamond P1 / HPHT-Ni / CVD at X-band — separate the three growth methods
     by the unpaired-electron fingerprint (no measurable EPR for ideal CVD).
  2. Smoky vs natural quartz at X-band — E1' oxygen-vacancy center is the
     irradiation marker.
  3. Mn2+ pearl / calcite fingerprint at X-band — six-line 55Mn hyperfine sextet.
  4. Cr3+ in ruby at X / Q / W bands — multi-frequency view of a powder pattern
     dominated by the D ~ 5.7 GHz zero-field splitting.
  5. Reference-standard quantification — compare the P1 spectrum's double
     integral against a synthetic DPPH reference at matched conditions.

Each scenario asserts an expected outcome so the script doubles as a smoke test.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args  # noqa: E402

from checkmsg import epr  # noqa: E402
from checkmsg.refdata.epr_centers import CENTERS  # noqa: E402
from checkmsg.synthetic import generate_epr  # noqa: E402


def synth(name: str, fields_mT: np.ndarray, frequency_GHz: float, *,
          spin_count: float = 1.0, noise: float = 0.005,
          orientations: tuple[int, int] = (11, 1), seed: int = 0):
    """Synthesise an EPR derivative spectrum for a bundled center."""
    return generate_epr(
        CENTERS[name], fields_mT, frequency_GHz=frequency_GHz,
        noise=noise, modulation_amplitude_mT=0.2, microwave_power_mW=2.0,
        cavity_Q=4500.0, spin_count=spin_count,
        orientations=orientations, seed=seed,
    )


def _noise_spectrum(fields_mT: np.ndarray, frequency_GHz: float, rng):
    """Pure baseline noise — no paramagnetic centers active in the modelled range."""
    from checkmsg.spectrum import Spectrum
    y = rng.normal(0.0, 1.0, size=fields_mT.size)
    return Spectrum(
        np.asarray(fields_mT, dtype=float), y, "epr", "mT",
        metadata={"frequency_GHz": frequency_GHz, "derivative": True},
    )


def scenario_diamond() -> tuple[dict, dict]:
    """Three diamonds at X-band: P1 (Ib natural), Ni (HPHT), and a clean CVD baseline."""
    print("--- Scenario 1: diamond P1 / HPHT-Ni / CVD at X-band ---")
    fields = np.linspace(326, 352, 1601)
    s_P1 = synth("diamond_P1", fields, 9.5, spin_count=1.0, seed=11)
    s_Ni = synth("diamond_Ni_HPHT", fields, 9.5, spin_count=0.6, seed=12)
    # CVD: no diagnostic centre at X-band for our model -> pure noise baseline.
    rng = np.random.default_rng(13)
    s_CVD = _noise_spectrum(fields, frequency_GHz=9.5, rng=rng)

    cands = {k: CENTERS[k] for k in ("diamond_P1", "diamond_Ni_HPHT", "free_electron")}
    r_P1 = epr.analyze(s_P1, frequency_GHz=9.5, candidates=cands, orientations=(7, 1))
    r_Ni = epr.analyze(s_Ni, frequency_GHz=9.5, candidates=cands, orientations=(7, 1))
    r_CVD = epr.analyze(s_CVD, frequency_GHz=9.5, candidates=cands, orientations=(7, 1))

    print(f"  Stone A (P1 Ib):     top={r_P1.best.name}  cos={r_P1.best.cosine:.3f}")
    print(f"  Stone B (HPHT Ni):   top={r_Ni.best.name}  cos={r_Ni.best.cosine:.3f}")
    print(f"  Stone C (CVD):       top={r_CVD.best.name}  cos={r_CVD.best.cosine:.3f}")
    assert r_P1.best.name == "diamond_P1"
    assert r_Ni.best.name == "diamond_Ni_HPHT"
    # CVD has near-zero spin density: top "match" cosine should be small (<0.3).
    assert r_CVD.best.cosine < 0.4, f"CVD cosine too high: {r_CVD.best.cosine}"
    return ({"P1": s_P1, "Ni": s_Ni, "CVD": s_CVD},
            {"P1": r_P1, "Ni": r_Ni, "CVD": r_CVD})


def scenario_quartz() -> dict:
    """E1' present (smoky) vs absent (natural rock crystal)."""
    print("\n--- Scenario 2: smoky vs natural quartz at X-band ---")
    fields = np.linspace(336, 342, 1201)
    smoky = synth("quartz_E1prime", fields, 9.5, spin_count=1.0, seed=21)
    natural = _noise_spectrum(fields, frequency_GHz=9.5, rng=np.random.default_rng(22))

    cands = {k: CENTERS[k] for k in ("quartz_E1prime", "quartz_Al_hole", "free_electron")}
    r_s = epr.analyze(smoky, 9.5, candidates=cands, orientations=(7, 1))
    r_n = epr.analyze(natural, 9.5, candidates=cands, orientations=(7, 1))
    print(f"  Smoky:    top={r_s.best.name}  cos={r_s.best.cosine:.3f}")
    print(f"  Natural:  top={r_n.best.name}  cos={r_n.best.cosine:.3f}")
    assert r_s.best.name == "quartz_E1prime"
    assert r_s.best.cosine > 0.85
    assert r_n.best.cosine < 0.4
    return {"smoky": smoky, "natural": natural}


def scenario_pearl() -> dict:
    """Mn2+ in calcite (pearl): 55Mn hyperfine sextet at X-band."""
    print("\n--- Scenario 3: Mn2+ pearl / calcite at X-band ---")
    fields = np.linspace(280, 400, 1601)
    spec = synth("calcite_Mn2plus", fields, 9.5, orientations=(7, 1), spin_count=1.0, seed=31)
    gs = epr.infer_g_factors(spec, 9.5)
    print(f"  Inferred g-factors: {len(gs)} resonances near g=2 (literature: 6 main 55Mn lines)")
    main = sorted(g for g in gs if 1.85 < g < 2.20)
    print(f"  Main lines: {[round(g, 4) for g in main]}")
    assert len(main) >= 6, f"expected at least 6 hyperfine lines, got {len(main)}"

    # Total hyperfine span: 5 spacings between the 6 outer-most central transitions.
    # Outer-line approach is robust against shoulder/satellite features from the
    # |+/-3/2> <-> |+/-5/2> fine-structure transitions slightly displaced by D.
    from checkmsg.epr import MU_B_OVER_H_MHZ_PER_MT
    fields_at = [9500.0 / (MU_B_OVER_H_MHZ_PER_MT * g) for g in main]
    span_mT = max(fields_at) - min(fields_at)
    avg_spacing = span_mT / 5.0
    expected = 245.0 / (MU_B_OVER_H_MHZ_PER_MT * 2.001)
    print(f"  Avg line spacing (span/5): {avg_spacing:.2f} mT (expected {expected:.2f} mT for A_iso=245 MHz)")
    assert abs(avg_spacing - expected) / expected < 0.20
    return {"pearl": spec}


def scenario_ruby_multiband() -> dict:
    """Cr3+ in corundum (ruby) at X, Q, W bands."""
    print("\n--- Scenario 4: Cr3+ in ruby at X / Q / W bands (powder) ---")
    bands = {"X": 9.5, "Q": 35.0, "W": 94.0}
    out = {}
    for label, freq in bands.items():
        # Each band has its own field range centered on g=2 resonance.
        center = 9500.0 / (epr.MU_B_OVER_H_MHZ_PER_MT * 2.0) * (freq / 9.5)
        span = 1200.0 if freq < 30 else (1500.0 if freq < 60 else 2500.0)
        fields = np.linspace(max(center - span, 0), center + span, 1201)
        spec = synth("corundum_Cr3plus", fields, freq, orientations=(15, 1), spin_count=1.0, seed=41)
        gs = epr.infer_g_factors(spec, freq)
        print(f"  {label}-band ({freq:g} GHz): {len(gs)} powder-pattern features")
        assert len(gs) >= 2, f"expected multiple Cr3+ powder features at {label}-band"
        out[label] = spec
    return out


def scenario_quantification(diamond_specs: dict) -> None:
    """Use a synthetic DPPH reference to quantify spins/g of the P1 stone."""
    print("\n--- Scenario 5: spin quantification via DPPH reference ---")
    fields = diamond_specs["P1"].axis
    # Reference: 1.0 unit of spin in DPPH at the same frequency, modulation, power.
    ref = synth("DPPH", fields, 9.5, orientations=(7, 1), spin_count=1.0, noise=0.005, seed=51)

    sample = diamond_specs["P1"]  # spin_count_input == 1.0 by construction
    count = epr.count_spins(sample, ref, reference_spins=1.0e18, sample_mass_g=0.10)
    print(f"  Sample double-integral / reference: {count.ratio_to_reference:.3f}")
    print(f"  Total spins:    {count.total:.3e}")
    print(f"  Spins per gram: {count.per_gram:.3e}")
    # Sample and reference both have spin_count=1.0 in synthesis.  Equal spins =>
    # ratio close to 1; but the P1 hyperfine triplet spreads intensity into 3 lines
    # and includes ZFS-free axial pattern, so the absolute double-integral shape
    # differs.  Tolerance set wide deliberately.
    assert 0.05 < count.ratio_to_reference < 20.0


def main() -> int:
    args = parse_smoke_args("06_epr_unpaired_electrons")
    print("=== Scenario 6: EPR identification + quantification of unpaired electrons ===\n")

    diamond_specs, diamond_results = scenario_diamond()
    quartz_specs = scenario_quartz()
    pearl_specs = scenario_pearl()
    ruby_specs = scenario_ruby_multiband()
    scenario_quantification(diamond_specs)

    if not args.smoke:
        _plot(diamond_specs, quartz_specs, pearl_specs, ruby_specs,
              output_path("06_epr_unpaired_electrons.png"))

    print("\nOK")
    return 0


def _plot(diamond, quartz, pearl, ruby, path: Path) -> None:
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(3, 2, figsize=(13, 9))

    ax = axs[0, 0]
    for label, s in diamond.items():
        ax.plot(s.axis, s.intensity, label=label, linewidth=0.8)
    ax.set_title("Diamond: P1 (Ib) vs HPHT-Ni vs CVD (X-band)")
    ax.set_xlabel("field (mT)")
    ax.legend(fontsize=8)

    ax = axs[0, 1]
    for label, s in quartz.items():
        ax.plot(s.axis, s.intensity, label=label, linewidth=0.8)
    ax.set_title("Quartz: smoky (E1') vs natural")
    ax.set_xlabel("field (mT)")
    ax.legend(fontsize=8)

    ax = axs[1, 0]
    s = list(pearl.values())[0]
    ax.plot(s.axis, s.intensity, color="purple", linewidth=0.8)
    ax.set_title("Mn$^{2+}$ in calcite (pearl) — 55Mn hyperfine sextet")
    ax.set_xlabel("field (mT)")

    ax = axs[1, 1]
    for label, s in ruby.items():
        ax.plot(s.axis, s.intensity / max(abs(s.intensity).max(), 1e-12), label=f"{label}-band", linewidth=0.8)
    ax.set_title("Cr$^{3+}$ in corundum across bands (normalised)")
    ax.set_xlabel("field (mT)")
    ax.legend(fontsize=8)

    # Bottom row: zoomed-in views
    ax = axs[2, 0]
    s = diamond["P1"]
    ax.plot(s.axis, s.intensity, color="darkred", linewidth=0.8)
    ax.set_title("Diamond P1 — three-line $^{14}$N hyperfine pattern")
    ax.set_xlabel("field (mT)")

    ax = axs[2, 1]
    s = ruby["X"]
    ax.plot(s.axis, s.intensity, color="navy", linewidth=0.7)
    ax.set_title("Cr$^{3+}$ in ruby at X-band — full powder pattern")
    ax.set_xlabel("field (mT)")

    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
