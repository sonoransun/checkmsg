"""Natural vs synthetic ruby — Raman alone cannot separate them; XRF/LIBS can.

Story: three corundum gems, all matching the corundum Raman fingerprint.
Discriminate origin via trace-element signatures:

  - Natural ruby (e.g. Mogok): Cr (color), Fe + Ti + V + Ga (host), Mg (charge balance).
  - Flame-fusion (Verneuil) synthetic: Cr only — feed powders are very pure;
    little to no Fe/Ti/V/Ga.
  - Flux-grown (Chatham, Kashan) synthetic: Cr + Pt or Mo flux residues.

Sources: Hughes 2017 (Ruby & Sapphire); Saminpanya & Sutherland 2011; GIA gem reference.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (  # noqa: E402
    CORUNDUM_RAMAN,
    install_raman_references,
    output_path,
    parse_smoke_args,
)

from checkmsg import libs as libs_mod  # noqa: E402
from checkmsg import raman  # noqa: E402
from checkmsg import xrf as xrf_mod  # noqa: E402
from checkmsg.synthetic import PeakSpec, generate  # noqa: E402

# Element abundances roughly mimicking observed XRF signal heights (relative units).
# Natural Mogok ruby: high Fe, moderate Ti and Ga, modest V, with Cr.
NATURAL_XRF = {"Al": 100.0, "Cr": 4.0, "Fe": 1.6, "Ti": 0.6, "V": 0.25, "Ga": 0.18, "Mg": 0.30}
# Verneuil flame-fusion: just Al + Cr; feed is ultra-pure aluminum oxide.
VERNEUIL_XRF = {"Al": 100.0, "Cr": 3.5}
# Flux-grown: Cr + flux residues. Pt flux from Kashan/Chatham, or PbO/MoO3.
FLUX_XRF = {"Al": 100.0, "Cr": 3.8, "Pt": 0.9, "Mo": 0.3}

# LIBS strong lines for the same elements (nm) — sourced from NIST ASD.
NATURAL_LIBS = [
    PeakSpec(396.152, intensity=1.00, sigma=0.10, gamma=0.05),  # Al
    PeakSpec(394.401, intensity=0.95, sigma=0.10, gamma=0.05),  # Al
    PeakSpec(425.435, intensity=0.55, sigma=0.10, gamma=0.05),  # Cr
    PeakSpec(427.480, intensity=0.50, sigma=0.10, gamma=0.05),  # Cr
    PeakSpec(371.994, intensity=0.40, sigma=0.10, gamma=0.05),  # Fe
    PeakSpec(334.941, intensity=0.20, sigma=0.10, gamma=0.05),  # Ti II
    PeakSpec(310.230, intensity=0.10, sigma=0.10, gamma=0.05),  # V II
    PeakSpec(287.424, intensity=0.08, sigma=0.10, gamma=0.05),  # Ga
    PeakSpec(279.553, intensity=0.07, sigma=0.10, gamma=0.05),  # Mg II
]
VERNEUIL_LIBS = [
    PeakSpec(396.152, intensity=1.00, sigma=0.10, gamma=0.05),
    PeakSpec(394.401, intensity=0.95, sigma=0.10, gamma=0.05),
    PeakSpec(425.435, intensity=0.50, sigma=0.10, gamma=0.05),
    PeakSpec(427.480, intensity=0.45, sigma=0.10, gamma=0.05),
]
FLUX_LIBS = [
    PeakSpec(396.152, intensity=1.00, sigma=0.10, gamma=0.05),
    PeakSpec(394.401, intensity=0.95, sigma=0.10, gamma=0.05),
    PeakSpec(425.435, intensity=0.55, sigma=0.10, gamma=0.05),
    PeakSpec(427.480, intensity=0.50, sigma=0.10, gamma=0.05),
    PeakSpec(265.945, intensity=0.45, sigma=0.10, gamma=0.05),  # Pt
    PeakSpec(379.825, intensity=0.20, sigma=0.10, gamma=0.05),  # Mo
]


def synth_xrf(abundances: dict[str, float], seed: int):
    from checkmsg.refdata.nist_xray import lines_for
    axis = np.linspace(0.5, 15.0, 4096)
    peaks: list[PeakSpec] = []
    for el, abund in abundances.items():
        for line in lines_for(el):
            peaks.append(PeakSpec(
                position=line.energy_keV,
                intensity=abund * line.relative_intensity,
                sigma=0.045,
                gamma=0.025,
            ))
    return generate(peaks, axis, technique="xrf", units="keV", noise=0.0015, seed=seed)


def synth_libs(peaks: list[PeakSpec], seed: int):
    axis = np.linspace(200.0, 800.0, 6000)
    return generate(peaks, axis, technique="libs", units="nm", noise=0.005, seed=seed)


def synth_raman(seed: int):
    axis = np.linspace(100.0, 1100.0, 1001)
    return generate(CORUNDUM_RAMAN, axis, technique="raman", units="cm-1", noise=0.01, seed=seed)


def classify(libs_result, xrf_result) -> tuple[str, list[str]]:
    notes: list[str] = []
    has_pt = libs_result.has("Pt") or any(e.element == "Pt" for e in xrf_result.elements)
    has_mo = libs_result.has("Mo") or any(e.element == "Mo" for e in xrf_result.elements)
    has_fe = libs_result.has("Fe") or any(e.element == "Fe" and e.confidence > 0.4 for e in xrf_result.elements)
    has_ti = libs_result.has("Ti") or any(e.element == "Ti" and e.confidence > 0.3 for e in xrf_result.elements)

    if has_pt or has_mo:
        notes.append("Pt and/or Mo detected — flux residues from a hydrothermal/flux growth process.")
        return "flux-grown synthetic", notes
    if not has_fe and not has_ti:
        notes.append("Cr present but no Fe/Ti/V/Ga: feedstock-pure aluminum oxide => flame-fusion (Verneuil).")
        return "Verneuil flame-fusion synthetic", notes
    notes.append("Trace Fe/Ti (and likely V/Ga) detected with Cr: natural geological host.")
    return "natural ruby", notes


def main() -> int:
    args = parse_smoke_args("02_natural_vs_synthetic_ruby")
    install_raman_references()

    samples = {
        "Stone 1 (Mogok-style natural)": (synth_raman(7), synth_xrf(NATURAL_XRF, 17), synth_libs(NATURAL_LIBS, 27)),
        "Stone 2 (Verneuil synthetic)": (synth_raman(8), synth_xrf(VERNEUIL_XRF, 18), synth_libs(VERNEUIL_LIBS, 28)),
        "Stone 3 (flux-grown synthetic)": (synth_raman(9), synth_xrf(FLUX_XRF, 19), synth_libs(FLUX_LIBS, 29)),
    }

    print("=== Scenario 2: Natural vs synthetic ruby (XRF + LIBS trace fingerprint) ===\n")
    verdicts: dict[str, str] = {}
    for label, (rspec, xspec, lspec) in samples.items():
        rresult = raman.analyze(rspec, candidates=["corundum", "ruby", "diamond", "beryl", "moissanite"])
        xresult = xrf_mod.identify_elements(xspec, tolerance_keV=0.05, min_snr=8.0)
        lresult = libs_mod.identify(lspec, tolerance_nm=0.4)

        print(f"{label}:")
        print(f"  Raman top: {rresult.best.mineral} (cosine={rresult.best.cosine:.2f})  "
              f"=> host is corundum")
        x_quant = xrf_mod.relative_quant(xresult)
        x_summary = ", ".join(
            f"{e.element}={x_quant.get(e.element, 0):.3f}"
            for e in xresult.elements[:8]
        )
        print(f"  XRF: {x_summary}")
        libs_summary = ", ".join(sorted(lresult.elements.keys()))
        print(f"  LIBS lines for: {libs_summary}")

        verdict, notes = classify(lresult, xresult)
        for n in notes:
            print(f"  - {n}")
        print(f"  => verdict: {verdict}\n")
        verdicts[label] = verdict

    assert verdicts["Stone 1 (Mogok-style natural)"] == "natural ruby"
    assert verdicts["Stone 2 (Verneuil synthetic)"] == "Verneuil flame-fusion synthetic"
    assert verdicts["Stone 3 (flux-grown synthetic)"] == "flux-grown synthetic"

    if not args.smoke:
        _plot(samples, output_path("02_natural_vs_synthetic_ruby.png"))
    print("OK")
    return 0


def _plot(samples: dict, path: Path) -> None:
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=False)
    for label, (_r, x_spec, _l) in samples.items():
        axs[0].plot(x_spec.axis, x_spec.intensity, label=label, linewidth=0.7)
    axs[0].set_xlabel("XRF energy (keV)")
    axs[0].set_xlim(0, 15)
    axs[0].legend(fontsize=8)
    axs[0].set_title("XRF — full ruby comparison")

    for label, (_r, _x, libs_spec) in samples.items():
        axs[1].plot(libs_spec.axis, libs_spec.intensity, label=label, linewidth=0.6)
    axs[1].set_xlabel("LIBS wavelength (nm)")
    axs[1].set_title("LIBS — trace fingerprint")
    axs[1].legend(fontsize=8)

    for label, (raman_spec, _x, _l) in samples.items():
        axs[2].plot(raman_spec.axis, raman_spec.intensity, label=label, linewidth=0.7)
    axs[2].set_xlabel("Raman shift (cm$^{-1}$)")
    axs[2].set_title("Raman — common corundum host")
    axs[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
