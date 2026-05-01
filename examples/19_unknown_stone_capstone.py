"""Capstone — diagnosing an unknown specimen with the full pipeline.

A buyer presents a green stone of unknown identity. The lab measures it on
all four core techniques (Raman + UV-VIS + XRF + LIBS) and feeds the results
to the unified `diagnose()` pipeline.

This example shows the full reasoning trace — every piece of evidence, what
it favours and what it rules out — so a learner can read the chain end to
end. It is the worked example for the catalog + pipeline architecture.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args  # noqa: E402

from checkmsg import minerals  # noqa: E402
from checkmsg.diagnose import diagnose  # noqa: E402


def main() -> int:
    args = parse_smoke_args("19_unknown_stone_capstone")
    print("=== Scenario 19: capstone — diagnose an unknown green stone ===\n")
    print("A jeweler hands the lab a translucent green cabochon. ")
    print("It could be tsavorite, demantoid, peridot, jadeite, or emerald. ")
    print("We measure four techniques and feed everything to diagnose().\n")

    # The "unknown" is actually tsavorite (V³⁺-grossular).
    truth = "tsavorite"
    profile = minerals.get(truth)
    print(f"(ground truth, kept for didactic comparison: {truth})\n")

    spectra = []
    spectra.append(minerals.synthesize_raman(profile, noise=0.005))
    spectra.append(minerals.synthesize_uvvis(profile, noise=0.01))
    spectra.append(minerals.synthesize_xrf(profile, noise=0.001))
    spectra.append(minerals.synthesize_libs(profile, noise=0.003))

    report = diagnose(spectra, frequency_GHz=9.5)
    print(report.render())

    assert report.verdict == truth, f"expected {truth}, got {report.verdict}"
    print(f"\n[capstone] correctly identified as {report.verdict} "
          f"with confidence {report.confidence:.2f}")

    if not args.smoke:
        _plot(spectra, output_path("19_unknown_stone_capstone.png"))
    print("\nOK")
    return 0


def _plot(spectra, path):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(13, 8))
    by_tech = {s.technique: s for s in spectra}
    placement = {"raman": axs[0, 0], "uvvis": axs[0, 1], "xrf": axs[1, 0], "libs": axs[1, 1]}
    titles = {"raman": "Raman (cm$^{-1}$)", "uvvis": "UV-VIS absorbance (nm)",
              "xrf": "XRF (keV)", "libs": "LIBS (nm)"}
    for tech, ax in placement.items():
        if tech in by_tech:
            s = by_tech[tech]
            ax.plot(s.axis, s.intensity, linewidth=0.8, color="navy")
            ax.set_xlabel(titles[tech])
            ax.set_title(f"unknown stone — {tech}")
    fig.suptitle("Capstone: full-pipeline diagnosis of an unknown specimen", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
