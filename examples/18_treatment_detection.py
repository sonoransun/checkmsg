"""Treatment-detection workflow — cascading techniques to catch enhanced gems.

Five treated stones present cosmetically as their natural counterparts but
carry diagnostic signatures of enhancement:

  heat-treated sapphire  silk dissolution; Fe/Ti redistribution
  glass-filled ruby      Pb-glass cavities; Pb XRF + Raman amorphous halo
  oiled emerald          oil/resin in fissures; CH stretches (3000 cm⁻¹)
  HPHT-treated diamond   catalyst metals (Fe/Co/Ni) at ppm levels via LA-ICP-MS
  Be-diffused sapphire   surface Be enrichment via LIBS depth profile

Pedagogy: treatments sit at the *intersection* of techniques. No single
spectrum reveals everything — the diagnostic comes from chains of evidence.
This example demonstrates how `diagnose()` reports a verdict for the *host*
mineral but uses follow-up recommendations to flag potential treatments.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args  # noqa: E402

from checkmsg import minerals  # noqa: E402
from checkmsg.diagnose import diagnose_profile  # noqa: E402


def main() -> int:
    args = parse_smoke_args("18_treatment_detection")
    print("=== Scenario 18: treatment-detection workflow ===\n")

    # All five treated stones share the host mineral identification path; the
    # treatment-specific diagnostic comes from ad-hoc inspection of evidence.
    treatment_cases = [
        ("Heat-treated blue sapphire", "sapphire_blue",
         "Heat treatment (1600+ °C) dissolves silk inclusions; XRF Fe/Ti unchanged. "
         "EPR Fe³⁺ centres survive but linewidth narrows."),
        ("Glass-filled ruby", "ruby",
         "Pb-glass cavity filling: high Pb in XRF + amorphous Raman halo over the cavity. "
         "Detectable via LA-ICP-MS Pb concentration spikes during depth profile."),
        ("Oiled emerald", "red_beryl",  # Stand-in for emerald (uses beryl Raman)
         "Cedar/synthetic oils show CH-stretch Raman bands ~2900-3000 cm⁻¹. "
         "GIA SSEF-style oiling clarity grade depends on fissure-fill volume."),
        ("HPHT-treated type IIa diamond", "diamond",
         "Catalyst residues (Fe + Co + Ni) at ppm levels via LA-ICP-MS distinguish "
         "synthetic-grown stones; HPHT-treated naturals show similar trace levels."),
        ("Be-diffused yellow sapphire", "sapphire_blue",  # corundum host
         "Be diffusion at high T pushes Be into the lattice up to ~10 µm depth. "
         "LA-ICP-MS depth profile is the only reliable detection method."),
    ]

    print(f"  {'specimen':<32} {'host verdict':<22} confidence")
    print("  " + "-" * 65)
    correct_host = 0
    for label, host_name, _ in treatment_cases:
        profile = minerals.get(host_name)
        report = diagnose_profile(profile)
        verdict = report.verdict or "?"
        match = verdict == host_name
        if match:
            correct_host += 1
        print(f"  {label:<32} {verdict:<22} {report.confidence:>5.2f}")
    print(f"\n  host identification accuracy: {correct_host}/{len(treatment_cases)}")
    assert correct_host >= len(treatment_cases) - 1

    print("\n  Treatment notes (follow-up techniques required for confirmation):")
    for label, _, note in treatment_cases:
        print(f"    {label}:")
        print(f"      {note}")

    if not args.smoke:
        _plot(treatment_cases, output_path("18_treatment_detection.png"))
    print("\nOK")
    return 0


def _plot(cases, path):
    import matplotlib.pyplot as plt

    from checkmsg import minerals as _m
    fig, ax = plt.subplots(figsize=(11, 6))
    for i, (label, host, _note) in enumerate(cases):
        profile = _m.get(host)
        spec = _m.synthesize_raman(profile, noise=0.005)
        ax.plot(spec.axis, spec.intensity + (len(cases) - 1 - i) * 1.05,
                linewidth=0.8, label=label)
    ax.set_xlabel("Raman shift (cm$^{-1}$)")
    ax.set_yticks([])
    ax.set_title("Treatment-detection workflow — host mineral Raman remains diagnostic")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
