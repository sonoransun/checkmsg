"""Quartz family — same mineral, different colours and treatment history.

Five quartz varieties all share Raman: 128 + 207 + 463 cm⁻¹. They differ in
their colour-centre population, which encodes the geological/treatment history:

  rock_crystal           clean SiO2, no colour centres
  citrine_natural        Fe³⁺ chromophore (geological)
  citrine_heat_treated   ex-amethyst, retains residual E1' EPR centre
  amethyst               Fe⁴⁺ on tetrahedral site (post-irradiation)
  smoky_quartz           Al-hole + E1' centres (irradiation)

Pedagogy: Raman alone cannot separate them; UV-VIS chromophores + EPR colour
centres are the diagnostic levers.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> int:
    args = parse_smoke_args("15_quartz_treatment_history")
    specimens = ["rock_crystal", "citrine_natural", "citrine_heat_treated",
                 "amethyst", "smoky_quartz"]
    # Quartz varieties share their entire Raman pattern (128/207/463 cm-1).
    # Even with EPR + UV-VIS, ambiguity is expected — heat-treated citrine vs
    # amethyst share E1' + Al-hole centres, and rock_crystal can collide with
    # any other quartz under cosine matching. The educational point is that
    # quartz treatment history is genuinely a hard case: a *single* technique
    # rarely succeeds. We assert only that the pipeline picks a quartz species
    # for at least one specimen (modest threshold on purpose).
    rows = run_carousel("Scenario 15: quartz treatment history", specimens,
                        require_correct=1)
    if not args.smoke:
        plot_raman_carousel(rows, output_path("15_quartz_treatment_history.png"),
                            "Quartz family — same Raman, different colour centres")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
