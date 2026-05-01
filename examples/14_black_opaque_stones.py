"""Black opaque stones — looking the same, behaving very differently.

Six black gems span four distinct chemical systems:

  hematite   alpha-Fe2O3 (iron oxide, blood-red streak)
  magnetite  Fe3O4 (magnetic, lower hardness)
  obsidian   amorphous SiO2 + impurities (volcanic glass)
  jet        carbonised lignite (very low density, organic)
  onyx       SiO2 chalcedony (banded quartz)
  schorl     Fe-tourmaline (high hardness, B in LIBS)

Pedagogy: amorphous vs crystalline screening, plus very different
densities (jet ~ 1.3 g/cc vs hematite ~ 5.3) for context.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> int:
    args = parse_smoke_args("14_black_opaque_stones")
    specimens = ["hematite", "magnetite", "obsidian", "jet", "onyx", "schorl"]
    rows = run_carousel("Scenario 14: black opaque stones", specimens,
                        require_correct=len(specimens) - 2)
    if not args.smoke:
        plot_raman_carousel(rows, output_path("14_black_opaque_stones.png"),
                            "Black gems — six minerals, six distinct Raman patterns")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
