"""Biomineral disambiguation — pearls, coral and ivory.

Five biogenic gem materials with overlapping appearance:

  pearl_natural_saltwater  aragonite, low Mn, marine Pb signature
  pearl_akoya              aragonite + freshwater-mussel bead nucleus
  pearl_freshwater         aragonite + high Mn (>500 ppm) from feedwater
  coral                    aragonite/calcite + carotenoid pigment
  ivory                    hydroxyapatite Ca5(PO4)3OH (NOT aragonite)

Raman is the first cut: ivory's 962 cm⁻¹ phosphate ν1 mode is unmistakable
against the aragonite carbonate 1086 cm⁻¹. Pearls and coral are then
discriminated via EPR Mn²⁺ + LA-ICP-MS Pb isotopes.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> int:
    args = parse_smoke_args("17_biomineral_disambiguation")
    specimens = ["pearl_natural_saltwater", "pearl_akoya",
                 "pearl_freshwater", "coral", "ivory"]
    # Three pearl varieties share aragonite Raman; Raman alone can't separate
    # them. Without LA-ICP-MS or full EPR, pearls all score equally.
    rows = run_carousel("Scenario 17: biomineral disambiguation", specimens,
                        require_correct=2)
    if not args.smoke:
        plot_raman_carousel(rows, output_path("17_biomineral_disambiguation.png"),
                            "Biominerals — pearls vs coral vs ivory")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
