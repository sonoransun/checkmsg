"""Red gems carousel — beyond ruby.

Red colour comes from many distinct chromophores. Without lab analysis a
buyer could mistake any of the following for ruby:

  ruby           Cr³⁺ in corundum
  red_spinel     Cr³⁺ in spinel (the "Black Prince's Ruby" was a spinel)
  red_beryl      Mn³⁺ in beryl (extremely rare, Wah Wah Mountains)
  rubellite      Mn²⁺ in elbaite tourmaline
  pyrope         Mg-rich garnet
  almandine      Fe-rich garnet
  rhodochrosite  Mn carbonate

The pipeline picks the right host mineral via Raman, then confirms the chromophore
species via UV-VIS.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> int:
    args = parse_smoke_args("13_red_gems_carousel")
    specimens = ["ruby", "red_spinel", "red_beryl", "rubellite",
                 "pyrope", "almandine", "rhodochrosite"]
    rows = run_carousel("Scenario 13: red gems beyond ruby", specimens,
                        require_correct=len(specimens) - 2)
    if not args.smoke:
        plot_raman_carousel(rows, output_path("13_red_gems_carousel.png"),
                            "Red gems carousel — many chromophores, many hosts")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
