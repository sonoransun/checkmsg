"""Blue-stone disambiguation — six chemically unrelated gems sharing a colour.

Sapphire, tanzanite, iolite, aquamarine, blue topaz and blue zircon all look
similarly blue to the unaided eye but belong to wholly different mineral
species — corundum, zoisite, cordierite, beryl, topaz and zircon respectively.

The combined Raman + UV-VIS pipeline resolves all six. The chromophore tells
us **why** the gem is blue (Fe-Ti IVCT in sapphire vs V³⁺ in tanzanite vs
Fe²⁺ in cordierite vs irradiation centre in topaz, etc.).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> int:
    args = parse_smoke_args("09_blue_stone_disambiguation")
    specimens = ["sapphire_blue", "tanzanite", "iolite",
                 "aquamarine", "blue_topaz", "blue_zircon"]
    rows = run_carousel("Scenario 9: blue-stone disambiguation", specimens,
                        require_correct=len(specimens) - 1)
    if not args.smoke:
        plot_raman_carousel(rows, output_path("09_blue_stone_disambiguation.png"),
                            "Six blue gems — distinct Raman fingerprints")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
