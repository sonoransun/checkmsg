"""Diamond simulant carousel — separating eight clear-stone "diamonds".

A jeweler offers eight colourless brilliants claimed to be diamond. Only one is.
Here every supposed diamond is run through the unified `diagnose()` pipeline,
which returns the canonical mineral name with a reasoning trace.

Specimens: diamond, moissanite, cubic_zirconia, GGG, YAG, white_sapphire,
white_topaz, white_spinel, glass_paste.

Pedagogy:
  - Raman is the workhorse — every mineral has a unique fingerprint
  - amorphous silicate (glass) is flagged by lack of sharp peaks
  - density (RI) data printed for context but not used in identification

Note: the spectra are synthesised from the bundled mineral catalog for
reproducible, offline didactic use. Real instrument data would feed in equally.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> int:
    args = parse_smoke_args("08_diamond_simulant_carousel")
    specimens = [
        "diamond", "moissanite", "cubic_zirconia", "GGG", "YAG",
        "white_sapphire", "white_topaz", "white_spinel", "glass_paste",
    ]
    rows = run_carousel("Scenario 8: diamond simulant carousel", specimens,
                        require_correct=len(specimens) - 1)
    if not args.smoke:
        plot_raman_carousel(rows, output_path("08_diamond_simulant_carousel.png"),
                            "Diamond simulants — Raman fingerprints")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
