"""Tourmaline species — same skeleton, four different cation flavours.

The tourmaline group shares the framework (BO3)3Si6O18 structure but the X-site
and Y-site occupants vary across species:

  elbaite       Na(Li,Al)3Al6...   (Li-rich, multicolor including watermelon)
  schorl        NaFe3Al6...        (Fe-rich, black)
  dravite       NaMg3Al6...        (Mg-rich, brown)
  liddicoatite  Ca(Li,Al)3Al6...   (Ca-Li, famous triangular zoning)

Pedagogy: Raman patterns are nearly identical. LIBS chemistry (Li/B/Mg/Fe/Ca
levels) is the discriminator.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> int:
    args = parse_smoke_args("12_tourmaline_species")
    specimens = ["elbaite", "schorl", "dravite", "liddicoatite"]
    rows = run_carousel("Scenario 12: tourmaline species", specimens,
                        require_correct=len(specimens) - 2)
    if not args.smoke:
        plot_raman_carousel(rows, output_path("12_tourmaline_species.png"),
                            "Tourmaline group — Raman common, chemistry diverges")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
