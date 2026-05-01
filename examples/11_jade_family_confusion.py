"""Jade family confusion — four very different minerals all sold as "jade".

  jadeite          NaAlSi2O6 (pyroxene, "imperial jade")
  nephrite         Ca2(Mg,Fe)5Si8O22(OH)2 (amphibole, traditional Asian jade)
  serpentine       (Mg,Fe)3Si2O5(OH)4 ("new jade", much softer)
  aventurine quartz SiO2 + fuchsite ("Indian jade", actually quartzite)
  prehnite         Ca2Al(AlSi3O10)(OH)2 (recently-marketed "olive jade")

Hardness alone (jadeite 6.5–7 vs serpentine 2.5–5.5) is enough to separate
some of these in jewellery practice. Raman gives a clean ID without scratching
the stone.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> int:
    args = parse_smoke_args("11_jade_family_confusion")
    specimens = ["jadeite", "nephrite", "serpentine", "aventurine_quartz", "prehnite"]
    rows = run_carousel("Scenario 11: jade family confusion", specimens,
                        require_correct=len(specimens) - 1)
    if not args.smoke:
        plot_raman_carousel(rows, output_path("11_jade_family_confusion.png"),
                            "Jade family — Raman fingerprints expose the imposters")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
