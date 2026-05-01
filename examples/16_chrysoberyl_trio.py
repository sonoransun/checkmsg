"""Chrysoberyl trio — same host, three optical phenomena.

  alexandrite        Cr³⁺-coloured chrysoberyl, exhibits "alexandrite effect"
                     (green in daylight, red under tungsten light)
  cymophane          chrysoberyl with rutile silk inclusions; cat's-eye effect
  chrysoberyl_yellow Fe-coloured ordinary chrysoberyl

All share Raman pattern 354/411/463/798/935 cm⁻¹. Cr³⁺ chromophore distinguishes
alexandrite via UV-VIS; Ti from rutile silk distinguishes cymophane via XRF.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> int:
    args = parse_smoke_args("16_chrysoberyl_trio")
    specimens = ["alexandrite", "cymophane", "chrysoberyl_yellow"]
    rows = run_carousel("Scenario 16: chrysoberyl trio", specimens,
                        require_correct=len(specimens) - 1)
    if not args.smoke:
        plot_raman_carousel(rows, output_path("16_chrysoberyl_trio.png"),
                            "Chrysoberyl trio — same Raman, different chromophores")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
