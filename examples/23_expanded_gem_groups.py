"""Expanded catalogue carousels — feldspar, copper, beryl, and accessory gems.

Four `diagnose()` carousels over the gem groups added in the catalogue
expansion. Each routes every specimen through the unified pipeline and asserts a
minimum-correct threshold (lowered for the genuinely confusable groups, in
keeping with the project's honesty about hard cases — e.g. colourless beryl
cannot be separated from other beryls by Raman + XRF alone).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> None:
    args = parse_smoke_args("23_expanded_gem_groups")

    run_carousel("Scenario 1: feldspar group",
                 ["orthoclase", "labradorite", "amazonite", "sunstone"],
                 require_correct=3)
    run_carousel("Scenario 2: copper minerals",
                 ["malachite", "azurite", "turquoise", "chrysocolla"],
                 require_correct=3)
    beryls = run_carousel("Scenario 3: beryl varieties (colourless is a hard case)",
                          ["emerald", "morganite", "heliodor", "goshenite"],
                          require_correct=2)
    run_carousel("Scenario 4: accessory & borosilicate gems",
                 ["sphene", "apatite", "danburite", "uvarovite", "kornerupine"],
                 require_correct=4)

    if not args.smoke:
        plot_raman_carousel(beryls, output_path("23_expanded_gem_groups.png"),
                            "Beryl varieties — Raman fingerprints")
        print(f"saved plot to {output_path('23_expanded_gem_groups.png')}")

    print("OK")


if __name__ == "__main__":
    main()
