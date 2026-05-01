"""Garnet species suite — five end-members of one mineral group.

All garnets share the cubic crystal structure and the SiO4 tetrahedron, so
their Raman patterns look closely similar. The chemical end-members are
discriminated by their dominant cation:

  pyrope     Mg3Al2(SiO4)3   (Mg-rich, red)
  almandine  Fe3Al2(SiO4)3   (Fe-rich, deep red)
  spessartine Mn3Al2(SiO4)3  (Mn-rich, orange)
  grossular  Ca3Al2(SiO4)3   (Ca-rich, green/orange)
  andradite  Ca3Fe2(SiO4)3   (Ca + Fe³⁺, green)

Plus rhodolite (a Mg-Fe pyralspite mix).

Pedagogy: same Raman, different XRF/LIBS chemistry — chemistry-led ID.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args, plot_raman_carousel, run_carousel  # noqa: E402


def main() -> int:
    args = parse_smoke_args("10_garnet_species_suite")
    specimens = ["pyrope", "almandine", "spessartine", "grossular", "andradite", "rhodolite"]
    rows = run_carousel("Scenario 10: garnet species suite", specimens,
                        require_correct=len(specimens) - 2)
    if not args.smoke:
        plot_raman_carousel(rows, output_path("10_garnet_species_suite.png"),
                            "Garnet group — same Raman, different chemistry")
    print("OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
