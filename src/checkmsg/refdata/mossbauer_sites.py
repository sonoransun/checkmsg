"""⁵⁷Fe Mössbauer reference site library.

Mössbauer spectroscopy resolves iron oxidation state (Fe²⁺ vs Fe³⁺) and site
geometry from the isomer shift (δ, the doublet centroid relative to α-Fe) and
the quadrupole splitting (ΔEQ, the doublet separation), both in mm/s at ~295 K.
Magnetically ordered oxides give a six-line pattern instead (a sextet).

Sources:
  - Amthauer, Annersten & Hafner 1976, *Z. Kristallogr.* 143:14 — garnet Fe²⁺/Fe³⁺ sites.
  - Goldman, Rossman & Parkin 1978, *Phys. Chem. Minerals* 3:225 — Fe in beryl.
  - Burns 1993, *Mineralogical Applications of Crystal Field Theory* (CUP) — Fe site systematics.
  - Dunlop & Özdemir 1997, *Rock Magnetism* (CUP) — hematite/magnetite hyperfine fields.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class FeSite:
    """A reference Fe environment: isomer shift + quadrupole splitting (mm/s)."""

    name: str
    valence: str                  # "Fe2+" | "Fe3+"
    isomer_shift_mm_s: float
    quadrupole_splitting_mm_s: float
    delta_tol: float = 0.10       # δ tolerance
    qs_tol: float = 0.30          # ΔEQ tolerance
    coordination: str = ""        # "octahedral" | "tetrahedral" | "dodecahedral"
    magnetic_sextet: bool = False  # ordered oxide → six-line pattern, not a doublet
    hyperfine_field_T: float = 0.0  # for sextets
    typical_in: tuple[str, ...] = ()
    notes: str = ""
    references: tuple[str, ...] = field(default_factory=tuple)


SITES: dict[str, FeSite] = {
    "corundum_Fe3plus": FeSite(
        "Fe3+ in corundum", "Fe3+", 0.37, 0.65, coordination="octahedral",
        typical_in=("sapphire_blue",),
        notes="Substitutional Fe3+ on the Al site in corundum.",
        references=("Burns 1993",)),
    "almandine_Fe2plus": FeSite(
        "Fe2+ in almandine (X-site)", "Fe2+", 1.28, 3.50, coordination="dodecahedral",
        typical_in=("almandine", "rhodolite"),
        notes="Large dodecahedral Fe2+ site; very large ΔEQ.",
        references=("Amthauer et al. 1976",)),
    "andradite_Fe3plus": FeSite(
        "Fe3+ in andradite (Y-site)", "Fe3+", 0.39, 0.55, coordination="octahedral",
        typical_in=("andradite", "demantoid"),
        notes="Octahedral Fe3+ in the andradite Y-site.",
        references=("Amthauer et al. 1976",)),
    "beryl_Fe2plus": FeSite(
        "Fe2+ in beryl", "Fe2+", 1.10, 2.50, coordination="octahedral",
        typical_in=("aquamarine",),
        notes="Octahedral Fe2+ in beryl; pairs with Fe3+ for the blue/green balance.",
        references=("Goldman et al. 1978",)),
    "beryl_Fe3plus": FeSite(
        "Fe3+ in beryl", "Fe3+", 0.38, 0.60, coordination="octahedral",
        typical_in=("aquamarine", "heliodor"),
        notes="Octahedral/channel Fe3+ in beryl.",
        references=("Goldman et al. 1978",)),
    "tourmaline_Fe2plus": FeSite(
        "Fe2+ in tourmaline", "Fe2+", 1.07, 2.40, coordination="octahedral",
        typical_in=("schorl", "dravite"),
        notes="Octahedral Fe2+ on the tourmaline Y/Z sites.",
        references=("Burns 1993",)),
    "olivine_Fe2plus": FeSite(
        "Fe2+ in olivine (M1+M2)", "Fe2+", 1.15, 2.95, coordination="octahedral",
        typical_in=("peridot",),
        notes="Fe2+ on the olivine M1/M2 octahedral sites.",
        references=("Burns 1993",)),
    "hematite_sextet": FeSite(
        "Fe3+ hematite sextet", "Fe3+", 0.37, -0.20, coordination="octahedral",
        magnetic_sextet=True, hyperfine_field_T=51.7,
        typical_in=("hematite",),
        notes="Magnetically split Fe3+ sextet, B_hf ≈ 51.7 T at 295 K.",
        references=("Dunlop & Özdemir 1997",)),
    "magnetite_sextet": FeSite(
        "Fe magnetite sextets", "Fe3+", 0.40, 0.0, coordination="octahedral",
        magnetic_sextet=True, hyperfine_field_T=49.0,
        typical_in=("magnetite",),
        notes="Inverse-spinel A+B site sextets (modelled as one mean sextet).",
        references=("Dunlop & Özdemir 1997",)),
}


def by_valence(valence: str) -> dict[str, FeSite]:
    return {k: v for k, v in SITES.items() if v.valence == valence}
