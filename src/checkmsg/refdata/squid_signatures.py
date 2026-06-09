"""Bundled SQUID magnetometry reference signatures for gemological identification.

Each `MagneticMineral` record captures the bulk-magnetism fingerprint that a
SQUID magnetometer would observe for a canonical mineral specimen. The values
combine room-temperature volume susceptibility, magnetic ordering type
(diamagnetic / paramagnetic / ferromagnetic / ferrimagnetic /
antiferromagnetic / canted-AFM), and the relevant Curie / Néel temperatures
plus saturation moment and coercivity for ordered phases.

Sources:
  - Dunlop & Özdemir 1997, *Rock Magnetism: Fundamentals and Frontiers* (CUP)
  - O'Reilly 1984, *Rock and Mineral Magnetism* (Blackie)
  - Hunt, Moskowitz, Banerjee 1995, *Rock Physics and Phase Relations* (AGU Ref. Shelf 3)
  - Morin 1950, *Phys. Rev.* 78:819 (hematite Morin transition)
  - Néel 1948, *Ann. Phys.* 3:137 (ferrimagnetism in spinel ferrites)
  - Wohlfarth 1980, *Ferromagnetic Materials* vol. 1 (North-Holland)

All susceptibilities are SI volume susceptibility (dimensionless). Saturation
moments are in emu/g (= A·m²/kg / 10³). Curie / Néel temperatures are in
kelvin. Coercivities are in mT.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

Ordering = Literal[
    "diamagnetic",
    "paramagnetic",
    "ferromagnetic",
    "ferrimagnetic",
    "antiferromagnetic",
    "canted-afm",
]

ORDERINGS: tuple[Ordering, ...] = (
    "diamagnetic",
    "paramagnetic",
    "ferromagnetic",
    "ferrimagnetic",
    "antiferromagnetic",
    "canted-afm",
)


@dataclass(frozen=True)
class MagneticMineral:
    """Canonical magnetic fingerprint for one mineral species or mineral class."""

    name: str
    ordering: Ordering
    curie_K: float = 0.0          # ferromagnet / ferrimagnet ordering temperature
    neel_K: float = 0.0            # antiferromagnet / canted-AFM ordering temperature
    weiss_K: float = 0.0           # paramagnetic Curie-Weiss intercept (negative for AFM coupling)
    saturation_emu_g: float = 0.0  # bulk saturation magnetisation at T<<Tc
    susceptibility_si: tuple[float, float] = (0.0, 0.0)  # T=295 K range (low, high)
    coercivity_mT: float = 0.0     # B-field at M=0 on hysteresis loop
    morin_K: float = 0.0           # spin-flop transition (hematite-class only)
    notes: str = ""
    references: tuple[str, ...] = field(default_factory=tuple)


# Room-temperature SI volume susceptibilities are wide ranges in the literature
# because grain-size and impurity content matter. Bracket values follow Hunt et
# al. 1995 Table 1 unless otherwise noted.

SIGNATURES: dict[str, MagneticMineral] = {
    # ---------- Strongly ordered phases ----------
    "magnetite": MagneticMineral(
        name="magnetite",
        ordering="ferrimagnetic",
        curie_K=858.0,
        saturation_emu_g=92.0,
        susceptibility_si=(1.0e-3, 5.7e-3),
        coercivity_mT=20.0,
        notes="Inverse spinel Fe3O4. Strongest natural ferrimagnet at room T.",
        references=("Dunlop & Özdemir 1997, Ch. 3", "O'Reilly 1984"),
    ),
    "hematite": MagneticMineral(
        name="hematite",
        ordering="canted-afm",
        neel_K=948.0,
        saturation_emu_g=0.4,
        susceptibility_si=(5.0e-4, 4.0e-3),
        coercivity_mT=300.0,
        morin_K=263.0,
        notes="α-Fe2O3 — canted AFM with weak parasitic ferromagnetism above the Morin transition (263 K).",
        references=("Morin 1950", "Dunlop & Özdemir 1997, Ch. 3"),
    ),
    "pyrrhotite_4c": MagneticMineral(
        name="pyrrhotite_4c",
        ordering="ferrimagnetic",
        curie_K=593.0,
        saturation_emu_g=20.0,
        susceptibility_si=(1.0e-3, 1.0e-2),
        coercivity_mT=50.0,
        notes="Fe7S8 monoclinic ('4C') pyrrhotite — vacancy-ordered ferrimagnet.",
        references=("Hunt, Moskowitz, Banerjee 1995",),
    ),
    "goethite": MagneticMineral(
        name="goethite",
        ordering="antiferromagnetic",
        neel_K=393.0,
        saturation_emu_g=0.05,
        susceptibility_si=(2.5e-5, 1.5e-4),
        coercivity_mT=2000.0,
        notes="α-FeOOH antiferromagnet with very high coercivity from defect spins.",
        references=("Hunt, Moskowitz, Banerjee 1995",),
    ),
    "ilmenite": MagneticMineral(
        name="ilmenite",
        ordering="antiferromagnetic",
        neel_K=55.0,
        saturation_emu_g=0.1,
        susceptibility_si=(1.0e-4, 1.5e-3),
        coercivity_mT=5.0,
        notes="FeTiO3 — paramagnetic at room T (T_N=55 K well below 295 K).",
        references=("O'Reilly 1984",),
    ),
    "feNiCo_catalyst": MagneticMineral(
        name="feNiCo_catalyst",
        ordering="ferromagnetic",
        curie_K=900.0,
        saturation_emu_g=180.0,
        susceptibility_si=(5.0e-1, 5.0e0),
        coercivity_mT=2.0,
        notes="Synthetic-diamond HPHT growth catalyst (Fe-Ni-Co metallic flux) — small "
              "ferromagnetic inclusion gives detectable moment in otherwise diamagnetic diamond.",
        references=("Sumiya & Satoh 1996",),
    ),

    # ---------- Paramagnetic chromophores ----------
    "Cr3plus_paramagnet": MagneticMineral(
        name="Cr3plus_paramagnet",
        ordering="paramagnetic",
        weiss_K=-2.0,
        susceptibility_si=(1.0e-7, 1.5e-5),
        notes="Dilute Cr3+ (S=3/2) substituting Al3+ in corundum/spinel/chrysoberyl. "
              "Curie-Weiss above ~10 K with small AFM Weiss constant.",
        references=("Manenkov & Prokhorov 1956",),
    ),
    "Mn2plus_paramagnet": MagneticMineral(
        name="Mn2plus_paramagnet",
        ordering="paramagnetic",
        weiss_K=-1.0,
        susceptibility_si=(1.0e-7, 5.0e-5),
        notes="Mn2+ (S=5/2 high-spin) in calcite/aragonite biominerals or carbonates. "
              "AC susceptibility scales linearly with [Mn2+] — non-destructive freshwater pearl screen.",
        references=("Bernstein 1979",),
    ),
    "Fe2plus_paramagnet": MagneticMineral(
        name="Fe2plus_paramagnet",
        ordering="paramagnetic",
        weiss_K=-5.0,
        susceptibility_si=(5.0e-6, 3.0e-4),
        notes="Octahedral Fe2+ in silicates (almandine, iolite, peridot, nephrite). "
              "Higher moment than Cr3+/Mn2+ at equal concentration.",
        references=("Burns 1993",),
    ),
    "Fe3plus_paramagnet": MagneticMineral(
        name="Fe3plus_paramagnet",
        ordering="paramagnetic",
        weiss_K=-3.0,
        susceptibility_si=(5.0e-6, 3.0e-4),
        notes="Fe3+ (S=5/2) in andradite, schorl tourmaline, citrine. "
              "High-spin like Mn2+ but on different site geometries.",
        references=("Burns 1993",),
    ),
    "Cu2plus_paramagnet": MagneticMineral(
        name="Cu2plus_paramagnet",
        ordering="paramagnetic",
        weiss_K=-5.0,
        susceptibility_si=(1.0e-5, 2.0e-4),
        notes="Cu2+ (S=1/2, d9) in malachite/azurite/turquoise. AFM superexchange "
              "between adjacent Cu gives a small negative Weiss constant.",
        references=("Burns 1993", "Dunlop & Özdemir 1997"),
    ),
    "Gd3plus_paramagnet": MagneticMineral(
        name="Gd3plus_paramagnet",
        ordering="paramagnetic",
        weiss_K=-1.0,
        susceptibility_si=(5.0e-3, 2.0e-2),
        notes="Gd3+ (S=7/2, half-filled 4f) in GGG — large paramagnetic moment.",
        references=("Brixner 1964",),
    ),
    "Mn3plus_paramagnet": MagneticMineral(
        name="Mn3plus_paramagnet",
        ordering="paramagnetic",
        weiss_K=-2.0,
        susceptibility_si=(1.0e-6, 8.0e-5),
        notes="Mn3+ (S=2, Jahn-Teller) in kunzite/sugilite.",
        references=("Burns 1993",),
    ),

    # ---------- Diamagnetic baselines ----------
    "diamond_diamagnetic": MagneticMineral(
        name="diamond_diamagnetic",
        ordering="diamagnetic",
        susceptibility_si=(-2.2e-5, -2.0e-5),
        notes="Pure type IIa diamond — diamagnetic from filled C-C bonding orbitals. "
              "Any positive moment indicates ferromagnetic catalyst inclusion (HPHT screening).",
        references=("Hudgens 1953",),
    ),
    "quartz_diamagnetic": MagneticMineral(
        name="quartz_diamagnetic",
        ordering="diamagnetic",
        susceptibility_si=(-1.5e-5, -1.3e-5),
        notes="Pure α-quartz / chalcedony — diamagnetic baseline. "
              "Any positive moment in coloured quartz indicates Fe trapped impurities.",
        references=("Hunt, Moskowitz, Banerjee 1995",),
    ),
    "calcite_diamagnetic": MagneticMineral(
        name="calcite_diamagnetic",
        ordering="diamagnetic",
        susceptibility_si=(-1.4e-5, -1.0e-5),
        notes="Pure CaCO3 (calcite/aragonite without Mn). "
              "Used to subtract host susceptibility before reading off Mn2+ in pearls.",
        references=("Hunt, Moskowitz, Banerjee 1995",),
    ),
}


# ---- Convenience accessors ----


def by_ordering(ordering: Ordering) -> dict[str, MagneticMineral]:
    """Filter the SIGNATURES registry by ordering type."""
    return {k: v for k, v in SIGNATURES.items() if v.ordering == ordering}


def diamagnetic_baselines() -> dict[str, MagneticMineral]:
    """Reference diamagnetic minerals used as host-subtraction baselines."""
    return by_ordering("diamagnetic")


def diagnostic_for(scenario: str) -> tuple[str, ...]:
    """Names of canonical signatures relevant to a given gemological scenario."""
    table: dict[str, tuple[str, ...]] = {
        "black-opaque-discrimination": ("magnetite", "hematite", "Fe3plus_paramagnet"),
        "hpht-diamond-screening": ("diamond_diamagnetic", "feNiCo_catalyst"),
        "pearl-mn-quantitation": ("Mn2plus_paramagnet", "calcite_diamagnetic"),
        "schorl-vs-elbaite": ("Fe3plus_paramagnet", "quartz_diamagnetic"),
    }
    return table.get(scenario, ())
