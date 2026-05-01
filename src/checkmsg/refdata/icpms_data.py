"""Bundled reference data for LA-ICP-MS analysis.

Contents:
  - ISOTOPES — natural abundances (mole fractions) for gem-relevant isotopes
  - NIST_SRM_612 — Pearce et al. 1997 preferred concentrations (ppm) for ~40 elements
  - NIST_SRM_610 — high-concentration variant (Pearce et al. 1997)
  - CHONDRITE_REE_PPM — McDonough & Sun 1995 CI-chondrite REE values
  - STACEY_KRAMERS_PB — present-day terrestrial Pb composition (Stacey & Kramers 1975)
  - INTERFERENCES — common polyatomic mass interferences in Ar-plasma ICP-MS
  - U_PB_DECAY_CONSTANTS — Steiger & Jäger 1977 / IUPAC

Sources:
  - IUPAC 2021 (Atomic Weights of the Elements) for natural abundances
  - Pearce, Perkins, Westgate, Gorton, Jackson, Neal & Chenery 1997,
    *Geostandards Newsletter* 21:115-144 (NIST SRM 612 / 610 preferred values)
  - McDonough & Sun 1995, *Chem. Geol.* 120:223-253 (CI chondrite REE)
  - Stacey & Kramers 1975, *Earth Planet. Sci. Lett.* 26:207-221 (terrestrial Pb model)
  - Steiger & Jäger 1977, *Earth Planet. Sci. Lett.* 36:359-362 (U-Pb decay constants)
  - Longerich, Jackson & Günther 1996, *J. Anal. At. Spectrom.* 11:899 (LA-ICP-MS quantitation)
"""

from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Decay constants (1/yr) from Steiger & Jäger 1977 (recommended by IUPAC).
# ---------------------------------------------------------------------------
LAMBDA_238 = 1.55125e-10
LAMBDA_235 = 9.8485e-10
U238_OVER_U235 = 137.818  # natural ratio (Hiess et al. 2012 reaffirmed)

# ---------------------------------------------------------------------------
# Natural-abundance isotope table.
# Keyed by "<Element><Mass>" (e.g. "Pb208"). Abundance is the molar fraction.
# Limited to gem-relevant elements; not the full periodic table.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class IsotopeRecord:
    element: str
    mass: int
    natural_abundance: float


def _i(el: str, mass: int, abund: float) -> tuple[str, IsotopeRecord]:
    return (f"{el}{mass}", IsotopeRecord(element=el, mass=mass, natural_abundance=abund))


ISOTOPES: dict[str, IsotopeRecord] = dict([
    # Light elements
    _i("H", 1, 0.99989),
    _i("B", 10, 0.199), _i("B", 11, 0.801),
    _i("C", 12, 0.9893), _i("C", 13, 0.0107),
    _i("N", 14, 0.99636), _i("N", 15, 0.00364),
    _i("O", 16, 0.99757), _i("O", 17, 0.00038), _i("O", 18, 0.00205),
    # Major rock-forming
    _i("Na", 23, 1.0),
    _i("Mg", 24, 0.7899), _i("Mg", 25, 0.10), _i("Mg", 26, 0.1101),
    _i("Al", 27, 1.0),
    _i("Si", 28, 0.9223), _i("Si", 29, 0.0467), _i("Si", 30, 0.0310),
    _i("P", 31, 1.0),
    _i("S", 32, 0.9499), _i("S", 33, 0.0075), _i("S", 34, 0.0425), _i("S", 36, 0.0001),
    _i("Cl", 35, 0.7576), _i("Cl", 37, 0.2424),
    _i("K", 39, 0.9326),  _i("K", 41, 0.0673),
    _i("Ca", 40, 0.9694), _i("Ca", 42, 0.00647), _i("Ca", 43, 0.00135),
    _i("Ca", 44, 0.02086), _i("Ca", 46, 0.00004), _i("Ca", 48, 0.00187),
    # First-row transition metals and friends
    _i("Sc", 45, 1.0),
    _i("Ti", 46, 0.0825), _i("Ti", 47, 0.0744), _i("Ti", 48, 0.7372),
    _i("Ti", 49, 0.0541), _i("Ti", 50, 0.0518),
    _i("V", 50, 0.00250), _i("V", 51, 0.99750),
    _i("Cr", 50, 0.04345), _i("Cr", 52, 0.83789), _i("Cr", 53, 0.09501), _i("Cr", 54, 0.02365),
    _i("Mn", 55, 1.0),
    _i("Fe", 54, 0.05845), _i("Fe", 56, 0.91754), _i("Fe", 57, 0.02119), _i("Fe", 58, 0.00282),
    _i("Co", 59, 1.0),
    _i("Ni", 58, 0.6808), _i("Ni", 60, 0.2622), _i("Ni", 61, 0.0114),
    _i("Ni", 62, 0.0363), _i("Ni", 64, 0.0093),
    _i("Cu", 63, 0.6917), _i("Cu", 65, 0.3083),
    _i("Zn", 64, 0.4863), _i("Zn", 66, 0.2790), _i("Zn", 67, 0.0410),
    _i("Zn", 68, 0.1875), _i("Zn", 70, 0.0062),
    _i("Ga", 69, 0.601), _i("Ga", 71, 0.399),
    # Heavier elements relevant to gems
    _i("As", 75, 1.0),
    _i("Br", 79, 0.5069), _i("Br", 81, 0.4931),
    _i("Rb", 85, 0.7217), _i("Rb", 87, 0.2783),
    _i("Sr", 84, 0.0056), _i("Sr", 86, 0.0986), _i("Sr", 87, 0.0700), _i("Sr", 88, 0.8258),
    _i("Y", 89, 1.0),
    _i("Zr", 90, 0.5145), _i("Zr", 91, 0.1122), _i("Zr", 92, 0.1715),
    _i("Zr", 94, 0.1738), _i("Zr", 96, 0.0280),
    _i("Nb", 93, 1.0),
    _i("Mo", 92, 0.1453), _i("Mo", 95, 0.1584), _i("Mo", 98, 0.2439),
    _i("Cd", 111, 0.1280), _i("Cd", 114, 0.2873),
    _i("Sn", 118, 0.2422), _i("Sn", 120, 0.3258),
    _i("Cs", 133, 1.0),
    _i("Ba", 137, 0.1123), _i("Ba", 138, 0.7170),
    # REE — at least one major isotope each
    _i("La", 139, 0.99910),
    _i("Ce", 140, 0.8845), _i("Ce", 142, 0.1114),
    _i("Pr", 141, 1.0),
    _i("Nd", 142, 0.272), _i("Nd", 144, 0.238), _i("Nd", 146, 0.172),
    _i("Sm", 147, 0.150), _i("Sm", 152, 0.267),
    _i("Eu", 151, 0.4781), _i("Eu", 153, 0.5219),
    _i("Gd", 157, 0.1565), _i("Gd", 158, 0.2484),
    _i("Tb", 159, 1.0),
    _i("Dy", 161, 0.1891), _i("Dy", 162, 0.2551), _i("Dy", 163, 0.2490),
    _i("Ho", 165, 1.0),
    _i("Er", 166, 0.3350), _i("Er", 167, 0.2287),
    _i("Tm", 169, 1.0),
    _i("Yb", 172, 0.2168), _i("Yb", 173, 0.1610), _i("Yb", 174, 0.3217),
    _i("Lu", 175, 0.97401),
    _i("Hf", 178, 0.2728), _i("Hf", 179, 0.1362), _i("Hf", 180, 0.3508),
    _i("Ta", 181, 0.99988),
    _i("W", 184, 0.3064), _i("W", 186, 0.2843),
    _i("Pt", 195, 0.3378),
    _i("Au", 197, 1.0),
    _i("Hg", 202, 0.2986),
    _i("Tl", 205, 0.7048),
    _i("Pb", 204, 0.014), _i("Pb", 206, 0.241), _i("Pb", 207, 0.221), _i("Pb", 208, 0.524),
    _i("Bi", 209, 1.0),
    _i("Th", 232, 1.0),
    _i("U", 234, 0.000054), _i("U", 235, 0.007204), _i("U", 238, 0.992742),
])


def isotope(key: str) -> IsotopeRecord:
    """Look up an isotope record by 'Pb208'-style key."""
    if key not in ISOTOPES:
        raise KeyError(f"unknown isotope {key!r}; example: 'Pb208', 'Sr88'")
    return ISOTOPES[key]


def isotopes_of(element: str) -> dict[str, IsotopeRecord]:
    """Return all bundled isotopes of a given element."""
    return {k: v for k, v in ISOTOPES.items() if v.element == element}


# ---------------------------------------------------------------------------
# NIST SRM 612 preferred concentrations (ppm) — Pearce et al. 1997, GSN 21:115.
# Subset chosen to cover all gem-relevant elements plus the matrix components.
# Major-element concentrations are quoted in wt% (multiply by 10000 for ppm).
# ---------------------------------------------------------------------------
NIST_SRM_612: dict[str, float] = {
    # Matrix (major elements)
    "Si": 336900.0, "Na": 103400.0, "Ca": 85000.0, "Al": 11200.0,
    # Trace elements (~38 ppm each by design; Pearce preferred values)
    "B": 35.0, "Mg": 77.0, "Sc": 41.0, "Ti": 50.1, "V": 39.0, "Cr": 36.4, "Mn": 38.7, "Fe": 51.0,
    "Co": 35.5, "Ni": 38.8, "Cu": 36.7, "Zn": 39.1, "Ga": 36.3,
    "Rb": 31.4, "Sr": 78.4, "Y": 38.0, "Zr": 37.9, "Nb": 38.9, "Mo": 38.0,
    "Cs": 41.6, "Ba": 39.3,
    "La": 35.8, "Ce": 38.4, "Pr": 37.2, "Nd": 35.5, "Sm": 36.7, "Eu": 35.6,
    "Gd": 36.7, "Tb": 36.0, "Dy": 35.5, "Ho": 37.7, "Er": 37.4, "Tm": 36.8,
    "Yb": 39.2, "Lu": 37.0,
    "Hf": 34.7, "Ta": 39.8, "W": 38.0, "Pt": 3.10, "Au": 5.10,
    "Pb": 38.6, "Th": 37.8, "U": 37.4,
}

# NIST SRM 610 — same matrix but ~10x higher trace concentrations.
NIST_SRM_610: dict[str, float] = {
    "Si": 327200.0, "Na": 99410.0, "Ca": 81600.0, "Al": 10750.0,
    "B": 351.0, "Mg": 432.0, "Ti": 434.0, "V": 442.0, "Cr": 405.0, "Mn": 433.0, "Fe": 458.0,
    "Co": 405.0, "Ni": 458.0, "Cu": 430.0, "Zn": 456.0, "Ga": 433.0,
    "Rb": 425.7, "Sr": 515.5, "Y": 462.0, "Zr": 448.0, "Nb": 419.0,
    "Cs": 360.0, "Ba": 435.0,
    "La": 457.0, "Ce": 453.0, "Nd": 431.0, "Sm": 451.0, "Eu": 461.0,
    "Gd": 444.0, "Dy": 437.0, "Er": 454.0, "Yb": 449.0,
    "Hf": 437.0, "Ta": 444.0, "W": 444.0,
    "Pb": 426.0, "Th": 457.2, "U": 461.5,
}

NIST_REFERENCE_VALUES = {"NIST612": NIST_SRM_612, "NIST610": NIST_SRM_610}

# ---------------------------------------------------------------------------
# CI chondrite REE concentrations (ppm) — McDonough & Sun 1995 Table 5.
# ---------------------------------------------------------------------------
CHONDRITE_REE_PPM: dict[str, float] = {
    "La": 0.237, "Ce": 0.612, "Pr": 0.0928, "Nd": 0.457, "Sm": 0.148,
    "Eu": 0.0563, "Gd": 0.199, "Tb": 0.0361, "Dy": 0.246, "Ho": 0.0546,
    "Er": 0.160, "Tm": 0.0247, "Yb": 0.161, "Lu": 0.0246,
}

REE_ELEMENTS: tuple[str, ...] = (
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
)


# ---------------------------------------------------------------------------
# Stacey & Kramers 1975 — present-day terrestrial Pb composition (model).
# Used as a reference common-Pb composition for U-Pb work and pearl provenance.
# ---------------------------------------------------------------------------
STACEY_KRAMERS_PB: dict[str, float] = {
    "206/204": 18.700,
    "207/204": 15.628,
    "208/204": 38.630,
    "207/206": 15.628 / 18.700,
    "208/206": 38.630 / 18.700,
}


# ---------------------------------------------------------------------------
# Common Ar-plasma polyatomic interferences. Used to flag risky channels;
# correction is not applied automatically (real workflows handle these via
# collision-cell or interference-equation corrections).
# ---------------------------------------------------------------------------
INTERFERENCES: dict[int, tuple[str, ...]] = {
    40: ("40Ar (plasma) on 40Ca",),
    45: ("29Si16O on 45Sc",),
    51: ("35Cl16O on 51V",),
    52: ("40Ar12C on 52Cr", "36Ar16O on 52Cr"),
    56: ("40Ar16O on 56Fe",),
    75: ("40Ar35Cl on 75As",),
    80: ("40Ar40Ar on 80Se",),
}
