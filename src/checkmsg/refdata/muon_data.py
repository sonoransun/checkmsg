"""Bundled reference data for muon imaging (transmission, scattering, muonic X-ray).

Constants and tables are kept small and literature-cited. Sources:

  - Particle Data Group 2024, *Phys. Rev. D* 110:030001 (Bethe-Bloch K, m_µ, m_e, R_y).
  - Tsai 1974, *Rev. Mod. Phys.* 46:815 (radiation lengths X_0 — used for the bundled
    materials below; values are ~5 % lower than the Geant4 default tables).
  - Engfer, Schneuwly, Vuilleumier, Walter & Zehnder 1974, *At. Data Nucl. Data Tables*
    14:509 (muonic K_α energies; tabulated for Z=1..100, we cull to gem-relevant Z).
  - Daniel 1979, *At. Data Nucl. Data Tables* 24:165 (corrections to muonic K_α
    energies for finite nuclear size).

The bundled `MATERIALS` registry covers the gem hosts and structural metals that
appear in `examples/20_muon_tomography.py`. Add new entries by computing X_0 from
Tsai 1974 Table I or via the closed-form fit ``X_0 = 716.4 * A / [Z(Z+1) ln(287/sqrt(Z))]``.
"""

from __future__ import annotations

from dataclasses import dataclass

# --- Physics constants (PDG 2024) -------------------------------------------------

M_MU_MeV: float = 105.6583755   # muon rest energy
M_E_MeV: float = 0.51099895     # electron rest energy
R_Y_eV: float = 13.605693       # Rydberg energy
K_BB: float = 0.30707           # Bethe-Bloch K = 4 pi N_A r_e^2 m_e c^2  (MeV cm² mol⁻¹)
ALPHA_FS: float = 7.2973525693e-3  # fine-structure constant


# --- Materials --------------------------------------------------------------------


@dataclass(frozen=True)
class Material:
    """A material with the parameters needed for Bethe-Bloch + Highland.

    Z_eff and A_eff are effective atomic and mass numbers (composition-weighted
    averages for compounds). I_eV is the mean excitation potential. X_0_g_cm2 is
    the radiation length expressed as a mass thickness — the natural variable for
    Highland multiple-scattering since the resulting RMS angle depends on x/X_0.
    """

    name: str
    Z_eff: float
    A_eff: float
    density_g_cc: float
    X_0_g_cm2: float
    I_eV: float


# Bundled common materials (gemological + structural).  X_0 from Tsai 1974,
# I from PDG; Z_eff / A_eff are stoichiometric averages for compounds.
MATERIALS: dict[str, Material] = {
    "vacuum": Material("vacuum", 0.0, 0.0, 0.0, float("inf"), 0.0),
    "air":      Material("air",      7.31,  14.36,  0.001205, 36.62, 85.7),
    "water":    Material("water",    7.42,  14.32,  1.000,    36.08, 75.0),
    "polymer":  Material("polymer (CH2)", 5.43, 10.39, 1.06,   44.77, 60.0),
    "quartz":   Material("quartz",   10.80, 21.38,  2.65,     27.05, 139.2),
    "calcite":  Material("calcite",  15.66, 31.36,  2.71,     19.51, 136.4),
    "corundum": Material("corundum", 11.31, 22.64,  4.00,     27.94, 145.2),
    "beryl":    Material("beryl",    9.55,  19.10,  2.75,     31.56, 132.3),
    "olivine":  Material("olivine (forsterite)", 12.41, 25.13, 3.27, 24.91, 143.6),
    "diamond":  Material("diamond",  6.0,   12.01,  3.51,     42.70, 78.0),
    "iron":     Material("iron",     26.0,  55.85,  7.874,    13.84, 286.0),
    "Fe-Ni":    Material("iron-nickel (taenite)", 26.74, 56.91, 7.9,  13.69, 290.0),
    "copper":   Material("copper",   29.0,  63.55,  8.96,     12.86, 322.0),
    "silver":   Material("silver",   47.0, 107.87, 10.49,      8.97, 470.0),
    "platinum": Material("platinum", 78.0, 195.08, 21.45,      6.54, 790.0),
    "gold":     Material("gold",     79.0, 196.97, 19.32,      6.46, 790.0),
    "lead":     Material("lead",     82.0, 207.2,  11.35,      6.37, 823.0),
    "uranium":  Material("uranium",  92.0, 238.03, 18.95,      6.00, 890.0),
}


def get_material(name: str) -> Material:
    """Look up a bundled `Material` by name (case-insensitive)."""
    key = name.lower()
    for canonical in MATERIALS:
        if canonical.lower() == key:
            return MATERIALS[canonical]
    raise KeyError(f"unknown material {name!r}; known: {sorted(MATERIALS)}")


# --- Muonic K_α energies ----------------------------------------------------------
# Energies in keV.  Engfer et al. 1974, Table I; values include reduced-mass and
# vacuum-polarisation corrections as tabulated.  For Z=82 the finite-nucleus
# correction lowers the K_α from the point-nucleus value by ~25 keV.

MUONIC_KALPHA_keV: dict[str, float] = {
    "H":   1.66,    "He":  6.31,    "Li":  18.7,    "Be":  33.4,    "B":   52.2,
    "C":   75.3,    "N":  102.4,    "O":  133.5,    "F":  168.5,    "Ne": 207.6,
    "Na": 250.3,    "Mg": 296.9,    "Al": 346.8,    "Si": 400.2,    "P":  457.0,
    "S":  517.0,    "Cl": 580.5,    "Ar": 647.4,    "K":  716.7,    "Ca": 783.7,
    "Ti": 931.0,    "V":  1009.3,   "Cr": 1090.4,   "Mn": 1174.0,   "Fe": 1257.6,
    "Co": 1344.3,   "Ni": 1432.3,   "Cu": 1521.8,   "Zn": 1613.3,   "Ga": 1707.9,
    "Sr": 2415.0,   "Y":  2515.0,   "Zr": 2618.5,   "Nb": 2723.9,   "Mo": 2829.6,
    "Ag": 3485.0,   "Sn": 3737.0,   "Sb": 3865.0,   "Ba": 4453.0,
    "W":  5645.0,   "Pt": 5963.0,   "Au": 6019.0,   "Hg": 6076.0,   "Pb": 5778.0,
    "Bi": 5839.0,   "Th": 6438.0,   "U":  6519.0,
}


def muonic_kalpha_keV(element: str) -> float:
    """Return the muonic K_α (2p → 1s) transition energy in keV for the given element.

    Raises KeyError if the element isn't in the bundled subset.
    """
    if element not in MUONIC_KALPHA_keV:
        raise KeyError(
            f"no muonic K_α tabulated for {element!r}; "
            f"known elements: {sorted(MUONIC_KALPHA_keV)}"
        )
    return MUONIC_KALPHA_keV[element]
