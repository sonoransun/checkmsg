"""Mineral / gemstone reference catalog.

`MineralProfile` is the structured per-mineral record that drives the curriculum
example scripts and the unified `diagnose()` pipeline. Each profile holds enough
diagnostic data to (a) synthesize realistic spectra for any of the six bundled
techniques and (b) answer "is this an X?" with explicit reasoning.

Entries are sourced from primary gemological literature:
  - Liddicoat 1989 (Handbook of Gem Identification)
  - Webster 1994 (Gems: Their Sources, Descriptions and Identification, 5th ed.)
  - RRUFF Raman reference database (rruff.info)
  - GIA technical bulletins / Gems & Gemology archive
  - Burns 1993 (Mineralogical Applications of Crystal Field Theory)

Catalog content covers ~50 gems across ten thematic groups: diamond simulants,
blue stones, garnet group end-members, jade family, red stones, black opaques,
tourmaline species, quartz colour-treatments, chrysoberyl, and biominerals.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

import numpy as np

from checkmsg.spectrum import Spectrum
from checkmsg.synthetic import PeakSpec, generate

# ---------------------------------------------------------------------------
# Profile dataclass
# ---------------------------------------------------------------------------

ElementLevel = str  # "major" | "minor" | "trace" | "absent"


@dataclass(frozen=True)
class MineralProfile:
    """A single gemstone or mineral with its diagnostic fingerprints."""

    name: str
    species: str
    aliases: tuple[str, ...] = ()
    chemical_formula: str = ""
    crystal_system: str = ""
    mohs_hardness: tuple[float, float] = (0.0, 0.0)
    density_g_cc: tuple[float, float] = (0.0, 0.0)
    refractive_index: tuple[float, float] = (0.0, 0.0)
    common_colors: tuple[str, ...] = ()
    raman_peaks_cm: tuple[tuple[float, float], ...] = ()
    uvvis_bands_nm: tuple[float, ...] = ()
    chromophores: tuple[str, ...] = ()
    xrf_signature: dict[str, ElementLevel] = field(default_factory=dict)
    libs_signature: dict[str, ElementLevel] = field(default_factory=dict)
    epr_centers: tuple[str, ...] = ()
    icpms_diagnostic_isotopes: tuple[str, ...] = ()
    confusables: tuple[str, ...] = ()
    diagnostic_features: tuple[str, ...] = ()
    references: tuple[str, ...] = ()
    is_amorphous: bool = False

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("MineralProfile.name must be non-empty")


# ---------------------------------------------------------------------------
# Catalog
# ---------------------------------------------------------------------------


def _p(*peaks: tuple[float, float]) -> tuple[tuple[float, float], ...]:
    return tuple(peaks)


CATALOG: dict[str, MineralProfile] = {
    # -------------------- Diamond simulants --------------------
    "diamond": MineralProfile(
        name="diamond", species="diamond", chemical_formula="C",
        crystal_system="cubic", mohs_hardness=(10.0, 10.0),
        density_g_cc=(3.51, 3.53), refractive_index=(2.417, 2.419),
        common_colors=("colorless", "yellow", "brown", "pink", "blue"),
        raman_peaks_cm=_p((1332.5, 1.0)),
        uvvis_bands_nm=(),
        xrf_signature={"C": "major"},
        confusables=("moissanite", "cubic_zirconia", "white_sapphire"),
        diagnostic_features=("razor-sharp Raman line at 1332.5 cm-1",
                             "transparent below 225 nm",
                             "highest hardness of any natural material"),
        references=("Solin & Ramdas 1970", "Knight & White 1992"),
    ),
    "moissanite": MineralProfile(
        name="moissanite", species="silicon carbide", chemical_formula="SiC",
        aliases=("synthetic moissanite", "6H-SiC"),
        crystal_system="hexagonal", mohs_hardness=(9.25, 9.5),
        density_g_cc=(3.21, 3.22), refractive_index=(2.65, 2.69),
        common_colors=("colorless", "yellow-green"),
        raman_peaks_cm=_p((149.0, 0.5), (767.0, 0.9), (789.0, 1.0), (965.0, 0.4)),
        uvvis_bands_nm=(425.0,), chromophores=("Moissanite UV cutoff",),
        xrf_signature={"Si": "major", "C": "major"},
        confusables=("diamond", "cubic_zirconia"),
        diagnostic_features=("folded LO/TO doublet 767 + 789 cm-1",
                             "UV cutoff at ~425 nm",
                             "double refraction visible under loupe"),
        references=("Nakashima & Harima 1997",),
    ),
    "cubic_zirconia": MineralProfile(
        name="cubic_zirconia", species="zirconium dioxide",
        aliases=("CZ", "zirconia"),
        chemical_formula="ZrO2 (Y-stabilised)",
        crystal_system="cubic", mohs_hardness=(8.0, 8.5),
        density_g_cc=(5.6, 6.0), refractive_index=(2.15, 2.18),
        common_colors=("colorless", "various dyed"),
        raman_peaks_cm=_p((269.0, 0.6), (471.0, 1.0), (641.0, 0.8)),
        xrf_signature={"Zr": "major", "Y": "minor"},
        confusables=("diamond", "moissanite", "GGG"),
        diagnostic_features=("only broad Raman bands, no sharp lines",
                             "high density (~5.7 g/cc) vs diamond's 3.52"),
        references=("Pemberton 1999", "Cai et al. 2003"),
    ),
    "GGG": MineralProfile(
        name="GGG", species="gadolinium gallium garnet",
        aliases=("Gd3Ga5O12", "gadolinium-gallium garnet"),
        chemical_formula="Gd3Ga5O12",
        crystal_system="cubic", mohs_hardness=(7.0, 7.5),
        density_g_cc=(7.05, 7.10), refractive_index=(1.97, 1.98),
        common_colors=("colorless",),
        raman_peaks_cm=_p((738.0, 1.0), (350.0, 0.4), (272.0, 0.3)),
        xrf_signature={"Gd": "major", "Ga": "major"},
        confusables=("diamond", "cubic_zirconia"),
        diagnostic_features=("Raman dominated by 738 cm-1 T2g mode",
                             "very high density 7.05 g/cc"),
        references=("Hurrell et al. 1968",),
    ),
    "YAG": MineralProfile(
        name="YAG", species="yttrium aluminium garnet",
        aliases=("Y3Al5O12",),
        chemical_formula="Y3Al5O12",
        crystal_system="cubic", mohs_hardness=(8.25, 8.5),
        density_g_cc=(4.55, 4.57), refractive_index=(1.83, 1.84),
        common_colors=("colorless", "various dyed"),
        raman_peaks_cm=_p((263.0, 0.4), (374.0, 0.5), (561.0, 0.7), (783.0, 1.0)),
        xrf_signature={"Y": "major", "Al": "major"},
        confusables=("diamond", "cubic_zirconia", "GGG"),
        diagnostic_features=("Raman 783 cm-1 dominant",
                             "RI 1.83 distinguishes from CZ (2.18) and diamond (2.42)"),
        references=("Hurrell et al. 1968",),
    ),
    "strontium_titanate": MineralProfile(
        name="strontium_titanate", species="strontium titanate",
        aliases=("fabulite", "diagem"),
        chemical_formula="SrTiO3",
        crystal_system="cubic", mohs_hardness=(5.5, 6.0),
        density_g_cc=(5.10, 5.13), refractive_index=(2.41, 2.42),
        common_colors=("colorless",),
        raman_peaks_cm=_p((100.0, 0.5), (255.0, 0.6), (540.0, 0.8), (810.0, 1.0)),
        xrf_signature={"Sr": "major", "Ti": "major"},
        confusables=("diamond", "moissanite"),
        diagnostic_features=("low Mohs hardness ~5.5 (poor wear)",
                             "very high dispersion creates 'rainbow' fire"),
        references=("Nilsen & Skinner 1968",),
    ),
    "white_sapphire": MineralProfile(
        name="white_sapphire", species="corundum",
        aliases=("colourless sapphire", "leucosapphire"),
        chemical_formula="Al2O3",
        crystal_system="trigonal", mohs_hardness=(9.0, 9.0),
        density_g_cc=(3.98, 4.05), refractive_index=(1.762, 1.770),
        common_colors=("colorless",),
        raman_peaks_cm=_p((378.0, 0.55), (417.0, 1.0), (430.0, 0.3),
                          (450.0, 0.2), (577.0, 0.45), (645.0, 0.5), (750.0, 0.2)),
        xrf_signature={"Al": "major"},
        confusables=("diamond", "white_topaz", "white_spinel"),
        diagnostic_features=("strong corundum 417 cm-1 mode",
                             "no Cr3+ chromophore distinguishes from ruby"),
        references=("Porto & Krishnan 1967",),
    ),
    "white_topaz": MineralProfile(
        name="white_topaz", species="topaz",
        aliases=("colourless topaz",),
        chemical_formula="Al2SiO4(F,OH)2",
        crystal_system="orthorhombic", mohs_hardness=(8.0, 8.0),
        density_g_cc=(3.49, 3.57), refractive_index=(1.609, 1.643),
        common_colors=("colorless", "blue", "yellow", "pink"),
        raman_peaks_cm=_p((269.0, 0.4), (404.0, 0.4), (928.0, 1.0), (1156.0, 0.3)),
        xrf_signature={"Al": "major", "Si": "major", "F": "minor"},
        confusables=("diamond", "white_sapphire", "white_spinel"),
        diagnostic_features=("928 cm-1 Si-O stretch",
                             "fluorine signature in XRF"),
        references=("Pinheiro et al. 2002",),
    ),
    "white_spinel": MineralProfile(
        name="white_spinel", species="spinel",
        aliases=("colourless spinel",),
        chemical_formula="MgAl2O4",
        crystal_system="cubic", mohs_hardness=(7.5, 8.0),
        density_g_cc=(3.55, 3.70), refractive_index=(1.715, 1.730),
        common_colors=("colorless",),
        raman_peaks_cm=_p((313.0, 0.3), (408.0, 0.5), (666.0, 0.7), (766.0, 1.0)),
        xrf_signature={"Mg": "major", "Al": "major"},
        confusables=("diamond", "white_sapphire", "white_topaz"),
        diagnostic_features=("Raman 666 + 766 cm-1 doublet diagnostic of spinel",
                             "isotropic (single refractive index)"),
        references=("Cynn et al. 1992",),
    ),
    "glass_paste": MineralProfile(
        name="glass_paste", species="silicate glass",
        aliases=("paste", "rhinestone glass"),
        chemical_formula="SiO2 (amorphous, with PbO/K2O fluxes)",
        crystal_system="amorphous", mohs_hardness=(5.0, 6.0),
        density_g_cc=(2.4, 4.5), refractive_index=(1.50, 1.70),
        common_colors=("colorless", "any dyed"),
        raman_peaks_cm=_p((460.0, 1.0), (800.0, 0.3), (1080.0, 0.25)),
        xrf_signature={"Si": "major", "Pb": "minor", "K": "minor"},
        confusables=("diamond", "cubic_zirconia"),
        diagnostic_features=("broad amorphous Raman envelope, no sharp peaks",
                             "Pb-glass paste has high density"),
        references=("McMillan 1984",),
        is_amorphous=True,
    ),

    # -------------------- Blue stones --------------------
    "sapphire_blue": MineralProfile(
        name="sapphire_blue", species="corundum",
        aliases=("blue sapphire",),
        chemical_formula="Al2O3:Fe,Ti",
        crystal_system="trigonal", mohs_hardness=(9.0, 9.0),
        density_g_cc=(3.95, 4.05), refractive_index=(1.762, 1.770),
        common_colors=("blue",),
        raman_peaks_cm=_p((378.0, 0.55), (417.0, 1.0), (430.0, 0.3),
                          (577.0, 0.45), (645.0, 0.5), (750.0, 0.2)),
        uvvis_bands_nm=(580.0,),
        chromophores=("Fe2+/Ti4+ IVCT (blue sapphire)",),
        xrf_signature={"Al": "major", "Fe": "trace", "Ti": "trace"},
        libs_signature={"Al": "major", "Fe": "trace", "Ti": "trace"},
        epr_centers=("corundum_Fe3plus",),
        confusables=("tanzanite", "iolite", "blue_topaz", "blue_zircon"),
        diagnostic_features=("417 cm-1 corundum Raman",
                             "580 nm Fe2+/Ti4+ IVCT chromophore"),
        references=("Ferguson & Fielding 1971",),
    ),
    "tanzanite": MineralProfile(
        name="tanzanite", species="zoisite",
        aliases=("blue zoisite",),
        chemical_formula="Ca2Al3(SiO4)(Si2O7)O(OH):V",
        crystal_system="orthorhombic", mohs_hardness=(6.0, 7.0),
        density_g_cc=(3.32, 3.38), refractive_index=(1.691, 1.700),
        common_colors=("blue", "violet"),
        raman_peaks_cm=_p((882.0, 1.0), (916.0, 0.5), (1030.0, 0.4)),
        uvvis_bands_nm=(595.0, 730.0),
        chromophores=("V3+ d-d (tsavorite, V-emerald)",),
        xrf_signature={"Ca": "major", "Al": "major", "Si": "major", "V": "trace"},
        libs_signature={"Ca": "major", "V": "trace"},
        confusables=("sapphire_blue", "iolite"),
        diagnostic_features=("Raman 882 cm-1 Si2O7 mode",
                             "strong V3+ chromophore at 595 nm",
                             "pleochroic blue/violet/burgundy"),
        references=("Liebscher 2004",),
    ),
    "iolite": MineralProfile(
        name="iolite", species="cordierite",
        aliases=("water sapphire", "dichroite"),
        chemical_formula="Mg2Al4Si5O18",
        crystal_system="orthorhombic", mohs_hardness=(7.0, 7.5),
        density_g_cc=(2.55, 2.66), refractive_index=(1.532, 1.547),
        common_colors=("blue", "violet"),
        raman_peaks_cm=_p((426.0, 0.7), (580.0, 0.9), (970.0, 1.0), (1180.0, 0.4)),
        uvvis_bands_nm=(426.0, 596.0),
        chromophores=("Fe2+ d-d (peridot)",),
        xrf_signature={"Mg": "major", "Al": "major", "Si": "major", "Fe": "trace"},
        confusables=("sapphire_blue", "tanzanite"),
        diagnostic_features=("strong pleochroism (blue/violet/yellow)",
                             "Raman 970 cm-1 Si-O stretch"),
        references=("Geiger et al. 2000",),
    ),
    "aquamarine": MineralProfile(
        name="aquamarine", species="beryl",
        aliases=("blue beryl",),
        chemical_formula="Be3Al2Si6O18:Fe2+",
        crystal_system="hexagonal", mohs_hardness=(7.5, 8.0),
        density_g_cc=(2.66, 2.80), refractive_index=(1.564, 1.595),
        common_colors=("blue", "blue-green"),
        raman_peaks_cm=_p((322.0, 0.4), (398.0, 0.6), (685.0, 1.0),
                          (1010.0, 0.25), (1067.0, 0.35)),
        uvvis_bands_nm=(372.0, 829.0),
        xrf_signature={"Be": "major", "Al": "major", "Si": "major", "Fe": "trace"},
        libs_signature={"Be": "major", "Al": "major", "Si": "major", "Fe": "trace"},
        confusables=("sapphire_blue", "blue_topaz", "iolite"),
        diagnostic_features=("Raman 685 cm-1 beryl ring breathing",
                             "Fe2+ in channel sites",
                             "Be diagnostic via LIBS"),
        references=("Wood & Nassau 1968", "Hagemann et al. 1990"),
    ),
    "blue_topaz": MineralProfile(
        name="blue_topaz", species="topaz",
        chemical_formula="Al2SiO4(F,OH)2 (irradiated)",
        crystal_system="orthorhombic", mohs_hardness=(8.0, 8.0),
        density_g_cc=(3.49, 3.57), refractive_index=(1.609, 1.643),
        common_colors=("sky blue", "Swiss blue", "London blue"),
        raman_peaks_cm=_p((269.0, 0.4), (404.0, 0.4), (928.0, 1.0), (1156.0, 0.3)),
        uvvis_bands_nm=(620.0,),
        xrf_signature={"Al": "major", "Si": "major", "F": "minor"},
        libs_signature={"Al": "major", "Si": "major"},
        confusables=("aquamarine", "sapphire_blue", "blue_zircon"),
        diagnostic_features=("928 cm-1 topaz Si-O stretch (vs 685 cm-1 beryl)",
                             "irradiation-induced colour center at 620 nm",
                             "no Be in LIBS (rules out beryl)"),
        references=("Pinheiro et al. 2002",),
    ),
    "blue_zircon": MineralProfile(
        name="blue_zircon", species="zircon",
        aliases=("starlite",),
        chemical_formula="ZrSiO4",
        crystal_system="tetragonal", mohs_hardness=(7.0, 7.5),
        density_g_cc=(4.6, 4.7), refractive_index=(1.92, 1.98),
        common_colors=("blue", "yellow", "red", "brown"),
        raman_peaks_cm=_p((357.0, 0.5), (438.0, 0.7), (974.0, 1.0), (1008.0, 0.5)),
        uvvis_bands_nm=(635.0, 685.0),
        xrf_signature={"Zr": "major", "Si": "major", "Hf": "trace", "U": "trace"},
        icpms_diagnostic_isotopes=("U238", "Pb206", "Pb207", "Hf178"),
        confusables=("aquamarine", "blue_topaz", "sapphire_blue"),
        diagnostic_features=("strong birefringence (~0.04)",
                             "974 cm-1 Si-O zircon mode",
                             "U-Pb dateable"),
        references=("Nasdala et al. 2002",),
    ),

    # -------------------- Garnet group --------------------
    "pyrope": MineralProfile(
        name="pyrope", species="garnet",
        chemical_formula="Mg3Al2(SiO4)3",
        crystal_system="cubic", mohs_hardness=(7.0, 7.5),
        density_g_cc=(3.62, 3.87), refractive_index=(1.714, 1.742),
        common_colors=("red", "purple-red", "pink"),
        raman_peaks_cm=_p((365.0, 0.5), (561.0, 0.8), (909.0, 0.7), (1063.0, 1.0)),
        xrf_signature={"Mg": "major", "Al": "major", "Si": "major", "Fe": "minor"},
        libs_signature={"Mg": "major", "Al": "major", "Si": "major", "Fe": "minor"},
        confusables=("almandine", "rhodolite", "ruby", "red_spinel"),
        diagnostic_features=("garnet-class Raman at 1063 cm-1",
                             "Mg-dominant XRF chemistry"),
        references=("Kolesov & Geiger 1998",),
    ),
    "almandine": MineralProfile(
        name="almandine", species="garnet",
        chemical_formula="Fe3Al2(SiO4)3",
        crystal_system="cubic", mohs_hardness=(7.0, 7.5),
        density_g_cc=(4.05, 4.32), refractive_index=(1.770, 1.820),
        common_colors=("dark red", "brownish red"),
        raman_peaks_cm=_p((350.0, 0.5), (555.0, 0.8), (916.0, 0.7), (1043.0, 1.0)),
        xrf_signature={"Fe": "major", "Al": "major", "Si": "major"},
        libs_signature={"Fe": "major", "Al": "major", "Si": "major"},
        confusables=("pyrope", "rhodolite", "ruby"),
        diagnostic_features=("Fe-dominant XRF chemistry",
                             "highest density of common garnets"),
        references=("Kolesov & Geiger 1998",),
    ),
    "spessartine": MineralProfile(
        name="spessartine", species="garnet",
        aliases=("Mandarin garnet",),
        chemical_formula="Mn3Al2(SiO4)3",
        crystal_system="cubic", mohs_hardness=(6.5, 7.5),
        density_g_cc=(4.12, 4.18), refractive_index=(1.795, 1.815),
        common_colors=("orange", "red-orange"),
        raman_peaks_cm=_p((350.0, 0.5), (552.0, 0.8), (906.0, 0.7), (1027.0, 1.0)),
        xrf_signature={"Mn": "major", "Al": "major", "Si": "major"},
        libs_signature={"Mn": "major", "Al": "major", "Si": "major"},
        confusables=("almandine", "pyrope"),
        diagnostic_features=("Mn-dominant chemistry — orange colour",
                             "intense fluorescence under longwave UV"),
        references=("Kolesov & Geiger 1998",),
    ),
    "grossular": MineralProfile(
        name="grossular", species="garnet",
        aliases=("hessonite",),
        chemical_formula="Ca3Al2(SiO4)3",
        crystal_system="cubic", mohs_hardness=(6.5, 7.5),
        density_g_cc=(3.50, 3.75), refractive_index=(1.730, 1.760),
        common_colors=("green", "orange", "colorless", "pink"),
        raman_peaks_cm=_p((372.0, 0.5), (549.0, 0.8), (879.0, 0.7), (1006.0, 1.0)),
        xrf_signature={"Ca": "major", "Al": "major", "Si": "major"},
        libs_signature={"Ca": "major", "Al": "major", "Si": "major"},
        confusables=("andradite", "tsavorite"),
        diagnostic_features=("Ca-dominant chemistry",
                             "wide colour range from V/Cr/Fe substitutions"),
        references=("Kolesov & Geiger 1998",),
    ),
    "andradite": MineralProfile(
        name="andradite", species="garnet",
        aliases=("demantoid",),
        chemical_formula="Ca3Fe2(SiO4)3",
        crystal_system="cubic", mohs_hardness=(6.5, 7.0),
        density_g_cc=(3.75, 3.85), refractive_index=(1.880, 1.890),
        common_colors=("green", "yellow", "brown", "black"),
        raman_peaks_cm=_p((371.0, 0.5), (510.0, 0.8), (875.0, 0.7), (994.0, 1.0)),
        xrf_signature={"Ca": "major", "Fe": "major", "Si": "major"},
        libs_signature={"Ca": "major", "Fe": "major", "Si": "major"},
        confusables=("grossular", "almandine"),
        diagnostic_features=("Ca + Fe3+ chemistry",
                             "high RI 1.89",
                             "demantoid green from Cr3+"),
        references=("Kolesov & Geiger 1998",),
    ),
    "rhodolite": MineralProfile(
        name="rhodolite", species="garnet",
        aliases=("rhodolite garnet",),
        chemical_formula="(Mg,Fe)3Al2(SiO4)3 (~ pyrope60-almandine40)",
        crystal_system="cubic", mohs_hardness=(7.0, 7.5),
        density_g_cc=(3.83, 3.95), refractive_index=(1.745, 1.760),
        common_colors=("purple-red", "raspberry"),
        raman_peaks_cm=_p((355.0, 0.5), (557.0, 0.8), (912.0, 0.7), (1052.0, 1.0)),
        xrf_signature={"Mg": "major", "Fe": "major", "Al": "major", "Si": "major"},
        confusables=("pyrope", "almandine"),
        diagnostic_features=("Mg + Fe mixed chemistry",
                             "raspberry purple distinguishes from pure pyrope/almandine"),
        references=("Kolesov & Geiger 1998",),
    ),
    "tsavorite": MineralProfile(
        name="tsavorite", species="garnet",
        aliases=("tsavorite garnet", "V-grossular"),
        chemical_formula="Ca3Al2(SiO4)3:V",
        crystal_system="cubic", mohs_hardness=(7.0, 7.5),
        density_g_cc=(3.57, 3.75), refractive_index=(1.730, 1.760),
        common_colors=("green",),
        raman_peaks_cm=_p((372.0, 0.5), (549.0, 0.8), (879.0, 0.7), (1006.0, 1.0)),
        uvvis_bands_nm=(430.0, 605.0),
        chromophores=("V3+ d-d (tsavorite, V-emerald)",),
        xrf_signature={"Ca": "major", "Al": "major", "Si": "major", "V": "trace"},
        libs_signature={"Ca": "major", "Al": "major", "V": "trace"},
        confusables=("grossular", "andradite"),
        diagnostic_features=("V3+ chromophore in grossular host",
                             "pure grossular Raman pattern"),
        references=("Schmetzer & Bank 1979",),
    ),

    # -------------------- Jade family --------------------
    "jadeite": MineralProfile(
        name="jadeite", species="pyroxene",
        aliases=("imperial jade",),
        chemical_formula="NaAlSi2O6",
        crystal_system="monoclinic", mohs_hardness=(6.5, 7.0),
        density_g_cc=(3.30, 3.38), refractive_index=(1.652, 1.688),
        common_colors=("green", "white", "lavender"),
        raman_peaks_cm=_p((374.0, 0.5), (695.0, 0.8), (990.0, 0.7), (1040.0, 1.0)),
        xrf_signature={"Na": "major", "Al": "major", "Si": "major"},
        confusables=("nephrite", "serpentine", "aventurine_quartz"),
        diagnostic_features=("Raman 990 + 1040 cm-1 doublet",
                             "Na pyroxene chemistry"),
        references=("Mao et al. 2007",),
    ),
    "nephrite": MineralProfile(
        name="nephrite", species="amphibole",
        chemical_formula="Ca2(Mg,Fe)5Si8O22(OH)2 (tremolite-actinolite)",
        crystal_system="monoclinic", mohs_hardness=(6.0, 6.5),
        density_g_cc=(2.95, 3.05), refractive_index=(1.605, 1.632),
        common_colors=("green", "cream", "black"),
        raman_peaks_cm=_p((224.0, 0.4), (369.0, 0.6), (678.0, 0.9), (1024.0, 1.0)),
        xrf_signature={"Ca": "major", "Mg": "major", "Si": "major", "Fe": "minor"},
        confusables=("jadeite", "serpentine"),
        diagnostic_features=("224 cm-1 amphibole low-frequency mode",
                             "lower density than jadeite"),
        references=("Mao et al. 2007",),
    ),
    "serpentine": MineralProfile(
        name="serpentine", species="serpentine group",
        aliases=("new jade", "olive jade", "antigorite"),
        chemical_formula="(Mg,Fe)3Si2O5(OH)4",
        crystal_system="monoclinic", mohs_hardness=(2.5, 5.5),
        density_g_cc=(2.45, 2.65), refractive_index=(1.555, 1.575),
        common_colors=("green", "yellow-green"),
        raman_peaks_cm=_p((234.0, 0.6), (385.0, 0.9), (685.0, 1.0)),
        xrf_signature={"Mg": "major", "Si": "major", "Fe": "minor"},
        confusables=("jadeite", "nephrite"),
        diagnostic_features=("low Mohs hardness 2.5–5.5 (scratched easily)",
                             "385 + 685 cm-1 serpentine pattern"),
        references=("Rinaudo et al. 2003",),
    ),
    "aventurine_quartz": MineralProfile(
        name="aventurine_quartz", species="quartz",
        aliases=("aventurine", "Indian jade"),
        chemical_formula="SiO2 + fuchsite mica inclusions",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.64, 2.69), refractive_index=(1.544, 1.553),
        common_colors=("green", "blue", "brown"),
        raman_peaks_cm=_p((128.0, 0.5), (207.0, 0.5), (463.0, 1.0)),
        xrf_signature={"Si": "major", "Cr": "trace"},
        confusables=("jadeite", "nephrite", "serpentine"),
        diagnostic_features=("classic quartz 463 cm-1 Raman",
                             "metallic-glittering inclusions = aventurescence"),
        references=("Etchepare et al. 1974",),
    ),
    "prehnite": MineralProfile(
        name="prehnite", species="phyllosilicate",
        chemical_formula="Ca2Al(AlSi3O10)(OH)2",
        crystal_system="orthorhombic", mohs_hardness=(6.0, 6.5),
        density_g_cc=(2.85, 2.95), refractive_index=(1.611, 1.669),
        common_colors=("green", "yellow"),
        raman_peaks_cm=_p((535.0, 0.6), (953.0, 1.0), (1072.0, 0.5)),
        xrf_signature={"Ca": "major", "Al": "major", "Si": "major"},
        confusables=("jadeite", "serpentine"),
        diagnostic_features=("953 cm-1 Si-O dominant",
                             "yellow-green colour with cleavage"),
        references=("Frost et al. 2007",),
    ),

    # -------------------- Red stones (beyond ruby in 02) --------------------
    "ruby": MineralProfile(
        name="ruby", species="corundum",
        chemical_formula="Al2O3:Cr",
        crystal_system="trigonal", mohs_hardness=(9.0, 9.0),
        density_g_cc=(3.97, 4.05), refractive_index=(1.762, 1.770),
        common_colors=("red", "pink-red"),
        raman_peaks_cm=_p((378.0, 0.55), (417.0, 1.0), (430.0, 0.3),
                          (577.0, 0.45), (645.0, 0.5), (750.0, 0.2)),
        uvvis_bands_nm=(405.0, 555.0),
        chromophores=("Cr3+ d-d (ruby/spinel)",),
        xrf_signature={"Al": "major", "Cr": "trace", "Fe": "trace"},
        libs_signature={"Al": "major", "Cr": "trace"},
        epr_centers=("corundum_Cr3plus",),
        confusables=("red_spinel", "rubellite", "almandine", "pyrope"),
        diagnostic_features=("corundum Raman + Cr3+ chromophore",
                             "fluorescent red R-line emission near 694 nm",
                             "Cr3+ EPR fine-structure pattern"),
        references=("Manenkov & Prokhorov 1956",),
    ),
    "red_spinel": MineralProfile(
        name="red_spinel", species="spinel",
        chemical_formula="MgAl2O4:Cr",
        crystal_system="cubic", mohs_hardness=(7.5, 8.0),
        density_g_cc=(3.55, 3.70), refractive_index=(1.715, 1.730),
        common_colors=("red", "pink"),
        raman_peaks_cm=_p((313.0, 0.3), (408.0, 0.5), (666.0, 0.7), (766.0, 1.0)),
        uvvis_bands_nm=(405.0, 540.0),
        chromophores=("Cr3+ d-d (ruby/spinel)",),
        xrf_signature={"Mg": "major", "Al": "major", "Cr": "trace"},
        libs_signature={"Mg": "major", "Al": "major", "Cr": "trace"},
        confusables=("ruby", "rubellite", "pyrope"),
        diagnostic_features=("Raman 666 + 766 cm-1 doublet (spinel)",
                             "Cr3+ chromophore but cubic crystal (not trigonal corundum)"),
        references=("Cynn et al. 1992", "Wood & Nassau 1968"),
    ),
    "red_beryl": MineralProfile(
        name="red_beryl", species="beryl",
        aliases=("bixbite",),
        chemical_formula="Be3Al2Si6O18:Mn",
        crystal_system="hexagonal", mohs_hardness=(7.5, 8.0),
        density_g_cc=(2.66, 2.80), refractive_index=(1.564, 1.595),
        common_colors=("red", "raspberry"),
        raman_peaks_cm=_p((322.0, 0.4), (398.0, 0.6), (685.0, 1.0),
                          (1010.0, 0.25), (1067.0, 0.35)),
        uvvis_bands_nm=(480.0,),
        xrf_signature={"Be": "major", "Al": "major", "Si": "major", "Mn": "trace"},
        confusables=("ruby", "rubellite", "almandine"),
        diagnostic_features=("beryl Raman pattern with Mn3+ at 480 nm",
                             "extremely rare gem (Utah's Wah Wah Mountains)"),
        references=("Shigley et al. 2003",),
    ),
    "rubellite": MineralProfile(
        name="rubellite", species="tourmaline",
        aliases=("Mn-elbaite", "red tourmaline"),
        chemical_formula="Na(Li,Al)3Al6(BO3)3Si6O18(OH)4:Mn",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.5),
        density_g_cc=(3.00, 3.10), refractive_index=(1.624, 1.644),
        common_colors=("red", "pink", "raspberry"),
        raman_peaks_cm=_p((220.0, 0.5), (372.0, 0.7), (707.0, 1.0), (1060.0, 0.6)),
        uvvis_bands_nm=(515.0,),
        xrf_signature={"Na": "major", "Al": "major", "Si": "major",
                       "Li": "minor", "Mn": "trace", "B": "minor"},
        libs_signature={"B": "minor", "Li": "minor", "Mn": "trace"},
        confusables=("ruby", "red_spinel", "red_beryl"),
        diagnostic_features=("707 cm-1 tourmaline ring mode",
                             "B + Li signature"),
        references=("Hawthorne & Henry 1999",),
    ),
    "rhodochrosite": MineralProfile(
        name="rhodochrosite", species="carbonate",
        chemical_formula="MnCO3",
        crystal_system="trigonal", mohs_hardness=(3.5, 4.0),
        density_g_cc=(3.50, 3.70), refractive_index=(1.578, 1.820),
        common_colors=("pink", "rose-red"),
        raman_peaks_cm=_p((184.0, 0.5), (290.0, 0.4), (720.0, 0.7), (1086.0, 1.0)),
        uvvis_bands_nm=(410.0, 540.0),
        chromophores=("Mn2+ d-d (rhodochrosite, rhodonite)",),
        xrf_signature={"Mn": "major", "C": "major"},
        confusables=("rubellite", "red_beryl"),
        diagnostic_features=("1086 cm-1 carbonate ν1 stretch",
                             "low Mohs hardness 3.5-4 (carbonate)"),
        references=("Edwards et al. 2005",),
    ),

    # -------------------- Black opaque stones --------------------
    "hematite": MineralProfile(
        name="hematite", species="iron oxide",
        chemical_formula="alpha-Fe2O3",
        crystal_system="trigonal", mohs_hardness=(5.5, 6.5),
        density_g_cc=(5.20, 5.30), refractive_index=(2.94, 3.22),
        common_colors=("black", "metallic gray", "blood-red streak"),
        raman_peaks_cm=_p((226.0, 0.5), (245.0, 0.4), (292.0, 1.0),
                          (411.0, 0.7), (497.0, 0.4), (612.0, 0.6)),
        xrf_signature={"Fe": "major"},
        confusables=("magnetite", "obsidian", "schorl"),
        diagnostic_features=("Raman 292 cm-1 + 411 cm-1 hematite signature",
                             "blood-red streak diagnostic vs magnetite"),
        references=("de Faria et al. 1997",),
    ),
    "magnetite": MineralProfile(
        name="magnetite", species="iron oxide",
        chemical_formula="Fe3O4",
        crystal_system="cubic", mohs_hardness=(5.5, 6.5),
        density_g_cc=(5.10, 5.20), refractive_index=(2.42, 2.42),
        common_colors=("black", "metallic"),
        raman_peaks_cm=_p((308.0, 0.4), (540.0, 0.6), (670.0, 1.0)),
        xrf_signature={"Fe": "major"},
        confusables=("hematite", "obsidian"),
        diagnostic_features=("magnetic — picks up steel pin",
                             "broad 670 cm-1 Raman A1g mode"),
        references=("Shebanova & Lazor 2003",),
    ),
    "obsidian": MineralProfile(
        name="obsidian", species="volcanic glass",
        chemical_formula="amorphous SiO2 + Fe + alkali oxides",
        crystal_system="amorphous", mohs_hardness=(5.0, 5.5),
        density_g_cc=(2.30, 2.60), refractive_index=(1.45, 1.55),
        common_colors=("black", "brown", "rainbow"),
        raman_peaks_cm=_p((460.0, 1.0), (800.0, 0.3), (1080.0, 0.25)),
        xrf_signature={"Si": "major", "K": "minor", "Fe": "minor", "Na": "minor"},
        confusables=("hematite", "jet", "onyx"),
        diagnostic_features=("amorphous Raman envelope (no sharp peaks)",
                             "concoidal fracture, lustrous"),
        references=("McMillan 1984",),
        is_amorphous=True,
    ),
    "jet": MineralProfile(
        name="jet", species="lignite",
        aliases=("Whitby jet",),
        chemical_formula="C (compacted lignite, ~ C25H30O3)",
        crystal_system="amorphous", mohs_hardness=(2.5, 4.0),
        density_g_cc=(1.30, 1.34), refractive_index=(1.640, 1.680),
        common_colors=("black", "dark brown"),
        raman_peaks_cm=_p((1340.0, 0.8), (1580.0, 1.0)),
        xrf_signature={"C": "major"},
        confusables=("obsidian", "hematite", "onyx"),
        diagnostic_features=("D + G carbon Raman bands at 1340 + 1580",
                             "very low density 1.3 g/cc",
                             "warm to touch (organic)"),
        references=("Tuinstra & Koenig 1970",),
        is_amorphous=True,
    ),
    "onyx": MineralProfile(
        name="onyx", species="quartz",
        aliases=("black chalcedony", "black agate"),
        chemical_formula="SiO2 (chalcedony)",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.60, 2.65), refractive_index=(1.535, 1.539),
        common_colors=("black", "banded white-black"),
        raman_peaks_cm=_p((128.0, 0.5), (207.0, 0.5), (463.0, 1.0)),
        xrf_signature={"Si": "major"},
        confusables=("obsidian", "jet"),
        diagnostic_features=("classic quartz 463 cm-1 Raman",
                             "Mohs 7 distinguishes from jet (2.5-4)"),
        references=("Etchepare et al. 1974",),
    ),
    "schorl": MineralProfile(
        name="schorl", species="tourmaline",
        aliases=("black tourmaline",),
        chemical_formula="NaFe3Al6(BO3)3Si6O18(OH)4",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.5),
        density_g_cc=(3.10, 3.25), refractive_index=(1.625, 1.655),
        common_colors=("black", "dark brown"),
        raman_peaks_cm=_p((220.0, 0.5), (372.0, 0.7), (707.0, 1.0), (1060.0, 0.6)),
        xrf_signature={"Na": "major", "Fe": "major", "Al": "major", "Si": "major"},
        libs_signature={"Fe": "major", "B": "minor"},
        confusables=("hematite", "obsidian"),
        diagnostic_features=("707 cm-1 tourmaline ring mode",
                             "Fe3 chemistry (vs Li-elbaite)"),
        references=("Hawthorne & Henry 1999",),
    ),

    # -------------------- Tourmaline species --------------------
    "elbaite": MineralProfile(
        name="elbaite", species="tourmaline",
        aliases=("Li-tourmaline",),
        chemical_formula="Na(Li,Al)3Al6(BO3)3Si6O18(OH)4",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.5),
        density_g_cc=(3.00, 3.10), refractive_index=(1.620, 1.652),
        common_colors=("pink", "green", "blue", "watermelon"),
        raman_peaks_cm=_p((220.0, 0.5), (372.0, 0.7), (707.0, 1.0), (1060.0, 0.6)),
        xrf_signature={"Na": "major", "Al": "major", "Si": "major", "Li": "minor", "B": "minor"},
        libs_signature={"Li": "minor", "B": "minor", "Al": "major"},
        confusables=("schorl", "dravite", "rubellite"),
        diagnostic_features=("Li chemistry distinguishes from schorl/dravite",
                             "color-zoned watermelon variety"),
        references=("Hawthorne & Henry 1999",),
    ),
    "dravite": MineralProfile(
        name="dravite", species="tourmaline",
        aliases=("Mg-tourmaline",),
        chemical_formula="NaMg3Al6(BO3)3Si6O18(OH)4",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.5),
        density_g_cc=(3.03, 3.18), refractive_index=(1.610, 1.640),
        common_colors=("brown", "yellow", "green"),
        raman_peaks_cm=_p((220.0, 0.5), (372.0, 0.7), (707.0, 1.0), (1060.0, 0.6)),
        xrf_signature={"Na": "major", "Mg": "major", "Al": "major", "Si": "major", "B": "minor"},
        libs_signature={"Mg": "major", "B": "minor"},
        confusables=("schorl", "elbaite"),
        diagnostic_features=("Mg-dominant tourmaline; brown coloration"),
        references=("Hawthorne & Henry 1999",),
    ),
    "liddicoatite": MineralProfile(
        name="liddicoatite", species="tourmaline",
        aliases=("Ca-tourmaline",),
        chemical_formula="Ca(Li,Al)3Al6(BO3)3Si6O18(OH)4",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.5),
        density_g_cc=(3.02, 3.10), refractive_index=(1.620, 1.652),
        common_colors=("multicolor zoned",),
        raman_peaks_cm=_p((220.0, 0.5), (372.0, 0.7), (707.0, 1.0), (1060.0, 0.6)),
        xrf_signature={"Ca": "major", "Al": "major", "Si": "major", "Li": "minor", "B": "minor"},
        libs_signature={"Ca": "major", "Li": "minor", "B": "minor"},
        confusables=("elbaite",),
        diagnostic_features=("Ca-dominant elbaite analog",
                             "famous for triangular cross-section colour zoning"),
        references=("Dunn et al. 1977",),
    ),

    # -------------------- Quartz colour treatments --------------------
    "rock_crystal": MineralProfile(
        name="rock_crystal", species="quartz",
        aliases=("clear quartz",),
        chemical_formula="SiO2",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.65, 2.66), refractive_index=(1.544, 1.553),
        common_colors=("colorless",),
        raman_peaks_cm=_p((128.0, 0.5), (207.0, 0.5), (463.0, 1.0)),
        xrf_signature={"Si": "major"},
        confusables=("citrine_natural", "amethyst"),
        diagnostic_features=("clean quartz 463 cm-1 with no colour-centre evidence",
                             "no EPR signal at room T"),
        references=("Etchepare et al. 1974",),
    ),
    "citrine_natural": MineralProfile(
        name="citrine_natural", species="quartz",
        chemical_formula="SiO2:Fe3+ (geological yellow)",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.65, 2.66), refractive_index=(1.544, 1.553),
        common_colors=("yellow", "golden"),
        raman_peaks_cm=_p((128.0, 0.5), (207.0, 0.5), (463.0, 1.0)),
        uvvis_bands_nm=(380.0, 450.0),
        chromophores=("Fe3+ spin-forbidden (yellow sapphire)",),
        xrf_signature={"Si": "major", "Fe": "trace"},
        confusables=("citrine_heat_treated", "amethyst"),
        diagnostic_features=("Fe3+ chromophore in geological context",
                             "weak Fe3+ EPR at low temperature"),
        references=("Cohen 1985",),
    ),
    "citrine_heat_treated": MineralProfile(
        name="citrine_heat_treated", species="quartz",
        aliases=("heated amethyst", "burnt amethyst"),
        chemical_formula="SiO2 (Fe3+ from heat-converted amethyst)",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.65, 2.66), refractive_index=(1.544, 1.553),
        common_colors=("yellow", "orange-yellow"),
        raman_peaks_cm=_p((128.0, 0.5), (207.0, 0.5), (463.0, 1.0)),
        uvvis_bands_nm=(380.0, 450.0),
        chromophores=("Fe3+ spin-forbidden (yellow sapphire)",),
        xrf_signature={"Si": "major", "Fe": "trace"},
        epr_centers=("quartz_E1prime", "quartz_Al_hole"),
        confusables=("citrine_natural", "amethyst"),
        diagnostic_features=("residual radiation-induced E1' centre from prior amethyst stage",
                             "Al-hole survivor in EPR despite heat treatment"),
        references=("Mackey & Sander 1972", "Lameiras et al. 2008"),
    ),
    "amethyst": MineralProfile(
        name="amethyst", species="quartz",
        chemical_formula="SiO2:Fe4+ (irradiated)",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.65, 2.66), refractive_index=(1.544, 1.553),
        common_colors=("violet", "purple"),
        raman_peaks_cm=_p((128.0, 0.5), (207.0, 0.5), (463.0, 1.0)),
        uvvis_bands_nm=(545.0,),
        xrf_signature={"Si": "major", "Fe": "trace"},
        epr_centers=("quartz_E1prime", "quartz_Al_hole"),
        confusables=("rock_crystal", "smoky_quartz"),
        diagnostic_features=("Fe4+ on tetrahedral site post-irradiation",
                             "violet 545 nm absorption band"),
        references=("Cohen 1985",),
    ),
    "smoky_quartz": MineralProfile(
        name="smoky_quartz", species="quartz",
        chemical_formula="SiO2 (Al-hole, irradiated)",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.65, 2.66), refractive_index=(1.544, 1.553),
        common_colors=("brown", "smoky", "near-black"),
        raman_peaks_cm=_p((128.0, 0.5), (207.0, 0.5), (463.0, 1.0)),
        uvvis_bands_nm=(450.0, 600.0),
        xrf_signature={"Si": "major", "Al": "trace"},
        epr_centers=("quartz_E1prime", "quartz_Al_hole"),
        confusables=("rock_crystal", "amethyst"),
        diagnostic_features=("Al-hole colour centre from irradiation",
                             "E1' EPR + brown body colour"),
        references=("Weil 1984",),
    ),

    # -------------------- Chrysoberyl --------------------
    "alexandrite": MineralProfile(
        name="alexandrite", species="chrysoberyl",
        chemical_formula="BeAl2O4:Cr",
        crystal_system="orthorhombic", mohs_hardness=(8.5, 8.5),
        density_g_cc=(3.70, 3.78), refractive_index=(1.745, 1.759),
        common_colors=("color-change green-to-red",),
        raman_peaks_cm=_p((354.0, 0.5), (411.0, 0.7), (463.0, 0.4),
                          (798.0, 0.6), (935.0, 1.0)),
        uvvis_bands_nm=(410.0, 580.0),
        chromophores=("Cr3+ d-d (emerald/alexandrite)",),
        xrf_signature={"Be": "major", "Al": "major", "Cr": "trace"},
        libs_signature={"Be": "major", "Al": "major", "Cr": "trace"},
        confusables=("chrysoberyl_yellow", "cymophane"),
        diagnostic_features=("935 cm-1 chrysoberyl Raman",
                             "Cr3+ chromophore producing the green-red 'alexandrite effect'"),
        references=("Schmetzer 2006",),
    ),
    "chrysoberyl_yellow": MineralProfile(
        name="chrysoberyl_yellow", species="chrysoberyl",
        chemical_formula="BeAl2O4 (Fe-yellow)",
        crystal_system="orthorhombic", mohs_hardness=(8.5, 8.5),
        density_g_cc=(3.70, 3.78), refractive_index=(1.745, 1.759),
        common_colors=("yellow", "yellow-green", "honey"),
        raman_peaks_cm=_p((354.0, 0.5), (411.0, 0.7), (463.0, 0.4),
                          (798.0, 0.6), (935.0, 1.0)),
        xrf_signature={"Be": "major", "Al": "major", "Fe": "trace"},
        libs_signature={"Be": "major", "Al": "major", "Fe": "trace"},
        confusables=("alexandrite", "cymophane"),
        diagnostic_features=("chrysoberyl Raman without Cr3+ chromophore",
                             "Fe-coloured yellow",),
        references=("Schmetzer 2006",),
    ),
    "cymophane": MineralProfile(
        name="cymophane", species="chrysoberyl",
        aliases=("cat's-eye chrysoberyl",),
        chemical_formula="BeAl2O4 + rutile inclusions",
        crystal_system="orthorhombic", mohs_hardness=(8.5, 8.5),
        density_g_cc=(3.70, 3.78), refractive_index=(1.745, 1.759),
        common_colors=("yellow", "honey", "brown"),
        raman_peaks_cm=_p((354.0, 0.5), (411.0, 0.7), (463.0, 0.4),
                          (798.0, 0.6), (935.0, 1.0)),
        xrf_signature={"Be": "major", "Al": "major", "Ti": "trace"},
        confusables=("alexandrite", "chrysoberyl_yellow"),
        diagnostic_features=("rutile silk inclusions create cat's-eye chatoyancy",
                             "Ti trace from rutile"),
        references=("Schmetzer 2006",),
    ),

    # -------------------- Biominerals --------------------
    "pearl_natural_saltwater": MineralProfile(
        name="pearl_natural_saltwater", species="aragonite",
        aliases=("Persian pearl", "natural pearl"),
        chemical_formula="CaCO3 (aragonite, biogenic)",
        crystal_system="orthorhombic", mohs_hardness=(2.5, 4.5),
        density_g_cc=(2.60, 2.85), refractive_index=(1.530, 1.685),
        common_colors=("white", "cream", "pink"),
        raman_peaks_cm=_p((153.0, 0.4), (206.0, 0.5), (705.0, 0.5), (1086.0, 1.0)),
        xrf_signature={"Ca": "major", "C": "major", "Sr": "trace", "Mn": "trace"},
        libs_signature={"Ca": "major", "Sr": "trace", "Mn": "trace"},
        epr_centers=("calcite_Mn2plus",),
        icpms_diagnostic_isotopes=("Sr88", "Mn55", "Pb206", "Pb204"),
        confusables=("pearl_akoya", "pearl_freshwater", "coral", "ivory"),
        diagnostic_features=("aragonite Raman pattern (1086 cm-1 ν1)",
                             "low Mn (~25 ppm), saltwater Pb isotope signature"),
        references=("Urmos et al. 1991",),
    ),
    "pearl_akoya": MineralProfile(
        name="pearl_akoya", species="aragonite",
        aliases=("akoya cultured pearl",),
        chemical_formula="CaCO3 + freshwater-mussel bead nucleus",
        crystal_system="orthorhombic", mohs_hardness=(2.5, 4.5),
        density_g_cc=(2.60, 2.85), refractive_index=(1.530, 1.685),
        common_colors=("white", "cream"),
        raman_peaks_cm=_p((153.0, 0.4), (206.0, 0.5), (705.0, 0.5), (1086.0, 1.0)),
        xrf_signature={"Ca": "major", "Sr": "trace", "Mn": "trace"},
        epr_centers=("calcite_Mn2plus",),
        icpms_diagnostic_isotopes=("Sr88", "Pb206", "Pb204"),
        confusables=("pearl_natural_saltwater", "pearl_freshwater"),
        diagnostic_features=("freshwater-shell bead Pb isotope signature",
                             "low ²⁰⁶Pb/²⁰⁴Pb (~16-17) from old continental Pb"),
        references=("Bolzicco et al. 2017",),
    ),
    "pearl_freshwater": MineralProfile(
        name="pearl_freshwater", species="aragonite",
        aliases=("freshwater cultured pearl",),
        chemical_formula="CaCO3 (aragonite)",
        crystal_system="orthorhombic", mohs_hardness=(2.5, 4.5),
        density_g_cc=(2.60, 2.85), refractive_index=(1.530, 1.685),
        common_colors=("white", "cream", "pink", "purple"),
        raman_peaks_cm=_p((153.0, 0.4), (206.0, 0.5), (705.0, 0.5), (1086.0, 1.0)),
        xrf_signature={"Ca": "major", "Mn": "minor", "Sr": "trace"},
        epr_centers=("calcite_Mn2plus",),
        icpms_diagnostic_isotopes=("Mn55", "Sr88", "Pb206"),
        confusables=("pearl_natural_saltwater", "pearl_akoya"),
        diagnostic_features=("high Mn (>500 ppm) from freshwater environment",
                             "low Sr compared to saltwater pearls"),
        references=("Bolzicco et al. 2017",),
    ),
    "coral": MineralProfile(
        name="coral", species="aragonite/calcite",
        aliases=("red coral", "precious coral"),
        chemical_formula="CaCO3 (aragonite or Mg-calcite)",
        crystal_system="trigonal/orthorhombic", mohs_hardness=(3.0, 4.0),
        density_g_cc=(2.60, 2.70), refractive_index=(1.486, 1.658),
        common_colors=("red", "pink", "white", "black"),
        raman_peaks_cm=_p((280.0, 0.5), (712.0, 0.6), (1086.0, 1.0)),
        xrf_signature={"Ca": "major", "C": "major", "Mg": "minor"},
        confusables=("pearl_natural_saltwater", "ivory"),
        diagnostic_features=("Mg-calcite chemistry (vs aragonite of pearls)",
                             "carotenoid Raman fluorescence in red varieties"),
        references=("Urmos et al. 1991",),
    ),
    "ivory": MineralProfile(
        name="ivory", species="hydroxyapatite",
        aliases=("elephant ivory", "tusk"),
        chemical_formula="Ca5(PO4)3(OH) + collagen",
        crystal_system="hexagonal", mohs_hardness=(2.5, 2.75),
        density_g_cc=(1.70, 1.95), refractive_index=(1.535, 1.555),
        common_colors=("cream", "ivory", "white"),
        raman_peaks_cm=_p((430.0, 0.4), (588.0, 0.4), (962.0, 1.0), (1450.0, 0.3)),
        xrf_signature={"Ca": "major", "P": "major"},
        confusables=("pearl_natural_saltwater", "coral"),
        diagnostic_features=("962 cm-1 PO4 ν1 hydroxyapatite",
                             "Schreger lines visible under microscope"),
        references=("Penel et al. 1998",),
    ),
}


# ---------------------------------------------------------------------------
# Catalog API
# ---------------------------------------------------------------------------

# Build a case-insensitive alias index for `get`.
_ALIAS_INDEX: dict[str, str] = {}
for _name, _prof in CATALOG.items():
    _ALIAS_INDEX[_name.lower()] = _name
    for _a in _prof.aliases:
        _ALIAS_INDEX.setdefault(_a.lower(), _name)


def get(name: str) -> MineralProfile:
    """Look up a mineral by name or alias (case-insensitive)."""
    canonical = _ALIAS_INDEX.get(name.lower())
    if canonical is None:
        raise KeyError(f"unknown mineral {name!r}; not in catalog or aliases")
    return CATALOG[canonical]


def by_species(species: str) -> dict[str, MineralProfile]:
    s = species.lower()
    return {k: v for k, v in CATALOG.items() if v.species.lower() == s}


def by_color(color: str) -> dict[str, MineralProfile]:
    c = color.lower()
    return {k: v for k, v in CATALOG.items()
            if any(c in col.lower() for col in v.common_colors)}


def by_confusable(name: str) -> dict[str, MineralProfile]:
    """Minerals that list the given name in their `confusables`."""
    n = name.lower()
    return {k: v for k, v in CATALOG.items()
            if any(c.lower() == n for c in v.confusables)}


def resolve_confusables(name: str) -> dict[str, MineralProfile]:
    """For a given mineral, return the catalog entries listed in its `confusables`."""
    profile = get(name)
    out: dict[str, MineralProfile] = {}
    for c in profile.confusables:
        try:
            out[c] = get(c)
        except KeyError:
            continue
    return out


def names() -> list[str]:
    """Sorted list of canonical names in the catalog."""
    return sorted(CATALOG.keys())


# ---------------------------------------------------------------------------
# Synthesis helpers
# ---------------------------------------------------------------------------

def _peak_specs(peaks: Iterable[tuple[float, float]],
                sigma: float = 2.5, gamma: float = 0.8) -> list[PeakSpec]:
    return [PeakSpec(position=pos, intensity=rel, sigma=sigma, gamma=gamma)
            for pos, rel in peaks]


def synthesize_raman(profile: MineralProfile, *,
                     fields_cm: np.ndarray | None = None,
                     noise: float = 0.01,
                     laser_nm: float | None = None,
                     temperature_K: float = 295.0,
                     seed: int | None = 0,
                     amorphous_fwhm: float = 60.0) -> Spectrum:
    """Synthesize a Raman spectrum from a `MineralProfile`'s peak table.

    Amorphous materials (`profile.is_amorphous=True`) get broad Voigts to
    mimic the noise-glass envelope.
    """
    if fields_cm is None:
        fields_cm = np.linspace(100.0, 1500.0, 1401)
    if not profile.raman_peaks_cm:
        raise ValueError(f"{profile.name} has no Raman peak data in catalog")
    if profile.is_amorphous:
        peaks = [PeakSpec(position=pos, intensity=rel,
                         sigma=amorphous_fwhm / 2.355, gamma=amorphous_fwhm / 4)
                 for pos, rel in profile.raman_peaks_cm]
    else:
        peaks = _peak_specs(profile.raman_peaks_cm)
    return generate(peaks, fields_cm, technique="raman", units="cm-1",
                    noise=noise, seed=seed, laser_nm=laser_nm,
                    temperature_K=temperature_K, fluorescence_amplitude=0.3)


def synthesize_uvvis(profile: MineralProfile, *,
                     fields_nm: np.ndarray | None = None,
                     noise: float = 0.005, seed: int | None = 0,
                     band_fwhm_nm: float = 50.0) -> Spectrum:
    """Synthesize a UV-VIS absorbance spectrum from a profile's chromophore bands."""
    if fields_nm is None:
        fields_nm = np.linspace(380.0, 800.0, 421)
    if not profile.uvvis_bands_nm:
        # No chromophores → flat baseline
        bands = []
    else:
        bands = [PeakSpec(position=b, intensity=1.0,
                         sigma=band_fwhm_nm / 2.355, gamma=band_fwhm_nm / 4)
                 for b in profile.uvvis_bands_nm]
    return generate(bands, fields_nm, technique="uvvis", units="nm",
                    noise=noise, seed=seed)


def synthesize_xrf(profile: MineralProfile, *,
                   fields_keV: np.ndarray | None = None,
                   noise: float = 0.001, seed: int | None = 0) -> Spectrum:
    """Synthesize an XRF spectrum from the profile's element-level signature."""
    from checkmsg.refdata.nist_xray import lines_for
    if fields_keV is None:
        fields_keV = np.linspace(0.5, 15.0, 4096)
    level_to_intensity = {"major": 100.0, "minor": 5.0, "trace": 0.5, "absent": 0.0}
    peaks: list[PeakSpec] = []
    for el, level in profile.xrf_signature.items():
        amp = level_to_intensity.get(level, 0.0)
        if amp <= 0:
            continue
        for line in lines_for(el):
            peaks.append(PeakSpec(position=line.energy_keV,
                                 intensity=amp * line.relative_intensity,
                                 sigma=0.045, gamma=0.025))
    return generate(peaks, fields_keV, technique="xrf", units="keV",
                    noise=noise, seed=seed)


_EPR_CACHE: dict[tuple[str, float], Spectrum] = {}


def synthesize_epr(profile: MineralProfile, *,
                   frequency_GHz: float = 9.5,
                   fields_mT: np.ndarray | None = None,
                   noise: float = 0.005, seed: int | None = 0) -> Spectrum | None:
    """Synthesise an EPR spectrum using the first bundled center listed by the profile.

    Returns None if the profile lists no EPR-active centers. Results are cached
    by (center, frequency) — the curriculum scripts run many specimens that share
    centres, and re-simulating the 36-dim Mn²⁺ Hamiltonian per call is the dominant
    test-runtime cost.
    """
    if not profile.epr_centers:
        return None
    from checkmsg.epr import simulate_field_sweep
    from checkmsg.refdata.epr_centers import CENTERS
    center_key = profile.epr_centers[0]
    if center_key not in CENTERS:
        return None
    cache_key = (center_key, frequency_GHz)
    if cache_key not in _EPR_CACHE:
        spin_system = CENTERS[center_key]
        if fields_mT is None:
            center = frequency_GHz * 1000.0 / (13.99624 * 2.0)
            # Tighter range + fewer fields than the full epr.simulate default.
            sweep = np.linspace(center - 50.0, center + 50.0, 401)
        else:
            sweep = fields_mT
        spec = simulate_field_sweep(spin_system, frequency_GHz=frequency_GHz,
                                    fields_mT=sweep, orientations=(5, 1))
        _EPR_CACHE[cache_key] = spec
    spec = _EPR_CACHE[cache_key]
    rng = np.random.default_rng(seed)
    peak = float(np.max(np.abs(spec.intensity)) or 1.0)
    noisy = spec.intensity + rng.normal(0.0, noise * peak, size=spec.intensity.size)
    return Spectrum(spec.axis.copy(), noisy, "epr", spec.units, dict(spec.metadata))


def synthesize_libs(profile: MineralProfile, *,
                    fields_nm: np.ndarray | None = None,
                    noise: float = 0.003, seed: int | None = 0) -> Spectrum:
    """Synthesize a LIBS emission spectrum from the profile's libs_signature."""
    from checkmsg.refdata.nist_asd import lines_for
    if fields_nm is None:
        fields_nm = np.linspace(200.0, 800.0, 6000)
    level_to_intensity = {"major": 1.0, "minor": 0.4, "trace": 0.1, "absent": 0.0}
    peaks: list[PeakSpec] = []
    for el, level in profile.libs_signature.items():
        amp = level_to_intensity.get(level, 0.0)
        if amp <= 0:
            continue
        for line in lines_for(el):
            peaks.append(PeakSpec(position=line.wavelength_nm,
                                 intensity=amp * line.relative_intensity,
                                 sigma=0.10, gamma=0.05))
    return generate(peaks, fields_nm, technique="libs", units="nm",
                    noise=noise, seed=seed)
