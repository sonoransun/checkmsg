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

Catalog content covers ~96 gems across thematic groups: diamond simulants
(incl. lab-grown CVD/HPHT diamond), blue stones, garnet group end-members, jade
family, red stones, black opaques, tourmaline species, quartz colour-treatments,
chrysoberyl, biominerals, feldspar group, beryl varieties, spodumene, copper
minerals, lapis/sodalite, and assorted silicates/halides.

Additional source families for the expanded catalog: Freeman et al. 2008
(feldspar Raman), Hofmeister & Rossman 1985 (amazonite), Frost et al. 2002
(malachite/azurite), Osticioli et al. 2009 (lazurite S3-), Kolesov & Geiger 2004
(olivine), Krishnan 1947 (fluorite), Nilsen 1969 (sphalerite), Smallwood et al.
1997 (opal), Wang et al. 2012 / D'Haenens-Johansson et al. 2015 (lab-grown
diamond), Edgar & Hutton 1978 (Cr3+ in beryl).
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
    # SQUID magnetometry signature (unset = not characterised).
    # `squid_ordering` is one of: "diamagnetic", "paramagnetic", "ferromagnetic",
    # "ferrimagnetic", "antiferromagnetic", "canted-afm", or "" if unknown.
    squid_ordering: str = ""
    squid_curie_K: float = 0.0
    squid_neel_K: float = 0.0
    squid_weiss_K: float = 0.0
    squid_saturation_emu_g: float = 0.0
    squid_susceptibility_si: tuple[float, float] = (0.0, 0.0)
    squid_coercivity_mT: float = 0.0
    squid_morin_K: float = 0.0
    # Modern-lab technique signatures (unset = not characterised).
    pl_centers: tuple[float, ...] = ()          # PL emission ZPL wavelengths (nm)
    ftir_bands: tuple[float, ...] = ()          # diagnostic FTIR band centres (cm-1)
    diamond_type: str = ""                      # FTIR IR type: "Ia"|"Ib"|"IIa"|"IIb"
    mossbauer_sites: tuple[str, ...] = ()       # keys into refdata.mossbauer_sites.SITES
    cl_bands: tuple[float, ...] = ()            # CL emission band centres (nm)

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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-2.2e-5, -2.0e-5),
        pl_centers=(415.2, 503.2),  # N3 + H3 (natural cape/aggregated)
        diamond_type="Ia",
        cl_bands=(440.0,),  # band-A
    ),
    "diamond_cvd": MineralProfile(
        name="diamond_cvd", species="diamond",
        aliases=("CVD diamond", "CVD synthetic diamond"),
        chemical_formula="C",
        crystal_system="cubic", mohs_hardness=(10.0, 10.0),
        density_g_cc=(3.51, 3.53), refractive_index=(2.417, 2.419),
        common_colors=("colorless", "brown", "pink"),
        raman_peaks_cm=_p((1332.5, 1.0)),
        xrf_signature={"C": "major"},
        confusables=("diamond", "diamond_hpht", "moissanite"),
        diagnostic_features=("Raman 1332.5 identical to natural — needs PL",
                             "Si-V centre at ~737 nm is the CVD growth marker",
                             "typically low-nitrogen Type IIa"),
        references=("Wang et al. 2012", "Eaton-Magaña & Shigley 2016"),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-2.2e-5, -2.0e-5),
        pl_centers=(575.0, 637.0, 736.6),  # NV0, NV-, Si-V (synthetic markers)
        diamond_type="IIa",
    ),
    "diamond_hpht": MineralProfile(
        name="diamond_hpht", species="diamond",
        aliases=("HPHT diamond", "HPHT synthetic diamond"),
        chemical_formula="C",
        crystal_system="cubic", mohs_hardness=(10.0, 10.0),
        density_g_cc=(3.51, 3.53), refractive_index=(2.417, 2.419),
        common_colors=("yellow", "colorless", "blue"),
        raman_peaks_cm=_p((1332.5, 1.0)),
        xrf_signature={"C": "major"},
        epr_centers=("diamond_Ni_HPHT",),
        confusables=("diamond", "diamond_cvd", "moissanite"),
        diagnostic_features=("Raman 1332.5 identical to natural — needs PL/EPR",
                             "Ni-related EPR centre from metallic growth catalyst",
                             "ferromagnetic Fe-Ni-Co flux inclusion → SQUID moment",
                             "often single-substitutional-N Type Ib"),
        references=("Sumiya & Satoh 1996", "D'Haenens-Johansson et al. 2015"),
        squid_ordering="ferromagnetic",
        squid_curie_K=900.0,
        squid_saturation_emu_g=2.0,
        squid_susceptibility_si=(1.0e-3, 1.0e-2),
        squid_coercivity_mT=2.0,
        pl_centers=(575.0, 637.0),  # NV pair
        diamond_type="Ib",
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.1e-5, -9.0e-6),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-1.0,
        squid_susceptibility_si=(5.0e-3, 2.0e-2),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.1e-5, -9.0e-6),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.1e-5, -9.0e-6),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.1e-5, -9.0e-6),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
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
        mossbauer_sites=("corundum_Fe3plus",),
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
        mossbauer_sites=("beryl_Fe2plus", "beryl_Fe3plus"),
        ftir_bands=(3600.0, 3700.0),  # structural water / OH
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
        cl_bands=(480.0, 575.0),  # REE3+ + broad zircon luminescence
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-1.0,
        squid_susceptibility_si=(1.0e-6, 5.0e-6),
        mossbauer_sites=("almandine_Fe2plus",),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-5.0,
        squid_susceptibility_si=(5.0e-5, 1.5e-4),
        mossbauer_sites=("almandine_Fe2plus",),
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
        aliases=(),
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
        aliases=(),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-3.0,
        squid_susceptibility_si=(2.0e-5, 1.0e-4),
        mossbauer_sites=("andradite_Fe3plus",),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-3.0,
        squid_susceptibility_si=(2.0e-5, 7.0e-5),
        mossbauer_sites=("almandine_Fe2plus",),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-5.0,
        squid_susceptibility_si=(1.0e-5, 5.0e-5),
        ftir_bands=(3675.0,),  # tremolite-actinolite OH stretch
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-3.0,
        squid_susceptibility_si=(5.0e-6, 3.0e-5),
        ftir_bands=(3675.0,),  # serpentine OH stretch
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.5e-5, -1.3e-5),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
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
        cl_bands=(694.0,),  # Cr3+ R-line red luminescence
        squid_ordering="paramagnetic",
        squid_weiss_K=-2.0,
        squid_susceptibility_si=(1.0e-7, 1.5e-5),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-2.0,
        squid_susceptibility_si=(1.0e-7, 5.0e-6),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-2.0,
        squid_susceptibility_si=(5.0e-6, 5.0e-5),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-2.0,
        squid_susceptibility_si=(5.0e-6, 3.0e-5),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-3.0,
        squid_susceptibility_si=(5.0e-5, 5.0e-4),
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
        squid_ordering="canted-afm",
        squid_neel_K=948.0,
        squid_saturation_emu_g=0.4,
        squid_susceptibility_si=(5.0e-4, 4.0e-3),
        squid_coercivity_mT=300.0,
        squid_morin_K=263.0,
        mossbauer_sites=("hematite_sextet",),
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
        squid_ordering="ferrimagnetic",
        squid_curie_K=858.0,
        squid_saturation_emu_g=92.0,
        squid_susceptibility_si=(1.0e-3, 5.7e-3),
        squid_coercivity_mT=20.0,
        mossbauer_sites=("magnetite_sextet",),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-5.0,
        squid_susceptibility_si=(5.0e-6, 3.0e-5),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.5e-5, -1.3e-5),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-3.0,
        squid_susceptibility_si=(5.0e-5, 3.0e-4),
        mossbauer_sites=("tourmaline_Fe2plus",),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.5e-5, -1.3e-5),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-1.0,
        squid_susceptibility_si=(1.0e-7, 5.0e-6),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-1.0,
        squid_susceptibility_si=(1.0e-7, 5.0e-6),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.5e-5, -1.3e-5),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.5e-5, -1.3e-5),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-2.0,
        squid_susceptibility_si=(1.0e-7, 1.5e-5),
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-1.0,
        squid_susceptibility_si=(1.0e-7, 5.0e-6),
        cl_bands=(610.0,),  # Mn2+-activated carbonate luminescence
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-1.0,
        squid_susceptibility_si=(1.0e-7, 5.0e-6),
        cl_bands=(610.0,),  # Mn2+-activated carbonate luminescence
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
        squid_ordering="paramagnetic",
        squid_weiss_K=-1.0,
        squid_susceptibility_si=(1.0e-5, 1.0e-4),
        cl_bands=(610.0,),  # Mn2+-activated carbonate luminescence (freshwater: stronger)
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.4e-5, -1.0e-5),
        cl_bands=(610.0,),  # Mn2+ in Mg-calcite
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
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
    ),

    # -------------------- Feldspar group --------------------
    "orthoclase": MineralProfile(
        name="orthoclase", species="feldspar",
        aliases=("moonstone", "adularia"),
        chemical_formula="KAlSi3O8",
        crystal_system="monoclinic", mohs_hardness=(6.0, 6.5),
        density_g_cc=(2.55, 2.63), refractive_index=(1.518, 1.526),
        common_colors=("colorless", "white", "blue sheen"),
        raman_peaks_cm=_p((290.0, 0.4), (476.0, 1.0), (513.0, 0.9)),
        xrf_signature={"K": "major", "Al": "major", "Si": "major"},
        libs_signature={"K": "major", "Al": "major", "Si": "major"},
        confusables=("rock_crystal", "labradorite", "amazonite"),
        diagnostic_features=("feldspar Raman doublet 476 + 513 cm-1",
                             "adularescent sheen in moonstone"),
        references=("Freeman et al. 2008",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-9.0e-6, -7.0e-6),
    ),
    "labradorite": MineralProfile(
        name="labradorite", species="feldspar",
        chemical_formula="(Ca,Na)(Al,Si)4O8",
        crystal_system="triclinic", mohs_hardness=(6.0, 6.5),
        density_g_cc=(2.68, 2.72), refractive_index=(1.559, 1.573),
        common_colors=("gray", "labradorescent blue-green"),
        raman_peaks_cm=_p((480.0, 1.0), (508.0, 0.9)),
        xrf_signature={"Na": "major", "Ca": "major", "Al": "major", "Si": "major"},
        libs_signature={"Na": "major", "Ca": "major", "Al": "major", "Si": "major"},
        confusables=("orthoclase", "amazonite", "sunstone"),
        diagnostic_features=("plagioclase Raman 480 + 508 cm-1",
                             "labradorescence from lamellar intergrowth"),
        references=("Freeman et al. 2008",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-9.0e-6, -7.0e-6),
    ),
    "amazonite": MineralProfile(
        name="amazonite", species="feldspar",
        aliases=("amazonstone",),
        chemical_formula="KAlSi3O8 (Pb)",
        crystal_system="triclinic", mohs_hardness=(6.0, 6.5),
        density_g_cc=(2.55, 2.63), refractive_index=(1.522, 1.530),
        common_colors=("blue-green",),
        raman_peaks_cm=_p((476.0, 1.0), (513.0, 0.9)),
        uvvis_bands_nm=(625.0,),
        xrf_signature={"K": "major", "Al": "major", "Si": "major", "Pb": "trace"},
        libs_signature={"K": "major", "Al": "major", "Si": "major"},
        confusables=("orthoclase", "turquoise", "labradorite"),
        diagnostic_features=("microcline Raman + Pb/H2O blue-green colour",
                             "625 nm absorption from structural Pb"),
        references=("Hofmeister & Rossman 1985",),
    ),
    "sunstone": MineralProfile(
        name="sunstone", species="feldspar",
        aliases=("heliolite",),
        chemical_formula="(Na,Ca)AlSi3O8 + Cu",
        crystal_system="triclinic", mohs_hardness=(6.0, 6.5),
        density_g_cc=(2.62, 2.67), refractive_index=(1.537, 1.547),
        common_colors=("orange", "red", "green"),
        raman_peaks_cm=_p((480.0, 1.0), (508.0, 0.8)),
        xrf_signature={"Na": "major", "Ca": "major", "Al": "major", "Si": "major", "Cu": "trace"},
        libs_signature={"Na": "major", "Al": "major", "Si": "major", "Cu": "trace"},
        confusables=("labradorite", "aventurine_quartz"),
        diagnostic_features=("plagioclase Raman + native-Cu platelet aventurescence",
                             "Cu detectable by LIBS"),
        references=("Hofmann et al. 2017",),
    ),

    # -------------------- Beryl varieties --------------------
    "emerald": MineralProfile(
        name="emerald", species="beryl",
        aliases=("green beryl",),
        chemical_formula="Be3Al2Si6O18:Cr,V",
        crystal_system="hexagonal", mohs_hardness=(7.5, 8.0),
        density_g_cc=(2.67, 2.78), refractive_index=(1.565, 1.602),
        common_colors=("green",),
        raman_peaks_cm=_p((322.0, 0.4), (398.0, 0.6), (685.0, 1.0), (1068.0, 0.4)),
        uvvis_bands_nm=(430.0, 605.0),
        chromophores=("Cr3+ d-d (emerald/alexandrite)",),
        xrf_signature={"Be": "major", "Al": "major", "Si": "major", "Cr": "trace"},
        libs_signature={"Be": "major", "Al": "major", "Si": "major", "Cr": "trace"},
        epr_centers=("beryl_Cr3plus",),
        ftir_bands=(3600.0, 3700.0),
        confusables=("tsavorite", "hiddenite", "peridot", "aquamarine"),
        diagnostic_features=("beryl 685 cm-1 + Cr3+ chromophore",
                             "Be via LIBS distinguishes from green garnet",
                             "FTIR water bands indicate growth method"),
        references=("Wood & Nassau 1968", "Edgar & Hutton 1978"),
    ),
    "morganite": MineralProfile(
        name="morganite", species="beryl",
        aliases=("pink beryl",),
        chemical_formula="Be3Al2Si6O18:Mn2+",
        crystal_system="hexagonal", mohs_hardness=(7.5, 8.0),
        density_g_cc=(2.71, 2.90), refractive_index=(1.572, 1.600),
        common_colors=("pink", "peach"),
        raman_peaks_cm=_p((322.0, 0.4), (398.0, 0.6), (685.0, 1.0)),
        uvvis_bands_nm=(410.0, 540.0),
        chromophores=("Mn2+ d-d (rhodochrosite, rhodonite)",),
        xrf_signature={"Be": "major", "Al": "major", "Si": "major", "Mn": "trace"},
        libs_signature={"Be": "major", "Al": "major", "Si": "major", "Mn": "trace"},
        confusables=("kunzite", "rose_quartz", "rubellite"),
        diagnostic_features=("beryl 685 cm-1 + Mn2+ pink chromophore",
                             "Be via LIBS"),
        references=("Wood & Nassau 1968",),
    ),
    "heliodor": MineralProfile(
        name="heliodor", species="beryl",
        aliases=("golden beryl",),
        chemical_formula="Be3Al2Si6O18:Fe3+",
        crystal_system="hexagonal", mohs_hardness=(7.5, 8.0),
        density_g_cc=(2.66, 2.80), refractive_index=(1.564, 1.595),
        common_colors=("yellow", "golden"),
        raman_peaks_cm=_p((322.0, 0.4), (398.0, 0.6), (685.0, 1.0)),
        uvvis_bands_nm=(380.0, 450.0),
        chromophores=("Fe3+ spin-forbidden (yellow sapphire)",),
        xrf_signature={"Be": "major", "Al": "major", "Si": "major", "Fe": "trace"},
        libs_signature={"Be": "major", "Al": "major", "Si": "major", "Fe": "trace"},
        mossbauer_sites=("beryl_Fe3plus",),
        confusables=("citrine_natural", "chrysoberyl_yellow", "goshenite"),
        diagnostic_features=("beryl 685 cm-1 + Fe3+ yellow",
                             "Mössbauer Fe3+ doublet"),
        references=("Wood & Nassau 1968",),
        squid_ordering="paramagnetic",
        squid_weiss_K=-3.0,
        squid_susceptibility_si=(5.0e-6, 3.0e-5),
    ),
    "goshenite": MineralProfile(
        name="goshenite", species="beryl",
        aliases=("colourless beryl",),
        chemical_formula="Be3Al2Si6O18",
        crystal_system="hexagonal", mohs_hardness=(7.5, 8.0),
        density_g_cc=(2.66, 2.80), refractive_index=(1.564, 1.595),
        common_colors=("colorless",),
        raman_peaks_cm=_p((322.0, 0.4), (398.0, 0.6), (685.0, 1.0)),
        xrf_signature={"Be": "major", "Al": "major", "Si": "major"},
        libs_signature={"Be": "major", "Al": "major", "Si": "major"},
        confusables=("rock_crystal", "white_topaz", "white_sapphire"),
        diagnostic_features=("beryl 685 cm-1, colourless (no chromophore)",
                             "Be via LIBS"),
        references=("Hagemann et al. 1990",),
    ),

    # -------------------- Spodumene --------------------
    "kunzite": MineralProfile(
        name="kunzite", species="spodumene",
        aliases=("pink spodumene",),
        chemical_formula="LiAlSi2O6:Mn3+",
        crystal_system="monoclinic", mohs_hardness=(6.5, 7.0),
        density_g_cc=(3.15, 3.21), refractive_index=(1.660, 1.676),
        common_colors=("pink", "violet"),
        raman_peaks_cm=_p((260.0, 0.4), (705.0, 1.0), (1060.0, 0.5)),
        uvvis_bands_nm=(545.0,),
        chromophores=("Mn3+ d-d (kunzite)",),
        xrf_signature={"Al": "major", "Si": "major", "Mn": "trace"},
        libs_signature={"Li": "major", "Al": "major", "Si": "major", "Mn": "trace"},
        confusables=("morganite", "rose_quartz", "rubellite"),
        diagnostic_features=("spodumene Raman 705 cm-1 + Mn3+ pink",
                             "Li seen by LIBS but invisible to XRF"),
        references=("Burns 1993", "Litvin et al. 1973"),
        squid_ordering="paramagnetic",
        squid_weiss_K=-2.0,
        squid_susceptibility_si=(1.0e-6, 8.0e-5),
    ),
    "hiddenite": MineralProfile(
        name="hiddenite", species="spodumene",
        aliases=("green spodumene",),
        chemical_formula="LiAlSi2O6:Cr",
        crystal_system="monoclinic", mohs_hardness=(6.5, 7.0),
        density_g_cc=(3.15, 3.21), refractive_index=(1.655, 1.681),
        common_colors=("green",),
        raman_peaks_cm=_p((260.0, 0.4), (705.0, 1.0), (1060.0, 0.5)),
        uvvis_bands_nm=(430.0, 605.0),
        chromophores=("Cr3+ d-d (emerald/alexandrite)",),
        xrf_signature={"Al": "major", "Si": "major", "Cr": "trace"},
        libs_signature={"Li": "major", "Al": "major", "Si": "major", "Cr": "trace"},
        confusables=("emerald", "tsavorite", "peridot"),
        diagnostic_features=("spodumene Raman 705 cm-1 + Cr3+ green",
                             "Li via LIBS"),
        references=("Burns 1993",),
    ),

    # -------------------- Quartz varieties (colour beyond treatments) --------------------
    "rose_quartz": MineralProfile(
        name="rose_quartz", species="quartz",
        chemical_formula="SiO2 (dumortierite micro-inclusions)",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.64, 2.66), refractive_index=(1.544, 1.553),
        common_colors=("pink",),
        raman_peaks_cm=_p((128.0, 0.4), (206.0, 0.5), (463.0, 1.0)),
        uvvis_bands_nm=(500.0,),
        xrf_signature={"Si": "major", "Ti": "trace"},
        confusables=("morganite", "rock_crystal", "amethyst"),
        diagnostic_features=("quartz 463 cm-1 + broad ~500 nm pink",
                             "colour from borosilicate nanofibres"),
        references=("Goreva, Ma & Rossman 2001",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.5e-5, -1.3e-5),
    ),
    "prasiolite": MineralProfile(
        name="prasiolite", species="quartz",
        aliases=("green amethyst", "vermarine"),
        chemical_formula="SiO2:Fe (heated amethyst)",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.64, 2.66), refractive_index=(1.544, 1.553),
        common_colors=("green",),
        raman_peaks_cm=_p((128.0, 0.4), (206.0, 0.5), (463.0, 1.0)),
        xrf_signature={"Si": "major", "Fe": "trace"},
        epr_centers=("quartz_E1prime",),
        confusables=("peridot", "rock_crystal", "amethyst"),
        diagnostic_features=("quartz 463 cm-1, green from heat-treated Fe",
                             "E1' radiation-history EPR centre"),
        references=("Nunes & Lameiras 2017",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.5e-5, -1.3e-5),
    ),

    # -------------------- Copper minerals & ornamentals --------------------
    "malachite": MineralProfile(
        name="malachite", species="copper carbonate",
        chemical_formula="Cu2CO3(OH)2",
        crystal_system="monoclinic", mohs_hardness=(3.5, 4.0),
        density_g_cc=(3.6, 4.05), refractive_index=(1.655, 1.909),
        common_colors=("green",),
        raman_peaks_cm=_p((433.0, 0.6), (539.0, 0.7), (1060.0, 1.0), (1492.0, 0.5)),
        uvvis_bands_nm=(700.0, 780.0),
        chromophores=("Cu2+ d-d (malachite/turquoise)",),
        xrf_signature={"Cu": "major"},
        libs_signature={"Cu": "major"},
        confusables=("azurite", "turquoise", "chrysocolla"),
        diagnostic_features=("carbonate 1060 cm-1 + Cu2+ green",
                             "banded green colour, low hardness"),
        references=("Frost et al. 2002",),
        squid_ordering="paramagnetic",
        squid_weiss_K=-5.0,
        squid_susceptibility_si=(1.0e-5, 2.0e-4),
    ),
    "azurite": MineralProfile(
        name="azurite", species="copper carbonate",
        chemical_formula="Cu3(CO3)2(OH)2",
        crystal_system="monoclinic", mohs_hardness=(3.5, 4.0),
        density_g_cc=(3.7, 3.9), refractive_index=(1.730, 1.838),
        common_colors=("deep blue",),
        raman_peaks_cm=_p((403.0, 0.7), (545.0, 0.4), (838.0, 0.6), (1095.0, 1.0)),
        uvvis_bands_nm=(700.0, 780.0),
        chromophores=("Cu2+ d-d (malachite/turquoise)",),
        xrf_signature={"Cu": "major"},
        libs_signature={"Cu": "major"},
        confusables=("malachite", "lapis_lazuli", "turquoise"),
        diagnostic_features=("carbonate 1095 cm-1 + Cu2+ deep blue",
                             "403 cm-1 azurite lattice mode"),
        references=("Frost et al. 2002",),
        squid_ordering="paramagnetic",
        squid_weiss_K=-5.0,
        squid_susceptibility_si=(1.0e-5, 2.0e-4),
    ),
    "turquoise": MineralProfile(
        name="turquoise", species="copper aluminium phosphate",
        chemical_formula="CuAl6(PO4)4(OH)8·4H2O",
        crystal_system="triclinic", mohs_hardness=(5.0, 6.0),
        density_g_cc=(2.6, 2.9), refractive_index=(1.610, 1.650),
        common_colors=("blue", "blue-green"),
        raman_peaks_cm=_p((430.0, 0.5), (820.0, 0.4), (1040.0, 1.0)),
        uvvis_bands_nm=(430.0, 720.0),
        xrf_signature={"Cu": "major", "Al": "major", "P": "major"},
        libs_signature={"Cu": "major", "Al": "major", "P": "trace"},
        confusables=("amazonite", "chrysocolla", "lapis_lazuli"),
        diagnostic_features=("phosphate 1040 cm-1 + Cu blue-green",
                             "Cu + P chemistry distinguishes from amazonite"),
        references=("Čejka et al. 2015",),
        squid_ordering="paramagnetic",
        squid_weiss_K=-5.0,
        squid_susceptibility_si=(1.0e-5, 1.0e-4),
    ),
    "chrysocolla": MineralProfile(
        name="chrysocolla", species="copper silicate",
        chemical_formula="(Cu,Al)2H2Si2O5(OH)4·nH2O",
        crystal_system="amorphous", mohs_hardness=(2.5, 3.5),
        density_g_cc=(1.9, 2.4), refractive_index=(1.575, 1.635),
        common_colors=("blue", "blue-green"),
        raman_peaks_cm=_p((430.0, 1.0), (670.0, 0.4)),
        uvvis_bands_nm=(700.0, 780.0),
        chromophores=("Cu2+ d-d (malachite/turquoise)",),
        xrf_signature={"Cu": "major", "Si": "major", "Al": "minor"},
        libs_signature={"Cu": "major", "Si": "major"},
        confusables=("turquoise", "malachite", "amazonite"),
        diagnostic_features=("amorphous Cu-silicate + Cu2+ blue",
                             "low hardness, porous"),
        references=("Frost & Xi 2013",),
        is_amorphous=True,
        squid_ordering="paramagnetic",
        squid_weiss_K=-5.0,
        squid_susceptibility_si=(1.0e-5, 1.0e-4),
    ),

    # -------------------- Lapis & sodalite group --------------------
    "lapis_lazuli": MineralProfile(
        name="lapis_lazuli", species="lazurite",
        aliases=("lapis",),
        chemical_formula="(Na,Ca)8(AlSiO4)6(S,SO4,Cl)2",
        crystal_system="cubic", mohs_hardness=(5.0, 5.5),
        density_g_cc=(2.7, 2.9), refractive_index=(1.500, 1.520),
        common_colors=("blue",),
        raman_peaks_cm=_p((258.0, 0.3), (548.0, 1.0), (1096.0, 0.4)),
        uvvis_bands_nm=(600.0,),
        chromophores=("S3- radical (lapis/lazurite)",),
        xrf_signature={"Na": "major", "Al": "major", "Si": "major", "S": "minor", "Ca": "minor"},
        libs_signature={"Na": "major", "Al": "major", "Si": "major"},
        confusables=("sodalite", "azurite", "turquoise"),
        diagnostic_features=("548 cm-1 S3- radical Raman (diagnostic)",
                             "pyrite flecks + S3- blue"),
        references=("Osticioli et al. 2009",),
    ),
    "sodalite": MineralProfile(
        name="sodalite", species="sodalite",
        chemical_formula="Na8Al6Si6O24Cl2",
        crystal_system="cubic", mohs_hardness=(5.5, 6.0),
        density_g_cc=(2.27, 2.33), refractive_index=(1.483, 1.487),
        common_colors=("blue",),
        raman_peaks_cm=_p((258.0, 0.5), (466.0, 0.7), (989.0, 1.0)),
        xrf_signature={"Na": "major", "Al": "major", "Si": "major", "Cl": "minor"},
        libs_signature={"Na": "major", "Al": "major", "Si": "major"},
        confusables=("lapis_lazuli", "amazonite"),
        diagnostic_features=("989 cm-1 framework mode + Cl chemistry",
                             "lower density than lapis"),
        references=("Hassan & Grundy 1984",),
    ),

    # -------------------- Olivine, manganese & other gems --------------------
    "peridot": MineralProfile(
        name="peridot", species="olivine",
        aliases=("chrysolite", "olivine"),
        chemical_formula="(Mg,Fe)2SiO4",
        crystal_system="orthorhombic", mohs_hardness=(6.5, 7.0),
        density_g_cc=(3.27, 3.48), refractive_index=(1.650, 1.703),
        common_colors=("green", "yellow-green"),
        raman_peaks_cm=_p((545.0, 0.3), (823.0, 0.9), (856.0, 1.0), (920.0, 0.4)),
        xrf_signature={"Mg": "major", "Si": "major", "Fe": "minor"},
        libs_signature={"Mg": "major", "Si": "major", "Fe": "minor"},
        mossbauer_sites=("olivine_Fe2plus",),
        confusables=("prasiolite", "tsavorite", "hiddenite"),
        diagnostic_features=("olivine 823 + 856 cm-1 Raman doublet (diagnostic)",
                             "Mössbauer Fe2+ M1/M2 doublet"),
        references=("Kolesov & Geiger 2004", "Burns 1993"),
        squid_ordering="paramagnetic",
        squid_weiss_K=-5.0,
        squid_susceptibility_si=(5.0e-5, 3.0e-4),
    ),
    "rhodonite": MineralProfile(
        name="rhodonite", species="pyroxenoid",
        chemical_formula="MnSiO3",
        crystal_system="triclinic", mohs_hardness=(5.5, 6.5),
        density_g_cc=(3.57, 3.76), refractive_index=(1.711, 1.752),
        common_colors=("pink", "red"),
        raman_peaks_cm=_p((558.0, 0.4), (660.0, 1.0), (1000.0, 0.6)),
        uvvis_bands_nm=(410.0, 540.0),
        chromophores=("Mn2+ d-d (rhodochrosite, rhodonite)",),
        xrf_signature={"Mn": "major", "Si": "major", "Ca": "minor"},
        libs_signature={"Mn": "major", "Si": "major", "Ca": "minor"},
        confusables=("rhodochrosite", "rubellite"),
        diagnostic_features=("pyroxenoid 660 cm-1 + Mn2+ pink",
                             "Mn silicate (vs Mn carbonate rhodochrosite)"),
        references=("Burns 1993",),
        squid_ordering="paramagnetic",
        squid_weiss_K=-2.0,
        squid_susceptibility_si=(1.0e-5, 8.0e-5),
    ),
    "fluorite": MineralProfile(
        name="fluorite", species="fluorite",
        aliases=("fluorspar",),
        chemical_formula="CaF2",
        crystal_system="cubic", mohs_hardness=(4.0, 4.0),
        density_g_cc=(3.0, 3.25), refractive_index=(1.433, 1.435),
        common_colors=("purple", "green", "blue", "colorless"),
        raman_peaks_cm=_p((322.0, 1.0)),
        uvvis_bands_nm=(580.0,),
        xrf_signature={"Ca": "major"},
        libs_signature={"Ca": "major"},
        confusables=("rock_crystal", "cubic_zirconia", "amethyst"),
        diagnostic_features=("single 322 cm-1 T2g Raman line (diagnostic)",
                             "low hardness 4, colour-centre absorption"),
        references=("Krishnan 1947", "Bill & Calas 1978"),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
    ),
    "sphalerite": MineralProfile(
        name="sphalerite", species="sphalerite",
        aliases=("zinc blende",),
        chemical_formula="ZnS",
        crystal_system="cubic", mohs_hardness=(3.5, 4.0),
        density_g_cc=(3.9, 4.1), refractive_index=(2.369, 2.371),
        common_colors=("yellow", "orange", "red", "brown"),
        raman_peaks_cm=_p((275.0, 0.6), (350.0, 1.0)),
        xrf_signature={"Zn": "major", "S": "major", "Fe": "trace"},
        libs_signature={"Zn": "major", "Fe": "trace"},
        confusables=("andradite", "cubic_zirconia"),
        diagnostic_features=("ZnS LO/TO Raman 350/275 cm-1",
                             "extreme dispersion (fire), low hardness"),
        references=("Nilsen 1969",),
    ),
    "opal_precious": MineralProfile(
        name="opal_precious", species="opal",
        aliases=("precious opal",),
        chemical_formula="SiO2·nH2O",
        crystal_system="amorphous", mohs_hardness=(5.5, 6.5),
        density_g_cc=(1.98, 2.25), refractive_index=(1.37, 1.47),
        common_colors=("white", "play-of-colour"),
        raman_peaks_cm=_p((400.0, 1.0), (790.0, 0.3), (970.0, 0.2)),
        xrf_signature={"Si": "major"},
        confusables=("glass_paste", "obsidian"),
        diagnostic_features=("broad amorphous-silica Raman envelope",
                             "play-of-colour from silica-sphere diffraction"),
        references=("Smallwood et al. 1997",),
        is_amorphous=True,
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.5e-5, -1.2e-5),
    ),
    "fire_opal": MineralProfile(
        name="fire_opal", species="opal",
        chemical_formula="SiO2·nH2O:Fe",
        crystal_system="amorphous", mohs_hardness=(5.5, 6.5),
        density_g_cc=(1.98, 2.20), refractive_index=(1.37, 1.47),
        common_colors=("orange", "red", "yellow"),
        raman_peaks_cm=_p((400.0, 1.0), (790.0, 0.3)),
        uvvis_bands_nm=(430.0, 560.0),
        xrf_signature={"Si": "major", "Fe": "trace"},
        confusables=("opal_precious", "glass_paste"),
        diagnostic_features=("amorphous-silica Raman + Fe orange colour",
                             "play-of-colour weaker than precious opal"),
        references=("Fritsch et al. 1999",),
        is_amorphous=True,
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.4e-5, -1.1e-5),
    ),

    # -------------------- Quartz colour varieties (extension) --------------------
    "ametrine": MineralProfile(
        name="ametrine", species="quartz",
        aliases=("amethyst-citrine",),
        chemical_formula="SiO2:Fe (amethyst+citrine zoned)",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.64, 2.66), refractive_index=(1.544, 1.553),
        common_colors=("purple-yellow zoned",),
        raman_peaks_cm=_p((128.0, 0.4), (206.0, 0.5), (463.0, 1.0)),
        uvvis_bands_nm=(450.0, 545.0),
        xrf_signature={"Si": "major", "Fe": "trace"},
        epr_centers=("quartz_E1prime",),
        confusables=("amethyst", "citrine_natural", "rock_crystal"),
        diagnostic_features=("quartz 463 cm-1 + zoned amethyst (545) + citrine (450) colour",
                             "natural bicolour from differential Fe oxidation"),
        references=("Vasconcelos et al. 1994",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.5e-5, -1.3e-5),
    ),
    "milky_quartz": MineralProfile(
        name="milky_quartz", species="quartz",
        chemical_formula="SiO2 (fluid inclusions)",
        crystal_system="trigonal", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.64, 2.66), refractive_index=(1.544, 1.553),
        common_colors=("white", "milky"),
        raman_peaks_cm=_p((128.0, 0.4), (206.0, 0.5), (463.0, 1.0)),
        xrf_signature={"Si": "major"},
        confusables=("rock_crystal", "rose_quartz", "opal_precious"),
        diagnostic_features=("quartz 463 cm-1; cloudiness from fluid/gas inclusions"),
        references=("Etchepare et al. 1974",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.5e-5, -1.3e-5),
    ),

    # -------------------- Garnets (extension + promotions) --------------------
    "uvarovite": MineralProfile(
        name="uvarovite", species="garnet",
        chemical_formula="Ca3Cr2(SiO4)3",
        crystal_system="cubic", mohs_hardness=(6.5, 7.5),
        density_g_cc=(3.71, 3.81), refractive_index=(1.860, 1.870),
        common_colors=("emerald green",),
        raman_peaks_cm=_p((370.0, 0.5), (550.0, 0.8), (880.0, 0.7), (1010.0, 1.0)),
        uvvis_bands_nm=(430.0, 605.0),
        chromophores=("Cr3+ d-d (emerald/alexandrite)",),
        xrf_signature={"Ca": "major", "Cr": "major", "Si": "major"},
        libs_signature={"Ca": "major", "Cr": "major", "Si": "major"},
        cl_bands=(694.0,),
        confusables=("tsavorite", "demantoid", "grossular"),
        diagnostic_features=("Cr-major garnet — intense green",
                             "Cr3+ chromophore + CL red"),
        references=("Kolesov & Geiger 1998",),
        squid_ordering="paramagnetic",
        squid_weiss_K=-2.0,
        squid_susceptibility_si=(1.0e-5, 8.0e-5),
    ),
    "mali_garnet": MineralProfile(
        name="mali_garnet", species="garnet",
        aliases=("grandite",),
        chemical_formula="Ca3(Al,Fe)2(SiO4)3 (grossular-andradite)",
        crystal_system="cubic", mohs_hardness=(6.5, 7.0),
        density_g_cc=(3.65, 3.85), refractive_index=(1.760, 1.820),
        common_colors=("yellow-green", "brown"),
        raman_peaks_cm=_p((371.0, 0.5), (530.0, 0.8), (875.0, 0.7), (1000.0, 1.0)),
        xrf_signature={"Ca": "major", "Al": "major", "Fe": "major", "Si": "major"},
        libs_signature={"Ca": "major", "Al": "major", "Fe": "major", "Si": "major"},
        mossbauer_sites=("andradite_Fe3plus",),
        confusables=("grossular", "andradite", "hessonite"),
        diagnostic_features=("grossular-andradite solid solution (Ca+Al+Fe)",
                             "high dispersion, vivid green-yellow"),
        references=("Johnson et al. 1995",),
        squid_ordering="paramagnetic",
        squid_weiss_K=-3.0,
        squid_susceptibility_si=(2.0e-5, 1.0e-4),
    ),
    "hessonite": MineralProfile(
        name="hessonite", species="garnet",
        aliases=("cinnamon stone",),
        chemical_formula="Ca3Al2(SiO4)3 (Fe/Mn)",
        crystal_system="cubic", mohs_hardness=(6.5, 7.5),
        density_g_cc=(3.55, 3.68), refractive_index=(1.738, 1.745),
        common_colors=("orange", "cinnamon brown"),
        raman_peaks_cm=_p((372.0, 0.5), (549.0, 0.8), (879.0, 0.7), (1006.0, 1.0)),
        uvvis_bands_nm=(380.0, 450.0),
        chromophores=("Fe3+ spin-forbidden (yellow sapphire)",),
        xrf_signature={"Ca": "major", "Al": "major", "Si": "major", "Fe": "minor", "Mn": "trace"},
        libs_signature={"Ca": "major", "Al": "major", "Si": "major"},
        mossbauer_sites=("andradite_Fe3plus",),
        confusables=("grossular", "spessartine", "andradite"),
        diagnostic_features=("Fe/Mn grossular — cinnamon orange",
                             "roiled 'heat-wave' inclusions; Fe3+ absorption"),
        references=("Amthauer et al. 1976", "Kolesov & Geiger 1998"),
        squid_ordering="paramagnetic",
        squid_weiss_K=-3.0,
        squid_susceptibility_si=(1.0e-5, 5.0e-5),
    ),
    "demantoid": MineralProfile(
        name="demantoid", species="garnet",
        chemical_formula="Ca3Fe2(SiO4)3:Cr (andradite)",
        crystal_system="cubic", mohs_hardness=(6.5, 7.0),
        density_g_cc=(3.81, 3.87), refractive_index=(1.880, 1.889),
        common_colors=("green", "yellow-green"),
        raman_peaks_cm=_p((371.0, 0.5), (510.0, 0.8), (875.0, 0.7), (994.0, 1.0)),
        uvvis_bands_nm=(430.0, 605.0),
        chromophores=("Cr3+ d-d (emerald/alexandrite)",),
        xrf_signature={"Ca": "major", "Fe": "major", "Si": "major", "Cr": "trace"},
        libs_signature={"Ca": "major", "Fe": "major", "Si": "major", "Cr": "trace"},
        mossbauer_sites=("andradite_Fe3plus",),
        cl_bands=(694.0,),
        confusables=("andradite", "tsavorite", "uvarovite", "peridot"),
        diagnostic_features=("Cr-bearing andradite — vivid green, high dispersion (fire)",
                             "horsetail byssolite inclusions; Cr3+ chromophore"),
        references=("Schmetzer 2006", "Kolesov & Geiger 1998"),
        squid_ordering="paramagnetic",
        squid_weiss_K=-3.0,
        squid_susceptibility_si=(2.0e-5, 1.0e-4),
    ),

    # -------------------- Borosilicates --------------------
    "danburite": MineralProfile(
        name="danburite", species="danburite",
        chemical_formula="CaB2Si2O8",
        crystal_system="orthorhombic", mohs_hardness=(7.0, 7.0),
        density_g_cc=(2.97, 3.03), refractive_index=(1.630, 1.636),
        common_colors=("colorless", "yellow", "pink"),
        raman_peaks_cm=_p((620.0, 0.6), (785.0, 0.5), (1080.0, 1.0)),
        xrf_signature={"Ca": "major", "Si": "major", "B": "major"},
        libs_signature={"Ca": "major", "Si": "major", "B": "minor"},
        confusables=("goshenite", "white_topaz", "apatite"),
        diagnostic_features=("Ca-borosilicate; B seen by LIBS (invisible to XRF)",
                             "1080 cm-1 framework mode"),
        references=("Best et al. 1994",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
    ),
    "kornerupine": MineralProfile(
        name="kornerupine", species="kornerupine",
        aliases=("prismatine",),
        chemical_formula="(Mg,Fe)4Al6(Si,Al,B)5O21(OH)",
        crystal_system="orthorhombic", mohs_hardness=(6.5, 7.0),
        density_g_cc=(3.27, 3.45), refractive_index=(1.660, 1.680),
        common_colors=("green", "brown", "yellow"),
        raman_peaks_cm=_p((620.0, 0.6), (910.0, 0.7), (980.0, 1.0)),
        xrf_signature={"Mg": "major", "Al": "major", "Si": "major", "B": "minor"},
        libs_signature={"Mg": "major", "Al": "major", "Si": "major", "B": "minor"},
        confusables=("peridot", "tsavorite", "danburite"),
        diagnostic_features=("Mg-Al borosilicate; strong pleochroism",
                             "B via LIBS"),
        references=("Grew et al. 1996",),
    ),

    # -------------------- Accessory & phosphate gems --------------------
    "sphene": MineralProfile(
        name="sphene", species="titanite",
        aliases=("titanite",),
        chemical_formula="CaTiSiO5",
        crystal_system="monoclinic", mohs_hardness=(5.0, 5.5),
        density_g_cc=(3.48, 3.60), refractive_index=(1.843, 2.110),
        common_colors=("yellow", "green", "brown"),
        raman_peaks_cm=_p((547.0, 0.7), (605.0, 1.0), (855.0, 0.4)),
        uvvis_bands_nm=(430.0,),
        xrf_signature={"Ca": "major", "Ti": "major", "Si": "major"},
        libs_signature={"Ca": "major", "Ti": "major", "Si": "major"},
        confusables=("demantoid", "sphalerite", "zircon_high"),
        diagnostic_features=("547+605 cm-1 titanite doublet",
                             "extreme dispersion (0.051) — strong fire"),
        references=("Su et al. 2018",),
    ),
    "apatite": MineralProfile(
        name="apatite", species="apatite",
        chemical_formula="Ca5(PO4)3(F,OH,Cl)",
        crystal_system="hexagonal", mohs_hardness=(5.0, 5.0),
        density_g_cc=(3.16, 3.23), refractive_index=(1.628, 1.650),
        common_colors=("blue", "green", "yellow", "violet"),
        raman_peaks_cm=_p((430.0, 0.4), (590.0, 0.5), (962.0, 1.0), (1040.0, 0.4)),
        xrf_signature={"Ca": "major", "P": "major"},
        libs_signature={"Ca": "major", "P": "minor"},
        icpms_diagnostic_isotopes=("Sr88", "U238", "Pb206"),
        confusables=("danburite", "ivory", "blue_zircon"),
        diagnostic_features=("962 cm-1 PO4 ν1 (geological, vs biogenic ivory)",
                             "REE + Sr by ICP-MS; low hardness 5"),
        references=("Penel et al. 1998",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
    ),
    "zircon_high": MineralProfile(
        name="zircon_high", species="zircon",
        aliases=("high zircon",),
        chemical_formula="ZrSiO4 (crystalline)",
        crystal_system="tetragonal", mohs_hardness=(7.0, 7.5),
        density_g_cc=(4.60, 4.70), refractive_index=(1.92, 1.98),
        common_colors=("colorless", "yellow", "brown"),
        raman_peaks_cm=_p((357.0, 0.5), (438.0, 0.7), (974.0, 1.0), (1008.0, 0.5)),
        xrf_signature={"Zr": "major", "Si": "major", "Hf": "trace", "U": "trace"},
        icpms_diagnostic_isotopes=("U238", "Pb206", "Pb207", "Hf178"),
        cl_bands=(480.0, 575.0),
        confusables=("blue_zircon", "zircon_low"),
        diagnostic_features=("sharp 1008 cm-1 zircon mode (undamaged lattice)",
                             "U-Pb dateable; REE CL"),
        references=("Nasdala et al. 2002",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
    ),
    "zircon_low": MineralProfile(
        name="zircon_low", species="zircon",
        aliases=("low zircon", "metamict zircon"),
        chemical_formula="ZrSiO4 (radiation-damaged)",
        crystal_system="tetragonal", mohs_hardness=(6.0, 6.5),
        density_g_cc=(3.90, 4.10), refractive_index=(1.78, 1.84),
        common_colors=("green", "brown", "orange"),
        raman_peaks_cm=_p((974.0, 1.0),),
        xrf_signature={"Zr": "major", "Si": "major", "U": "minor"},
        icpms_diagnostic_isotopes=("U238", "Pb206", "Pb207"),
        confusables=("zircon_high", "blue_zircon"),
        diagnostic_features=("broadened/weak 974 cm-1 from metamictisation",
                             "lower RI/density than high zircon; higher U"),
        references=("Nasdala et al. 2002",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
    ),

    # -------------------- Manganese & rare ornamentals --------------------
    "sugilite": MineralProfile(
        name="sugilite", species="sugilite",
        aliases=("luvulite",),
        chemical_formula="KNa2(Fe,Mn,Al)2Li3Si12O30",
        crystal_system="hexagonal", mohs_hardness=(5.5, 6.5),
        density_g_cc=(2.74, 2.80), refractive_index=(1.607, 1.610),
        common_colors=("violet", "purple", "magenta"),
        raman_peaks_cm=_p((370.0, 0.5), (660.0, 0.8), (1010.0, 1.0)),
        uvvis_bands_nm=(550.0,),
        chromophores=("Mn3+ d-d (kunzite)",),
        xrf_signature={"K": "major", "Na": "major", "Si": "major", "Mn": "minor", "Fe": "minor"},
        libs_signature={"K": "major", "Na": "major", "Si": "major", "Mn": "minor"},
        confusables=("charoite", "rhodonite"),
        diagnostic_features=("Mn3+ violet cyclosilicate",
                             "Li (LIBS) + Mn chemistry"),
        references=("Shigley et al. 1987",),
        squid_ordering="paramagnetic",
        squid_weiss_K=-2.0,
        squid_susceptibility_si=(1.0e-5, 8.0e-5),
    ),
    "charoite": MineralProfile(
        name="charoite", species="charoite",
        chemical_formula="(K,Sr)(Ca,Na)2Si4O10(OH,F)·H2O",
        crystal_system="monoclinic", mohs_hardness=(5.0, 6.0),
        density_g_cc=(2.54, 2.78), refractive_index=(1.550, 1.559),
        common_colors=("lilac", "violet", "purple"),
        raman_peaks_cm=_p((540.0, 0.6), (1010.0, 1.0)),
        xrf_signature={"K": "major", "Ca": "major", "Na": "major", "Si": "major"},
        libs_signature={"K": "major", "Ca": "major", "Na": "major", "Si": "major"},
        confusables=("sugilite", "lapis_lazuli"),
        diagnostic_features=("fibrous swirling violet silicate",
                             "K-Ca-Na chemistry distinguishes from sugilite (Mn)"),
        references=("Rozhdestvenskaya et al. 2010",),
    ),

    # -------------------- Treated stones --------------------
    "jadeite_polymer": MineralProfile(
        name="jadeite_polymer", species="jadeite",
        aliases=("B-jade", "polymer-impregnated jadeite"),
        chemical_formula="NaAlSi2O6 + epoxy resin",
        crystal_system="monoclinic", mohs_hardness=(6.5, 7.0),
        density_g_cc=(3.25, 3.36), refractive_index=(1.652, 1.688),
        common_colors=("green", "lavender", "white"),
        raman_peaks_cm=_p((375.0, 0.5), (700.0, 0.8), (990.0, 1.0)),
        xrf_signature={"Na": "major", "Al": "major", "Si": "major"},
        ftir_bands=(2870.0, 2930.0),
        confusables=("jadeite", "nephrite", "serpentine"),
        diagnostic_features=("jadeite Raman + FTIR polymer C-H stretch (2870-2930)",
                             "bleached-and-impregnated B-jade treatment flag"),
        references=("Fritsch et al. 1992",),
        squid_ordering="diamagnetic",
        squid_susceptibility_si=(-1.0e-5, -8.0e-6),
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


def _profile_squid_susceptibility(profile: MineralProfile) -> float:
    """Return the midpoint of the profile's room-T volume susceptibility range."""
    lo, hi = profile.squid_susceptibility_si
    if lo == 0.0 and hi == 0.0:
        return 0.0
    return 0.5 * (lo + hi)


def synthesize_squid_mh(profile: MineralProfile, *,
                        fields_mT: np.ndarray | None = None,
                        temperature_K: float = 295.0,
                        noise: float = 0.005,
                        seed: int | None = 0):
    """Synthesise a dc-SQUID M(H) hysteresis loop for the profile.

    Returns ``None`` when the profile has no `squid_ordering` set.
    """
    if not profile.squid_ordering:
        return None
    from checkmsg.squid import simulate_mh
    chi = _profile_squid_susceptibility(profile)
    return simulate_mh(
        profile.squid_ordering,
        fields_mT=fields_mT,
        saturation_emu_g=profile.squid_saturation_emu_g,
        coercivity_mT=profile.squid_coercivity_mT,
        susceptibility_si=chi,
        temperature_K=temperature_K,
        noise=noise,
        seed=seed,
    )


def synthesize_squid_chi_T(profile: MineralProfile, *,
                           temperatures_K: np.ndarray | None = None,
                           applied_field_mT: float = 10.0,
                           noise: float = 0.005,
                           seed: int | None = 0):
    """Synthesise an rf-SQUID χ(T) sweep from a profile."""
    if not profile.squid_ordering:
        return None
    from checkmsg.squid import simulate_chi_T
    chi = _profile_squid_susceptibility(profile)
    return simulate_chi_T(
        profile.squid_ordering,
        temperatures_K=temperatures_K,
        curie_K=profile.squid_curie_K,
        neel_K=profile.squid_neel_K,
        weiss_K=profile.squid_weiss_K,
        susceptibility_si=chi,
        saturation_emu_g=profile.squid_saturation_emu_g,
        applied_field_mT=applied_field_mT,
        morin_K=profile.squid_morin_K,
        noise=noise,
        seed=seed,
    )


def synthesize_squid_chi_ac(profile: MineralProfile, *,
                            frequencies_Hz: np.ndarray | None = None,
                            temperature_K: float = 295.0,
                            blocking_temp_K: float = 100.0,
                            noise: float = 0.005,
                            seed: int | None = 0):
    """Synthesise an rf-SQUID AC susceptibility sweep from a profile.

    Uses the profile's room-T susceptibility as χ_T and 1% of it as χ_S
    (Casimir-du Pré high-frequency limit). Diamagnetic profiles produce a flat
    negative response without a relaxation peak.
    """
    if not profile.squid_ordering:
        return None
    from checkmsg.squid import simulate_chi_ac
    chi_T = abs(_profile_squid_susceptibility(profile)) or 1.0e-5
    chi_S = chi_T * 0.01
    return simulate_chi_ac(
        profile.squid_ordering,
        frequencies_Hz=frequencies_Hz,
        blocking_temp_K=blocking_temp_K,
        temperature_K=temperature_K,
        chi_T=chi_T,
        chi_S=chi_S,
        noise=noise,
        seed=seed,
    )


def synthesize_pl(profile: MineralProfile, *,
                  fields_nm: np.ndarray | None = None,
                  noise: float = 0.004, seed: int | None = 0,
                  line_fwhm_nm: float = 3.0) -> Spectrum | None:
    """Synthesise a PL emission spectrum from a profile's defect ZPLs.

    Returns ``None`` when the profile lists no PL centres.
    """
    if not profile.pl_centers:
        return None
    if fields_nm is None:
        fields_nm = np.linspace(400.0, 800.0, 601)
    peaks = [PeakSpec(position=c, intensity=1.0,
                      sigma=line_fwhm_nm / 2.355, gamma=line_fwhm_nm / 4)
             for c in profile.pl_centers]
    return generate(peaks, fields_nm, technique="pl", units="nm", noise=noise, seed=seed)


def synthesize_cl(profile: MineralProfile, *,
                  fields_nm: np.ndarray | None = None,
                  noise: float = 0.004, seed: int | None = 0,
                  band_fwhm_nm: float = 35.0) -> Spectrum | None:
    """Synthesise a CL emission spectrum from a profile's activator bands."""
    if not profile.cl_bands:
        return None
    if fields_nm is None:
        fields_nm = np.linspace(350.0, 750.0, 401)
    peaks = [PeakSpec(position=c, intensity=1.0,
                      sigma=band_fwhm_nm / 2.355, gamma=band_fwhm_nm / 4)
             for c in profile.cl_bands]
    return generate(peaks, fields_nm, technique="cl", units="nm", noise=noise, seed=seed)


def synthesize_ftir(profile: MineralProfile, *,
                    fields_cm: np.ndarray | None = None,
                    noise: float = 0.003, seed: int | None = 0,
                    band_fwhm_cm: float = 18.0) -> Spectrum | None:
    """Synthesise an FTIR absorption spectrum from a profile's bands / diamond type.

    Returns ``None`` when the profile has neither ``ftir_bands`` nor a
    ``diamond_type``.
    """
    from checkmsg.refdata.ftir_bands import DIAMOND_TYPE_BANDS
    centers: list[float] = list(profile.ftir_bands)
    if profile.diamond_type:
        centers.extend(DIAMOND_TYPE_BANDS.get(profile.diamond_type, ()))
    if not centers:
        return None
    if fields_cm is None:
        fields_cm = np.linspace(400.0, 4000.0, 1801)
    peaks = [PeakSpec(position=c, intensity=1.0,
                      sigma=band_fwhm_cm / 2.355, gamma=band_fwhm_cm / 4)
             for c in centers]
    return generate(peaks, fields_cm, technique="ftir", units="cm-1", noise=noise, seed=seed)


def synthesize_mossbauer(profile: MineralProfile, *,
                         noise: float = 0.004, seed: int | None = 0) -> Spectrum | None:
    """Synthesise a ⁵⁷Fe Mössbauer spectrum from the profile's dominant Fe site.

    Returns ``None`` when the profile lists no Mössbauer sites. Only the first
    (dominant) site is modelled so the doublet/sextet extraction is unambiguous.
    """
    if not profile.mossbauer_sites:
        return None
    from checkmsg.mossbauer import simulate_mossbauer
    from checkmsg.refdata.mossbauer_sites import SITES
    key = profile.mossbauer_sites[0]
    if key not in SITES:
        return None
    return simulate_mossbauer(SITES[key], noise=noise, seed=seed)
