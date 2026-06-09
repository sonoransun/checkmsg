"""FTIR (Fourier-transform infrared) absorption-band reference library.

FTIR is the standard for diamond nitrogen-aggregation type (Ia/Ib/IIa/IIb) and
for detecting structural water in beryl and polymer impregnation in jade.
Positions are absorption-band wavenumbers in cm-1.

Sources:
  - Field 1992, *The Properties of Natural and Synthetic Diamond* (Academic Press).
  - Zaitsev 2001, *Optical Properties of Diamond* (Springer).
  - Woods, van Wyk & Collins 1990, *Phil. Mag. B* 62:589 — single-substitutional N (Ib).
  - Collins & Williams 1971, *J. Phys. C* 4:1789 — boron acceptor (IIb).
  - Wood & Nassau 1968, *Am. Mineral.* 53:777 — beryl OH / water.
  - Fritsch, Wu, Moses et al. 1992, *Gems & Gemology* 28:176 — polymer-impregnated (B-) jade.
  - Farmer 1974, *The Infrared Spectra of Minerals* (Mineralogical Society).
"""

from __future__ import annotations

from checkmsg.refdata.line_bands import BandSet

# Diamond IR type → its diagnostic band multiplet. Type IIa is defined by the
# ABSENCE of nitrogen/boron bands, so it carries no positive bandset.
DIAMOND_TYPE_BANDS: dict[str, tuple[float, ...]] = {
    "Ia": (1282.0, 1175.0, 1370.0),  # A-aggregate + platelet + B'
    "Ib": (1130.0, 1344.0),          # single substitutional N
    "IIb": (2800.0,),                # boron acceptor continuum
}

FTIR_BANDS: tuple[BandSet, ...] = (
    BandSet("diamond Type Ia (aggregated N)", (1282.0, 1175.0, 1370.0), 12.0, ("diamond",),
            "A-aggregate (N-pair) 1282 + platelet 1175 + B' 1370; most natural diamond.",
            ("Field 1992", "Zaitsev 2001")),
    BandSet("diamond Type Ib (single N)", (1130.0, 1344.0), 10.0, ("diamond",),
            "Single substitutional nitrogen; rare in nature, common in HPHT synthetics.",
            ("Woods et al. 1990", "Zaitsev 2001")),
    BandSet("diamond Type IIb (boron)", (2800.0,), 40.0, ("diamond",),
            "Uncompensated boron acceptor continuum; semiconducting blue diamond.",
            ("Collins & Williams 1971",)),
    BandSet("beryl structural water", (3600.0, 3700.0), 30.0, ("aquamarine", "emerald"),
            "Type I/II H2O + OH in beryl channels; 5270 combination band in extended range.",
            ("Wood & Nassau 1968",)),
    BandSet("polymer impregnation (C-H)", (2870.0, 2930.0), 25.0, ("jadeite",),
            "Epoxy/resin C-H stretch 2870-2970; flags B-jade treatment.",
            ("Fritsch et al. 1992",)),
    BandSet("amphibole/serpentine OH", (3675.0,), 20.0, ("nephrite", "serpentine"),
            "Tremolite-actinolite / serpentine OH stretch.",
            ("Farmer 1974",)),
)
