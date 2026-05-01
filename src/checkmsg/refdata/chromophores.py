"""UV-VIS chromophore band table for gemstone identification.

Bands are characteristic absorption features in transmission/absorbance spectra.
Sources: Burns 1993 (Mineralogical Applications of Crystal Field Theory),
Fritsch & Rossman 1987-88 (Gems & Gemology), Nassau 2001.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Chromophore:
    name: str
    bands_nm: tuple[float, ...]
    tolerance_nm: float
    typical_in: tuple[str, ...]
    notes: str


# Centers are approximate; tolerance covers crystal-field shifts across hosts.
CHROMOPHORES: tuple[Chromophore, ...] = (
    Chromophore(
        name="Cr3+ d-d (ruby/spinel)",
        bands_nm=(405, 555),
        tolerance_nm=20,
        typical_in=("ruby", "spinel"),
        notes="Cr3+ in corundum host: 4A2 -> 4T1, 4T2 transitions. Narrow R-line near 694 nm in ruby.",
    ),
    Chromophore(
        name="Cr3+ d-d (emerald/alexandrite)",
        bands_nm=(430, 605),
        tolerance_nm=22,
        typical_in=("emerald", "alexandrite"),
        notes="Cr3+ in beryl/chrysoberyl: smaller crystal field shifts bands to longer wavelength than ruby.",
    ),
    Chromophore(
        name="Fe2+/Ti4+ IVCT (blue sapphire)",
        bands_nm=(580,),
        tolerance_nm=30,
        typical_in=("sapphire-blue",),
        notes="Intervalence charge transfer giving blue color in corundum.",
    ),
    Chromophore(
        name="Fe3+ spin-forbidden (yellow sapphire)",
        bands_nm=(380, 450),
        tolerance_nm=20,
        typical_in=("sapphire-yellow", "citrine"),
        notes="UV/violet edges; weak in absorbance.",
    ),
    Chromophore(
        name="Fe2+ d-d (peridot)",
        bands_nm=(1050,),
        tolerance_nm=80,
        typical_in=("peridot",),
        notes="Near-IR; needs extended UV-VIS-NIR range.",
    ),
    Chromophore(
        name="Co2+ tetrahedral (blue spinel/glass)",
        bands_nm=(540, 580, 625),
        tolerance_nm=15,
        typical_in=("spinel-cobalt", "glass-cobalt"),
        notes="Triplet of sharp d-d bands diagnostic of Co.",
    ),
    Chromophore(
        name="V3+ d-d (tsavorite, V-emerald)",
        bands_nm=(430, 605),
        tolerance_nm=20,
        typical_in=("tsavorite", "emerald-vanadium"),
        notes="Similar to Cr but slightly shifted; LIBS V signal disambiguates.",
    ),
    Chromophore(
        name="Mn2+ d-d (rhodochrosite, rhodonite)",
        bands_nm=(410, 540),
        tolerance_nm=20,
        typical_in=("rhodochrosite", "rhodonite"),
        notes="Spin-forbidden, weak.",
    ),
    Chromophore(
        name="Diamond N3 center",
        bands_nm=(415,),
        tolerance_nm=5,
        typical_in=("diamond-cape",),
        notes="Sharp UV line in cape-series diamonds.",
    ),
    Chromophore(
        name="Moissanite UV cutoff",
        bands_nm=(425,),
        tolerance_nm=10,
        typical_in=("moissanite",),
        notes="6H-SiC band-gap edge near 3.0 eV.",
    ),
)


def assign(positions_nm: list[float], tolerance_scale: float = 1.0) -> list[tuple[float, Chromophore]]:
    """Match observed bands against the chromophore table.

    A chromophore is reported only when **all** of its expected bands are present
    in `positions_nm` within tolerance — this prevents single-band coincidences
    from being labelled as a multi-band fingerprint (e.g. one peak at ~430 nm
    triggering both Cr3+ and V3+ matches).
    """
    out: list[tuple[float, Chromophore]] = []
    for c in CHROMOPHORES:
        tol = c.tolerance_nm * tolerance_scale
        matched_bands: list[float] = []
        for expected in c.bands_nm:
            hits = [p for p in positions_nm if abs(p - expected) <= tol]
            if not hits:
                matched_bands = []
                break
            matched_bands.append(min(hits, key=lambda p: abs(p - expected)))
        if matched_bands:
            for pos in matched_bands:
                out.append((pos, c))
    return out
