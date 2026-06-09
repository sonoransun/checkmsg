"""Cathodoluminescence (CL) emission-band reference library.

CL excites luminescence with an electron beam, revealing activator centres and
(in imaging mode) growth zoning that separates natural from synthetic growth.
Positions are emission-band maxima in nm. Growth-zoning discrimination genuinely
needs imaging, not the 1-D spectrum modelled here — captured in the caveat.

Sources:
  - Marshall 1988, *Cathodoluminescence of Geological Materials* (Unwin Hyman).
  - Gaft, Reisfeld & Panczer 2005, *Modern Luminescence Spectroscopy of Minerals* (Springer).
  - Dean 1965, *Phys. Rev.* 139:A588 — diamond band-A donor-acceptor luminescence.
  - Nasdala, Zhang, Kempe et al. 2003, *Rev. Mineral. Geochem.* 53:427 — zircon CL.
"""

from __future__ import annotations

from checkmsg.refdata.line_bands import BandSet

CL_BANDS: tuple[BandSet, ...] = (
    BandSet("diamond band-A", (440.0,), 30.0, ("diamond",),
            "Broad blue donor-acceptor-pair luminescence; zoning pattern (octahedral vs "
            "cuboctahedral vs striated) separates natural / HPHT / CVD in imaging mode.",
            ("Dean 1965", "Marshall 1988")),
    BandSet("zircon REE3+ / broad", (480.0, 575.0), 25.0, ("blue_zircon",),
            "Broad yellow-green band + Dy3+/Sm3+ sharp REE lines.",
            ("Nasdala et al. 2003", "Gaft et al. 2005")),
    BandSet("Mn2+ carbonate orange", (610.0,), 25.0, ("pearl_natural_saltwater", "coral"),
            "Mn2+-activated orange luminescence in calcite/aragonite biominerals.",
            ("Marshall 1988",)),
    BandSet("Cr3+ corundum red", (694.0,), 8.0, ("ruby", "sapphire_blue"),
            "Cr3+ R-line red luminescence in corundum.",
            ("Gaft et al. 2005",)),
)
