"""Photoluminescence (PL) defect-centre reference library.

PL probes sharp zero-phonon lines (ZPLs) of point defects in gem materials —
the modern frontier for separating natural, HPHT-synthetic, and CVD-synthetic
diamond. Positions are room-/low-temperature ZPL wavelengths in nm.

Sources:
  - Zaitsev 2001, *Optical Properties of Diamond* (Springer) — comprehensive ZPL catalogue.
  - Davies 1977, *Chem. Phys. Carbon* 13:1 — GR1, NV, H3/H4 centres.
  - Clark, Kanda, Kiflawi & Sittas 1995, *Phys. Rev. B* 51:16681 — Si-V centre (737 nm).
  - Wang, Moses, Linares et al. 2012, *Gems & Gemology* 48:80 — CVD-synthetic identification (NV/SiV).
"""

from __future__ import annotations

from checkmsg.refdata.line_bands import BandSet

# Si-V (≈737 nm) and the NV pair are the CVD/HPHT-synthetic markers.
SYNTHETIC_MARKER_NM: tuple[float, ...] = (736.6,)

PL_CENTERS: tuple[BandSet, ...] = (
    BandSet("N3 (cape series)", (415.2,), 3.0, ("diamond",),
            "N3 (three N + vacancy) ZPL; natural cape-series diamond marker.",
            ("Zaitsev 2001",)),
    BandSet("H4", (496.0,), 3.0, ("diamond",),
            "N4V2 aggregate centre; natural / annealed diamond.",
            ("Davies 1977", "Zaitsev 2001")),
    BandSet("H3 (N-V-N)", (503.2,), 3.0, ("diamond",),
            "Neutral N-V-N centre; annealing / treatment indicator.",
            ("Davies 1977", "Zaitsev 2001")),
    BandSet("3H", (503.5,), 3.0, ("diamond",),
            "Irradiation-induced interstitial centre.",
            ("Zaitsev 2001",)),
    BandSet("NV0", (575.0,), 3.0, ("diamond",),
            "Neutral nitrogen-vacancy ZPL; strong in synthetic diamond.",
            ("Davies 1977", "Zaitsev 2001")),
    BandSet("NV-", (637.0,), 3.0, ("diamond",),
            "Negative nitrogen-vacancy ZPL; dominant in HPHT/CVD synthetics.",
            ("Davies 1977", "Zaitsev 2001")),
    BandSet("SiV-", (736.6,), 3.0, ("diamond",),
            "Silicon-vacancy doublet ~736.6/736.9 nm; near-definitive CVD-synthetic marker.",
            ("Clark et al. 1995", "Wang et al. 2012")),
)
