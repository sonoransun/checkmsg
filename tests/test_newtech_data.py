"""Integrity tests for the new technique reference tables (PL/FTIR/CL/Mössbauer)."""

from __future__ import annotations

from checkmsg.refdata import cl_bands, ftir_bands, mossbauer_sites, pl_lines
from checkmsg.refdata.line_bands import BandSet, match_bands


def _check_bandsets(bandsets):
    assert bandsets
    for bs in bandsets:
        assert isinstance(bs, BandSet)
        assert bs.name and bs.centers and bs.tolerance > 0
        assert all(c > 0 for c in bs.centers)
        assert bs.references, f"{bs.name} missing a citation"


def test_pl_table():
    _check_bandsets(pl_lines.PL_CENTERS)
    names = {b.name for b in pl_lines.PL_CENTERS}
    assert {"NV-", "NV0", "SiV-", "N3 (cape series)"} <= names
    assert pl_lines.SYNTHETIC_MARKER_NM  # SiV marker present


def test_ftir_table():
    _check_bandsets(ftir_bands.FTIR_BANDS)
    assert set(ftir_bands.DIAMOND_TYPE_BANDS) == {"Ia", "Ib", "IIb"}


def test_cl_table():
    _check_bandsets(cl_bands.CL_BANDS)


def test_mossbauer_sites():
    sites = mossbauer_sites.SITES
    assert sites
    for site in sites.values():
        assert site.valence in ("Fe2+", "Fe3+")
        assert site.isomer_shift_mm_s > 0
        assert site.references
        if site.magnetic_sextet:
            assert site.hyperfine_field_T > 0
    assert mossbauer_sites.by_valence("Fe2+")
    assert mossbauer_sites.by_valence("Fe3+")


def test_match_bands_require_all_vs_any():
    bs = (BandSet("multi", (100.0, 200.0), 5.0, ()),)
    # require_all: both centres must be present
    assert match_bands([100.0], bs, require_all=True) == []
    assert len(match_bands([100.0, 199.0], bs, require_all=True)) == 2
    # require_all=False: any centre counts
    assert len(match_bands([100.0], bs, require_all=False)) == 1
