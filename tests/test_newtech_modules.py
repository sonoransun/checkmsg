"""Analyze round-trip tests for the PL, FTIR, CL, and Mössbauer modules."""

from __future__ import annotations

import numpy as np

from checkmsg import cl, ftir, mossbauer, pl
from checkmsg.refdata.mossbauer_sites import SITES
from checkmsg.synthetic import PeakSpec, generate


def test_pl_identifies_nv_and_siv():
    axis = np.linspace(400.0, 800.0, 601)
    spec = generate([PeakSpec(637.0, 1.0, 1.2, 0.4), PeakSpec(736.6, 0.7, 1.2, 0.4)],
                    axis, technique="pl", units="nm", noise=0.004, seed=1)
    res = pl.analyze(spec)
    names = {b.name for b in res.defects()}
    assert "NV-" in names and "SiV-" in names
    assert res.has_synthetic_marker()


def test_ftir_diamond_type_ia():
    axis = np.linspace(400.0, 4000.0, 1801)
    spec = generate([PeakSpec(p, 1.0, 8.0, 4.0) for p in (1282.0, 1175.0, 1370.0)],
                    axis, technique="ftir", units="cm-1", noise=0.003, seed=1)
    res = ftir.analyze(spec)
    assert res.diamond_type == "Ia"


def test_ftir_polymer_flag():
    axis = np.linspace(400.0, 4000.0, 1801)
    spec = generate([PeakSpec(2870.0, 1.0, 12.0, 6.0), PeakSpec(2930.0, 0.9, 12.0, 6.0)],
                    axis, technique="ftir", units="cm-1", noise=0.003, seed=2)
    res = ftir.analyze(spec)
    assert res.has_polymer()


def test_cl_band_a():
    axis = np.linspace(350.0, 750.0, 401)
    spec = generate([PeakSpec(440.0, 1.0, 18.0, 7.0)],
                    axis, technique="cl", units="nm", noise=0.004, seed=1)
    res = cl.analyze(spec)
    assert any("band-A" in b.name for b in res.emitters())


def test_mossbauer_fe2_doublet():
    spec = mossbauer.simulate_mossbauer(SITES["almandine_Fe2plus"], seed=1)
    res = mossbauer.analyze(spec)
    assert not res.extracted["is_sextet"]
    assert res.extracted["valence"] == "Fe2+"
    assert abs(res.extracted["delta"] - 1.28) < 0.15
    assert res.best.name == "almandine_Fe2plus"


def test_mossbauer_fe3_doublet_valence():
    spec = mossbauer.simulate_mossbauer(SITES["corundum_Fe3plus"], seed=1)
    res = mossbauer.analyze(spec)
    assert res.extracted["valence"] == "Fe3+"  # δ≈0.37 → Fe3+


def test_mossbauer_sextet():
    spec = mossbauer.simulate_mossbauer(SITES["hematite_sextet"], seed=1)
    res = mossbauer.analyze(spec)
    assert res.extracted["is_sextet"]
    assert res.best.name == "hematite_sextet"
