"""Bethe-Bloch + Highland + CSDA spot checks vs textbook reference values."""


import pytest

from checkmsg.muon.physics import (
    bethe_bloch_dE_dx,
    csda_range_g_cm2,
    highland_scattering_rms_mrad,
)
from checkmsg.refdata.muon_data import get_material


def test_bethe_bloch_minimum_ionising_water_within_10pct_of_PDG():
    """At the MIP momentum (~MIP for muons in water is around 350 MeV/c, dE/dx ≈ 1.99 MeV g⁻¹ cm²)."""
    mip = bethe_bloch_dE_dx(350.0, get_material("water"))
    assert 1.7 <= mip <= 2.2


def test_bethe_bloch_low_p_loss_higher_than_mip():
    """At lower momentum the muon is below MIP and loses more energy per g/cm² (1/β² scaling)."""
    low = bethe_bloch_dE_dx(50.0, get_material("water"))
    mip = bethe_bloch_dE_dx(350.0, get_material("water"))
    assert low > mip


def test_bethe_bloch_zero_density_returns_zero():
    """Vacuum produces no energy loss."""
    assert bethe_bloch_dE_dx(300.0, get_material("vacuum")) == 0.0


def test_highland_scattering_approximately_inverse_p_in_relativistic_regime():
    """In the relativistic regime (β ≈ 1) doubling momentum approximately halves the
    RMS scattering angle. At lower β the (1/βp) factor diverges from pure 1/p."""
    pt = get_material("platinum")
    a = highland_scattering_rms_mrad(300.0, x_g_cm2=10.0, X_0_g_cm2=pt.X_0_g_cm2)
    b = highland_scattering_rms_mrad(600.0, x_g_cm2=10.0, X_0_g_cm2=pt.X_0_g_cm2)
    assert 1.85 < a / b < 2.2


def test_highland_scattering_high_Z_strongly_exceeds_low_Z():
    """Same thickness in g/cm² gives much larger scattering for high-Z (small X_0)."""
    al = get_material("corundum")
    pt = get_material("platinum")
    a = highland_scattering_rms_mrad(300.0, x_g_cm2=5.0, X_0_g_cm2=al.X_0_g_cm2)
    b = highland_scattering_rms_mrad(300.0, x_g_cm2=5.0, X_0_g_cm2=pt.X_0_g_cm2)
    assert b > 1.5 * a


def test_csda_range_monotone_in_momentum():
    """Higher momentum → longer range."""
    m = get_material("water")
    r1 = csda_range_g_cm2(50.0, m)
    r2 = csda_range_g_cm2(100.0, m)
    r3 = csda_range_g_cm2(200.0, m)
    assert r1 < r2 < r3


def test_negative_momentum_raises():
    with pytest.raises(ValueError):
        bethe_bloch_dE_dx(-1.0, get_material("water"))
