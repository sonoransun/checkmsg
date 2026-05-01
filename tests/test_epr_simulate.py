import numpy as np
import pytest

from checkmsg.epr import (
    MU_B_OVER_H_MHZ_PER_MT,
    Hyperfine,
    SpinSystem,
    infer_g_factors,
    simulate_field_sweep,
)


@pytest.mark.parametrize("frequency_GHz", [1.1, 3.0, 5.0, 9.5, 24.0, 35.0, 94.0])
def test_free_electron_resonance_field_at_each_band(frequency_GHz):
    """A simple S=1/2 isotropic system resonates at hnu / (g mu_B)."""
    sys_ = SpinSystem(name="e", S=0.5, g=2.0023, linewidth_mT=0.3)
    expected = frequency_GHz * 1000.0 / (MU_B_OVER_H_MHZ_PER_MT * 2.0023)
    fields = np.linspace(expected - 5.0, expected + 5.0, 1001)
    spec = simulate_field_sweep(sys_, frequency_GHz=frequency_GHz, fields_mT=fields)
    gs = infer_g_factors(spec, frequency_GHz=frequency_GHz)
    assert gs, f"no resonance found at {frequency_GHz} GHz"
    assert abs(gs[0] - 2.0023) < 0.001


def test_hyperfine_sextet_spacing_matches_A():
    """55Mn-like (S=5/2, I=5/2) -> six allowed central transitions, equispaced by ~A/(g*mu_B)."""
    A_MHz = 250.0
    sys_ = SpinSystem(name="Mn", S=2.5, g=2.001, D_MHz=20.0,
                     hyperfine=(Hyperfine("55Mn", I=2.5, A_iso_MHz=A_MHz),),
                     linewidth_mT=0.4)
    fields = np.linspace(280, 400, 2001)
    spec = simulate_field_sweep(sys_, frequency_GHz=9.5, fields_mT=fields, orientations=(7, 1))
    gs = sorted(g for g in infer_g_factors(spec, 9.5) if 1.85 < g < 2.20)
    assert len(gs) >= 6
    fields_at = [9500.0 / (MU_B_OVER_H_MHZ_PER_MT * g) for g in gs]
    span = max(fields_at) - min(fields_at)
    avg_spacing_mT = span / (len(gs) - 1)
    expected = A_MHz / (MU_B_OVER_H_MHZ_PER_MT * 2.001)
    # Allow loose tolerance because shoulder satellites from D shift the spacing slightly.
    assert 0.7 * expected < avg_spacing_mT < 1.5 * expected


def test_cr3plus_corundum_powder_pattern_has_multiple_features():
    """Cr3+ in corundum (S=3/2, axial D=5.7 GHz) at X-band: powder pattern with several resonances."""
    sys_ = SpinSystem(name="Cr3+", S=1.5, g=(1.984, 1.984, 1.984),
                     D_MHz=5715.0, E_MHz=0.0, linewidth_mT=0.6)
    fields = np.linspace(50, 700, 1301)
    spec = simulate_field_sweep(sys_, frequency_GHz=9.5, fields_mT=fields, orientations=(15, 1))
    gs = infer_g_factors(spec, 9.5)
    # Expect multiple distinct features (parallel / perpendicular and inter-doublet transitions).
    assert len(gs) >= 3


def test_spectrum_metadata_records_frequency_and_system():
    sys_ = SpinSystem(name="probe", S=0.5, g=2.0)
    fields = np.linspace(335, 345, 401)
    spec = simulate_field_sweep(sys_, frequency_GHz=9.5, fields_mT=fields)
    assert spec.metadata["frequency_GHz"] == 9.5
    assert spec.metadata["spin_system"] == "probe"
    assert spec.technique == "epr"
    assert spec.units == "mT"


def test_derivative_flag_changes_output():
    sys_ = SpinSystem(name="probe", S=0.5, g=2.0)
    fields = np.linspace(335, 345, 401)
    deriv = simulate_field_sweep(sys_, 9.5, fields, derivative=True)
    absorp = simulate_field_sweep(sys_, 9.5, fields, derivative=False)
    # The two should differ; absorption is positive everywhere, derivative changes sign.
    assert np.min(absorp.intensity) >= -1e-3
    assert np.min(deriv.intensity) < 0 < np.max(deriv.intensity)


def test_rejects_non_epr_spectrum():
    from checkmsg.spectrum import Spectrum
    s = Spectrum(np.linspace(0, 1, 10), np.zeros(10), "raman", "cm-1")
    with pytest.raises(ValueError, match="expected epr"):
        infer_g_factors(s, 9.5)
