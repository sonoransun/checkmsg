"""Forward-model invariants for SQUID magnetometry simulators."""

from __future__ import annotations

import numpy as np
import pytest

from checkmsg import squid


def test_ferrimagnetic_hysteresis_has_remanence_and_coercivity():
    m = squid.simulate_mh("ferrimagnetic", saturation_emu_g=92.0, coercivity_mT=20.0,
                          noise=0.0)
    Hc = squid.extract_coercivity(m)
    Ms = squid.extract_saturation(m)
    Mr = squid.extract_remanence(m)
    assert abs(Hc - 20.0) < 2.0
    assert abs(Ms - 92.0) < 1.0
    # Remanence > 50% of saturation for soft-tanh hysteresis with Hc/Hk ~ 1.43
    assert Mr > 0.5 * Ms


def test_diamagnetic_loop_has_negative_slope_and_no_remanence():
    m = squid.simulate_mh("diamagnetic", susceptibility_si=-2.2e-5, noise=0.0)
    fields = m.axis
    sig = m.signal
    slope = np.polyfit(fields, sig, 1)[0]
    assert slope < 0.0
    Mr = squid.extract_remanence(m)
    assert Mr < 1e-3


def test_paramagnetic_curie_scaling_with_temperature():
    """Paramagnetic χ should scale ∝ 1/T, so M(H) at 100 K is ~3× M(H) at 295 K."""
    chi_si = 1.0e-4
    m_warm = squid.simulate_mh("paramagnetic", susceptibility_si=chi_si,
                                temperature_K=295.0, noise=0.0)
    m_cold = squid.simulate_mh("paramagnetic", susceptibility_si=chi_si,
                                temperature_K=100.0, noise=0.0)
    slope_warm = np.polyfit(m_warm.axis, m_warm.signal, 1)[0]
    slope_cold = np.polyfit(m_cold.axis, m_cold.signal, 1)[0]
    # Cold sample should have higher susceptibility by ~295/100 = 2.95
    assert slope_cold > 2.0 * slope_warm


def test_curie_weiss_obeyed_above_tc():
    """χ(T) = C / (T - θ) for a paramagnet — verify by linear fit of 1/χ vs T."""
    ct = squid.simulate_chi_T("paramagnetic", susceptibility_si=1.0e-4,
                              weiss_K=-3.0, noise=0.0)
    C, theta = squid.fit_curie_weiss(ct, T_min_K=50.0)
    assert C > 0.0
    assert abs(theta - (-3.0)) < 5.0  # generous tolerance for noise + grid


def test_chi_T_curie_temperature_recovered_for_magnetite():
    ct = squid.simulate_chi_T("ferrimagnetic", curie_K=858.0,
                              susceptibility_si=2.0e-3, saturation_emu_g=92.0,
                              noise=0.0)
    Tc = squid.extract_curie_temperature(ct)
    # Allow ±5% — the cusp is sharp but the dχ/dT minimum lives just above Tc.
    assert abs(Tc - 858.0) < 50.0


def test_chi_T_neel_temperature_recovered_for_hematite_above_morin():
    ct = squid.simulate_chi_T("canted-afm", neel_K=948.0,
                              susceptibility_si=2.0e-3, morin_K=263.0,
                              noise=0.0)
    TN = squid.extract_curie_temperature(ct)
    assert abs(TN - 948.0) < 50.0


def test_ac_susceptibility_loss_peak_at_1_over_tau():
    """χ''(ω) peaks at ωτ = 1, so loss-peak frequency f_peak = 1 / (2π·τ)."""
    tau_s = 1.0e-4
    ac = squid.simulate_chi_ac("paramagnetic", tau_s=tau_s, chi_T=1.0e-3, chi_S=1.0e-5,
                                noise=0.0)
    f_peak = squid.extract_loss_peak(ac)
    expected = 1.0 / (2.0 * np.pi * tau_s)
    assert 0.7 * expected < f_peak < 1.4 * expected


def test_diamagnetic_chi_T_is_flat_negative():
    ct = squid.simulate_chi_T("diamagnetic", susceptibility_si=-1.5e-5, noise=0.0)
    assert np.all(ct.signal < 0.0)
    assert np.std(ct.signal) < 1e-9


def test_chi_ac_signal_has_real_and_imag_packed():
    ac = squid.simulate_chi_ac("paramagnetic", tau_s=1e-3, chi_T=1e-3, chi_S=1e-5)
    assert ac.signal.size == 2 * ac.axis.size
    re, im = ac.chi_complex()
    assert re.shape == ac.axis.shape
    assert im.shape == ac.axis.shape


def test_simulate_mh_unknown_ordering_raises():
    with pytest.raises(ValueError):
        squid.simulate_mh("unicorn-magnet", saturation_emu_g=10.0)


def test_to_spectrum_round_trip_dc_mh():
    m = squid.simulate_mh("ferrimagnetic", saturation_emu_g=92.0, coercivity_mT=20.0)
    s = squid.to_spectrum(m)
    assert s.technique == "squid-mh"
    assert s.units == "mT"
    m2 = squid.from_spectrum(s)
    assert m2.mode == "dc-mh"
    np.testing.assert_array_equal(m.axis, m2.axis)
    np.testing.assert_array_equal(m.signal, m2.signal)


def test_to_spectrum_round_trip_chi_ac():
    ac = squid.simulate_chi_ac("paramagnetic", tau_s=1e-3, chi_T=1e-3, chi_S=1e-5)
    s = squid.to_spectrum(ac)
    assert s.technique == "squid-chi"
    assert s.units == "Hz"
    ac2 = squid.from_spectrum(s)
    assert ac2.mode == "rf-chi-ac"
    re1, im1 = ac.chi_complex()
    re2, im2 = ac2.chi_complex()
    np.testing.assert_array_equal(re1, re2)
    np.testing.assert_array_equal(im1, im2)


def test_chi_T_canted_afm_morin_drops_below_transition():
    """Hematite-style: spin-flop drops in-plane χ below Morin temperature."""
    ct = squid.simulate_chi_T("canted-afm", neel_K=948.0, susceptibility_si=2e-3,
                              morin_K=263.0, noise=0.0)
    T = ct.axis
    chi = ct.signal
    above_morin = chi[(T > 270) & (T < 290)].mean()
    below_morin = chi[(T > 230) & (T < 250)].mean()
    assert below_morin < 0.5 * above_morin


def test_invalid_signal_size_raises():
    with pytest.raises(ValueError):
        squid.SquidMeasurement(
            mode="dc-mh", axis=np.array([1.0, 2.0, 3.0]), signal=np.array([1.0, 2.0]),
        )
