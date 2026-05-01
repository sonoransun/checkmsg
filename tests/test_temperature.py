
import pytest

from checkmsg.temperature import (
    LN2_K,
    ROOM_K,
    bose_einstein,
    fwhm_factor,
    infer_temperature,
    phonon_shift,
    stokes_antistokes_ratio,
)


def test_bose_einstein_decreases_with_freq():
    assert bose_einstein(100, 295) > bose_einstein(1000, 295)


def test_bose_einstein_increases_with_temperature():
    assert bose_einstein(417, 600) > bose_einstein(417, 295) > bose_einstein(417, 77)


def test_bose_einstein_zero_temp_raises():
    with pytest.raises(ValueError):
        bose_einstein(100, 0)


def test_stokes_antistokes_ratio_known_values():
    # 417 cm-1 mode at room temp: AS/S ~ 0.13
    r = stokes_antistokes_ratio(417.0, ROOM_K)
    assert 0.10 < r < 0.16
    # 1332 cm-1 at room temp: AS/S ~ 1.5e-3 (very small)
    r2 = stokes_antistokes_ratio(1332.0, ROOM_K)
    assert r2 < 0.005


def test_thermometry_round_trip():
    omega = 417.0
    T_truth = 250.0
    ratio = stokes_antistokes_ratio(omega, T_truth)
    recovered = infer_temperature(omega, ratio)
    assert abs(recovered - T_truth) < 0.1


def test_infer_temperature_rejects_invalid_ratio():
    with pytest.raises(ValueError):
        infer_temperature(417.0, 0.0)
    with pytest.raises(ValueError):
        infer_temperature(417.0, 1.5)


def test_fwhm_factor_increases_with_temperature():
    # Cooler -> narrower peaks (smaller factor).
    assert fwhm_factor(LN2_K) < fwhm_factor(ROOM_K)
    assert fwhm_factor(ROOM_K) == pytest.approx(1.0, abs=1e-6)
    assert fwhm_factor(600) > 1.0


def test_fwhm_factor_floor():
    # Very low T should clamp.
    assert fwhm_factor(1.0) == pytest.approx(0.30, abs=1e-6)


def test_phonon_shift_blue_at_cold():
    # Cooling should produce a positive (blue) shift.
    s = phonon_shift(417.0, LN2_K)
    assert s > 0
    # Heating should produce a red shift.
    s_warm = phonon_shift(417.0, 600)
    assert s_warm < 0


def test_phonon_shift_zero_at_reference():
    assert phonon_shift(1332.0, ROOM_K) == pytest.approx(0.0, abs=1e-9)
