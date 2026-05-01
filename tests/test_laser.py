
import pytest

from checkmsg.laser import SUPPORTED_LASERS, LaserConfig, by_wavelength


def test_supported_lasers_sorted_uv_to_nir():
    assert list(SUPPORTED_LASERS) == sorted(SUPPORTED_LASERS)
    assert SUPPORTED_LASERS[0] == 275.0
    assert SUPPORTED_LASERS[-1] == 830.0


def test_regime_labels():
    assert LaserConfig(275).regime == "deep-UV"
    assert LaserConfig(325).regime == "UV"
    assert LaserConfig(405).regime == "violet"
    assert LaserConfig(457).regime == "blue"
    assert LaserConfig(488).regime == "green"
    assert LaserConfig(514).regime == "green"
    assert LaserConfig(633).regime == "red"
    assert LaserConfig(830).regime == "NIR"


def test_photon_energy_consistent():
    e_405 = LaserConfig(405).photon_eV
    assert abs(e_405 - 3.061) < 0.005


def test_lambda4_scale_decreases_with_wavelength():
    assert LaserConfig(275).lambda4_scale > LaserConfig(514).lambda4_scale
    assert LaserConfig(514).lambda4_scale > LaserConfig(830).lambda4_scale


def test_shift_to_wavelength_round_trip():
    laser = LaserConfig(514)
    scattered = laser.shift_to_wavelength(1332.0)
    recovered = laser.wavelength_to_shift(scattered)
    assert abs(recovered - 1332.0) < 0.01


def test_shift_to_wavelength_known_value():
    # Diamond 1332 cm-1 line with 532 nm laser scatters at 572.57 nm.
    laser = LaserConfig(532)
    assert abs(laser.shift_to_wavelength(1332.0) - 572.57) < 0.1


def test_resonance_enhancement_peaks_at_match():
    # Resonance with 555 nm Cr3+ band: 514 nm > 488 nm > 830 nm.
    e_514 = LaserConfig(514).resonance_enhancement(555.0, sigma_nm=30.0, max_enhancement=20.0)
    e_488 = LaserConfig(488).resonance_enhancement(555.0, sigma_nm=30.0, max_enhancement=20.0)
    e_830 = LaserConfig(830).resonance_enhancement(555.0, sigma_nm=30.0, max_enhancement=20.0)
    assert e_514 > e_488 > 1.0
    assert e_830 == pytest.approx(1.0, abs=1e-6)


def test_fluorescence_factor_table():
    # Visible (488/514) should be much larger than NIR (830) and deep-UV (275).
    assert LaserConfig(514).fluorescence_factor > 0.9
    assert LaserConfig(830).fluorescence_factor < 0.1
    assert LaserConfig(275).fluorescence_factor < 0.2


def test_fluorescence_factor_interpolates_between_table_entries():
    # 532 isn't in the table; should fall between 514 and 633.
    f = LaserConfig(532).fluorescence_factor
    assert 0.4 <= f <= 1.0


def test_by_wavelength_validates_against_catalogue():
    by_wavelength(514)
    with pytest.raises(ValueError, match="not a supported laser"):
        by_wavelength(532)


def test_above_bandgap():
    # Diamond bandgap ~5.45 eV. 275 nm = 4.51 eV -> below; 200 nm = 6.2 eV -> above.
    assert LaserConfig(275).above_bandgap(5.45) is False
    assert LaserConfig(200).above_bandgap(5.45) is True


def test_zero_or_negative_wavelength_raises():
    with pytest.raises(ValueError):
        LaserConfig(0)
    with pytest.raises(ValueError):
        LaserConfig(-100)
