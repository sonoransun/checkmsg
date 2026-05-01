import pytest

from checkmsg.microwave import (
    SUPPORTED_BANDS,
    MicrowaveBand,
    by_name,
    g_factor,
    resonance_field_mT,
)


def test_supported_bands_sorted_low_to_high():
    freqs = [b.frequency_GHz for b in SUPPORTED_BANDS]
    assert freqs == sorted(freqs)


def test_band_lookup_by_name():
    assert by_name("X").frequency_GHz == 9.5
    assert by_name("q").frequency_GHz == 35.0  # case-insensitive
    with pytest.raises(KeyError):
        by_name("Z")


def test_regime_labels():
    assert by_name("L").regime == "LF"
    assert by_name("X").regime == "X-class"
    assert by_name("W").regime == "high-field"


def test_resonance_field_known_value():
    # X-band, g=2 -> 339.4 mT
    x = by_name("X")
    assert abs(x.resonance_field_mT(2.0) - 339.38) < 0.05


def test_resonance_field_round_trip():
    band = by_name("Q")
    g_target = 2.0024
    field = band.resonance_field_mT(g_target)
    g_recovered = band.g_at(field)
    assert abs(g_recovered - g_target) < 1e-6


def test_standalone_helpers_match_class():
    assert abs(resonance_field_mT(9.5, 2.0) - by_name("X").resonance_field_mT(2.0)) < 1e-9
    assert abs(g_factor(339.0, 9.5) - by_name("X").g_at(339.0)) < 1e-9


def test_invalid_input_raises():
    with pytest.raises(ValueError):
        MicrowaveBand("bad", 0.0)
    with pytest.raises(ValueError):
        by_name("X").resonance_field_mT(0.0)
    with pytest.raises(ValueError):
        by_name("X").g_at(0.0)
