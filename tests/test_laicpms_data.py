import pytest

from checkmsg.refdata.icpms_data import (
    CHONDRITE_REE_PPM,
    ISOTOPES,
    LAMBDA_235,
    LAMBDA_238,
    NIST_SRM_612,
    REE_ELEMENTS,
    STACEY_KRAMERS_PB,
    U238_OVER_U235,
    isotope,
    isotopes_of,
)


def test_lookup_pb208():
    rec = isotope("Pb208")
    assert rec.element == "Pb"
    assert rec.mass == 208
    assert 0.5 < rec.natural_abundance < 0.6


def test_unknown_isotope_raises():
    with pytest.raises(KeyError):
        isotope("Pb999")


def test_isotopes_of_pb_returns_four_isotopes():
    pb_isotopes = isotopes_of("Pb")
    assert set(pb_isotopes) == {"Pb204", "Pb206", "Pb207", "Pb208"}
    assert sum(rec.natural_abundance for rec in pb_isotopes.values()) == pytest.approx(1.0, abs=0.01)


def test_uranium_isotope_ratio_consistent_with_constant():
    u235 = isotope("U235")
    u238 = isotope("U238")
    measured_ratio = u238.natural_abundance / u235.natural_abundance
    assert abs(measured_ratio - U238_OVER_U235) / U238_OVER_U235 < 0.05


def test_decay_constants_have_expected_magnitudes():
    # Steiger & Jäger 1977 values to 4 sig figs
    assert abs(LAMBDA_238 - 1.55125e-10) < 1e-15
    assert abs(LAMBDA_235 - 9.8485e-10) < 1e-14


def test_nist_612_has_pb_th_u_at_typical_levels():
    # Pearce 1997 preferred values: ~38 ppm for trace elements
    assert 30.0 < NIST_SRM_612["Pb"] < 50.0
    assert 30.0 < NIST_SRM_612["Th"] < 50.0
    assert 30.0 < NIST_SRM_612["U"] < 50.0


def test_chondrite_ree_keys_match_ree_elements():
    assert set(CHONDRITE_REE_PPM.keys()) == set(REE_ELEMENTS)
    # All values strictly positive
    assert all(v > 0 for v in CHONDRITE_REE_PPM.values())


def test_stacey_kramers_present_day_terrestrial_pb():
    # Stacey & Kramers 1975 model yields ²⁰⁶/²⁰⁴ ≈ 18.7 today
    assert 18.5 < STACEY_KRAMERS_PB["206/204"] < 18.9
    # Self-consistency between ratios
    assert abs(STACEY_KRAMERS_PB["207/206"] -
               STACEY_KRAMERS_PB["207/204"] / STACEY_KRAMERS_PB["206/204"]) < 1e-9


def test_isotope_table_size_reasonable():
    # We bundle ~130 gem-relevant isotopes, not the full periodic table.
    assert 100 < len(ISOTOPES) < 200
