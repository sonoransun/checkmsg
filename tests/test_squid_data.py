"""Integrity tests for the bundled SQUID magnetic-mineral signature library."""

from __future__ import annotations

from checkmsg.refdata.squid_signatures import (
    ORDERINGS,
    SIGNATURES,
    by_ordering,
    diagnostic_for,
    diamagnetic_baselines,
)


def test_signatures_has_required_minerals():
    expected = {"magnetite", "hematite", "pyrrhotite_4c", "goethite", "ilmenite",
                "feNiCo_catalyst", "Cr3plus_paramagnet", "Mn2plus_paramagnet",
                "Fe2plus_paramagnet", "Fe3plus_paramagnet",
                "diamond_diamagnetic", "quartz_diamagnetic", "calcite_diamagnetic"}
    assert expected.issubset(SIGNATURES.keys())


def test_all_orderings_in_allowed_set():
    for sig in SIGNATURES.values():
        assert sig.ordering in ORDERINGS


def test_curie_or_neel_set_for_ordered_phases():
    for sig in SIGNATURES.values():
        if sig.ordering in ("ferromagnetic", "ferrimagnetic"):
            assert sig.curie_K > 0.0, f"{sig.name} missing curie_K"
        if sig.ordering in ("antiferromagnetic", "canted-afm"):
            assert sig.neel_K > 0.0, f"{sig.name} missing neel_K"


def test_diamagnetic_susceptibilities_negative():
    for sig in diamagnetic_baselines().values():
        lo, hi = sig.susceptibility_si
        assert hi < 0.0


def test_paramagnetic_susceptibilities_positive():
    for sig in by_ordering("paramagnetic").values():
        lo, hi = sig.susceptibility_si
        assert lo > 0.0
        assert hi > lo


def test_ferrimagnetic_have_high_saturation():
    for sig in by_ordering("ferrimagnetic").values():
        assert sig.saturation_emu_g >= 10.0


def test_morin_transition_only_on_canted_afm():
    for sig in SIGNATURES.values():
        if sig.morin_K > 0.0:
            assert sig.ordering == "canted-afm"


def test_diagnostic_for_returns_known_signatures():
    for scenario in ("black-opaque-discrimination", "hpht-diamond-screening",
                     "pearl-mn-quantitation", "schorl-vs-elbaite"):
        keys = diagnostic_for(scenario)
        assert keys, f"{scenario} returned no signatures"
        for k in keys:
            assert k in SIGNATURES


def test_minerals_catalog_squid_orderings_valid():
    """Every populated MineralProfile.squid_ordering matches the allowed set."""
    from checkmsg.minerals import CATALOG

    for name, p in CATALOG.items():
        if not p.squid_ordering:
            continue
        assert p.squid_ordering in ORDERINGS, (
            f"{name} has invalid squid_ordering={p.squid_ordering!r}"
        )


def test_minerals_catalog_has_squid_baseline_minerals():
    """Diamond, hematite, magnetite must carry SQUID data — used in HPHT screening."""
    from checkmsg.minerals import CATALOG

    for name in ("diamond", "magnetite", "hematite"):
        assert CATALOG[name].squid_ordering, f"{name} missing squid_ordering"
