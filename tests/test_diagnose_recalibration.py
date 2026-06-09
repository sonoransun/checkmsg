"""Regression guards for the scoring recalibration.

The recalibration refactored the inline weight literals into ``scoring.WEIGHTS``
and added qualitative confidence bands + conflict detection. None of that may
change which verdict ``diagnose`` produces. This module freezes the catalog-wide
verdicts as a baseline so any future weight edit that shifts a verdict is caught.
"""

from __future__ import annotations

from checkmsg import minerals, scoring
from checkmsg.diagnose import diagnose_profile

# Frozen baseline: diagnose_profile(get(name)).verdict for every catalog entry.
# Entries that resolve to a *different* name (e.g. cymophane -> alexandrite,
# smoky_quartz -> amethyst) are the documented hard/confusable cases; the test
# asserts they remain STABLE, not that they are gemologically "correct".
BASELINE_VERDICTS = {
    "GGG": "GGG", "YAG": "YAG", "alexandrite": "alexandrite",
    "almandine": "almandine", "amazonite": "amazonite", "amethyst": "amethyst",
    "ametrine": "amethyst", "andradite": "andradite", "apatite": "apatite",
    "aquamarine": "aquamarine", "aventurine_quartz": "aventurine_quartz", "azurite": "azurite",
    "blue_topaz": "blue_topaz", "blue_zircon": "zircon_high", "charoite": "charoite",
    "chrysoberyl_yellow": "chrysoberyl_yellow", "chrysocolla": "chrysocolla", "citrine_heat_treated": "citrine_heat_treated",
    "citrine_natural": "citrine_natural", "coral": "coral", "cubic_zirconia": "cubic_zirconia",
    "cymophane": "cymophane", "danburite": "danburite", "demantoid": "demantoid",
    "diamond": "diamond", "diamond_cvd": "diamond_cvd", "diamond_hpht": "diamond_hpht",
    "dravite": "dravite", "elbaite": "elbaite", "emerald": "emerald",
    "fire_opal": "glass_paste", "fluorite": "fluorite", "glass_paste": "glass_paste",
    "goshenite": "morganite", "grossular": "grossular", "heliodor": "heliodor",
    "hematite": "hematite", "hessonite": "hessonite", "hiddenite": "hiddenite",
    "iolite": "iolite", "ivory": "ivory", "jadeite": "jadeite",
    "jadeite_polymer": "jadeite_polymer", "jet": "glass_paste", "kornerupine": "kornerupine",
    "kunzite": "kunzite", "labradorite": "labradorite", "lapis_lazuli": "lapis_lazuli",
    "liddicoatite": "liddicoatite", "magnetite": "magnetite", "malachite": "malachite",
    "mali_garnet": "mali_garnet", "milky_quartz": "aventurine_quartz", "moissanite": "moissanite",
    "morganite": "morganite", "nephrite": "nephrite", "obsidian": "obsidian",
    "onyx": "aventurine_quartz", "opal_precious": "glass_paste", "orthoclase": "orthoclase",
    "pearl_akoya": "pearl_natural_saltwater", "pearl_freshwater": "pearl_natural_saltwater", "pearl_natural_saltwater": "pearl_natural_saltwater",
    "peridot": "peridot", "prasiolite": "amethyst", "prehnite": "prehnite",
    "pyrope": "pyrope", "red_beryl": "red_beryl", "red_spinel": "red_spinel",
    "rhodochrosite": "rhodochrosite", "rhodolite": "rhodolite", "rhodonite": "rhodonite",
    "rock_crystal": "aventurine_quartz", "rose_quartz": "aventurine_quartz", "rubellite": "rubellite",
    "ruby": "ruby", "sapphire_blue": "sapphire_blue", "schorl": "schorl",
    "serpentine": "serpentine", "smoky_quartz": "amethyst", "sodalite": "sodalite",
    "spessartine": "spessartine", "sphalerite": "sphalerite", "sphene": "sphene",
    "strontium_titanate": "strontium_titanate", "sugilite": "sugilite", "sunstone": "sunstone",
    "tanzanite": "tanzanite", "tsavorite": "tsavorite", "turquoise": "turquoise",
    "uvarovite": "uvarovite", "white_sapphire": "white_sapphire", "white_spinel": "white_spinel",
    "white_topaz": "white_topaz", "zircon_high": "zircon_high", "zircon_low": "zircon_low",
}


def test_catalog_verdicts_are_stable():
    """Every catalog entry round-trips to its baseline verdict (no drift)."""
    drift = {}
    for name in minerals.names():
        report = diagnose_profile(minerals.get(name))
        expected = BASELINE_VERDICTS.get(name)
        if report.verdict != expected:
            drift[name] = (expected, report.verdict)
    assert not drift, f"verdict drift detected: {drift}"


def test_baseline_covers_whole_catalog():
    assert set(BASELINE_VERDICTS) == set(minerals.names())


def test_reports_carry_a_valid_band_and_caveats():
    for name in ("ruby", "diamond", "magnetite"):
        report = diagnose_profile(minerals.get(name))
        assert report.confidence_band in scoring.BANDS
        assert report.verdict is not None
        assert report.caveats, f"{name}: expected non-empty caveats"
        assert report.evidence_agreement >= 1
