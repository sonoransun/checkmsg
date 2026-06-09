"""Tests for the plain-language glossary."""

from __future__ import annotations

import pytest

from checkmsg import glossary as g

# Shorthand strings that genuinely appear in diagnose() output and therefore
# MUST resolve to a glossary entry (the core "no unexplained jargon" guarantee).
JARGON_IN_OUTPUT = [
    "Raman", "XRF", "LIBS", "UV-VIS", "EPR", "SQUID",
    "cm-1", "keV", "nm", "mT", "emu/g", "ppm",
    "Cr3+", "Fe2+", "Fe3+", "Ti4+", "V3+", "Mn2+", "Co2+",
    "IVCT", "d-d", "FWHM", "chromophore", "cosine",
    "diamagnetic", "paramagnetic", "ferromagnetic", "ferrimagnetic",
    "antiferromagnetic", "canted-afm",
]


@pytest.mark.parametrize("term", JARGON_IN_OUTPUT)
def test_known_jargon_is_defined(term):
    gt = g.define(term)
    assert gt is not None, f"{term!r} not in glossary"
    assert gt.description


def test_alias_resolution():
    assert g.define("epr").canonical == "EPR"
    assert g.define("uvvis").canonical == "UV-VIS"
    assert g.define("laicpms").canonical == "LA-ICP-MS"
    # EPR-centre key resolves to its descriptive name.
    assert g.define("diamond_P1") is not None
    assert g.define("diamond_P1").category == "center"


def test_case_insensitive_canonical_fallback():
    assert g.define("epr") is not None
    assert g.define("FWHM") is not None
    assert g.define("fwhm") is not None


def test_expand_inline_form():
    assert g.expand("EPR") == "EPR (electron paramagnetic resonance)"
    # Unknown terms pass through unchanged.
    assert g.expand("zzznope") == "zzznope"


def test_harvested_layers_present():
    cats = {t.category for t in g.all_terms()}
    assert {"technique", "unit", "ion", "concept", "magnetic",
            "chromophore", "center", "mineral"} <= cats
    # A specific harvested chromophore (description pulled from refdata notes).
    chro = g.define("Cr3+ d-d (ruby/spinel)")
    assert chro is not None and "corundum" in chro.description
    # A harvested mineral.
    assert g.define("ruby") is not None and g.define("ruby").category == "mineral"


def test_find_terms_and_expand_first_use():
    text = "dominant peak at 1332.5 cm-1 (FWHM 110.5)"
    found = {t.canonical for t in g.find_terms(text)}
    assert "cm-1" in found and "FWHM" in found
    seen: set[str] = set()
    out = g.expand_first_use(text, seen)
    assert "cm-1 (wavenumbers)" in out
    # Second pass with the same `seen` does not expand again.
    out2 = g.expand_first_use("another cm-1 value", seen)
    assert "(wavenumbers)" not in out2


def test_glossary_for_report_duck_typed():
    class _Ev:
        def __init__(self, obs):
            self.observation = obs

    class _Step:
        finding = "magnetic ordering: ferrimagnetic"
        implication = ""

    class _Rep:
        verdict = "magnetite"
        evidence = [_Ev("saturation moment 92.0 emu/g matches magnetite")]
        reasoning_trace = [_Step()]

    terms = {t.canonical for t in g.glossary_for(_Rep())}
    assert "ferrimagnetic" in terms
    assert "emu/g" in terms
