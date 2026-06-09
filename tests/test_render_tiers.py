"""Tests for the three-tier DiagnosticReport rendering."""

from __future__ import annotations

import pytest

from checkmsg import glossary as g
from checkmsg import minerals
from checkmsg.diagnose import Tier, diagnose_profile


@pytest.fixture(scope="module")
def ruby_report():
    return diagnose_profile(minerals.get("ruby"))


@pytest.fixture(scope="module")
def magnetite_report():
    return diagnose_profile(minerals.get("magnetite"))


def test_novice_is_plain_and_short(ruby_report):
    text = ruby_report.render("novice")
    assert text.startswith("Most likely:")
    assert ruby_report.confidence_band in text
    # Novice tier must not leak expert scaffolding.
    assert "Reasoning trace" not in text
    assert "candidate scores" not in text.lower()
    assert "Disclaimer" not in text
    assert len(text.splitlines()) <= 4


def test_practitioner_explains_and_glossaries(ruby_report):
    text = ruby_report.render("practitioner")
    assert "Verdict:" in text
    assert "(corundum" in text  # verdict enriched with species
    assert "Why:" in text
    assert "Glossary:" in text
    # The raw chromophore key must be translated, not dumped verbatim as a label.
    assert "colour analysis" in text.lower()


def test_expert_is_backward_compatible(ruby_report):
    text = ruby_report.render("expert")
    assert "Verdict" in text
    assert "Reasoning trace" in text
    assert "Top-5 candidate scores" in text
    assert "Disclaimer" in text


def test_default_render_is_expert(ruby_report):
    assert ruby_report.render() == ruby_report.render("expert")
    assert ruby_report.render(Tier.EXPERT) == ruby_report.render("expert")


def test_every_emitted_term_resolves_in_glossary(ruby_report, magnetite_report):
    """No shorthand surfaced in any report goes undefined."""
    for report in (ruby_report, magnetite_report):
        for term in g.glossary_for(report):
            assert g.define(term.canonical) is not None
        # And the rendered practitioner glossary block only lists defined terms.
        lines = report.render("practitioner").splitlines()
        in_block = False
        for line in lines:
            if line.strip() == "Glossary:":
                in_block = True
                continue
            if in_block and " — " in line:
                head = line.split(" — ", 1)[0].strip()
                assert g.define(head) is not None, head


def test_magnetite_novice_mentions_band(magnetite_report):
    text = magnetite_report.render("novice")
    assert "Magnetite" in text
    assert magnetite_report.confidence_band in text


def test_invalid_tier_raises(ruby_report):
    with pytest.raises(ValueError):
        ruby_report.render("wizard")
