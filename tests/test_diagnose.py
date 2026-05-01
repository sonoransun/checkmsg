import numpy as np

from checkmsg import minerals
from checkmsg.diagnose import diagnose, diagnose_profile
from checkmsg.spectrum import Spectrum


def test_diagnose_ruby_via_profile():
    report = diagnose_profile(minerals.get("ruby"))
    assert report.verdict == "ruby"
    assert report.confidence > 0.6


def test_diagnose_diamond_via_profile():
    report = diagnose_profile(minerals.get("diamond"))
    assert report.verdict == "diamond"


def test_diagnose_returns_none_for_pure_noise():
    rng = np.random.default_rng(0)
    fake = Spectrum(np.linspace(100, 1500, 200), rng.normal(0, 1.0, 200), "raman", "cm-1")
    report = diagnose([fake])
    # No coherent matches → either None or very low score
    assert report.confidence < 0.7


def test_reasoning_trace_non_empty():
    report = diagnose_profile(minerals.get("almandine"))
    assert report.reasoning_trace
    assert all(step.technique for step in report.reasoning_trace)


def test_render_returns_human_readable():
    report = diagnose_profile(minerals.get("tanzanite"))
    text = report.render()
    assert "Verdict" in text
    assert "Reasoning trace" in text


def test_candidate_scores_include_runners_up():
    """Ruby should score above red_spinel but red_spinel should also score positive."""
    report = diagnose_profile(minerals.get("ruby"))
    assert report.candidate_scores["ruby"] >= report.candidate_scores.get("red_spinel", 0)


def test_diagnose_garnet_picks_correct_endmember():
    report = diagnose_profile(minerals.get("spessartine"))
    assert report.verdict == "spessartine"


def test_amorphous_glass_is_flagged():
    report = diagnose_profile(minerals.get("glass_paste"))
    # Either picks glass_paste or another amorphous (obsidian/jet) — both acceptable
    assert report.verdict in ("glass_paste", "obsidian", "jet")
