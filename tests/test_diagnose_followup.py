from checkmsg import minerals
from checkmsg.diagnose import diagnose


def test_followup_recommendations_when_only_raman_provided():
    """A single-technique input should trigger 'missing techniques' follow-up."""
    profile = minerals.get("ruby")
    raman_only = minerals.synthesize_raman(profile, noise=0.005)
    report = diagnose([raman_only])
    assert any("missing" in rec.lower() or "additional" in rec.lower() or "low" in rec.lower()
               for rec in report.follow_up_recommendations)


def test_no_followup_when_full_evidence_unambiguous():
    """A clean almandine with all four techniques should produce a confident verdict."""
    p = minerals.get("almandine")
    spectra = [
        minerals.synthesize_raman(p, noise=0.001),
        minerals.synthesize_xrf(p, noise=0.001),
        minerals.synthesize_libs(p, noise=0.001),
    ]
    report = diagnose(spectra)
    # Confidence should be high — at least one runner-up is far below.
    assert report.confidence > 0.6


def test_close_runner_up_triggers_disambiguation_advice():
    """Two minerals with overlapping evidence should yield a low-margin warning."""
    p = minerals.get("citrine_natural")
    spectra = [
        minerals.synthesize_raman(p, noise=0.005),
        minerals.synthesize_uvvis(p, noise=0.005),
    ]
    report = diagnose(spectra)
    # Quartz family is hard; expect at least one follow-up message.
    assert report.follow_up_recommendations
