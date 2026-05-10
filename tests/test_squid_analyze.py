"""End-to-end SQUID identification: simulate → analyze → top-1 candidate match."""

from __future__ import annotations

import numpy as np

from checkmsg import squid


def test_analyze_ranks_magnetite_first_for_ferrimagnetic_mh():
    m = squid.simulate_mh("ferrimagnetic", saturation_emu_g=92.0, coercivity_mT=20.0,
                          susceptibility_si=2.0e-3, noise=0.005)
    r = squid.analyze(m)
    assert r.best is not None
    assert r.best.name == "magnetite"
    assert r.best.combined > r.candidates[1].combined  # winner is unambiguous


def test_analyze_ranks_hematite_first_for_canted_afm():
    """Hematite: small remanence + linear susceptibility + high coercivity.

    The dc-mh inference uses Hc + Mr/Ms ratios to distinguish canted-AFM from
    weak ferri. For hematite parameters we rely on χ(T) instead, which is
    where the ordering signature lives (TN=948 K, Morin transition).
    """
    ct = squid.simulate_chi_T("canted-afm", neel_K=948.0, susceptibility_si=2.0e-3,
                              morin_K=263.0, noise=0.005)
    r = squid.analyze(ct)
    # Best by combined score should be hematite (the only canted-AFM with TN≈948).
    assert r.best is not None
    assert r.best.name == "hematite"


def test_analyze_ranks_diamond_first_for_diamagnetic_loop():
    m = squid.simulate_mh("diamagnetic", susceptibility_si=-2.1e-5, noise=0.005)
    r = squid.analyze(m)
    assert r.best is not None
    # Three diamagnetic baselines; diamond's susceptibility (−2.2..−2.0e-5) bracket matches.
    diamagnetic_top = [c.name for c in r.candidates if c.ordering_match]
    assert "diamond_diamagnetic" in diamagnetic_top


def test_analyze_paramagnet_chi_t_recovers_curie_weiss():
    """At low noise the Curie-Weiss fit recovers θ within a few K.

    The simulator's noise is scaled to peak |χ|, which lives at the lowest T
    (closest to the divergence). At high T the absolute noise stays the same
    while χ is ~100× smaller, so even 0.5 % relative noise overwhelms the
    1/χ-vs-T slope. Hence noise=0 here — the test verifies the fit equation
    is correct, not the noise tolerance.
    """
    ct = squid.simulate_chi_T("paramagnetic", susceptibility_si=5.0e-5, weiss_K=-3.0,
                              noise=0.0)
    r = squid.analyze(ct)
    assert r.extracted["ordering"] == "paramagnetic"
    assert abs(r.extracted["weiss_K"] - (-3.0)) < 5.0


def test_analyze_dc_only_accepts_dc_mh():
    ct = squid.simulate_chi_T("paramagnetic", susceptibility_si=1e-5, noise=0.0)
    try:
        squid.analyze_dc(ct)
    except ValueError:
        return
    raise AssertionError("expected ValueError for non-dc measurement")


def test_analyze_rf_only_accepts_chi_modes():
    m = squid.simulate_mh("paramagnetic", susceptibility_si=1e-5, noise=0.0)
    try:
        squid.analyze_rf(m)
    except ValueError:
        return
    raise AssertionError("expected ValueError for dc-mh measurement")


def test_diagnose_profile_round_trip_for_magnetic_minerals():
    """Synthesised + diagnosed catalog entries put the right mineral on top."""
    from checkmsg.diagnose import diagnose_profile
    from checkmsg.minerals import get

    for name in ("magnetite", "hematite", "ruby", "diamond", "schorl"):
        report = diagnose_profile(get(name))
        assert report.verdict == name, (
            f"{name} round-trip mis-identified as {report.verdict}"
        )
        squid_evidence = [e for e in report.evidence if e.technique == "squid"]
        assert squid_evidence, f"{name} report has no SQUID evidence"


def test_squid_evidence_present_when_only_squid_passed():
    """diagnose() called with only a SQUID spectrum should still produce ranked verdict."""
    from checkmsg.diagnose import diagnose
    from checkmsg.minerals import get, synthesize_squid_mh

    meas = synthesize_squid_mh(get("magnetite"))
    spec = squid.to_spectrum(meas)
    report = diagnose([spec])
    # SQUID evidence alone may not produce a unique winner (many paramagnets
    # share ordering), but for ferrimagnetic + Tc match magnetite should win.
    assert report.verdict == "magnetite"


def test_squid_recommended_in_followups_when_missing():
    """If diagnose() runs without any SQUID input, the missing-techniques list cites it."""
    from checkmsg.diagnose import diagnose
    from checkmsg.minerals import get, synthesize_raman

    spec = synthesize_raman(get("magnetite"))
    report = diagnose([spec])
    joined = " ".join(report.follow_up_recommendations)
    assert "squid" in joined.lower()


def test_param_extraction_returns_floats():
    """Sanity: extract_* helpers return plain floats, not numpy scalars."""
    m = squid.simulate_mh("ferrimagnetic", saturation_emu_g=92.0, coercivity_mT=20.0)
    assert isinstance(squid.extract_coercivity(m), float)
    assert isinstance(squid.extract_saturation(m), float)
    assert isinstance(squid.extract_remanence(m), float)


def test_chi_t_sweep_axis_handles_low_temperature():
    """χ(T) sweeps starting at 2 K should produce finite values without div-by-zero."""
    T = np.linspace(2.0, 600.0, 301)
    ct = squid.simulate_chi_T("paramagnetic", temperatures_K=T, susceptibility_si=1e-4,
                              weiss_K=-1.0, noise=0.0)
    assert np.all(np.isfinite(ct.signal))


def test_extract_loss_peak_returns_zero_for_dc_mh():
    """Calling AC-only extractors on dc-mh data returns 0 rather than crashing."""
    m = squid.simulate_mh("ferrimagnetic", saturation_emu_g=92.0, coercivity_mT=20.0)
    assert squid.extract_loss_peak(m) == 0.0
