import numpy as np

from checkmsg import minerals


def test_synthesize_raman_dominant_peak_matches_catalog():
    """The synthesized Raman spectrum's tallest peak should be the catalog's #1 peak."""
    profile = minerals.get("ruby")
    spec = minerals.synthesize_raman(profile, noise=0.001)
    expected_top = max(profile.raman_peaks_cm, key=lambda pr: pr[1])[0]
    actual_top = float(spec.axis[int(np.argmax(spec.intensity))])
    assert abs(actual_top - expected_top) < 5.0


def test_synthesize_uvvis_bands_present():
    profile = minerals.get("ruby")
    spec = minerals.synthesize_uvvis(profile, noise=0.0)
    # Each catalog band should be near a local maximum
    for band in profile.uvvis_bands_nm:
        local = spec.intensity[np.argmin(np.abs(spec.axis - band))]
        assert local > 0.3


def test_synthesize_xrf_includes_majors():
    from checkmsg import xrf as xrf_mod
    profile = minerals.get("almandine")  # Fe-major garnet
    spec = minerals.synthesize_xrf(profile, noise=0.0)
    res = xrf_mod.identify_elements(spec, tolerance_keV=0.05, min_snr=8.0)
    detected = {e.element for e in res.elements}
    assert "Fe" in detected
    assert "Al" in detected
    assert "Si" in detected


def test_synthesize_libs_finds_diagnostic_elements():
    from checkmsg import libs as libs_mod
    profile = minerals.get("aquamarine")
    spec = minerals.synthesize_libs(profile, noise=0.0)
    res = libs_mod.identify(spec, tolerance_nm=0.4, min_snr=5.0)
    assert "Be" in res.elements


def test_synthesize_epr_returns_none_when_no_centers():
    profile = minerals.get("diamond")  # no EPR centers in profile
    out = minerals.synthesize_epr(profile)
    assert out is None


def test_synthesize_epr_produces_spectrum_when_centers_exist():
    profile = minerals.get("smoky_quartz")  # has E1' + Al-hole
    out = minerals.synthesize_epr(profile)
    assert out is not None
    assert out.technique == "epr"
    assert out.metadata.get("frequency_GHz") == 9.5


def test_synthesize_uvvis_for_chromophore_free_mineral_is_flat():
    profile = minerals.get("rock_crystal")  # no UV-VIS bands
    spec = minerals.synthesize_uvvis(profile, noise=0.0)
    # Just baseline and noise, max < 0.1 by construction (no peaks added)
    assert float(np.max(np.abs(spec.intensity))) < 0.5
