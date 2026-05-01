import numpy as np

from checkmsg import laicpms
from checkmsg.synthetic import generate_laicpms_run


def test_blank_subtraction_removes_baseline():
    blank = {"Pb208": 100.0}
    sample = {"Pb208": 1100.0}
    run = generate_laicpms_run(sample, blank_cps=blank, sample_duration_s=15,
                               blank_duration_s=5, noise_factor=0.01, seed=1)
    cleaned = run.blank_subtracted_cps()
    assert 950 < cleaned["Pb208"] < 1050


def test_segment_detection_finds_two_plateaus():
    # Manually build a run with a clear surface->bulk step.
    t = np.linspace(0.0, 50.0, 5000)
    iso = laicpms.Isotope.lookup("Pb208")
    rng = np.random.default_rng(0)
    y = np.full_like(t, 1.0)
    y[(t >= 12) & (t < 22)] = 5000.0
    y[(t >= 22) & (t < 45)] = 100.0
    y = y + rng.normal(0.0, 0.05 * np.maximum(y, 1.0))
    y = np.clip(y, 0, None)
    transient = laicpms.IcpmsTransient(isotope=iso, time_s=t, intensity_cps=y)
    run = laicpms.IcpmsRun(
        transients={"Pb208": transient},
        blank_window_s=(0.5, 10.0),
        sample_window_s=(12.5, 44.5),
    )
    segs = run.detect_segments("Pb208", threshold_factor=6.0, min_gap_s=1.0,
                               merge_fraction=0.30)
    assert 2 <= len(segs) <= 4
    # First segment should be the surface plateau.
    s0_mean = run.integrate("Pb208", segs[0])
    s_last_mean = run.integrate("Pb208", segs[-1])
    assert s0_mean > 5 * s_last_mean


def test_to_spectrum_orders_by_mass():
    sample = {"Pb208": 100.0, "Pb204": 5.0, "Sr88": 50.0}
    run = generate_laicpms_run(sample, blank_cps={k: 0.0 for k in sample},
                               noise_factor=0.0, seed=0)
    spec = run.to_spectrum()
    assert list(spec.axis) == sorted(spec.axis)
    assert spec.technique == "laicpms"
    assert spec.units == "m/z"


def test_run_from_spectrum_round_trip():
    sample = {"Pb208": 100.0, "Sr88": 50.0}
    run1 = generate_laicpms_run(sample, blank_cps={k: 0.0 for k in sample},
                                noise_factor=0.0, seed=0)
    spec = run1.to_spectrum()
    run2 = laicpms.run_from_spectrum(spec)
    cps2 = run2.blank_subtracted_cps()
    assert "Pb208" in cps2 and "Sr88" in cps2
    # Values within 5% (small jitter from noise + integration windowing)
    assert abs(cps2["Pb208"] - 100.0) / 100.0 < 0.05


def test_integrate_with_explicit_window():
    t = np.linspace(0, 60, 6000)
    y = np.where(t < 30, 100.0, 200.0)
    iso = laicpms.Isotope.lookup("Pb208")
    transient = laicpms.IcpmsTransient(isotope=iso, time_s=t, intensity_cps=y)
    run = laicpms.IcpmsRun(transients={"Pb208": transient},
                          blank_window_s=(0, 5), sample_window_s=(40, 55))
    assert abs(run.integrate("Pb208", (0, 25)) - 100.0) < 1.0
    assert abs(run.integrate("Pb208", (35, 55)) - 200.0) < 1.0
