"""Tests for laser/temperature wiring in synthetic.generate."""

import numpy as np

from checkmsg import raman
from checkmsg.peaks import fit_voigt
from checkmsg.synthetic import PeakSpec, generate
from checkmsg.temperature import LN2_K, ROOM_K


def _ruby_peak(resonant: bool = True):
    band = 555.0 if resonant else None
    return [PeakSpec(417.0, intensity=1.0, sigma=2.5, gamma=0.8,
                    resonance_band_nm=band, resonance_max=20.0)]


def test_metadata_records_laser_and_temperature():
    s = generate(_ruby_peak(False), np.linspace(100, 1100, 1001),
                 "raman", "cm-1", noise=0.0, seed=0,
                 laser_nm=514, temperature_K=LN2_K)
    assert s.metadata["laser_nm"] == 514.0
    assert s.metadata["temperature_K"] == LN2_K


def test_resonance_enhances_peak_at_514_vs_830():
    axis = np.linspace(300, 700, 800)
    s_514 = generate(_ruby_peak(True), axis, "raman", "cm-1", noise=0.0, seed=0,
                     laser_nm=514, temperature_K=ROOM_K, fluorescence_amplitude=0.0)
    s_830 = generate(_ruby_peak(True), axis, "raman", "cm-1", noise=0.0, seed=0,
                     laser_nm=830, temperature_K=ROOM_K, fluorescence_amplitude=0.0)
    h_514 = float(np.max(s_514.intensity))
    h_830 = float(np.max(s_830.intensity))
    # Resonance + lambda^4 should make 514 dominate.
    assert h_514 > 5 * h_830


def test_lambda4_scaling_without_resonance():
    axis = np.linspace(100, 1100, 1001)
    s_275 = generate(_ruby_peak(False), axis, "raman", "cm-1", noise=0.0, seed=0,
                     laser_nm=275, temperature_K=ROOM_K, fluorescence_amplitude=0.0)
    s_830 = generate(_ruby_peak(False), axis, "raman", "cm-1", noise=0.0, seed=0,
                     laser_nm=830, temperature_K=ROOM_K, fluorescence_amplitude=0.0)
    # Expect (830/275)^4 ~ 82x
    ratio = float(np.max(s_275.intensity) / np.max(s_830.intensity))
    assert 60 < ratio < 110


def test_fluorescence_baseline_higher_at_488_than_830():
    axis = np.linspace(100, 1100, 1001)
    # No peaks: only fluorescence baseline.
    s_488 = generate([], axis, "raman", "cm-1", noise=0.0, seed=0,
                     laser_nm=488, temperature_K=ROOM_K, fluorescence_amplitude=1.0)
    s_830 = generate([], axis, "raman", "cm-1", noise=0.0, seed=0,
                     laser_nm=830, temperature_K=ROOM_K, fluorescence_amplitude=1.0)
    assert s_488.intensity.mean() > 5 * s_830.intensity.mean()


def test_lower_temperature_narrows_peaks():
    axis = np.linspace(380, 460, 800)
    s_rt = generate(_ruby_peak(False), axis, "raman", "cm-1", noise=0.0005, seed=0,
                    temperature_K=ROOM_K)
    s_ln = generate(_ruby_peak(False), axis, "raman", "cm-1", noise=0.0005, seed=0,
                    temperature_K=LN2_K)
    fit_rt = fit_voigt(s_rt, around=420.0, window=20.0)
    fit_ln = fit_voigt(s_ln, around=420.0, window=20.0)
    assert fit_ln.width < fit_rt.width


def test_lower_temperature_blue_shifts_peaks():
    axis = np.linspace(380, 460, 800)
    s_rt = generate(_ruby_peak(False), axis, "raman", "cm-1", noise=0.0005, seed=0,
                    temperature_K=ROOM_K)
    s_ln = generate(_ruby_peak(False), axis, "raman", "cm-1", noise=0.0005, seed=0,
                    temperature_K=LN2_K)
    fit_rt = fit_voigt(s_rt, around=420.0, window=20.0)
    fit_ln = fit_voigt(s_ln, around=420.0, window=20.0)
    assert fit_ln.position > fit_rt.position


def test_fingerprint_position_invariant_across_lasers():
    """Raman shift is laser-invariant; the strongest peak lands at the same cm-1."""
    axis = np.linspace(380, 460, 800)
    positions = []
    for laser in (275, 405, 514, 830):
        s = generate(_ruby_peak(False), axis, "raman", "cm-1", noise=0.0005, seed=0,
                     laser_nm=laser, temperature_K=ROOM_K, fluorescence_amplitude=0.0)
        fit = fit_voigt(s, around=420.0, window=20.0)
        positions.append(fit.position)
    positions = np.asarray(positions)
    assert positions.max() - positions.min() < 0.5


def test_antistokes_present_when_requested():
    axis = np.linspace(-700, 700, 1601)
    s = generate(_ruby_peak(False), axis, "raman", "cm-1", noise=0.0005, seed=0,
                 temperature_K=ROOM_K, include_antistokes=True)
    # Should have a Stokes peak near +417 and a (smaller) AS peak near -417.
    pos_part = s.intensity[s.axis > 0]
    neg_part = s.intensity[s.axis < 0]
    assert pos_part.max() > 5 * neg_part.max()  # Stokes dominant
    assert neg_part.max() > 5 * neg_part.mean() + 1e-6  # AS clearly above background


def test_thermometry_recovers_room_temperature():
    axis = np.linspace(-700, 700, 1601)
    s = generate(_ruby_peak(False), axis, "raman", "cm-1", noise=0.0005, seed=0,
                 temperature_K=ROOM_K, fluorescence_amplitude=0.0,
                 include_antistokes=True)
    T = raman.infer_temperature(s, mode_cm=417.0, window_cm=25.0)
    assert abs(T - ROOM_K) < 30


def test_thermometry_recovers_ln2_temperature():
    axis = np.linspace(-700, 700, 1601)
    s = generate(_ruby_peak(False), axis, "raman", "cm-1", noise=0.0001, seed=0,
                 temperature_K=LN2_K, fluorescence_amplitude=0.0,
                 include_antistokes=True)
    T = raman.infer_temperature(s, mode_cm=417.0, window_cm=25.0)
    assert abs(T - LN2_K) < 30
