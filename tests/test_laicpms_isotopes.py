import math

from checkmsg import laicpms
from checkmsg.synthetic import generate_laicpms_run


def _pb_run(r206_204, r207_204, r208_204, total_pb_cps=10000.0, seed=0):
    """Build a Pb-only run with prescribed isotope ratios."""
    a204_frac = 1.0 / (1.0 + r206_204 + r207_204 + r208_204)
    signals = {
        "Pb204": total_pb_cps * a204_frac,
        "Pb206": total_pb_cps * a204_frac * r206_204,
        "Pb207": total_pb_cps * a204_frac * r207_204,
        "Pb208": total_pb_cps * a204_frac * r208_204,
    }
    return generate_laicpms_run(signals, blank_cps={k: 0.5 for k in signals},
                                sample_duration_s=15, blank_duration_s=5,
                                noise_factor=0.02, seed=seed)


def test_pb_ratios_recover_known_values():
    run = _pb_run(r206_204=18.70, r207_204=15.63, r208_204=38.65, seed=1)
    ratios = laicpms.pb_ratios(run)
    assert abs(ratios["206/204"] - 18.70) / 18.70 < 0.02
    assert abs(ratios["207/204"] - 15.63) / 15.63 < 0.02
    assert abs(ratios["208/204"] - 38.65) / 38.65 < 0.02
    assert abs(ratios["207/206"] - 15.63 / 18.70) / (15.63 / 18.70) < 0.02


def test_pb_ratios_empty_when_no_pb_channels():
    run = generate_laicpms_run({"Mn55": 1000.0}, seed=2)
    assert laicpms.pb_ratios(run) == {}


def test_sr_ratio_with_mass_bias_correction():
    # Build Sr signals at exactly natural ratios but with a small mass bias.
    # 87/86 ≈ 0.7099 (modern seawater); 86/88 ≈ 0.1194 (Steiger & Jäger).
    target_87_86 = 0.7099
    sr88 = 10000.0
    sr86 = sr88 * 0.1194
    sr87 = sr86 * target_87_86
    signals = {"Sr86": sr86, "Sr87": sr87, "Sr88": sr88}
    run = generate_laicpms_run(signals, blank_cps={k: 0.5 for k in signals},
                               sample_duration_s=15, blank_duration_s=5,
                               noise_factor=0.005, seed=3)
    ratio = laicpms.sr_ratio(run, mass_bias_correct=True)
    assert abs(ratio - target_87_86) / target_87_86 < 0.005


def test_sr_ratio_returns_nan_when_channels_missing():
    run = generate_laicpms_run({"Pb208": 1000.0}, seed=4)
    ratio = laicpms.sr_ratio(run)
    assert math.isnan(ratio)
