import pytest

from checkmsg import laicpms
from checkmsg.refdata.icpms_data import NIST_SRM_612
from checkmsg.synthetic import generate_laicpms_run

SENS = 1000.0  # cps per ppm — uniform across channels for these tests


def _signals(true_ppm):
    """Build per-isotope cps from element concentrations using uniform sensitivity."""
    out = {}
    for el, ppm in true_ppm.items():
        for key, rec in laicpms.ISOTOPES.items():
            if rec.element == el:
                out[key] = ppm * rec.natural_abundance * SENS
    return out


def _calibration_run(seed=0):
    sample = _signals(NIST_SRM_612)
    return generate_laicpms_run(sample, sample_duration_s=15, blank_duration_s=5,
                                noise_factor=0.02, seed=seed)


def test_quant_recovers_known_ppm_within_5_percent():
    truth = {"Mn": 50.0, "Sr": 200.0, "Pb": 25.0}
    sample = generate_laicpms_run(_signals(truth), sample_duration_s=15, blank_duration_s=5,
                                  noise_factor=0.02, seed=1)
    cal = _calibration_run(seed=2)
    out = laicpms.quantify(sample, cal, internal_standard_element=None,
                           internal_standard_ppm=None, reference="NIST612")
    for el, expected in truth.items():
        c = out[el]
        rel_err = abs(c.ppm - expected) / expected
        assert rel_err < 0.05, f"{el}: got {c.ppm:.2f}, expected {expected}"


def test_internal_standard_correction_recovers_within_2_percent():
    # Sample matrix has 800,000 ppm Ca (calcite); calibration is silicate (NIST 612).
    # With IS=Ca, Longerich correction scales sample concentrations to match the IS.
    truth = {"Ca": 400000.0, "Mn": 50.0, "Sr": 200.0}
    sample = generate_laicpms_run(_signals(truth), sample_duration_s=15, blank_duration_s=5,
                                  noise_factor=0.02, seed=3)
    cal = _calibration_run(seed=4)
    out = laicpms.quantify(sample, cal, internal_standard_element="Ca",
                           internal_standard_ppm=400000.0, reference="NIST612")
    assert abs(out["Mn"].ppm - 50.0) / 50.0 < 0.02
    assert abs(out["Sr"].ppm - 200.0) / 200.0 < 0.02


def test_unknown_reference_raises():
    sample = generate_laicpms_run({"Pb208": 1000.0}, seed=5)
    with pytest.raises(KeyError):
        laicpms.quantify(sample, sample, reference="bogus")


def test_concentration_includes_detection_limit():
    sample = generate_laicpms_run(_signals({"Pb": 25.0}), seed=6)
    cal = _calibration_run(seed=7)
    out = laicpms.quantify(sample, cal)
    assert out["Pb"].detection_limit_ppm >= 0.0
