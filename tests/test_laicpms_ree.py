import math

from checkmsg import laicpms
from checkmsg.refdata.icpms_data import CHONDRITE_REE_PPM, REE_ELEMENTS


def _conc(el: str, ppm: float) -> laicpms.Concentration:
    return laicpms.Concentration(element=el, ppm=ppm)


def test_ree_pattern_round_trip_against_chondrite():
    # Sample with concentrations exactly equal to chondrite => normalised pattern of 1.0.
    sample = {el: _conc(el, ppm) for el, ppm in CHONDRITE_REE_PPM.items()}
    pattern = laicpms.ree_pattern(sample)
    for el in REE_ELEMENTS:
        assert abs(pattern[el] - 1.0) < 1e-9, f"{el} normalisation failed"


def test_ree_pattern_skips_missing_elements():
    sample = {"La": _conc("La", 0.5), "Yb": _conc("Yb", 1.0)}
    pattern = laicpms.ree_pattern(sample)
    assert set(pattern) == {"La", "Yb"}


def test_eu_anomaly_detects_negative():
    # Igneous zircon: strong negative Eu anomaly. Set Eu deliberately low vs Sm/Gd.
    sample = {
        "Sm": _conc("Sm", 1.0 * CHONDRITE_REE_PPM["Sm"]),
        "Eu": _conc("Eu", 0.3 * CHONDRITE_REE_PPM["Eu"]),
        "Gd": _conc("Gd", 1.0 * CHONDRITE_REE_PPM["Gd"]),
    }
    pattern = laicpms.ree_pattern(sample)
    assert laicpms.eu_anomaly(pattern) < 0.5


def test_eu_anomaly_returns_unity_for_no_anomaly():
    sample = {el: _conc(el, ppm) for el, ppm in CHONDRITE_REE_PPM.items()}
    pattern = laicpms.ree_pattern(sample)
    assert abs(laicpms.eu_anomaly(pattern) - 1.0) < 1e-9


def test_ce_anomaly_detects_positive():
    # Hydrothermal zircon shows positive Ce anomaly (Ce4+ enrichment).
    sample = {
        "La": _conc("La", 1.0 * CHONDRITE_REE_PPM["La"]),
        "Ce": _conc("Ce", 5.0 * CHONDRITE_REE_PPM["Ce"]),
        "Pr": _conc("Pr", 1.0 * CHONDRITE_REE_PPM["Pr"]),
    }
    pattern = laicpms.ree_pattern(sample)
    assert laicpms.ce_anomaly(pattern) > 3.0


def test_anomalies_return_nan_for_missing_data():
    pattern = {"La": 1.0}  # missing Ce, Pr, Sm, Eu, Gd
    assert math.isnan(laicpms.eu_anomaly(pattern))
    assert math.isnan(laicpms.ce_anomaly(pattern))
