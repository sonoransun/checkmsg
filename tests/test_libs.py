import numpy as np

from checkmsg import libs
from checkmsg.refdata.nist_asd import all_elements, lines_for
from checkmsg.synthetic import PeakSpec, generate


def test_table_has_iron_lines():
    fe = lines_for("Fe")
    assert fe and any(abs(line.wavelength_nm - 371.994) < 0.01 for line in fe)


def test_all_elements_includes_chromophores():
    elements = all_elements()
    assert {"Cr", "Fe", "Ti"}.issubset(elements)


def test_identify_picks_iron_and_titanium():
    axis = np.linspace(250.0, 800.0, 8000)
    s = generate([
        PeakSpec(371.994, 1.0, 0.10, 0.05),
        PeakSpec(334.941, 0.7, 0.10, 0.05),
    ], axis, "libs", "nm", noise=0.001, seed=2)
    res = libs.identify(s, tolerance_nm=0.4, min_snr=6.0)
    assert res.has("Fe")
    assert res.has("Ti")


def test_trace_ratios_basic():
    axis = np.linspace(250.0, 800.0, 8000)
    s = generate([
        PeakSpec(371.994, 2.0, 0.10, 0.05),  # Fe
        PeakSpec(334.941, 1.0, 0.10, 0.05),  # Ti
    ], axis, "libs", "nm", noise=0.001, seed=3)
    res = libs.identify(s, tolerance_nm=0.4)
    ratios = libs.trace_ratios(res, [("Fe", "Ti")])
    assert ratios["Fe/Ti"] > 1.0
