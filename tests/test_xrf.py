import numpy as np

from checkmsg import xrf
from checkmsg.refdata.nist_xray import candidates_at, lines_for
from checkmsg.synthetic import PeakSpec, generate


def test_nist_table_loads_and_has_iron():
    lines = lines_for("Fe")
    assert any(line.line == "Ka" and abs(line.energy_keV - 6.404) < 0.01 for line in lines)


def test_candidates_within_tolerance():
    cands = candidates_at(6.40, tolerance_keV=0.05)
    assert any(c.element == "Fe" for c in cands)


def test_identify_elements_picks_iron():
    axis = np.linspace(0.5, 12.0, 4096)
    s = generate([
        PeakSpec(6.404, 1.0, 0.045, 0.025),
        PeakSpec(7.058, 0.13, 0.045, 0.025),
    ], axis, "xrf", "keV", noise=0.001, seed=5)
    res = xrf.identify_elements(s, tolerance_keV=0.05, min_snr=8.0)
    assert any(e.element == "Fe" for e in res.elements)


def test_relative_quant_sums_to_one():
    axis = np.linspace(0.5, 12.0, 4096)
    s = generate([
        PeakSpec(1.740, 1.0, 0.045, 0.025),  # Si
        PeakSpec(8.048, 0.4, 0.045, 0.025),  # Cu
    ], axis, "xrf", "keV", noise=0.001, seed=6)
    res = xrf.identify_elements(s, tolerance_keV=0.05, min_snr=8.0)
    rq = xrf.relative_quant(res)
    assert abs(sum(rq.values()) - 1.0) < 1e-6


def test_xrf_rejects_wrong_technique():
    import pytest

    from checkmsg.spectrum import Spectrum
    with pytest.raises(ValueError, match="expected xrf"):
        xrf.preprocess_xrf(Spectrum(np.linspace(0, 1, 10), np.zeros(10), "raman", "cm-1"))
