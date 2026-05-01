import numpy as np
import pytest

from checkmsg.epr import infer_g_factors, simulate_field_sweep
from checkmsg.refdata.epr_centers import CENTERS, by_host, diagnostic_for, quant_standards


@pytest.mark.parametrize("name", list(CENTERS))
def test_every_center_simulates(name):
    """All bundled centers must produce a non-empty spectrum at X-band without raising."""
    sys_ = CENTERS[name]
    # Sweep wide enough to catch resonances spread by D / hyperfine.
    fields = np.linspace(50, 700, 651)
    spec = simulate_field_sweep(sys_, frequency_GHz=9.5, fields_mT=fields, orientations=(7, 1))
    assert spec.technique == "epr"
    # At least one detectable feature in the swept range.
    assert float(np.max(np.abs(spec.intensity))) > 0


def test_quant_standards_only_returns_standards():
    standards = quant_standards()
    assert "DPPH" in standards
    assert "free_electron" in standards
    assert "diamond_P1" not in standards


def test_by_host_filters_correctly():
    diamonds = by_host("diamond")
    assert "diamond_P1" in diamonds
    assert "diamond_Ni_HPHT" in diamonds
    assert "DPPH" not in diamonds


def test_diagnostic_for_known_scenarios():
    assert "diamond_P1" in diagnostic_for("diamond-natural-vs-synthetic")
    assert "calcite_Mn2plus" in diagnostic_for("pearl-fingerprint")
    assert diagnostic_for("nonexistent-scenario") == ()


def test_p1_produces_three_lines_in_powder():
    """P1 anisotropic 14N hyperfine produces a 3-line pattern."""
    sys_ = CENTERS["diamond_P1"]
    fields = np.linspace(330, 350, 1001)
    spec = simulate_field_sweep(sys_, frequency_GHz=9.5, fields_mT=fields, orientations=(11, 1))
    gs = infer_g_factors(spec, 9.5)
    assert 3 <= len(gs) <= 5  # central + two satellites; some powder shoulders OK


def test_e1prime_produces_single_sharp_line():
    sys_ = CENTERS["quartz_E1prime"]
    fields = np.linspace(336, 342, 1001)
    spec = simulate_field_sweep(sys_, frequency_GHz=9.5, fields_mT=fields)
    gs = infer_g_factors(spec, 9.5)
    assert len(gs) == 1
    assert abs(gs[0] - 2.0006) < 0.001
