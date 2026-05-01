import numpy as np

from checkmsg.identify import combined_report
from checkmsg.refdata import rruff
from checkmsg.synthetic import PeakSpec, generate


def _install_diamond_ref():
    axis = np.linspace(100, 1700, 1601)
    spec = generate([PeakSpec(1332, 1.0, 1.5, 0.4)], axis, "raman", "cm-1", noise=0.001, seed=0)
    rruff.install_synthetic_fallback("diamond", spec)


def test_combined_report_with_one_technique():
    _install_diamond_ref()
    raman_spec = generate([PeakSpec(1332, 1.0, 1.5, 0.5)],
                          np.linspace(100, 1700, 1601), "raman", "cm-1", noise=0.005, seed=0)
    result = combined_report([raman_spec])
    assert result.raman is not None
    assert result.raman.best.mineral == "diamond"
    assert result.xrf is None
    assert "Raman" in result.headline()


def test_combined_report_with_amorphous_note():
    glass = generate([PeakSpec(450, 1.0, 70.0, 30.0)],
                     np.linspace(100, 1700, 1601), "raman", "cm-1", noise=0.012, seed=0)
    result = combined_report([glass])
    assert any("amorphous" in n for n in result.notes)
    assert result.report()  # smoke


def test_combined_report_xrf_lists_iron():
    s = generate([PeakSpec(6.404, 1.0, 0.045, 0.025)],
                 np.linspace(0.5, 12, 4096), "xrf", "keV", noise=0.001, seed=0)
    result = combined_report([s])
    assert any(e.element == "Fe" for e in result.xrf.elements)
