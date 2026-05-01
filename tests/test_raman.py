import numpy as np

from checkmsg import raman
from checkmsg.refdata import rruff
from checkmsg.synthetic import PeakSpec, generate


def _diamond(seed=0):
    return generate([PeakSpec(1332, 1.0, 1.5, 0.5)],
                    np.linspace(100, 1700, 1601), "raman", "cm-1", noise=0.005, seed=seed)


def _moissanite(seed=0):
    return generate([PeakSpec(767, 0.9, 2.0, 1.0), PeakSpec(789, 1.0, 2.0, 1.0)],
                    np.linspace(100, 1700, 1601), "raman", "cm-1", noise=0.005, seed=seed)


def _amorphous_glass(seed=0):
    return generate([PeakSpec(450, 1.0, 70.0, 30.0)],
                    np.linspace(100, 1700, 1601), "raman", "cm-1", noise=0.012, seed=seed)


def _install_refs():
    axis = np.linspace(100, 1700, 1601)
    for name, peaks in {
        "diamond": [PeakSpec(1332, 1.0, 1.5, 0.4)],
        "moissanite": [PeakSpec(767, 0.9, 2.0, 1.0), PeakSpec(789, 1.0, 2.0, 1.0)],
        "corundum": [PeakSpec(417, 1.0, 2.5, 0.8), PeakSpec(645, 0.5, 2.5, 0.8)],
        "ruby": [PeakSpec(417, 1.0, 2.5, 0.8), PeakSpec(645, 0.5, 2.5, 0.8)],
        "beryl": [PeakSpec(685, 1.0, 2.5, 0.7), PeakSpec(1067, 0.35, 2.5, 0.7)],
    }.items():
        rruff.install_synthetic_fallback(name, generate(peaks, axis, "raman", "cm-1", noise=0.001, seed=0))


def test_analyze_diamond_top_match():
    _install_refs()
    res = raman.analyze(_diamond(seed=11),
                        candidates=["diamond", "moissanite", "corundum", "ruby", "beryl"])
    assert res.best.mineral == "diamond"
    assert res.best.cosine > 0.85


def test_analyze_moissanite_top_match():
    _install_refs()
    res = raman.analyze(_moissanite(seed=12),
                        candidates=["diamond", "moissanite", "corundum", "ruby", "beryl"])
    assert res.best.mineral == "moissanite"


def test_is_amorphous_true_for_broad_envelope():
    assert raman.is_amorphous(_amorphous_glass(seed=0)) is True


def test_is_amorphous_false_for_diamond():
    assert raman.is_amorphous(_diamond(seed=0)) is False


def test_analyze_rejects_non_raman():
    import pytest

    from checkmsg.spectrum import Spectrum
    with pytest.raises(ValueError, match="expected raman"):
        raman.preprocess_raman(Spectrum(np.linspace(0, 1, 10), np.zeros(10), "xrf", "keV"))
