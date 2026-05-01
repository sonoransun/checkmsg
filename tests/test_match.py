import numpy as np

from checkmsg.match import cosine, peak_list_match, rank
from checkmsg.peaks import Peak
from checkmsg.spectrum import Spectrum
from checkmsg.synthetic import PeakSpec, generate


def _spectrum(positions, seed=0):
    return generate(
        [PeakSpec(p, 1.0, 1.5, 0.5) for p in positions],
        np.linspace(100, 1700, 1601), "raman", "cm-1", noise=0.001, seed=seed,
    )


def test_cosine_self_is_one():
    s = _spectrum([700, 1100])
    assert cosine(s, s) > 0.999


def test_cosine_unrelated_is_low():
    a = _spectrum([700], seed=0)
    b = _spectrum([1300], seed=1)
    assert cosine(a, b) < 0.4


def test_cosine_technique_mismatch_raises():
    a = Spectrum(np.linspace(0, 1, 10), np.zeros(10), "raman", "cm-1")
    b = Spectrum(np.linspace(0, 1, 10), np.zeros(10), "xrf", "keV")
    import pytest
    with pytest.raises(ValueError, match="technique mismatch"):
        cosine(a, b)


def test_peak_list_match_perfect_overlap():
    pa = [Peak(700, 1.0, 5.0, 50.0), Peak(1100, 0.5, 5.0, 30.0)]
    pb = [Peak(701, 1.0, 5.0, 50.0), Peak(1099, 0.5, 5.0, 30.0)]
    res = peak_list_match(pa, pb, tolerance=5.0)
    assert res["score"] == 1.0
    assert len(res["matches"]) == 2


def test_peak_list_match_no_overlap():
    pa = [Peak(700, 1.0, 5.0, 50.0)]
    pb = [Peak(1300, 1.0, 5.0, 50.0)]
    res = peak_list_match(pa, pb, tolerance=5.0)
    assert res["score"] == 0.0


def test_rank_returns_top_candidate():
    target = _spectrum([800, 1200])
    candidates = {"correct": _spectrum([800, 1200], seed=2), "wrong": _spectrum([300, 600], seed=3)}
    out = rank(target, candidates, top=2)
    assert out[0].name == "correct"
