"""RRUFF fetcher tests — network is mocked; CHECKMSG_OFFLINE behavior is exercised."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from checkmsg.refdata import rruff
from checkmsg.refdata.cache import cache_root

SAMPLE_RRUFF_TXT = """\
##NAME=Diamond
##RRUFFID=R050204
##LASER_WAVELENGTH=532
##METADATA=mocked
100.0,0.05
500.0,0.10
1330.0,0.20
1332.0,1.00
1334.0,0.30
1700.0,0.05
"""


def test_fetch_writes_to_cache(monkeypatch):
    monkeypatch.delenv("CHECKMSG_OFFLINE", raising=False)
    mock_resp = MagicMock(text=SAMPLE_RRUFF_TXT, status_code=200)
    mock_resp.raise_for_status = MagicMock()
    monkeypatch.setattr(rruff.requests, "get", MagicMock(return_value=mock_resp))
    s = rruff.fetch("diamond")
    assert s.technique == "raman"
    assert s.units == "cm-1"
    assert (cache_root() / "rruff" / "diamond.txt").exists()


def test_fetch_returns_cached_on_second_call(monkeypatch):
    monkeypatch.delenv("CHECKMSG_OFFLINE", raising=False)
    mock_get = MagicMock(return_value=MagicMock(text=SAMPLE_RRUFF_TXT, status_code=200, raise_for_status=MagicMock()))
    monkeypatch.setattr(rruff.requests, "get", mock_get)
    rruff.fetch("diamond")
    rruff.fetch("diamond")
    assert mock_get.call_count == 1


def test_unknown_mineral_raises():
    with pytest.raises(KeyError):
        rruff.fetch("unobtainium")


def test_offline_mode_with_no_cache_fails(monkeypatch):
    monkeypatch.setenv("CHECKMSG_OFFLINE", "1")
    with pytest.raises(RuntimeError, match="cache miss"):
        rruff.fetch("diamond")


def test_install_synthetic_fallback_round_trip():
    import numpy as np

    from checkmsg.spectrum import Spectrum
    spec = Spectrum(np.array([100, 500, 1332.0, 1500.0]),
                    np.array([0.0, 0.1, 1.0, 0.2]), "raman", "cm-1")
    rruff.install_synthetic_fallback("diamond", spec)
    loaded = rruff.load_cached("diamond")
    assert loaded.technique == "raman"
    assert abs(loaded.intensity.max() - 1.0) < 1e-6
