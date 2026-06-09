"""Tests for the spectral-embedding similarity search."""

from __future__ import annotations

import numpy as np
import pytest

from checkmsg import minerals, similarity


def test_embed_is_unit_vector():
    spec = minerals.synthesize_raman(minerals.get("diamond"), seed=0)
    v = similarity.embed(spec)
    assert abs(float(np.linalg.norm(v)) - 1.0) < 1e-9


@pytest.mark.parametrize("name", ["diamond", "magnetite", "fluorite", "ruby"])
def test_self_mineral_is_top_match(name):
    spec = minerals.synthesize_raman(minerals.get(name), seed=0)
    hits = similarity.nearest(spec, k=3)
    assert hits[0].name == name
    assert 0.0 <= hits[0].score <= 1.0 + 1e-9


def test_unsupported_technique_returns_empty():
    # muon-xray has no canonical grid / catalog references.
    from checkmsg.spectrum import Spectrum
    spec = Spectrum(np.linspace(1, 100, 100), np.ones(100), "laicpms", "m/z")
    assert similarity.nearest(spec) == []


def test_index_is_cached():
    similarity._index.cache_clear()
    similarity.nearest(minerals.synthesize_raman(minerals.get("ruby"), seed=0))
    info = similarity._index.cache_info()
    similarity.nearest(minerals.synthesize_raman(minerals.get("emerald"), seed=0))
    assert similarity._index.cache_info().hits > info.hits
