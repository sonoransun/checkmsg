"""Reconstruction algorithm sanity tests: forward-then-back recovers a known phantom."""

import math

import numpy as np

from checkmsg.muon.reconstruct import (
    _forward_project,
    reconstruct_art,
    reconstruct_fbp,
    reconstruct_mlem,
    reconstruct_sart,
)


def _disk_phantom(size: int = 32) -> np.ndarray:
    """Simple disk phantom: 1.0 inside a circle of radius size/4, 0.0 outside."""
    coords = np.arange(size) - (size - 1) / 2.0
    xs, ys = np.meshgrid(coords, coords, indexing="ij")
    r2 = xs ** 2 + ys ** 2
    return (r2 < (size / 4) ** 2).astype(float)


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = a.flatten() - a.mean()
    b = b.flatten() - b.mean()
    norm = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(a @ b / norm) if norm > 0 else 0.0


def test_forward_then_fbp_recovers_phantom():
    phantom = _disk_phantom(32)
    angles = np.linspace(0, math.pi, 60, endpoint=False)
    sino = _forward_project(phantom, angles, 32)
    rec = reconstruct_fbp(sino, angles, 32)
    assert _pearson(rec, phantom) > 0.7


def test_sart_converges_below_target_rmse():
    phantom = _disk_phantom(24)
    angles = np.linspace(0, math.pi, 40, endpoint=False)
    sino = _forward_project(phantom, angles, 24)
    rec = reconstruct_sart(sino, angles, 24, iterations=8, relaxation=0.15)
    assert _pearson(rec, phantom) > 0.8


def test_mlem_produces_nonnegative_image():
    phantom = _disk_phantom(20)
    angles = np.linspace(0, math.pi, 30, endpoint=False)
    sino = _forward_project(phantom, angles, 20)
    rec = reconstruct_mlem(sino, angles, 20, iterations=4)
    assert (rec >= 0).all()
    assert _pearson(rec, phantom) > 0.5


def test_art_recovers_low_resolution_phantom():
    phantom = _disk_phantom(16)
    angles = np.linspace(0, math.pi, 24, endpoint=False)
    sino = _forward_project(phantom, angles, 16)
    rec = reconstruct_art(sino, angles, 16, iterations=4, relaxation=0.15)
    assert _pearson(rec, phantom) > 0.5
