from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import spsolve


def als_baseline(y: np.ndarray, lam: float = 1e5, p: float = 0.001, n_iter: int = 10) -> np.ndarray:
    """Asymmetric least squares baseline (Eilers & Boelens 2005)."""
    y = np.asarray(y, dtype=float)
    L = y.size
    D = diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(L - 2, L), dtype=float, format="csc")
    DtD = lam * (D.T @ D)
    w = np.ones(L)
    z = y
    for _ in range(n_iter):
        W = diags(w, 0, shape=(L, L), format="csc")
        z = spsolve(csc_matrix(W + DtD), w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def snip_baseline(y: np.ndarray, iterations: int = 40) -> np.ndarray:
    """Statistics-sensitive Non-linear Iterative Peak (SNIP) clipping baseline.

    Standard for X-ray spectra: monotonic clipping that follows broad continua without removing peaks.
    """
    y = np.asarray(y, dtype=float)
    # LLS transform for stability
    v = np.log(np.log(np.sqrt(np.abs(y) + 1.0) + 1.0) + 1.0)
    L = v.size
    for w in range(1, iterations + 1):
        v_new = v.copy()
        # vectorized: for indices where i-w and i+w are valid
        if L > 2 * w:
            i = np.arange(w, L - w)
            avg = 0.5 * (v[i - w] + v[i + w])
            v_new[i] = np.minimum(v[i], avg)
        v = v_new
    # invert LLS
    z = (np.exp(np.exp(v) - 1.0) - 1.0) ** 2 - 1.0
    return z


def savgol(y: np.ndarray, window: int = 11, order: int = 3) -> np.ndarray:
    if window % 2 == 0:
        window += 1
    window = max(window, order + 2 + ((order + 1) % 2))
    return savgol_filter(y, window, order)


def min_max_normalize(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    lo, hi = float(np.min(y)), float(np.max(y))
    if hi <= lo:
        return np.zeros_like(y)
    return (y - lo) / (hi - lo)


def area_normalize(y: np.ndarray, x: np.ndarray | None = None) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    a = np.trapezoid(np.abs(y), x) if x is not None else np.sum(np.abs(y))
    return y / a if a else y
