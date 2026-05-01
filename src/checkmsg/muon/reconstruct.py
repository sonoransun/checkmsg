"""Reconstruction algorithms for muon tomography sinograms.

  - `reconstruct_fbp`   — analytic filtered back projection (ramp filter via FFT).
  - `reconstruct_art`   — Algebraic Reconstruction Technique (Kaczmarz iteration).
  - `reconstruct_sart`  — Simultaneous ART (SART), the production-grade default.
  - `reconstruct_mlem`  — Maximum-Likelihood Expectation-Maximisation (Poisson model).

All four operate on **2-D** central slices: the input is a sinogram (n_proj × n_pix)
and the output is an n_pix × n_pix slice. For full 3-D reconstruction we slice the
volume into z-layers and reconstruct each layer independently — adequate for the
parallel-beam geometry produced by `simulate_transmission` / `simulate_scattering`.
"""

from __future__ import annotations

import math

import numpy as np


def _ramp_filter(n: int) -> np.ndarray:
    """Cropped ramp filter in the frequency domain (length n, even)."""
    freqs = np.fft.fftfreq(n)
    # Use Ram-Lak ramp |f| with a Hann taper for noise control.
    ramp = 2 * np.abs(freqs)
    hann = 0.5 * (1 + np.cos(2 * math.pi * freqs))
    return ramp * hann


def _backproject(sinogram: np.ndarray, angles_rad: np.ndarray,
                 image_size: int) -> np.ndarray:
    """Plain back projection of a (n_proj, n_pix) sinogram into a square image."""
    n_proj, n_pix = sinogram.shape
    image = np.zeros((image_size, image_size), dtype=float)
    centre = (image_size - 1) / 2.0
    pix_centre = (n_pix - 1) / 2.0
    coords = np.arange(image_size) - centre
    xs, ys = np.meshgrid(coords, coords, indexing="ij")
    for i, angle in enumerate(angles_rad):
        # For each pixel (x, y), the sinogram bin is x cos θ + y sin θ.
        proj = xs * math.cos(angle) + ys * math.sin(angle)
        bins = proj + pix_centre
        # Linear interpolation between adjacent sinogram bins.
        lo = np.clip(np.floor(bins).astype(int), 0, n_pix - 1)
        hi = np.clip(lo + 1, 0, n_pix - 1)
        frac = np.clip(bins - lo, 0.0, 1.0)
        contrib = (1 - frac) * sinogram[i, lo] + frac * sinogram[i, hi]
        image += contrib
    return image / n_proj


def _forward_project(image: np.ndarray, angles_rad: np.ndarray,
                     n_pixels: int) -> np.ndarray:
    """Plain Radon forward projection: image (Ny, Nx) → sinogram (n_proj, n_pixels)."""
    image_size = image.shape[0]
    centre = (image_size - 1) / 2.0
    pix_centre = (n_pixels - 1) / 2.0
    coords = np.arange(image_size) - centre
    xs, ys = np.meshgrid(coords, coords, indexing="ij")
    sino = np.zeros((len(angles_rad), n_pixels), dtype=float)
    for i, angle in enumerate(angles_rad):
        # For each image pixel, accumulate into the projection bin.
        bins = xs * math.cos(angle) + ys * math.sin(angle) + pix_centre
        lo = np.clip(np.floor(bins).astype(int), 0, n_pixels - 1)
        hi = np.clip(lo + 1, 0, n_pixels - 1)
        frac = np.clip(bins - lo, 0.0, 1.0)
        np.add.at(sino[i], lo, image * (1 - frac))
        np.add.at(sino[i], hi, image * frac)
    return sino


# --- Public reconstructors ---------------------------------------------------------


def reconstruct_fbp(sinogram: np.ndarray, angles_rad: np.ndarray,
                    image_size: int) -> np.ndarray:
    """Filtered back projection (ramp + Hann)."""
    filt = _ramp_filter(sinogram.shape[1])
    fft_sino = np.fft.fft(sinogram, axis=1)
    filtered = np.real(np.fft.ifft(fft_sino * filt[None, :], axis=1))
    return _backproject(filtered, angles_rad, image_size)


def _sart_update(residual_sino: np.ndarray, angles_rad: np.ndarray,
                  image_size: int, ones_image: np.ndarray) -> np.ndarray:
    """Properly normalised SART update: D^{-1} A^T M^{-1} residual.

    M^{-1} normalises each ray by its forward-projected length (so a ray that
    crosses N pixels contributes residual/N back to each touched pixel);
    D^{-1} normalises each pixel by the back-projection of ones (column sum).
    Without these the update overshoots and the clip-to-zero traps the iteration.
    """
    n_pixels = residual_sino.shape[1]
    row_weights = _forward_project(ones_image, angles_rad, n_pixels)
    row_weights = np.where(row_weights > 1e-9, row_weights, 1.0)
    weighted = residual_sino / row_weights
    col_weights = _backproject(np.ones_like(residual_sino), angles_rad, image_size)
    col_weights = np.where(col_weights > 1e-9, col_weights, 1.0)
    return _backproject(weighted, angles_rad, image_size) / col_weights


def reconstruct_art(sinogram: np.ndarray, angles_rad: np.ndarray,
                    image_size: int, *, iterations: int = 5,
                    relaxation: float = 0.2) -> np.ndarray:
    """Kaczmarz / ART — sequential per-ray updates."""
    image = np.zeros((image_size, image_size), dtype=float)
    ones = np.ones((image_size, image_size), dtype=float)
    for _ in range(iterations):
        for i, angle in enumerate(angles_rad):
            est = _forward_project(image, np.array([angle]), sinogram.shape[1])
            residual = sinogram[i:i + 1] - est
            update = _sart_update(residual, np.array([angle]), image_size, ones)
            image = image + relaxation * update
            image = np.clip(image, 0.0, None)
    return image


def reconstruct_sart(sinogram: np.ndarray, angles_rad: np.ndarray,
                     image_size: int, *, iterations: int = 5,
                     relaxation: float = 0.2) -> np.ndarray:
    """Simultaneous ART with proper M^{-1} A^T D^{-1} update normalisation.

    Far more stable than naive ART; converges in 5–10 iterations on parallel-beam data.
    """
    image = np.zeros((image_size, image_size), dtype=float)
    ones = np.ones((image_size, image_size), dtype=float)
    for _ in range(iterations):
        est = _forward_project(image, angles_rad, sinogram.shape[1])
        residual = sinogram - est
        update = _sart_update(residual, angles_rad, image_size, ones)
        image = image + relaxation * update
        image = np.clip(image, 0.0, None)
    return image


def reconstruct_mlem(sinogram: np.ndarray, angles_rad: np.ndarray,
                     image_size: int, *, iterations: int = 8) -> np.ndarray:
    """Maximum-Likelihood Expectation-Maximisation for Poisson-distributed counts.

    The sinogram should be non-negative (count-like). The reconstruction is
    multiplicative: image_{k+1} = image_k × (BP[sino / FP[image_k]] / sensitivity).
    """
    sino = np.maximum(sinogram, 0.0)
    image = np.ones((image_size, image_size), dtype=float)
    sensitivity = _backproject(np.ones_like(sino), angles_rad, image_size)
    sensitivity = np.where(sensitivity > 0, sensitivity, 1.0)
    for _ in range(iterations):
        est = _forward_project(image, angles_rad, sino.shape[1])
        ratio = sino / np.where(est > 1e-9, est, 1.0)
        correction = _backproject(ratio, angles_rad, image_size) / sensitivity
        image = image * correction
        image = np.clip(image, 0.0, None)
    return image
