"""Muonic-atom K_α spectroscopy.

When a negative muon is captured by an atom it cascades into the K shell and
emits a series of muonic X-rays. The K_α line (2p → 1s) carries the highest
intensity and has Z²-scaling, making it diagnostic of the absorbing element.

`simulate_muonic_xray` takes a list of stopping voxels (returned by
`forward.simulate_*`) and a `VoxelGrid`, accumulates K_α counts per element, and
returns a `Spectrum` (technique='muon-xray', units='keV') that can be fed to
existing peak-detection code.
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from checkmsg.muon.voxel import VoxelGrid
from checkmsg.refdata.muon_data import MUONIC_KALPHA_keV, muonic_kalpha_keV
from checkmsg.spectrum import Spectrum
from checkmsg.synthetic import voigt_pseudo


def _element_from_material(material) -> str | None:
    """Best-guess element symbol for a `Material`. Compounds return None."""
    elemental = {
        "iron": "Fe", "copper": "Cu", "silver": "Ag", "platinum": "Pt",
        "gold": "Au", "lead": "Pb", "uranium": "U", "diamond": "C",
    }
    return elemental.get(material.name.split()[0].lower())


def _composition_elements(material) -> list[str]:
    """Decompose a material into its elemental constituents (rough)."""
    table = {
        "corundum": ["Al", "O"],
        "beryl": ["Be", "Al", "Si", "O"],
        "quartz": ["Si", "O"],
        "calcite": ["Ca", "C", "O"],
        "olivine (forsterite)": ["Mg", "Si", "O"],
        "iron-nickel (taenite)": ["Fe", "Ni"],
        "polymer (CH2)": ["C", "H"],
        "water": ["H", "O"],
        "air": ["N", "O"],
    }
    if material.name in table:
        return table[material.name]
    elem = _element_from_material(material)
    return [elem] if elem else []


def simulate_muonic_xray(
    stopping_voxels: list[tuple[int, int, int]],
    grid: VoxelGrid,
    *,
    energies_keV: np.ndarray | None = None,
    fwhm_keV: float = 5.0,
    noise: float = 0.02,
    seed: int = 0,
) -> Spectrum:
    """Build a muonic-X-ray spectrum from a list of stopping voxels.

    Each stopped muon emits one K_α photon; element identity comes from the
    `VoxelGrid`'s material at that voxel. Compounds contribute Z-weighted counts
    proportionally across their elements.

    Returns a `Spectrum` with technique='muon-xray', units='keV'.
    """
    rng = np.random.default_rng(seed)
    if energies_keV is None:
        energies_keV = np.linspace(50.0, 7000.0, 2048)

    # Tally one count per stopping muon, distributed across the voxel material's
    # elements weighted by their stoichiometric Z fraction (rough proxy for muon
    # capture probability — Fermi-Teller Z law).
    counts = Counter[str]()
    for ix, iy, iz in stopping_voxels:
        if not (0 <= ix < grid.shape[0] and 0 <= iy < grid.shape[1] and 0 <= iz < grid.shape[2]):
            continue
        material = grid.materials[ix, iy, iz]
        elements = _composition_elements(material)
        # Weight by Z (Fermi-Teller capture probability scales ~Z; here we treat
        # equal-weight as a coarse proxy unless the material is monoelemental).
        weights = []
        for el in elements:
            try:
                weights.append(muonic_kalpha_keV(el))  # any non-zero weight; not the energy
            except KeyError:
                weights.append(0.0)
        total = sum(weights) or 1.0
        for el, w in zip(elements, weights, strict=False):
            if w > 0:
                counts[el] += w / total

    # Build the spectrum: one Voigt per element at its K_α energy.
    intensity = np.zeros_like(energies_keV)
    sigma = fwhm_keV / 2.355
    for element, n in counts.items():
        if n <= 0:
            continue
        energy_keV = MUONIC_KALPHA_keV[element]
        if energy_keV < energies_keV.min() or energy_keV > energies_keV.max():
            continue
        intensity = intensity + voigt_pseudo(
            energies_keV, center=energy_keV, sigma=sigma, gamma=0.5 * sigma,
            amplitude=float(n),
        )

    if noise > 0:
        peak = float(np.max(np.abs(intensity)) or 1.0)
        intensity = intensity + rng.normal(0.0, noise * peak, size=intensity.size)

    return Spectrum(
        energies_keV, intensity, "muon-xray", "keV",
        metadata={
            "fwhm_keV": fwhm_keV,
            "n_stopped_muons": len(stopping_voxels),
            "elements_detected": sorted(counts.keys()),
        },
    )
