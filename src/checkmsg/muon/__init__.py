"""Experimental muon imaging mode (transmission + scattering + muonic X-ray).

Public API:

    from checkmsg.muon import (
        MuonSource,
        Material, VoxelGrid,
        MuonTrack, simulate_transmission, simulate_scattering, trace_one,
        reconstruct_fbp, reconstruct_sart, reconstruct_art, reconstruct_mlem,
        MuonImage,
        muonic_kalpha_keV, simulate_muonic_xray,
        analyze,
    )

The package is documented in `docs/techniques.md` (Muon imaging section) and
worked through in `examples/20_muon_tomography.py`.
"""

from __future__ import annotations

from checkmsg.muon.analyze import analyze
from checkmsg.muon.forward import (
    MuonTrack,
    ScatteringSinogram,
    TransmissionSinogram,
    simulate_scattering,
    simulate_transmission,
    trace_one,
)
from checkmsg.muon.image import MuonImage
from checkmsg.muon.muonic_xray import muonic_kalpha_keV, simulate_muonic_xray
from checkmsg.muon.reconstruct import (
    reconstruct_art,
    reconstruct_fbp,
    reconstruct_mlem,
    reconstruct_sart,
)
from checkmsg.muon.source import MuonSource
from checkmsg.muon.voxel import Material, VoxelGrid

__all__ = [
    "MuonSource",
    "Material",
    "VoxelGrid",
    "MuonTrack",
    "TransmissionSinogram",
    "ScatteringSinogram",
    "simulate_transmission",
    "simulate_scattering",
    "trace_one",
    "reconstruct_fbp",
    "reconstruct_art",
    "reconstruct_sart",
    "reconstruct_mlem",
    "MuonImage",
    "muonic_kalpha_keV",
    "simulate_muonic_xray",
    "analyze",
]
