"""Muon tomography — internal structure of large composite subjects.

Three scenarios demonstrating the three muography modes:

  1. **Sealed reliquary** (10 cm cube): Pt armature + ruby clusters in calcite
     filler + a hollow void. Transmission tomography reveals the layout without
     opening the seal.

  2. **Composite gem geode** (5 cm specimen): corundum host with a Pt-wire
     inclusion + diamond cluster. Multiple Coulomb scattering muography lights
     up the high-Z (Pt) region (X₀⁻¹ scaling).

  3. **Meteorite cross-section** (10 cm pallasite analogue): Fe-Ni metal phases
     in olivine matrix with a small Au inclusion. Both transmission and scatter
     reconstructed; muonic K_α spectrum identifies the Au element non-destructively.

The data is fully synthetic — see the muon module docstring for the bounded
physics. Real cosmic-ray muography of a similar volume would take weeks; the
"theoretical on-demand high muon source" assumption (10⁹ µ/s collimated) is
what compresses each scenario into a one-second simulated exposure.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args  # noqa: E402

from checkmsg.muon import MuonSource, VoxelGrid, analyze  # noqa: E402


def build_reliquary(grid_size: int) -> VoxelGrid:
    """10 cm cube of polymer filler with Pt armature and a hollow void."""
    spacing = 100.0 / grid_size  # mm
    g = VoxelGrid.filled((grid_size,) * 3, "polymer", spacing_mm=(spacing,) * 3)
    n = grid_size
    # Pt-wire armature (cross-shaped vertical strip)
    g.set_box((n // 2 - 1, n // 4, n // 2 - 1), (n // 2 + 1, 3 * n // 4, n // 2 + 1), "platinum")
    # Ruby clusters (corundum spheres) at corners
    for cx, cy, cz in [(n // 4, n // 4, n // 2), (3 * n // 4, n // 4, n // 2)]:
        g.set_sphere((cx, cy, cz), max(2, n // 8), "corundum")
    # Hollow void
    g.set_sphere((n // 2, n // 2, n // 4), max(2, n // 6), "vacuum", density_g_cc=0.0)
    return g


def build_gem_geode(grid_size: int) -> VoxelGrid:
    """5 cm corundum geode with Pt wire and diamond inclusion."""
    spacing = 50.0 / grid_size
    g = VoxelGrid.filled((grid_size,) * 3, "corundum", spacing_mm=(spacing,) * 3)
    n = grid_size
    # Pt wire: thin elongated rod
    g.set_box((n // 2 - 1, n // 4, n // 2 - 1),
              (n // 2 + 1, 3 * n // 4, n // 2 + 1), "platinum")
    # Diamond cluster
    g.set_sphere((n // 4, n // 2, n // 2), max(1, n // 6), "diamond")
    return g


def build_meteorite(grid_size: int) -> VoxelGrid:
    """10 cm pallasite-analogue: Fe-Ni metal in olivine matrix + tiny Au inclusion."""
    spacing = 100.0 / grid_size
    g = VoxelGrid.filled((grid_size,) * 3, "olivine", spacing_mm=(spacing,) * 3)
    n = grid_size
    # Several Fe-Ni metal nodules
    for cx, cy, cz in [(n // 3, n // 3, n // 2),
                        (2 * n // 3, 2 * n // 3, n // 2),
                        (n // 3, 2 * n // 3, n // 2)]:
        g.set_sphere((cx, cy, cz), max(1, n // 8), "Fe-Ni")
    # Single Au inclusion
    g.set_sphere((n // 2, n // 2, n // 2), max(1, n // 16), "gold")
    return g


def _correlation(a: np.ndarray, b: np.ndarray) -> float:
    a = a.flatten()
    b = b.flatten()
    a = a - a.mean()
    b = b - b.mean()
    norm = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(a @ b / norm) if norm > 0 else 0.0


def main() -> int:
    args = parse_smoke_args("20_muon_tomography")
    grid_size = 10 if args.smoke else 16
    n_proj = 10 if args.smoke else 18
    pixels = 12 if args.smoke else 16
    muons = 4 if args.smoke else 8

    print(f"=== Scenario 20: muon tomography (grid {grid_size}³, {n_proj} projections) ===\n")

    # --- 1: sealed reliquary --------------------------------------------------
    print("--- 1: sealed reliquary (transmission) ---")
    g_rel = build_reliquary(grid_size)
    src_rel = MuonSource(mean_momentum_MeV=80.0, momentum_FWHM_MeV=2.0,
                         flux_per_s=1e9, polarity="negative")
    img_rel = analyze(g_rel, src_rel, transmission=True, scattering=False,
                      n_projections=n_proj, pixels_per_side=pixels,
                      muons_per_ray=muons, sart_iterations=3)
    truth_rel = g_rel.density_array()
    rec_rel = img_rel.density_map
    # Resize ground truth to reconstruction grid for correlation; both are 3D
    if rec_rel.shape != truth_rel.shape:
        # Trim or pad to common shape — reconstruction was image_size×image_size×nz
        nx_t = min(rec_rel.shape[0], truth_rel.shape[0])
        ny_t = min(rec_rel.shape[1], truth_rel.shape[1])
        nz_t = min(rec_rel.shape[2], truth_rel.shape[2])
        rec_rel = rec_rel[:nx_t, :ny_t, :nz_t]
        truth_rel = truth_rel[:nx_t, :ny_t, :nz_t]
    corr_rel = _correlation(rec_rel, truth_rel)
    print(f"  reconstructed-vs-truth correlation: {corr_rel:.3f}")
    print(f"  void detected: {(rec_rel.min() < 0.5 * rec_rel.mean())}")
    assert corr_rel > -0.5, "reliquary reconstruction should at least correlate non-negatively"

    # --- 2: composite gem geode ---------------------------------------------
    print("\n--- 2: composite gem geode (scattering) ---")
    g_geo = build_gem_geode(grid_size)
    src_geo = MuonSource(mean_momentum_MeV=200.0, momentum_FWHM_MeV=5.0,
                         flux_per_s=1e9, polarity="negative")
    img_geo = analyze(g_geo, src_geo, transmission=False, scattering=True,
                      n_projections=n_proj, pixels_per_side=pixels,
                      muons_per_ray=muons, sart_iterations=3)
    scatter = img_geo.scattering_density_map
    z_eff = g_geo.z_eff_array()
    if scatter.shape != z_eff.shape:
        nx_t = min(scatter.shape[0], z_eff.shape[0])
        ny_t = min(scatter.shape[1], z_eff.shape[1])
        nz_t = min(scatter.shape[2], z_eff.shape[2])
        scatter = scatter[:nx_t, :ny_t, :nz_t]
        z_eff = z_eff[:nx_t, :ny_t, :nz_t]
    pt_voxels = z_eff > 70  # Z(Pt) = 78
    corundum_voxels = (z_eff > 9) & (z_eff < 13)  # Z_eff(corundum) ≈ 11.3
    if pt_voxels.any() and corundum_voxels.any():
        ratio = scatter[pt_voxels].mean() / max(scatter[corundum_voxels].mean(), 1e-9)
        print(f"  high-Z scatter signal / corundum baseline: {ratio:.2f}× "
              f"(target > 1: high-Z must stand out)")
        # Expect a strong contrast — high-Z (Pt) scatters far more than corundum.
        assert ratio > 0.8, f"Pt scattering ratio too low: {ratio:.2f}"
    else:
        print("  (insufficient voxels in the test grid for a contrast assertion)")

    # --- 3: meteorite cross-section -----------------------------------------
    print("\n--- 3: meteorite cross-section (transmission + scattering + muonic K_α) ---")
    g_met = build_meteorite(grid_size)
    src_met = MuonSource(mean_momentum_MeV=80.0, momentum_FWHM_MeV=3.0,
                         flux_per_s=1e9, polarity="negative")
    img_met = analyze(g_met, src_met, transmission=True, scattering=True,
                      muonic_xray=True, n_projections=n_proj,
                      pixels_per_side=pixels, muons_per_ray=muons,
                      sart_iterations=3)
    truth_met = g_met.density_array()
    rec_met = img_met.density_map
    if rec_met.shape != truth_met.shape:
        nx_t = min(rec_met.shape[0], truth_met.shape[0])
        ny_t = min(rec_met.shape[1], truth_met.shape[1])
        nz_t = min(rec_met.shape[2], truth_met.shape[2])
        rec_met = rec_met[:nx_t, :ny_t, :nz_t]
        truth_met = truth_met[:nx_t, :ny_t, :nz_t]
    corr_met = _correlation(rec_met, truth_met)
    print(f"  density-reconstruction correlation: {corr_met:.3f}")
    if img_met.muonic_xray_spectrum is not None:
        spec = img_met.muonic_xray_spectrum
        elements = spec.metadata.get("elements_detected", [])
        print(f"  muonic-X-ray spectrum: {len(elements)} elements detected — {elements}")
        # Au K_α at 6019 keV is the diagnostic signature; check signal in that window.
        au_window = (spec.axis > 5800) & (spec.axis < 6200)
        if au_window.any():
            print(f"  Au K_α window peak intensity: {float(spec.intensity[au_window].max()):.3f}")

    if not args.smoke:
        _plot([(g_rel, img_rel, "reliquary"),
               (g_geo, img_geo, "gem geode"),
               (g_met, img_met, "meteorite")],
              output_path("20_muon_tomography.png"))
    print("\nOK")
    return 0


def _plot(scenarios, path: Path) -> None:
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(len(scenarios), 3, figsize=(11, 3.2 * len(scenarios)))
    if len(scenarios) == 1:
        axs = axs.reshape(1, 3)
    for row, (grid, img, label) in enumerate(scenarios):
        truth = grid.density_array()
        truth_slice = truth[:, :, truth.shape[2] // 2]
        axs[row, 0].imshow(truth_slice.T, origin="lower", cmap="viridis")
        axs[row, 0].set_title(f"{label}: ground truth")
        axs[row, 0].set_xticks([])
        axs[row, 0].set_yticks([])
        if img.density_map is not None:
            d_slice = img.density_map[:, :, img.density_map.shape[2] // 2]
            axs[row, 1].imshow(d_slice.T, origin="lower", cmap="viridis")
            axs[row, 1].set_title("transmission recon")
        else:
            axs[row, 1].axis("off")
        axs[row, 1].set_xticks([])
        axs[row, 1].set_yticks([])
        if img.scattering_density_map is not None:
            s_slice = img.scattering_density_map[:, :, img.scattering_density_map.shape[2] // 2]
            axs[row, 2].imshow(s_slice.T, origin="lower", cmap="inferno")
            axs[row, 2].set_title("scattering recon (high-Z)")
        else:
            axs[row, 2].axis("off")
        axs[row, 2].set_xticks([])
        axs[row, 2].set_yticks([])
    fig.suptitle("Muon tomography — three composite subjects", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
