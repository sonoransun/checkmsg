"""Generate the six instrument-physics schematic figures used in docs/techniques.md.

Each schematic is a stylised cartoon — not a CAD-accurate drawing — that
illustrates the physics flow from probe → sample → detector for one of the
toolkit's analytical techniques.

Style: white background, navy/orange/green palette, 12pt labels, no gridlines.
Every figure is reproducible (no RNG) and saves to `docs/figures/` at 150 dpi.

Usage:
    python tools/build_schematics.py [--output docs/figures]
    python tools/build_schematics.py --check        # verify committed PNGs match
"""

from __future__ import annotations

import argparse
import filecmp
import sys
import tempfile
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

NAVY = "#1f3a5f"
ORANGE = "#e07a3c"
GREEN = "#2f7a3a"
GREY = "#666666"


def _box(ax, x, y, w, h, label, *, fill=NAVY, text="white"):
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04",
                                  linewidth=1.4, edgecolor=NAVY, facecolor=fill, alpha=0.9)
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
            fontsize=10.5, color=text, weight="bold")


def _arrow(ax, x0, y0, x1, y1, *, color=NAVY, label="", label_above=True):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops={"arrowstyle": "-|>", "color": color, "lw": 1.8, "shrinkA": 4, "shrinkB": 4})
    if label:
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        ax.text(mx, my + (0.06 if label_above else -0.18), label,
                ha="center", va="bottom" if label_above else "top",
                fontsize=9, color=color, style="italic")


def _setup_ax(ax, title):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=13, weight="bold", color=NAVY, loc="left", pad=10)


def draw_raman(ax):
    _setup_ax(ax, "Raman scattering — vibrational fingerprint")
    _box(ax, 0.3, 2.3, 1.6, 1.0, "532 nm\nlaser", fill=GREEN)
    _box(ax, 4.0, 2.3, 1.6, 1.0, "sample", fill=GREY)
    _box(ax, 7.6, 2.3, 2.0, 1.0, "Czerny-Turner\nspectrometer", fill=NAVY)
    _arrow(ax, 1.9, 2.8, 4.0, 2.8, color=GREEN, label="ν₀")
    _arrow(ax, 5.6, 2.8, 7.6, 2.8, color=ORANGE, label="ν₀ ± νphonon")
    # Energy-level inset
    ax.plot([2.5, 4.0], [4.6, 4.6], color=GREY, lw=1.2)
    ax.plot([2.5, 4.0], [4.95, 4.95], color=GREY, lw=1.2)
    ax.plot([2.5, 4.0], [5.5, 5.5], color=GREY, lw=1.2, linestyle="--")  # virtual
    ax.text(2.3, 4.6, "v=0", fontsize=8, ha="right", va="center", color=GREY)
    ax.text(2.3, 4.95, "v=1", fontsize=8, ha="right", va="center", color=GREY)
    ax.text(2.3, 5.5, "virtual", fontsize=8, ha="right", va="center", color=GREY)
    _arrow(ax, 2.8, 4.6, 2.8, 5.5, color=GREEN)
    _arrow(ax, 3.5, 5.5, 3.5, 4.95, color=ORANGE)
    ax.text(5.0, 5.0, "Stokes shift\n= mode frequency", fontsize=9, color=NAVY, va="center")
    ax.text(5.0, 0.7, "→ Raman shift in cm⁻¹ identifies the mineral",
            fontsize=10, color=NAVY, weight="bold", ha="center", va="center")


def draw_xrf(ax):
    _setup_ax(ax, "X-ray fluorescence — elemental fingerprint")
    _box(ax, 0.3, 2.3, 1.7, 1.0, "X-ray\ntube (Rh)", fill=NAVY)
    _box(ax, 4.0, 2.3, 1.6, 1.0, "sample", fill=GREY)
    _box(ax, 7.6, 2.3, 2.0, 1.0, "SDD\ndetector", fill=NAVY)
    _arrow(ax, 2.0, 2.8, 4.0, 2.8, color=NAVY, label="primary X-ray")
    _arrow(ax, 5.6, 2.8, 7.6, 2.8, color=ORANGE, label="K_α / K_β fluorescence")
    # Atom: nucleus + shells
    cx, cy = 4.8, 4.8
    ax.add_patch(patches.Circle((cx, cy), 0.06, color=NAVY))
    for r in (0.4, 0.7):
        ax.add_patch(patches.Circle((cx, cy), r, fill=False, edgecolor=GREY, lw=1.0))
    ax.text(cx + 0.55, cy + 0.55, "K", fontsize=9, color=GREY)
    ax.text(cx + 0.85, cy + 0.85, "L", fontsize=9, color=GREY)
    _arrow(ax, cx - 0.35, cy + 0.05, cx - 0.85, cy + 0.05, color=GREEN, label_above=False)
    ax.text(cx - 0.6, cy - 0.3, "ejected e⁻", fontsize=8, color=GREEN, ha="center")
    _arrow(ax, cx + 0.7, cy + 0.4, cx + 0.4, cy + 0.0, color=ORANGE)
    ax.text(cx + 0.95, cy + 0.65, "L → K\nphoton", fontsize=8, color=ORANGE)
    ax.text(5.0, 0.7, "→ peak energy in keV identifies the element",
            fontsize=10, color=NAVY, weight="bold", ha="center", va="center")


def draw_libs(ax):
    _setup_ax(ax, "LIBS — laser-induced plasma emission")
    _box(ax, 0.3, 2.3, 2.0, 1.0, "pulsed Nd:YAG\n(1064 nm, ns)", fill=GREEN)
    _box(ax, 4.0, 2.3, 1.6, 1.0, "sample", fill=GREY)
    _box(ax, 7.6, 2.3, 2.0, 1.0, "echelle\nspectrometer", fill=NAVY)
    _arrow(ax, 2.3, 2.8, 4.0, 2.8, color=GREEN, label="laser pulse")
    # Plasma plume — concentric ellipses
    for r, alpha in [(0.7, 0.15), (0.5, 0.25), (0.3, 0.4)]:
        ellipse = patches.Ellipse((4.95, 3.7), r * 1.5, r, color=ORANGE, alpha=alpha)
        ax.add_patch(ellipse)
    ax.text(4.95, 4.4, "plasma\nT~10⁴ K", fontsize=9, ha="center", va="center",
            color=ORANGE, weight="bold")
    _arrow(ax, 5.6, 2.8, 7.6, 2.8, color=ORANGE, label="atomic line emission")
    # Time gate inset
    t = np.linspace(0, 5, 100)
    sig = 0.8 * np.exp(-((t - 2.0) ** 2) / 0.4)
    ax.plot(0.5 + t * 0.5, 4.7 + sig * 0.5, color=NAVY, lw=1.2)
    ax.text(0.5, 5.5, "delay-gated detection", fontsize=8, color=NAVY)
    ax.text(5.0, 0.7, "→ emission wavelengths in nm identify elements",
            fontsize=10, color=NAVY, weight="bold", ha="center", va="center")


def draw_uvvis(ax):
    _setup_ax(ax, "UV-VIS — chromophore absorption")
    _box(ax, 0.3, 3.0, 1.5, 0.8, "tungsten\n+ deuterium", fill=GREEN)
    _box(ax, 2.3, 3.0, 1.6, 0.8, "monochromator", fill=NAVY)
    _box(ax, 4.5, 2.6, 1.5, 1.5, "sample\ncuvette", fill=GREY)
    _box(ax, 6.7, 3.0, 1.5, 0.8, "PMT\ndetector", fill=NAVY)
    _arrow(ax, 1.8, 3.4, 2.3, 3.4, color=GREEN)
    _arrow(ax, 3.9, 3.4, 4.5, 3.4, color=GREEN, label="λ-selected light")
    _arrow(ax, 6.0, 3.4, 6.7, 3.4, color=ORANGE, label="transmitted")
    # Absorption spectrum inset
    wl = np.linspace(0, 6, 200)
    spec = np.exp(-((wl - 2.0) ** 2) / 0.6) + 0.5 * np.exp(-((wl - 4.0) ** 2) / 0.8)
    ax.plot(1.5 + wl * 0.7, 1.0 + spec * 0.5, color=ORANGE, lw=1.5)
    ax.text(1.5, 0.5, "λ", fontsize=9, color=GREY)
    ax.text(1.0, 1.5, "abs", fontsize=9, color=GREY, rotation=90)
    ax.text(5.5, 1.2, "Cr³⁺ d-d bands", fontsize=9, color=ORANGE, va="center")
    ax.text(5.0, 5.0, "→ absorption bands in nm identify chromophores",
            fontsize=10, color=NAVY, weight="bold", ha="center", va="center")


def draw_epr(ax):
    _setup_ax(ax, "EPR — electron spin resonance")
    _box(ax, 0.3, 2.3, 1.7, 1.0, "X-band\nklystron\n(9.5 GHz)", fill=GREEN)
    _box(ax, 3.5, 2.3, 2.5, 1.0, "resonant cavity\n+ sample\n+ field B", fill=GREY)
    _box(ax, 7.5, 2.3, 2.0, 1.0, "lock-in\ndetector", fill=NAVY)
    _arrow(ax, 2.0, 2.8, 3.5, 2.8, color=GREEN, label="μW")
    _arrow(ax, 6.0, 2.8, 7.5, 2.8, color=ORANGE, label="absorption")
    # Field-modulation coils + magnet
    ax.add_patch(patches.Rectangle((3.4, 1.3), 0.3, 1.0, color=NAVY))
    ax.add_patch(patches.Rectangle((5.8, 1.3), 0.3, 1.0, color=NAVY))
    ax.text(4.75, 1.6, "B₀", fontsize=11, ha="center", color=NAVY, weight="bold")
    # Zeeman split inset
    ax.plot([2.0, 3.0], [4.4, 4.4], color=GREY)
    ax.plot([3.5, 4.5], [4.7, 4.7], color=GREY)
    ax.plot([3.5, 4.5], [4.1, 4.1], color=GREY)
    ax.text(2.5, 4.5, "B=0", fontsize=8, ha="center", color=GREY)
    ax.text(4.0, 4.85, "+½", fontsize=8, ha="center", color=GREY)
    ax.text(4.0, 3.95, "−½", fontsize=8, ha="center", color=GREY)
    _arrow(ax, 4.2, 4.1, 4.2, 4.7, color=GREEN)
    ax.text(4.5, 4.4, "hν=gμ_B B", fontsize=9, color=GREEN, va="center")
    ax.text(5.0, 5.5, "→ resonance field gives g-factor; hyperfine reveals nuclei",
            fontsize=10, color=NAVY, weight="bold", ha="center", va="center")


def draw_laicpms(ax):
    _setup_ax(ax, "LA-ICP-MS — ablation, ionisation, mass-resolved counting")
    _box(ax, 0.1, 2.3, 1.8, 1.0, "193 nm ArF\nexcimer", fill=GREEN)
    _box(ax, 2.3, 2.3, 1.4, 1.0, "ablation\ncell", fill=GREY)
    _box(ax, 4.0, 2.3, 1.4, 1.0, "Ar plasma\n(ICP)", fill=ORANGE)
    _box(ax, 5.7, 2.3, 1.5, 1.0, "quadrupole\nm/z filter", fill=NAVY)
    _box(ax, 7.5, 2.3, 2.0, 1.0, "electron\nmultiplier", fill=NAVY)
    _arrow(ax, 1.9, 2.8, 2.3, 2.8, color=GREEN)
    _arrow(ax, 3.7, 2.8, 4.0, 2.8, color=GREY, label="aerosol")
    _arrow(ax, 5.4, 2.8, 5.7, 2.8, color=ORANGE, label="ions")
    _arrow(ax, 7.2, 2.8, 7.5, 2.8, color=NAVY)
    # Time-resolved transient inset
    t = np.linspace(0, 6, 200)
    sig = np.zeros_like(t)
    sig[(t > 1) & (t < 2)] = 0.0  # blank
    sig[(t >= 2) & (t < 4.5)] = 1.0
    sig += 0.05 * np.random.RandomState(0).randn(len(t))
    ax.plot(1.0 + t * 0.8, 4.5 + sig * 0.5, color=NAVY, lw=1.2)
    ax.text(1.0, 4.2, "blank", fontsize=8, color=GREY)
    ax.text(3.5, 5.2, "sample window", fontsize=8, color=NAVY)
    ax.text(5.0, 0.7, "→ ppm concentrations + isotope ratios + depth profile",
            fontsize=10, color=NAVY, weight="bold", ha="center", va="center")


def draw_muon(ax):
    _setup_ax(ax, "Muon imaging — 3-D tomography of large composite subjects")
    _box(ax, 0.1, 2.3, 1.7, 1.0, "θ-source\n(theoretical\n10⁹ µ/s)", fill=GREEN)
    _box(ax, 2.0, 2.3, 1.0, 1.0, "collimator", fill=NAVY)
    _box(ax, 3.2, 2.3, 1.4, 1.0, "entry tracker\n(drift chambers)", fill=NAVY)
    _box(ax, 4.9, 1.8, 2.0, 2.0, "subject\n(decimeter\ncomposite)", fill=GREY)
    _box(ax, 7.2, 2.3, 1.4, 1.0, "exit tracker\n(drift chambers)", fill=NAVY)
    _box(ax, 8.8, 2.3, 1.0, 1.0, "DAQ", fill=NAVY)
    _arrow(ax, 1.8, 2.8, 2.0, 2.8, color=GREEN)
    _arrow(ax, 3.0, 2.8, 3.2, 2.8, color=GREEN, label="µ⁻")
    _arrow(ax, 4.6, 2.8, 4.9, 2.8, color=GREEN)
    _arrow(ax, 6.9, 2.8, 7.2, 2.8, color=ORANGE, label="scattered µ⁻")
    _arrow(ax, 8.6, 2.8, 8.8, 2.8, color=NAVY)
    # Inset: scattering angle measurement geometry
    ax.plot([3.6, 5.5], [4.6, 4.6], color=GREY, lw=1.2, linestyle="--")
    ax.plot([5.5, 7.0], [4.6, 4.95], color=ORANGE, lw=1.4)
    ax.plot([5.5, 7.0], [4.6, 4.6], color=GREY, lw=1.0, linestyle=":")
    ax.text(3.5, 4.65, "incoming", fontsize=8, ha="right", va="center", color=GREY)
    ax.text(7.05, 5.0, "θ_scatter", fontsize=9, color=ORANGE)
    # Inset: voxel grid
    for i in range(4):
        for j in range(4):
            colour = GREY if (i, j) != (2, 2) else ORANGE
            ax.add_patch(patches.Rectangle((1.0 + 0.2 * i, 4.5 + 0.2 * j), 0.18, 0.18,
                                            edgecolor=NAVY, facecolor=colour, alpha=0.5))
    ax.text(0.95, 5.5, "voxel grid", fontsize=8, color=NAVY)
    ax.text(5.0, 0.7, "→ density (transmission), Z² (scattering), and muonic K_α (stopping)",
            fontsize=10, color=NAVY, weight="bold", ha="center", va="center")


SCHEMATIC_BUILDERS = {
    "raman": draw_raman,
    "xrf": draw_xrf,
    "libs": draw_libs,
    "uvvis": draw_uvvis,
    "epr": draw_epr,
    "laicpms": draw_laicpms,
    "muon": draw_muon,
}


def build_one(name: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.4))
    SCHEMATIC_BUILDERS[name](ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_all(output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for name in SCHEMATIC_BUILDERS:
        out = output_dir / f"{name}_schematic.png"
        build_one(name, out)
        written.append(out)
    return written


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser("build_schematics")
    p.add_argument("--output", type=Path, default=Path("docs/figures"))
    p.add_argument("--check", action="store_true",
                   help="render to a temp dir and verify byte-identical match against committed files")
    args = p.parse_args(argv)

    if args.check:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            tmp_files = build_all(tmp_dir)
            mismatched: list[str] = []
            for tmp_file in tmp_files:
                committed = args.output / tmp_file.name
                if not committed.exists():
                    mismatched.append(f"{tmp_file.name}: missing committed copy")
                    continue
                if not filecmp.cmp(tmp_file, committed, shallow=False):
                    mismatched.append(f"{tmp_file.name}: drift detected")
            if mismatched:
                for m in mismatched:
                    print(f"  drift: {m}", file=sys.stderr)
                return 1
            print(f"all {len(tmp_files)} schematics match committed files")
            return 0

    written = build_all(args.output)
    for f in written:
        print(f"wrote {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
