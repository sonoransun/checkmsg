"""Shared helpers for example scripts.

Each example is self-contained: it builds the unknown spectra it needs, installs
the synthetic reference data it needs (so the demos are reproducible offline,
without requiring an actual RRUFF download), and runs the pipeline.

Peak positions and trace ratios are taken from cited literature. Spectra are
synthesized — these are not real instrument files.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from checkmsg.refdata import rruff
from checkmsg.synthetic import PeakSpec, generate

# ---------- Raman peak tables (cm-1) ----------
# Diamond: single very sharp F2g peak. Source: Solin & Ramdas 1970.
DIAMOND_RAMAN = [PeakSpec(1332.5, intensity=1.00, sigma=1.5, gamma=0.4)]

# Moissanite (6H polytype, dominant in synthetic SiC gems): folded LO/TO doublet.
# Source: Nakashima & Harima 1997.
MOISSANITE_RAMAN = [
    PeakSpec(767.0, intensity=0.9, sigma=2.0, gamma=1.0),
    PeakSpec(789.0, intensity=1.0, sigma=2.0, gamma=1.0),
    PeakSpec(965.0, intensity=0.4, sigma=2.5, gamma=1.0),
    PeakSpec(149.0, intensity=0.5, sigma=2.5, gamma=1.0),
]

# Cubic zirconia (yttria-stabilized): broad bands from disordered cation sublattice.
# Source: Pemberton 1999, Cai et al. 2003.
CZ_RAMAN = [
    PeakSpec(269.0, intensity=0.6, sigma=18.0, gamma=8.0),
    PeakSpec(471.0, intensity=1.0, sigma=22.0, gamma=10.0),
    PeakSpec(641.0, intensity=0.8, sigma=20.0, gamma=8.0),
]

# Corundum (sapphire/ruby host): A1g + Eg modes. Source: Porto & Krishnan 1967.
CORUNDUM_RAMAN = [
    PeakSpec(378.0, intensity=0.55, sigma=2.5, gamma=0.8),
    PeakSpec(417.0, intensity=1.00, sigma=2.5, gamma=0.8),
    PeakSpec(430.0, intensity=0.30, sigma=2.5, gamma=0.8),
    PeakSpec(450.0, intensity=0.20, sigma=2.5, gamma=0.8),
    PeakSpec(577.0, intensity=0.45, sigma=3.0, gamma=1.0),
    PeakSpec(645.0, intensity=0.50, sigma=2.5, gamma=0.8),
    PeakSpec(750.0, intensity=0.20, sigma=3.0, gamma=1.0),
]

# Beryl (emerald/aquamarine host): ring breathing, Si-O modes. Source: Hagemann et al. 1990.
BERYL_RAMAN = [
    PeakSpec(322.0, intensity=0.40, sigma=2.5, gamma=0.7),
    PeakSpec(398.0, intensity=0.60, sigma=2.5, gamma=0.7),
    PeakSpec(685.0, intensity=1.00, sigma=2.5, gamma=0.7),
    PeakSpec(1010.0, intensity=0.25, sigma=2.5, gamma=0.7),
    PeakSpec(1067.0, intensity=0.35, sigma=2.5, gamma=0.7),
]

# Common silicate-glass amorphous envelope (Si-O-Si bending/stretching).
GLASS_RAMAN_ENVELOPE = [
    PeakSpec(450.0, intensity=0.7, sigma=70.0, gamma=30.0),
    PeakSpec(800.0, intensity=0.3, sigma=80.0, gamma=30.0),
    PeakSpec(1080.0, intensity=0.25, sigma=70.0, gamma=30.0),
]


def install_raman_references() -> None:
    """Populate the RRUFF cache with literature-derived synthetic references.

    Lets the examples run fully offline. In a real lab you'd let `rruff.fetch`
    pull genuine RRUFF spectra; the library API is identical either way.
    """
    axis = np.linspace(100, 1700, 1601)
    bundle = {
        "diamond": DIAMOND_RAMAN,
        "moissanite": MOISSANITE_RAMAN,
        "corundum": CORUNDUM_RAMAN,
        "ruby": CORUNDUM_RAMAN,
        "beryl": BERYL_RAMAN,
    }
    for mineral, peaks in bundle.items():
        spec = generate(peaks, axis, technique="raman", units="cm-1", noise=0.001, seed=0)
        rruff.install_synthetic_fallback(mineral, spec)


def parse_smoke_args(name: str) -> argparse.Namespace:
    p = argparse.ArgumentParser(name)
    p.add_argument("--smoke", action="store_true", help="Skip plotting; for CI smoke tests.")
    return p.parse_args()


def output_path(filename: str) -> Path:
    out = Path(__file__).resolve().parent / "output"
    out.mkdir(parents=True, exist_ok=True)
    return out / filename


def run_carousel(title: str, specimens: list[str], *,
                 require_correct: int | None = None,
                 epr_freq_GHz: float = 9.5,
                 seed_offset: int = 0):
    """Run `diagnose_profile` on each specimen and print a verdict table.

    Returns a list of (specimen_name, report) tuples for downstream plotting.
    """
    from checkmsg import minerals
    from checkmsg.diagnose import diagnose_profile

    print(f"=== {title} ===\n")
    print(f"  {'specimen':<28} {'verdict':<28} confidence  match")
    print("  " + "-" * 78)
    correct = 0
    rows: list = []
    for i, name in enumerate(specimens):
        profile = minerals.get(name)
        report = diagnose_profile(profile, epr_frequency_GHz=epr_freq_GHz, seed=seed_offset + i)
        verdict = report.verdict or "?"
        match = "yes" if verdict == name else "no"
        if verdict == name:
            correct += 1
        print(f"  {name:<28} {verdict:<28} {report.confidence:>5.2f}      {match}")
        rows.append((name, report, profile))
    print(f"\n  accuracy: {correct}/{len(specimens)}")
    threshold = require_correct if require_correct is not None else len(specimens) - 1
    assert correct >= threshold, f"too few correct verdicts ({correct}/{len(specimens)})"
    return rows


def plot_raman_carousel(rows, path, title: str):
    """Stack the synthesized Raman spectra of every specimen on one axis with offsets."""
    import matplotlib.pyplot as plt

    from checkmsg import minerals as _mins
    fig, ax = plt.subplots(figsize=(10, max(5, 0.7 * len(rows))))
    for i, (name, _report, profile) in enumerate(rows):
        try:
            spec = _mins.synthesize_raman(profile, noise=0.005)
        except Exception:
            continue
        offset = (len(rows) - 1 - i) * 1.05
        ax.plot(spec.axis, spec.intensity + offset, linewidth=0.8, label=name)
        ax.text(spec.axis[0] - 50, offset + 0.05, name, fontsize=8, ha="right")
    ax.set_xlabel("Raman shift (cm$^{-1}$)")
    ax.set_yticks([])
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")
