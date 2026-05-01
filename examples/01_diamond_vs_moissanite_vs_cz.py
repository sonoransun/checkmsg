"""Diamond vs moissanite vs cubic zirconia — the classic three-way Raman test.

Story: a buyer presents three colorless brilliants claimed to be diamond. One is
genuine, one is moissanite (6H-SiC), one is cubic zirconia. A 532 nm Raman
spectrometer separates them in seconds:

  - diamond: razor-sharp 1332 cm-1 line, nothing else
  - moissanite: folded-mode doublet at 767/789 cm-1 + 149 cm-1 acoustic + 965
  - cubic zirconia: only broad bands (~269, 471, 641) — disordered cation sublattice
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (  # noqa: E402
    CZ_RAMAN,
    DIAMOND_RAMAN,
    MOISSANITE_RAMAN,
    install_raman_references,
    output_path,
    parse_smoke_args,
)

from checkmsg import raman  # noqa: E402
from checkmsg.synthetic import PeakSpec, generate  # noqa: E402


def synth(peaks: list[PeakSpec], seed: int) -> Spectrum:  # noqa: F821
    axis = np.linspace(100, 1700, 1601)
    return generate(peaks, axis, technique="raman", units="cm-1", noise=0.012, seed=seed)


def main() -> int:
    args = parse_smoke_args("01_diamond_vs_moissanite_vs_cz")
    install_raman_references()

    samples = {
        "Stone A (claimed diamond)": synth(DIAMOND_RAMAN, seed=42),
        "Stone B (claimed diamond)": synth(MOISSANITE_RAMAN, seed=43),
        "Stone C (claimed diamond)": synth(CZ_RAMAN, seed=44),
    }

    print("=== Scenario 1: Diamond / Moissanite / CZ Raman discrimination ===\n")
    verdicts: dict[str, str] = {}
    for label, spec in samples.items():
        result = raman.analyze(spec, candidates=["diamond", "moissanite", "corundum", "beryl", "ruby"])
        print(f"{label}:")
        for c in result.candidates[:3]:
            print(f"   {c.mineral:<12} cosine={c.cosine:.3f}  peak_score={c.peak_score:.3f}  combined={c.combined:.3f}")
        # CZ test: dominant Raman feature is broad (FWHM > 25 cm-1) with poor RRUFF match.
        dominant = max(result.peaks, key=lambda p: p.height) if result.peaks else None
        cz_like = dominant is not None and dominant.width > 25.0 and (
            not result.best or result.best.combined < 0.40)
        if cz_like:
            verdicts[label] = "cubic_zirconia"
            print(f"   -> dominant peak FWHM={dominant.width:.1f} cm-1 (broad), no good crystalline match")
            print("      => cubic zirconia (disordered cation sublattice).")
        elif result.best:
            verdicts[label] = result.best.mineral
            print(f"   -> top match: {result.best.mineral}  (combined={result.best.combined:.2f})")
        else:
            verdicts[label] = "unknown"
        print()

    print("Verdicts:")
    for k, v in verdicts.items():
        print(f"  {k:<32} -> {v}")

    # Sanity assertions for smoke mode
    assert verdicts["Stone A (claimed diamond)"] == "diamond"
    assert verdicts["Stone B (claimed diamond)"] == "moissanite"
    assert verdicts["Stone C (claimed diamond)"] == "cubic_zirconia"

    if not args.smoke:
        _plot(samples, output_path("01_diamond_vs_moissanite_vs_cz.png"))
    print("OK")
    return 0


def _plot(samples: dict, path: Path) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5))
    for label, s in samples.items():
        ax.plot(s.axis, s.intensity, label=label, linewidth=0.9)
    ax.set_xlabel("Raman shift (cm$^{-1}$)")
    ax.set_ylabel("intensity (a.u.)")
    ax.set_title("Diamond vs Moissanite vs CZ — 532 nm Raman")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
