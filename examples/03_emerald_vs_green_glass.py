"""Emerald vs green glass — distinguishing crystalline beryl from amorphous imitation.

Story: a parcel of green stones is offered as 'emerald'. Some are real beryl
coloured by Cr3+; some are green glass passed off as gem material.

Two complementary signatures separate them cleanly:

  Raman:
    - Beryl: sharp ring breathing mode at ~685 cm-1 + Si-O at ~1067 cm-1.
    - Glass: broad amorphous Si-O envelope (~450 cm-1, ~800 cm-1) — no sharp peaks.

  UV-VIS:
    - Emerald: Cr3+ d-d bands (~430 nm, ~600 nm) with red transmission window.
    - Green glass: featureless or unrelated colourant absorption.

Sources: Hagemann et al. 1990 (beryl Raman), Wood & Nassau 1968 (Cr3+ in beryl).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (  # noqa: E402
    BERYL_RAMAN,
    GLASS_RAMAN_ENVELOPE,
    install_raman_references,
    output_path,
    parse_smoke_args,
)

from checkmsg import raman, uvvis  # noqa: E402
from checkmsg.synthetic import PeakSpec, generate  # noqa: E402

# Cr3+ in beryl produces two broad bands plus the red-window transmission.
EMERALD_UVVIS = [
    PeakSpec(430.0, intensity=0.85, sigma=18.0, gamma=8.0),
    PeakSpec(605.0, intensity=0.95, sigma=22.0, gamma=10.0),
]
# Green glass coloured by Fe2+/Fe3+ or Cu — broad single envelope, no Cr3+ doublet.
GREEN_GLASS_UVVIS = [
    PeakSpec(660.0, intensity=0.70, sigma=70.0, gamma=30.0),
]


def synth_raman_emerald(seed: int):
    axis = np.linspace(100.0, 1500.0, 1401)
    return generate(BERYL_RAMAN, axis, technique="raman", units="cm-1", noise=0.012, seed=seed)


def synth_raman_glass(seed: int):
    axis = np.linspace(100.0, 1500.0, 1401)
    return generate(GLASS_RAMAN_ENVELOPE, axis, technique="raman", units="cm-1", noise=0.012, seed=seed)


def synth_uvvis(peaks: list[PeakSpec], seed: int):
    axis = np.linspace(380.0, 800.0, 421)
    return generate(peaks, axis, technique="uvvis", units="nm", noise=0.01, seed=seed)


def main() -> int:
    args = parse_smoke_args("03_emerald_vs_green_glass")
    install_raman_references()

    samples = {
        "Stone X (claimed emerald)": (synth_raman_emerald(11), synth_uvvis(EMERALD_UVVIS, 21)),
        "Stone Y (claimed emerald)": (synth_raman_glass(12), synth_uvvis(GREEN_GLASS_UVVIS, 22)),
    }

    print("=== Scenario 3: Emerald vs green glass (Raman + UV-VIS) ===\n")
    verdicts: dict[str, str] = {}
    for label, (rspec, uspec) in samples.items():
        rresult = raman.analyze(rspec, candidates=["beryl", "diamond", "moissanite", "corundum", "ruby"])
        amorphous = raman.is_amorphous(rspec)
        uresult = uvvis.assign_bands(uspec, polarity="absorbance", min_snr=4.0)

        print(f"{label}:")
        print(f"  Raman top: {rresult.best.mineral} (cosine={rresult.best.cosine:.3f}, "
              f"peak_score={rresult.best.peak_score:.3f})")
        print(f"  Raman amorphous-test: {amorphous}")
        print(f"  UV-VIS bands at: {[round(b.position, 1) for b in uresult.bands]}")
        chromophores = uresult.chromophores()
        print(f"  UV-VIS chromophores: {[c.name for c in chromophores] or ['none recognised']}")

        is_beryl = (
            not amorphous
            and rresult.best.mineral == "beryl"
            and rresult.best.combined > 0.4
        )
        has_cr3_band = any("Cr3+" in c.name for c in chromophores)

        if amorphous:
            verdict = "imitation (amorphous — green glass or paste)"
        elif is_beryl and has_cr3_band:
            verdict = "natural emerald (beryl + Cr3+ chromophore)"
        elif is_beryl:
            verdict = "beryl (no Cr3+ chromophore — possibly aquamarine or pale stone)"
        else:
            verdict = "unknown (does not match beryl)"
        print(f"  => verdict: {verdict}\n")
        verdicts[label] = verdict

    assert "natural emerald" in verdicts["Stone X (claimed emerald)"]
    assert "imitation" in verdicts["Stone Y (claimed emerald)"]

    if not args.smoke:
        _plot(samples, output_path("03_emerald_vs_green_glass.png"))
    print("OK")
    return 0


def _plot(samples: dict, path: Path) -> None:
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for label, (r, _u) in samples.items():
        axs[0].plot(r.axis, r.intensity, label=label, linewidth=0.8)
    axs[0].set_xlabel("Raman shift (cm$^{-1}$)")
    axs[0].set_title("Raman: crystalline (beryl) vs amorphous (glass)")
    axs[0].legend(fontsize=8)

    for label, (_r, u) in samples.items():
        axs[1].plot(u.axis, u.intensity, label=label, linewidth=1.0)
    axs[1].set_xlabel("wavelength (nm)")
    axs[1].set_title("UV-VIS: Cr$^{3+}$ doublet vs broad glass colourant")
    axs[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
