"""Sapphire geographic origin from LIBS trace-element fingerprint.

Story: a parcel of blue sapphires is offered with an unverified origin claim.
Origin determines value: Kashmir > Burma >> Madagascar/Montana. The host
mineral (corundum) is the same, so Raman cannot help. Trace elements measured
by LIBS produce locality-specific patterns that classify probabilistically.

Per-origin trace-element centroids are loaded from
`examples/data/sapphire_origins.json`. We synthesize four "unknown" stones
(one from each locality), run LIBS, build feature vectors, and classify by
minimum Mahalanobis-like distance to each centroid.

Sources for centroid trends: Peucat et al. 2007 (Lithos); Sutherland et al. 2009
(Australian J. Earth Sciences); GIA Origin Determination references.

This is a teaching example, not a forensic origin determination — real
classifiers in commercial labs use additional isotopic + inclusion data.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import output_path, parse_smoke_args  # noqa: E402

from checkmsg import libs as libs_mod  # noqa: E402
from checkmsg.synthetic import PeakSpec, generate  # noqa: E402

# LIBS strong line per element (nm) — uses NIST ASD entries already bundled.
ELEMENT_LINE = {
    "Fe": 371.994,
    "Ti": 334.941,
    "Ga": 287.424,
    "Mg": 279.553,
    "V": 310.230,
    "Cr": 425.435,
}


def load_centroids() -> dict:
    p = Path(__file__).resolve().parent / "data" / "sapphire_origins.json"
    return json.loads(p.read_text())


def synth_libs_for_origin(centroid: list[float], elements: list[str], sigma: float, seed: int):
    """Build a LIBS spectrum where each element's emission line height ≈ centroid value."""
    rng = np.random.default_rng(seed)
    axis = np.linspace(250.0, 500.0, 5000)
    peaks: list[PeakSpec] = []
    for el, c in zip(elements, centroid, strict=True):
        # Add per-stone variability (geological scatter) on top of centroid.
        height = max(0.001, c + rng.normal(0.0, sigma))
        peaks.append(PeakSpec(
            position=ELEMENT_LINE[el],
            intensity=height,
            sigma=0.10,
            gamma=0.05,
        ))
    # Add Al host emission for realism (always strong in corundum).
    peaks.append(PeakSpec(396.152, intensity=2.5, sigma=0.10, gamma=0.05))
    peaks.append(PeakSpec(394.401, intensity=2.4, sigma=0.10, gamma=0.05))
    return generate(peaks, axis, technique="libs", units="nm", noise=0.002, seed=seed + 1000)


def classify(features: np.ndarray, centroids: dict[str, np.ndarray]) -> dict[str, float]:
    """Return distance from `features` to each named centroid (smaller = closer)."""
    out = {}
    for name, c in centroids.items():
        out[name] = float(np.linalg.norm(features - c))
    return out


def main() -> int:
    args = parse_smoke_args("04_sapphire_origin")
    cfg = load_centroids()
    elements = cfg["elements"]
    centroids = {k: np.asarray(v, dtype=float) for k, v in cfg["centroids"].items()}
    sigma = float(cfg["noise_sigma"])

    print("=== Scenario 4: Blue sapphire geographic origin via LIBS trace fingerprint ===\n")
    print(f"Elements used: {elements}\n")

    truth = {}
    samples = {}
    for i, origin in enumerate(centroids):
        seed = 100 + i
        spec = synth_libs_for_origin(centroids[origin], elements, sigma, seed)
        truth[f"unknown #{i+1}"] = origin
        samples[f"unknown #{i+1}"] = spec

    correct = 0
    for label, spec in samples.items():
        result = libs_mod.identify(spec, tolerance_nm=0.4)
        feature = np.array([result.height(el) for el in elements], dtype=float)
        # Renormalize to centroid scale (max non-Al height ~1.0)
        scale = max(np.max(feature), 1e-6)
        feature = feature / scale * max(np.max(c) for c in centroids.values())
        distances = classify(feature, centroids)
        ranking = sorted(distances.items(), key=lambda kv: kv[1])
        predicted = ranking[0][0]
        is_right = predicted == truth[label]
        correct += int(is_right)
        print(f"{label}  (truth = {truth[label]}):")
        print("  measured trace heights: " + ", ".join(f"{e}={v:.2f}" for e, v in zip(elements, feature, strict=True)))
        for name, d in ranking:
            mark = "<-- predicted" if name == predicted else ""
            print(f"    {name:<12} distance={d:.3f} {mark}")
        verdict = "MATCH" if is_right else "MISS"
        print(f"  => predicted: {predicted}  [{verdict}]\n")

    print(f"Accuracy on synthetic test set: {correct}/{len(samples)}")
    assert correct == len(samples), f"only {correct}/{len(samples)} correct"

    if not args.smoke:
        _plot(samples, output_path("04_sapphire_origin.png"))
    print("OK")
    return 0


def _plot(samples: dict, path: Path) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11, 5))
    for label, s in samples.items():
        ax.plot(s.axis, s.intensity, label=label, linewidth=0.6)
    ax.set_xlabel("LIBS wavelength (nm)")
    ax.set_ylabel("emission (a.u.)")
    ax.set_title("Sapphire LIBS trace-element fingerprint")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"plot saved: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
