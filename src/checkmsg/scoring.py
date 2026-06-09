"""Cross-cutting scoring registries, confidence bands, caveats, and conflict
detection — the single source of truth shared by the diagnosis pipeline, the
CLI, and the web service.

This module is intentionally dependency-light: it imports nothing from
`diagnose` (which would create an import cycle) and nothing from numpy/scipy.
It holds four things the rest of the toolkit composes around:

  * ``WEIGHTS`` — every additive scoring weight used by ``diagnose`` with a
    plain-language rationale and (where available) a primary-source citation.
    ``diagnose`` reads ``WEIGHTS[key].value`` instead of magic numbers, so the
    numeric behaviour is unchanged while the rationale becomes inspectable.
  * ``confidence_band`` — maps the numeric confidence onto a qualitative band.
    The band is the primary user-facing signal; the float is a *separation
    ratio*, not a calibrated probability of correctness.
  * ``CAVEATS`` / ``SERVICE_DISCLAIMER`` — the honesty markers scattered through
    the technique modules, surfaced once so every front-end can show them.
  * ``detect_conflicts`` / ``band_after_conflicts`` — flag when the evidence
    contradicts the chosen verdict and demote the qualitative band accordingly
    (without ever changing the numeric confidence or the verdict itself).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

# ---------------------------------------------------------------------------
# Confidence bands
# ---------------------------------------------------------------------------

# Qualitative bands, ordered weakest -> strongest. Used both for labelling and
# for the one-step downgrade applied when conflicts are present.
BANDS: tuple[str, ...] = ("inconclusive", "low", "medium", "high")

# Thresholds on the numeric confidence in [0, 1], checked high-to-low.
BAND_THRESHOLDS: tuple[tuple[float, str], ...] = (
    (0.85, "high"),
    (0.65, "medium"),
    (0.40, "low"),
)


def confidence_band(confidence: float) -> str:
    """Map a numeric confidence in [0, 1] to a qualitative band.

    Bands — not the raw float — are the primary user-facing signal. The float
    is a *separation ratio* (how cleanly the top candidate out-scores the
    runner-up, nudged up when several independent techniques agree). It is NOT
    a calibrated P(correct): there is no validation set, so we deliberately do
    not present it as a probability.

    Cutoffs (operational, not frequency-calibrated):
      * ``high``        >= 0.85 — winner dominates and/or many techniques agree.
      * ``medium``      >= 0.65 — clear winner with a non-trivial runner-up.
      * ``low``         >= 0.40 — near-tie or a single weak technique.
      * ``inconclusive`` < 0.40 — no verdict, or worse-than-coin-flip separation.
    """
    for cutoff, band in BAND_THRESHOLDS:
        if confidence >= cutoff:
            return band
    return "inconclusive"


def downgrade_band(band: str, steps: int = 1) -> str:
    """Lower a qualitative band by ``steps`` (clamped at ``inconclusive``)."""
    try:
        idx = BANDS.index(band)
    except ValueError:
        return band
    return BANDS[max(0, idx - steps)]


# ---------------------------------------------------------------------------
# Weight rationale registry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class WeightSpec:
    """One additive scoring weight with its rationale and optional citation."""

    value: float
    rationale: str
    citation: str = ""


# NOTE: the ``value`` of every entry is copied verbatim from the literals that
# previously lived inline in ``diagnose.py``. ``tests/test_scoring.py`` asserts
# this parity so the registry refactor stays behaviour-preserving.
WEIGHTS: dict[str, WeightSpec] = {
    "raman.dominant_peak": WeightSpec(
        1.0,
        "The dominant Raman peak localises the phase but does not by itself "
        "favour a particular catalog name; recorded as an informational anchor.",
    ),
    "raman.amorphous": WeightSpec(
        1.5,
        "Amorphous-vs-crystalline is a strong structural discriminator: a broad "
        "envelope with no sharp modes rules out every crystalline catalog entry.",
    ),
    "raman.match_unit": WeightSpec(
        1.0,
        "Per-mineral fingerprint credit, scaled by matched/total catalog peaks; "
        "a full fingerprint earns the full unit.",
        "RRUFF reference spectra",
    ),
    "raman.no_peaks": WeightSpec(
        0.5,
        "Absence of any Raman peak weakly rules out crystalline catalog entries.",
    ),
    "uvvis.chromophore": WeightSpec(
        0.6,
        "A crystal-field chromophore band is diagnostic of the transition-metal "
        "centre, but several minerals share a chromophore, so the weight is moderate.",
        "Burns 1993; Fritsch & Rossman 1987-88",
    ),
    "uvvis.no_chromophore": WeightSpec(
        0.3,
        "A featureless UV-VIS spectrum weakly favours colourless / empty-band profiles.",
    ),
    "xrf.detected": WeightSpec(
        0.4,
        "Generic 'elements detected' anchor; informational, not name-specific.",
    ),
    "xrf.major_set_unit": WeightSpec(
        0.4,
        "Per-major-element credit awarded when a profile's full major set is "
        "present (multiplied by the number of majors).",
        "NIST X-ray line table",
    ),
    "xrf.per_element": WeightSpec(
        0.3,
        "A non-matrix diagnostic element is present (ubiquitous Al/Si/Ca/Mg/Na/K "
        "are excluded); discriminates e.g. Cr (ruby) from Ti (blue sapphire).",
    ),
    "libs.detected": WeightSpec(
        0.3,
        "Generic LIBS emission anchor; informational.",
        "NIST Atomic Spectra Database",
    ),
    "libs.major_set_unit": WeightSpec(
        0.3,
        "Per-major credit for LIBS, which recovers light elements XRF cannot see "
        "(multiplied by the number of majors).",
    ),
    "libs.per_element": WeightSpec(
        0.25,
        "A non-matrix diagnostic LIBS line; weighted below the XRF equivalent "
        "because optical-emission line identification is noisier.",
    ),
    "epr.center_match": WeightSpec(
        0.7,
        "A paramagnetic centre matched with cosine > 0.5 is highly specific to a "
        "particular defect or transition-metal ion.",
        "Loubser & van Wyk 1978; Manenkov & Prokhorov 1956",
    ),
    "laicpms.isotope_set": WeightSpec(
        0.5,
        "The diagnostic isotope set for a mineral is present; strong but "
        "mode-dependent (see the LA-ICP-MS caveat).",
        "Pearce et al. 1997; IUPAC 2021",
    ),
    "squid.ordering": WeightSpec(
        0.7,
        "Magnetic ordering type is the single most discriminative SQUID signal: "
        "it favours every entry with the same ordering and rules out every entry "
        "with a different non-empty ordering.",
        "Dunlop & Ozdemir 1997",
    ),
    "squid.tc_tn": WeightSpec(
        0.5,
        "The observed Curie/Neel temperature is within 5% of a catalog value.",
        "Dunlop & Ozdemir 1997",
    ),
    "squid.saturation": WeightSpec(
        0.4,
        "The observed saturation moment is within 20% of a catalog value.",
    ),
    "pl.line_match": WeightSpec(
        0.7,
        "A sharp defect zero-phonon line (NV/SiV/N3) is highly specific to a "
        "colour centre and to synthetic vs natural diamond.",
        "Zaitsev 2001; Davies 1977",
    ),
    "pl.synthetic_marker": WeightSpec(
        0.5,
        "The Si-V centre (~737 nm) is a near-definitive CVD/HPHT synthetic-diamond "
        "growth marker.",
        "Wang et al. 2012",
    ),
    "ftir.diamond_type": WeightSpec(
        0.7,
        "Diamond nitrogen-aggregation type (Ia/Ib/IIa/IIb) is a primary "
        "classification; favours the matching type and rules out other types.",
        "Field 1992; Zaitsev 2001",
    ),
    "ftir.band_match": WeightSpec(
        0.5,
        "A diagnostic IR absorption band (beryl OH/H2O, nephrite/serpentine OH) "
        "supports the host chemistry.",
        "Wood & Nassau 1968; Farmer 1974",
    ),
    "ftir.polymer_flag": WeightSpec(
        0.5,
        "A C-H stretch near 2870-2970 cm-1 flags polymer/resin impregnation (a "
        "treatment, e.g. B-jade).",
        "Fritsch et al. 1992",
    ),
    "mossbauer.valence": WeightSpec(
        0.7,
        "Fe2+ vs Fe3+ from the isomer shift is the single strongest Mössbauer "
        "discriminator; rules out entries with the wrong dominant valence.",
        "Burns 1993; Amthauer et al. 1976",
    ),
    "mossbauer.site_match": WeightSpec(
        0.5,
        "The observed (isomer shift, quadrupole splitting) doublet matches a "
        "catalog Fe site within tolerance.",
        "Amthauer et al. 1976; Goldman et al. 1978",
    ),
    "cl.band_match": WeightSpec(
        0.5,
        "A cathodoluminescence emission band (band-A, Cr3+ red, Mn2+ orange, REE) "
        "corroborates the host but is shared across several minerals.",
        "Marshall 1988; Gaft et al. 2005",
    ),
}


# ---------------------------------------------------------------------------
# Caveats / provenance registry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Caveat:
    """A per-technique honesty marker surfaced to users at runtime."""

    technique: str
    severity: str  # "info" | "approximation" | "formula-only" | "pedagogical"
    measurement_accurate: bool
    text: str
    source_marker: str  # "file:line" pointer to the in-code honesty marker


CAVEATS: dict[str, Caveat] = {
    "raman": Caveat(
        "raman", "info", True,
        "Raman peak matching is reference-spectrum based and reliable for phase "
        "identification. Synthetic demo spectra are noise-clean by construction.",
        "diagnose.py:14",
    ),
    "uvvis": Caveat(
        "uvvis", "approximation", True,
        "Chromophore band centres are approximate; the matching tolerance covers "
        "crystal-field shifts across host minerals.",
        "refdata/chromophores.py:22",
    ),
    "xrf": Caveat(
        "xrf", "approximation", False,
        "XRF here is a relative indicator (peak-area fractions), not "
        "fundamental-parameter concentration; light elements (Z < 11) are invisible.",
        "xrf.py:61",
    ),
    "libs": Caveat(
        "libs", "approximation", False,
        "LIBS line intensities are relative; matrix effects and self-absorption "
        "are not modelled, so treat element presence, not absolute concentration, "
        "as reliable.",
        "libs.py:5",
    ),
    "epr": Caveat(
        "epr", "formula-only", False,
        "EPR absolute spin counts are formula-correct but not measurement "
        "calibrated; treat g-factor / centre identification, not concentration, "
        "as the reliable output.",
        "epr.py:500",
    ),
    "laicpms": Caveat(
        "laicpms", "approximation", False,
        "Without a matched calibration glass and internal standard, LA-ICP-MS "
        "concentrations are sensitivity-only estimates that degrade for "
        "off-matrix samples.",
        "laicpms.py:226",
    ),
    "squid": Caveat(
        "squid", "pedagogical", False,
        "SQUID forward models are deterministic, pedagogical approximations "
        "(including a toy mT->emu/g factor); the ordering type is meaningful, but "
        "absolute moments are illustrative.",
        "squid.py:31",
    ),
    "laser": Caveat(
        "laser", "approximation", False,
        "Laser excitation effects are physically-motivated approximations for "
        "synthesis, not first-principles models.",
        "laser.py:19",
    ),
    "temperature": Caveat(
        "temperature", "approximation", False,
        "Thermometric and peak-broadening helpers are heuristic but physically "
        "grounded.",
        "temperature.py:13",
    ),
    "muon": Caveat(
        "muon", "pedagogical", False,
        "Muon imaging uses textbook formulas valid to ~5-10%; hadronic "
        "interactions and decay-in-flight are out of scope, and it is not fed "
        "into diagnose().",
        "muon/physics.py:5",
    ),
    "pl": Caveat(
        "pl", "approximation", False,
        "PL line positions are tabulated zero-phonon-line values; synthetic "
        "spectra are noise-clean Voigt lines, not true ZPL + phonon-sideband "
        "shapes, and ignore temperature/strain shifts.",
        "pl.py:1",
    ),
    "ftir": Caveat(
        "ftir", "approximation", False,
        "FTIR diamond-type assignment uses textbook band positions; real spectra "
        "need quantitative one-phonon deconvolution for nitrogen concentration.",
        "ftir.py:1",
    ),
    "mossbauer": Caveat(
        "mossbauer", "pedagogical", False,
        "Doublet/sextet shapes are synthesised from tabulated isomer shift and "
        "quadrupole splitting with Lorentzian lines; no thickness, texture, or "
        "recoil-free-fraction correction.",
        "mossbauer.py:1",
    ),
    "cl": Caveat(
        "cl", "pedagogical", False,
        "CL band centres are representative; natural-vs-synthetic growth-zoning "
        "discrimination genuinely requires imaging, not the 1-D spectrum here.",
        "cl.py:1",
    ),
    "calibration": Caveat(
        "calibration", "pedagogical", False,
        "The calibrated probability is fit on noise-clean synthetic self-diagnoses "
        "via leave-mineral-out cross-validation; it estimates P(pipeline verdict "
        "equals the true catalog mineral on synthetic data), not a field guarantee.",
        "calibrate.py:1",
    ),
    "similarity": Caveat(
        "similarity", "pedagogical", False,
        "Spectral-embedding search retrieves the nearest SYNTHETIC catalog "
        "references; it is a look-alike aid and does not change the additive verdict.",
        "similarity.py:1",
    ),
}


SERVICE_DISCLAIMER: str = (
    "Check M.S.G. is a pedagogical toolkit. The confidence band reflects how "
    "cleanly the evidence separates the top candidate from the runners-up — it is "
    "NOT a probability that the verdict is correct. Reference data come from "
    "primary literature, but several modules are physically-motivated "
    "approximations rather than measurement-calibrated instruments. Confirm any "
    "high-value stone with an accredited gemological laboratory."
)


def caveats_for_techniques(techniques: set[str]) -> list[Caveat]:
    """Return the relevant caveats for the techniques that contributed evidence.

    The two SQUID spectrum literals (``squid-mh`` / ``squid-chi``) map onto the
    single ``squid`` caveat, mirroring the follow-up normalisation in
    ``diagnose._follow_ups``.
    """
    norm = set(techniques)
    if norm & {"squid-mh", "squid-chi"}:
        norm = (norm - {"squid-mh", "squid-chi"}) | {"squid"}
    return [CAVEATS[t] for t in sorted(norm) if t in CAVEATS]


# ---------------------------------------------------------------------------
# Conflict detection
# ---------------------------------------------------------------------------


@runtime_checkable
class EvidenceLike(Protocol):
    """Structural type for the ``Evidence`` records produced by ``diagnose``.

    Declared here (rather than imported) to keep ``scoring`` free of an import
    cycle with ``diagnose``.
    """

    technique: str
    weight: float
    favors: tuple[str, ...]
    rules_out: tuple[str, ...]


@dataclass(frozen=True)
class Conflict:
    """A transparency-relevant disagreement in the collected evidence."""

    kind: str  # "rules_out_verdict" | "favors_non_winner" | "technique_disagreement"
    technique: str
    detail: str
    weight: float
    against: str = ""  # the conflicting candidate name, when applicable


def detect_conflicts(
    evidence: list[EvidenceLike],
    verdict: str | None,
    scores: dict[str, float],
    *,
    strong_weight: float = 0.5,
    near_fraction: float = 0.75,
) -> list[Conflict]:
    """Detect conflicts between the evidence and the chosen verdict.

    Three read-only detectors over the existing ``Evidence.favors/rules_out``:

      1. ``rules_out_verdict`` — a piece of evidence with weight >=
         ``strong_weight`` whose ``rules_out`` contains the verdict (a direct
         contradiction; this is what triggers a band downgrade).
      2. ``favors_non_winner`` — a strong piece of evidence favours a non-verdict
         candidate whose total score is within ``near_fraction`` of the winner's
         (a strong technique pointing at a near-tie runner-up).
      3. ``technique_disagreement`` — one technique favours the verdict while a
         *different* technique rules it out.

    Returns an empty list when ``verdict`` is ``None``.
    """
    if not verdict:
        return []

    conflicts: list[Conflict] = []
    top_score = scores.get(verdict, 0.0)

    # Detector 1: evidence directly ruling out the verdict.
    for ev in evidence:
        if ev.weight >= strong_weight and verdict in ev.rules_out:
            conflicts.append(Conflict(
                kind="rules_out_verdict",
                technique=ev.technique,
                detail=f"{ev.technique} evidence rules out the verdict '{verdict}'",
                weight=ev.weight,
                against=verdict,
            ))

    # Detector 2: a strong technique favours a near-tie non-winner.
    if top_score > 0.0:
        seen_against: set[str] = set()
        for ev in evidence:
            if ev.weight < strong_weight:
                continue
            for name in ev.favors:
                if name == verdict or name in seen_against:
                    continue
                if scores.get(name, 0.0) >= near_fraction * top_score:
                    seen_against.add(name)
                    conflicts.append(Conflict(
                        kind="favors_non_winner",
                        technique=ev.technique,
                        detail=(f"{ev.technique} evidence strongly favours '{name}', "
                                f"a close runner-up to '{verdict}'"),
                        weight=ev.weight,
                        against=name,
                    ))

    # Detector 3: cross-technique disagreement on the verdict.
    favoring = {ev.technique for ev in evidence if verdict in ev.favors}
    ruling_out = {ev.technique for ev in evidence if verdict in ev.rules_out}
    for a in sorted(favoring):
        for b in sorted(ruling_out):
            if a == b:
                continue
            conflicts.append(Conflict(
                kind="technique_disagreement",
                technique=b,
                detail=(f"{a} supports '{verdict}' but {b} contradicts it"),
                weight=0.0,
                against=verdict,
            ))

    return conflicts


def band_after_conflicts(band: str, conflicts: list[Conflict]) -> str:
    """Downgrade the qualitative band one step if the verdict is contradicted.

    Only a direct ``rules_out_verdict`` conflict demotes the band. The numeric
    confidence and the verdict itself are never altered — this keeps the float
    formula stable while making the user-facing band honest.
    """
    if any(c.kind == "rules_out_verdict" for c in conflicts):
        return downgrade_band(band, 1)
    return band
