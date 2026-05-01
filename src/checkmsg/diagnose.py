"""Multi-technique diagnostic pipeline with explicit reasoning trace.

`diagnose(spectra)` accepts spectra from any of the six bundled techniques and
produces a `DiagnosticReport` containing:

  - a top-ranked verdict (canonical mineral name) with confidence score
  - a per-candidate score table
  - the evidence collected from each technique
  - a reasoning trace explaining how each piece of evidence narrowed the
    candidate list (which minerals it ruled in and which it ruled out)
  - follow-up technique recommendations when confidence is low

The pipeline scores every catalog entry against the observed evidence using
additive, transparent rules. It deliberately avoids learned classifiers — the
goal is pedagogical clarity, not maximum accuracy.

See `examples/19_unknown_stone_capstone.py` for a worked walk-through.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

from checkmsg import epr as epr_mod
from checkmsg import libs as libs_mod
from checkmsg import minerals
from checkmsg import raman as raman_mod
from checkmsg import uvvis as uvvis_mod
from checkmsg import xrf as xrf_mod
from checkmsg.minerals import CATALOG, MineralProfile
from checkmsg.refdata.epr_centers import CENTERS as EPR_CENTERS
from checkmsg.spectrum import Spectrum

# ---------- Result data primitives ----------


@dataclass
class Evidence:
    technique: str
    observation: str
    weight: float = 1.0
    favors: tuple[str, ...] = ()    # mineral names this evidence favours
    rules_out: tuple[str, ...] = ()  # mineral names this evidence rules out


@dataclass
class TraceStep:
    step: int
    technique: str
    finding: str
    implication: str
    confusables_ruled_out: list[str] = field(default_factory=list)


@dataclass
class DiagnosticReport:
    verdict: str | None
    confidence: float
    candidate_scores: dict[str, float]
    evidence: list[Evidence]
    reasoning_trace: list[TraceStep]
    follow_up_recommendations: list[str]

    def render(self) -> str:
        """Multi-paragraph human-readable diagnostic report."""
        lines: list[str] = ["=== Check M.S.G. Diagnostic Report ==="]
        if self.verdict:
            lines.append(f"Verdict: {self.verdict}  (confidence {self.confidence:.2f})")
        else:
            lines.append("Verdict: insufficient evidence")
        lines.append("")
        lines.append("Top-5 candidate scores:")
        ranked = sorted(self.candidate_scores.items(), key=lambda kv: kv[1], reverse=True)
        for name, score in ranked[:5]:
            mark = "  <-- verdict" if name == self.verdict else ""
            lines.append(f"  {name:<28} {score:6.2f}{mark}")
        lines.append("")
        lines.append("Reasoning trace:")
        for step in self.reasoning_trace:
            lines.append(f"  [{step.step}] {step.technique}: {step.finding}")
            lines.append(f"       -> {step.implication}")
            if step.confusables_ruled_out:
                lines.append(f"       ruled out: {', '.join(step.confusables_ruled_out)}")
        lines.append("")
        if self.follow_up_recommendations:
            lines.append("Recommended follow-up:")
            for r in self.follow_up_recommendations:
                lines.append(f"  - {r}")
        return "\n".join(lines)


# ---------- Per-technique evidence collection ----------


def _evidence_from_raman(spectrum: Spectrum) -> list[Evidence]:
    """Match the spectrum's dominant peaks against catalog Raman peaks."""
    cleaned = raman_mod.preprocess_raman(spectrum)
    detected = raman_mod.detect(cleaned, min_snr=8.0)
    if not detected:
        return [Evidence(technique="raman", observation="no Raman peaks detected", weight=0.5,
                        rules_out=tuple(n for n, p in CATALOG.items()
                                        if p.raman_peaks_cm and not p.is_amorphous))]
    out: list[Evidence] = []
    # Identify the dominant peak.
    biggest = max(detected, key=lambda p: p.height)
    out.append(Evidence(
        technique="raman",
        observation=f"dominant peak at {biggest.position:.1f} cm-1 (FWHM {biggest.width:.1f})",
        weight=1.0,
    ))
    # Test amorphous: dominant FWHM > 60 cm-1 + no narrow peak above 50% threshold.
    if biggest.width > 50.0 and not any(p.width < 15 and p.height > 0.5 * biggest.height for p in detected):
        out.append(Evidence(
            technique="raman",
            observation="amorphous-like envelope (no sharp peaks)",
            weight=1.5,
            favors=tuple(n for n, p in CATALOG.items() if p.is_amorphous),
            rules_out=tuple(n for n, p in CATALOG.items()
                          if p.raman_peaks_cm and not p.is_amorphous),
        ))
    # Catalog match: for each profile, count detected peaks within 10 cm-1 of catalog peaks.
    detected_positions = [p.position for p in detected]
    for name, profile in CATALOG.items():
        if not profile.raman_peaks_cm or profile.is_amorphous:
            continue
        matched = 0
        for pos, _rel in profile.raman_peaks_cm:
            if any(abs(pos - dp) < 10.0 for dp in detected_positions):
                matched += 1
        # Match threshold scales with peak count: 1-peak minerals need 1 match,
        # multi-peak minerals need at least half their peaks present.
        threshold = max(1, len(profile.raman_peaks_cm) // 2)
        if matched >= threshold:
            out.append(Evidence(
                technique="raman",
                observation=f"matched {matched}/{len(profile.raman_peaks_cm)} catalog peaks for {name}",
                weight=1.0 * matched / len(profile.raman_peaks_cm),
                favors=(name,),
            ))
    return out


def _evidence_from_uvvis(spectrum: Spectrum) -> list[Evidence]:
    """Identify chromophore bands and link to candidate minerals."""
    res = uvvis_mod.assign_bands(spectrum)
    out: list[Evidence] = []
    for ch in res.chromophores():
        favoured: list[str] = []
        for name, p in CATALOG.items():
            if ch.name in p.chromophores:
                favoured.append(name)
        out.append(Evidence(
            technique="uvvis",
            observation=f"chromophore {ch.name}",
            weight=0.6,
            favors=tuple(favoured),
        ))
    if not res.chromophores():
        out.append(Evidence(
            technique="uvvis",
            observation="no recognised chromophore bands",
            weight=0.3,
            favors=tuple(n for n, p in CATALOG.items() if not p.uvvis_bands_nm),
        ))
    return out


XRF_INVISIBLE: frozenset[str] = frozenset({"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"})


def _evidence_from_xrf(spectrum: Spectrum) -> list[Evidence]:
    """Compare detected XRF elements against each profile's xrf_signature.

    Light elements (Z < 11) are below the typical XRF detector window and aren't
    included in the bundled NIST X-ray line table. Profiles that list them as
    'major' (e.g. Be in aquamarine) get those entries silently excluded from
    the matching test so their other majors can still satisfy the check.
    """
    res = xrf_mod.identify_elements(spectrum, tolerance_keV=0.05, min_snr=8.0)
    detected = {e.element for e in res.elements}
    out: list[Evidence] = []
    if detected:
        out.append(Evidence(technique="xrf",
                            observation=f"elements detected: {', '.join(sorted(detected))}",
                            weight=0.4))
    for name, p in CATALOG.items():
        if not p.xrf_signature:
            continue
        majors_required = {el for el, lvl in p.xrf_signature.items()
                           if lvl == "major" and el not in XRF_INVISIBLE}
        if majors_required and majors_required.issubset(detected):
            out.append(Evidence(
                technique="xrf",
                observation=f"all major elements of {name} present",
                weight=0.4 * len(majors_required),
                favors=(name,),
            ))
    # Per-element diagnostic credit: each detected trace/minor element favors
    # every profile that lists it. This discriminates e.g. Cr (ruby) from
    # Ti (blue sapphire) when Al-major alone cannot.
    for el in detected:
        if el in ("Al", "Si", "Ca", "Mg", "Na", "K"):  # too common to be diagnostic
            continue
        favoured = [n for n, p in CATALOG.items()
                    if p.xrf_signature.get(el) in ("major", "minor", "trace")]
        if favoured and len(favoured) < 20:
            out.append(Evidence(
                technique="xrf",
                observation=f"diagnostic element {el} present",
                weight=0.3,
                favors=tuple(favoured),
            ))
    return out


def _evidence_from_libs(spectrum: Spectrum) -> list[Evidence]:
    """Compare detected LIBS elements against each profile's libs_signature."""
    res = libs_mod.identify(spectrum, tolerance_nm=0.4, min_snr=5.0)
    detected = set(res.elements.keys())
    out: list[Evidence] = []
    if detected:
        out.append(Evidence(technique="libs",
                            observation=f"emission lines: {', '.join(sorted(detected))}",
                            weight=0.3))
    for name, p in CATALOG.items():
        if not p.libs_signature:
            continue
        majors_required = {el for el, lvl in p.libs_signature.items() if lvl == "major"}
        if majors_required and majors_required.issubset(detected):
            out.append(Evidence(
                technique="libs",
                observation=f"LIBS supports {name} chemistry",
                weight=0.3 * len(majors_required),
                favors=(name,),
            ))
    # Per-element credit (excludes ubiquitous matrix elements).
    for el in detected:
        if el in ("Al", "Si", "Ca", "Mg", "Na", "K"):
            continue
        favoured = [n for n, p in CATALOG.items()
                    if p.libs_signature.get(el) in ("major", "minor", "trace")]
        if favoured and len(favoured) < 20:
            out.append(Evidence(
                technique="libs",
                observation=f"LIBS detected diagnostic {el}",
                weight=0.25,
                favors=tuple(favoured),
            ))
    return out


def _evidence_from_epr(spectrum: Spectrum, frequency_GHz: float | None) -> list[Evidence]:
    freq = frequency_GHz or spectrum.metadata.get("frequency_GHz")
    if freq is None:
        return []
    try:
        res = epr_mod.analyze(spectrum, frequency_GHz=freq, candidates=EPR_CENTERS)
    except Exception:
        return []
    out: list[Evidence] = []
    if res.best and res.best.cosine > 0.5:
        center_name = res.best.name
        favoured = [n for n, p in CATALOG.items() if center_name in p.epr_centers]
        out.append(Evidence(
            technique="epr",
            observation=f"top EPR centre: {center_name} (cosine {res.best.cosine:.2f})",
            weight=0.7,
            favors=tuple(favoured),
        ))
    return out


def _evidence_from_laicpms(spectrum: Spectrum) -> list[Evidence]:
    if spectrum.technique != "laicpms":
        return []
    isotope_keys = spectrum.metadata.get("isotope_keys", [])
    if not isotope_keys:
        return []
    out: list[Evidence] = []
    detected = set(isotope_keys)
    for name, p in CATALOG.items():
        if not p.icpms_diagnostic_isotopes:
            continue
        if set(p.icpms_diagnostic_isotopes).issubset(detected):
            out.append(Evidence(
                technique="laicpms",
                observation=f"diagnostic isotopes for {name} present",
                weight=0.5,
                favors=(name,),
            ))
    return out


# ---------- Scoring + reasoning ----------


def _aggregate_scores(evidence: list[Evidence]) -> dict[str, float]:
    scores: dict[str, float] = {n: 0.0 for n in CATALOG}
    for ev in evidence:
        for name in ev.favors:
            if name in scores:
                scores[name] += ev.weight
        for name in ev.rules_out:
            if name in scores:
                scores[name] -= ev.weight
    return scores


def _build_trace(evidence: list[Evidence], scores: dict[str, float]) -> list[TraceStep]:
    trace: list[TraceStep] = []
    for i, ev in enumerate(evidence, start=1):
        if ev.favors:
            implication = f"favors: {', '.join(ev.favors[:5])}"
            ruled = sorted(set(ev.rules_out))[:8]
        elif ev.rules_out:
            implication = f"rules out: {', '.join(ev.rules_out[:5])}"
            ruled = list(ev.rules_out)[:8]
        else:
            implication = "informational only"
            ruled = []
        trace.append(TraceStep(
            step=i,
            technique=ev.technique,
            finding=ev.observation,
            implication=implication,
            confusables_ruled_out=ruled,
        ))
    return trace


def _follow_ups(spectra_techniques: set[str], top_score: float, second_score: float) -> list[str]:
    """Recommend additional techniques if confidence is low or runner-up is close."""
    recs: list[str] = []
    margin = top_score - second_score
    if top_score < 1.0:
        recs.append("Confidence is low; provide additional spectra to narrow the diagnosis.")
    missing = {"raman", "xrf", "libs", "uvvis", "epr", "laicpms"} - spectra_techniques
    if missing:
        recs.append("Missing techniques: " + ", ".join(sorted(missing)))
    if margin < 0.5 and top_score > 0:
        recs.append("Top two candidates are close; prefer LA-ICP-MS or EPR for stable discrimination.")
    return recs


# ---------- Public API ----------


def diagnose(
    spectra: Iterable[Spectrum] | dict[str, Spectrum],
    *,
    candidates: list[str] | None = None,
    frequency_GHz: float | None = None,
) -> DiagnosticReport:
    """Run the full diagnostic pipeline on a set of spectra and return a report."""
    if isinstance(spectra, dict):
        spectra = list(spectra.values())
    spectra_list = list(spectra)
    spectra_techniques = {s.technique for s in spectra_list}

    evidence: list[Evidence] = []
    for s in spectra_list:
        if s.technique == "raman":
            evidence.extend(_evidence_from_raman(s))
        elif s.technique == "uvvis":
            evidence.extend(_evidence_from_uvvis(s))
        elif s.technique == "xrf":
            evidence.extend(_evidence_from_xrf(s))
        elif s.technique == "libs":
            evidence.extend(_evidence_from_libs(s))
        elif s.technique == "epr":
            evidence.extend(_evidence_from_epr(s, frequency_GHz))
        elif s.technique == "laicpms":
            evidence.extend(_evidence_from_laicpms(s))

    scores = _aggregate_scores(evidence)
    if candidates is not None:
        # Only consider scores for the named subset.
        scores = {n: scores.get(n, 0.0) for n in candidates}

    ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
    top_name, top_score = ranked[0] if ranked else (None, 0.0)
    second_score = ranked[1][1] if len(ranked) > 1 else 0.0

    # Verdict + confidence.
    if top_score <= 0.0:
        verdict = None
        confidence = 0.0
    else:
        verdict = top_name
        # Confidence = top score / (top + second), clipped to [0,1]; if no runner-up, use 1.
        denom = top_score + max(second_score, 0.0)
        confidence = float(top_score / denom) if denom > 0 else 0.0
        # Bias upward when many independent evidence pieces favour the verdict
        confidence = min(1.0, confidence + 0.1 * sum(
            1 for ev in evidence if verdict in ev.favors))

    trace = _build_trace(evidence, scores)
    follow_ups = _follow_ups(spectra_techniques, top_score, second_score)

    return DiagnosticReport(
        verdict=verdict,
        confidence=confidence,
        candidate_scores={k: float(v) for k, v in scores.items()},
        evidence=evidence,
        reasoning_trace=trace,
        follow_up_recommendations=follow_ups,
    )


def diagnose_profile(profile: MineralProfile, *,
                     laser_nm: float | None = None,
                     temperature_K: float = 295.0,
                     epr_frequency_GHz: float = 9.5,
                     seed: int = 0) -> DiagnosticReport:
    """Helper: synthesize all available technique spectra for a profile and diagnose.

    Useful for unit tests and for the curriculum's "synthesise + diagnose" round-trip.
    """
    spectra: list[Spectrum] = []
    if profile.raman_peaks_cm:
        spectra.append(minerals.synthesize_raman(
            profile, noise=0.005, laser_nm=laser_nm,
            temperature_K=temperature_K, seed=seed,
        ))
    if profile.uvvis_bands_nm:
        spectra.append(minerals.synthesize_uvvis(profile, seed=seed))
    if profile.xrf_signature:
        spectra.append(minerals.synthesize_xrf(profile, seed=seed))
    if profile.libs_signature:
        spectra.append(minerals.synthesize_libs(profile, seed=seed))
    if profile.epr_centers:
        epr_spec = minerals.synthesize_epr(profile, frequency_GHz=epr_frequency_GHz, seed=seed)
        if epr_spec is not None:
            spectra.append(epr_spec)
    return diagnose(spectra, frequency_GHz=epr_frequency_GHz)
