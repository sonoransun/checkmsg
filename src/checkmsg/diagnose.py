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
from enum import StrEnum

from checkmsg import epr as epr_mod
from checkmsg import glossary as glossary_mod
from checkmsg import libs as libs_mod
from checkmsg import minerals, scoring
from checkmsg import raman as raman_mod
from checkmsg import squid as squid_mod
from checkmsg import uvvis as uvvis_mod
from checkmsg import xrf as xrf_mod
from checkmsg.minerals import CATALOG, MineralProfile
from checkmsg.refdata.epr_centers import CENTERS as EPR_CENTERS
from checkmsg.scoring import WEIGHTS, Caveat, Conflict
from checkmsg.spectrum import Spectrum

# Shorthand for the scoring-weight registry (single source of truth in scoring.py).
W = WEIGHTS


class Tier(StrEnum):
    """Sophistication tier for rendering a :class:`DiagnosticReport`."""

    NOVICE = "novice"            # plain verdict + verbal confidence
    PRACTITIONER = "practitioner"  # explained evidence, jargon expanded
    EXPERT = "expert"            # full trace, raw scores, citations, caveats

# ---------- Result data primitives ----------


@dataclass
class Evidence:
    technique: str
    observation: str
    weight: float = 1.0
    favors: tuple[str, ...] = ()    # mineral names this evidence favours
    rules_out: tuple[str, ...] = ()  # mineral names this evidence rules out
    weight_key: str = ""             # key into scoring.WEIGHTS (for rationale lookup)


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
    # --- recalibration additions (all defaulted; back-compatible) ---
    conflicts: list[Conflict] = field(default_factory=list)
    confidence_band: str = "inconclusive"
    evidence_agreement: int = 0
    caveats: list[Caveat] = field(default_factory=list)
    # Opt-in calibrated probability (None unless diagnose(..., calibrated=...)).
    calibrated_confidence: float | None = None
    calibration_method: str | None = None

    def render(self, tier: Tier | str = Tier.EXPERT) -> str:
        """Human-readable report at one of three sophistication tiers.

        ``expert`` (the default) preserves the original detailed layout and is
        what the CLI and examples consume. ``practitioner`` explains the
        evidence in plain language with jargon expanded on first use, and
        ``novice`` returns a single plain-language verdict line.
        """
        tier = Tier(tier)
        if tier is Tier.NOVICE:
            return self._render_novice()
        if tier is Tier.PRACTITIONER:
            return self._render_practitioner()
        return self._render_expert()

    # ---- verdict / glossary helpers ----

    def _verdict_display(self, *, with_formula: bool = True) -> str:
        """Enrich the bare verdict name with species (and optionally formula)."""
        if not self.verdict:
            return ""
        label = self.verdict.replace("_", " ")
        try:
            p = minerals.get(self.verdict)
        except Exception:
            return label
        extra: list[str] = []
        if p.species and p.species != self.verdict:
            extra.append(p.species)
        if with_formula and p.chemical_formula:
            extra.append(p.chemical_formula)
        return f"{label} ({', '.join(extra)})" if extra else label

    def _glossary_block(self) -> list[str]:
        terms = glossary_mod.glossary_for(self)
        if not terms:
            return []
        out = ["Glossary:"]
        for t in terms:
            out.append(f"  {t.canonical} — {t.description}")
        return out

    # ---- tier renderers ----

    def _render_novice(self) -> str:
        if not self.verdict:
            return ("Most likely: no confident identification yet — "
                    "more measurements are needed.")
        label = self.verdict.replace("_", " ").capitalize()
        lines = [f"Most likely: {label} — {self.confidence_band} confidence."]
        try:
            p = minerals.get(self.verdict)
            if p.common_colors or p.species:
                colours = ", ".join(p.common_colors) if p.common_colors else "natural"
                lines.append(f"This looks like a {colours} {p.species or 'gem'}.")
        except Exception:
            pass
        if self.confidence_band in ("low", "inconclusive"):
            lines.append("This is tentative — confirm with further testing.")
        if self.conflicts:
            lines.append("Some measurements disagree, so treat this as provisional.")
        return "\n".join(lines)

    def _render_practitioner(self) -> str:
        seen: set[str] = set()
        lines: list[str] = []
        if self.verdict:
            lines.append(
                f"Verdict: {self._verdict_display()}  (confidence: {self.confidence_band})")
        else:
            lines.append("Verdict: no confident identification")
        why: list[str] = []
        for ev in self.evidence:
            if self.verdict and self.verdict in ev.favors:
                why.append(_practitioner_evidence_line(ev, seen))
        if why:
            lines.append("Why:")
            for w in why[:8]:
                lines.append(f"  - {w}")
        ruled: list[str] = []
        for st in self.reasoning_trace:
            for c in st.confusables_ruled_out:
                if c not in ruled and c != self.verdict:
                    ruled.append(c)
        if ruled:
            lines.append("Look-alikes ruled out: "
                         + ", ".join(r.replace("_", " ") for r in ruled[:8]))
        if self.conflicts:
            lines.append("Conflicts:")
            for c in self.conflicts:
                lines.append(f"  ! {glossary_mod.expand_first_use(c.detail, seen)}")
        if self.follow_up_recommendations:
            lines.append("Suggested next steps:")
            for r in self.follow_up_recommendations[:4]:
                lines.append(f"  - {r}")
        block = self._glossary_block()
        if block:
            lines.append("")
            lines.extend(block)
        return "\n".join(lines)

    def _render_expert(self) -> str:
        lines: list[str] = ["=== Check M.S.G. Diagnostic Report ==="]
        if self.verdict:
            lines.append(
                f"Verdict: {self.verdict}  (confidence: {self.confidence_band.upper()} "
                f"— separation ratio {self.confidence:.2f}, "
                f"{self.evidence_agreement} technique(s) agree)"
            )
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
        if self.conflicts:
            lines.append("Conflicts:")
            for c in self.conflicts:
                lines.append(f"  ! {c.detail}")
            lines.append("")
        if self.follow_up_recommendations:
            lines.append("Recommended follow-up:")
            for r in self.follow_up_recommendations:
                lines.append(f"  - {r}")
            lines.append("")
        block = self._glossary_block()
        if block:
            lines.extend(block)
            lines.append("")
        if self.caveats:
            lines.append("Caveats:")
            for cav in self.caveats:
                acc = "" if cav.measurement_accurate else " [not measurement-accurate]"
                lines.append(f"  ~ {cav.technique}: {cav.text}{acc}")
            lines.append("")
        if self.calibrated_confidence is not None:
            lines.append(
                f"Calibrated P(correct): {self.calibrated_confidence:.2f} "
                f"(method={self.calibration_method}; synthetic-trained — not a field guarantee)")
            lines.append("")
        lines.append(f"Disclaimer: {scoring.SERVICE_DISCLAIMER}")
        return "\n".join(lines)


def _practitioner_evidence_line(ev: Evidence, seen: set[str]) -> str:
    """Render one Evidence in plain language, expanding jargon on first use."""
    obs = ev.observation
    if obs.startswith("chromophore "):
        name = obs[len("chromophore "):]
        gt = glossary_mod.define(name)
        if gt and gt.description:
            return f"colour analysis — {gt.description}"
    return f"{ev.technique}: {glossary_mod.expand_first_use(obs, seen)}"


# ---------- Per-technique evidence collection ----------


def _evidence_from_raman(spectrum: Spectrum) -> list[Evidence]:
    """Match the spectrum's dominant peaks against catalog Raman peaks."""
    cleaned = raman_mod.preprocess_raman(spectrum)
    detected = raman_mod.detect(cleaned, min_snr=8.0)
    if not detected:
        return [Evidence(technique="raman", observation="no Raman peaks detected",
                        weight=W["raman.no_peaks"].value, weight_key="raman.no_peaks",
                        rules_out=tuple(n for n, p in CATALOG.items()
                                        if p.raman_peaks_cm and not p.is_amorphous))]
    out: list[Evidence] = []
    # Identify the dominant peak.
    biggest = max(detected, key=lambda p: p.height)
    out.append(Evidence(
        technique="raman",
        observation=f"dominant peak at {biggest.position:.1f} cm-1 (FWHM {biggest.width:.1f})",
        weight=W["raman.dominant_peak"].value, weight_key="raman.dominant_peak",
    ))
    # Test amorphous: dominant FWHM > 60 cm-1 + no narrow peak above 50% threshold.
    if biggest.width > 50.0 and not any(p.width < 15 and p.height > 0.5 * biggest.height for p in detected):
        out.append(Evidence(
            technique="raman",
            observation="amorphous-like envelope (no sharp peaks)",
            weight=W["raman.amorphous"].value, weight_key="raman.amorphous",
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
                weight=W["raman.match_unit"].value * matched / len(profile.raman_peaks_cm),
                weight_key="raman.match_unit",
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
            weight=W["uvvis.chromophore"].value, weight_key="uvvis.chromophore",
            favors=tuple(favoured),
        ))
    if not res.chromophores():
        out.append(Evidence(
            technique="uvvis",
            observation="no recognised chromophore bands",
            weight=W["uvvis.no_chromophore"].value, weight_key="uvvis.no_chromophore",
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
                            weight=W["xrf.detected"].value, weight_key="xrf.detected"))
    for name, p in CATALOG.items():
        if not p.xrf_signature:
            continue
        majors_required = {el for el, lvl in p.xrf_signature.items()
                           if lvl == "major" and el not in XRF_INVISIBLE}
        if majors_required and majors_required.issubset(detected):
            out.append(Evidence(
                technique="xrf",
                observation=f"all major elements of {name} present",
                weight=W["xrf.major_set_unit"].value * len(majors_required),
                weight_key="xrf.major_set_unit",
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
                weight=W["xrf.per_element"].value, weight_key="xrf.per_element",
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
                            weight=W["libs.detected"].value, weight_key="libs.detected"))
    for name, p in CATALOG.items():
        if not p.libs_signature:
            continue
        majors_required = {el for el, lvl in p.libs_signature.items() if lvl == "major"}
        if majors_required and majors_required.issubset(detected):
            out.append(Evidence(
                technique="libs",
                observation=f"LIBS supports {name} chemistry",
                weight=W["libs.major_set_unit"].value * len(majors_required),
                weight_key="libs.major_set_unit",
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
                weight=W["libs.per_element"].value, weight_key="libs.per_element",
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
            weight=W["epr.center_match"].value, weight_key="epr.center_match",
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
                weight=W["laicpms.isotope_set"].value, weight_key="laicpms.isotope_set",
                favors=(name,),
            ))
    return out


def _evidence_from_squid(spectrum: Spectrum) -> list[Evidence]:
    """Translate a SQUID magnetometry observation into Evidence entries.

    Ordering type (ferri/ferro/canted-AFM/AFM/paramagnetic/diamagnetic) is
    strongly discriminatory: ferrimagnetism alone separates magnetite from
    ~50 other catalog entries, and a measurable ferromagnetic moment in an
    otherwise diamagnetic host (e.g. type IIa diamond) is a near-definitive
    HPHT-treatment flag.

    Weights:
      * ordering match: +0.7 (favours every catalog entry whose
        `squid_ordering` matches; rules out every entry with a different
        non-empty ordering).
      * Tc/TN within 5 % of catalog: +0.5 (favours just that mineral).
      * saturation moment within 20 % of catalog: +0.4.
    """
    if spectrum.technique not in ("squid-mh", "squid-chi"):
        return []
    from checkmsg.squid import from_spectrum as _from_spectrum
    try:
        meas = _from_spectrum(spectrum)
        result = squid_mod.analyze(meas)
    except Exception:
        return []

    out: list[Evidence] = []
    extracted = result.extracted
    obs_ord = extracted.get("ordering", "")
    if not obs_ord:
        return out

    # Ordering match emits one Evidence: favours every profile whose
    # squid_ordering matches, rules out every profile with a *different*
    # non-empty ordering.
    favoured = [n for n, p in CATALOG.items() if p.squid_ordering == obs_ord]
    ruled_out = [n for n, p in CATALOG.items()
                 if p.squid_ordering and p.squid_ordering != obs_ord]
    out.append(Evidence(
        technique="squid",
        observation=f"magnetic ordering: {obs_ord}",
        weight=W["squid.ordering"].value, weight_key="squid.ordering",
        favors=tuple(favoured),
        rules_out=tuple(ruled_out),
    ))

    # Curie / Néel temperature match.
    tc_obs = float(extracted.get("curie_K", 0.0) or extracted.get("neel_K", 0.0))
    if tc_obs > 0.0:
        for name, p in CATALOG.items():
            tc_ref = p.squid_curie_K or p.squid_neel_K
            if tc_ref <= 0.0:
                continue
            if abs(tc_obs - tc_ref) <= 0.05 * tc_ref:
                out.append(Evidence(
                    technique="squid",
                    observation=f"ordering temperature {tc_obs:.0f} K matches {name}",
                    weight=W["squid.tc_tn"].value, weight_key="squid.tc_tn",
                    favors=(name,),
                ))

    # Saturation moment match (only meaningful for ordered phases).
    ms_obs = float(extracted.get("saturation_emu_g", 0.0))
    if ms_obs > 1.0:  # below 1 emu/g, the slope is paramagnetic noise
        for name, p in CATALOG.items():
            if p.squid_saturation_emu_g <= 0.0:
                continue
            if abs(ms_obs - p.squid_saturation_emu_g) <= 0.20 * p.squid_saturation_emu_g:
                out.append(Evidence(
                    technique="squid",
                    observation=f"saturation moment {ms_obs:.1f} emu/g matches {name}",
                    weight=W["squid.saturation"].value, weight_key="squid.saturation",
                    favors=(name,),
                ))

    return out


def _evidence_from_pl(spectrum: Spectrum) -> list[Evidence]:
    if spectrum.technique != "pl":
        return []
    from checkmsg import pl as pl_mod
    try:
        res = pl_mod.analyze(spectrum)
    except Exception:
        return []
    out: list[Evidence] = []
    seen: set[str] = set()
    for pos, bs in res.assignments:
        if bs.name in seen:
            continue
        seen.add(bs.name)
        favoured = [n for n, p in CATALOG.items()
                    if any(abs(pos - c) <= 4.0 for c in p.pl_centers)]
        if favoured:
            out.append(Evidence(
                technique="pl", observation=f"PL line {bs.name} at {pos:.0f} nm",
                weight=W["pl.line_match"].value, weight_key="pl.line_match",
                favors=tuple(favoured)))
    if res.has_synthetic_marker():
        favoured = [n for n, p in CATALOG.items()
                    if any(abs(736.6 - c) <= 4.0 for c in p.pl_centers)]
        out.append(Evidence(
            technique="pl", observation="Si-V centre present (CVD-synthetic marker)",
            weight=W["pl.synthetic_marker"].value, weight_key="pl.synthetic_marker",
            favors=tuple(favoured)))
    return out


def _evidence_from_ftir(spectrum: Spectrum) -> list[Evidence]:
    if spectrum.technique != "ftir":
        return []
    from checkmsg import ftir as ftir_mod
    try:
        res = ftir_mod.analyze(spectrum)
    except Exception:
        return []
    out: list[Evidence] = []
    if res.diamond_type:
        favoured = [n for n, p in CATALOG.items() if p.diamond_type == res.diamond_type]
        ruled = [n for n, p in CATALOG.items()
                 if p.diamond_type and p.diamond_type != res.diamond_type]
        out.append(Evidence(
            technique="ftir", observation=f"diamond IR Type {res.diamond_type}",
            weight=W["ftir.diamond_type"].value, weight_key="ftir.diamond_type",
            favors=tuple(favoured), rules_out=tuple(ruled)))
    seen: set[str] = set()
    for pos, bs in res.assignments:
        if bs.name.startswith("diamond Type") or bs.name in seen:
            continue
        seen.add(bs.name)
        favoured = [n for n, p in CATALOG.items()
                    if any(abs(pos - b) <= 20.0 for b in p.ftir_bands)]
        if not favoured:
            continue
        if "polymer" in bs.name:
            out.append(Evidence(
                technique="ftir", observation=f"polymer impregnation band ({bs.name})",
                weight=W["ftir.polymer_flag"].value, weight_key="ftir.polymer_flag",
                favors=tuple(favoured)))
        else:
            out.append(Evidence(
                technique="ftir", observation=f"FTIR band {bs.name}",
                weight=W["ftir.band_match"].value, weight_key="ftir.band_match",
                favors=tuple(favoured)))
    return out


def _evidence_from_cl(spectrum: Spectrum) -> list[Evidence]:
    if spectrum.technique != "cl":
        return []
    from checkmsg import cl as cl_mod
    try:
        res = cl_mod.analyze(spectrum)
    except Exception:
        return []
    out: list[Evidence] = []
    seen: set[str] = set()
    for pos, bs in res.assignments:
        if bs.name in seen:
            continue
        seen.add(bs.name)
        favoured = [n for n, p in CATALOG.items()
                    if any(abs(pos - b) <= bs.tolerance for b in p.cl_bands)]
        if favoured:
            out.append(Evidence(
                technique="cl", observation=f"CL band {bs.name}",
                weight=W["cl.band_match"].value, weight_key="cl.band_match",
                favors=tuple(favoured)))
    return out


def _evidence_from_mossbauer(spectrum: Spectrum) -> list[Evidence]:
    if spectrum.technique != "mossbauer":
        return []
    from checkmsg import mossbauer as moss_mod
    from checkmsg.refdata.mossbauer_sites import SITES
    try:
        res = moss_mod.analyze(spectrum)
    except Exception:
        return []
    out: list[Evidence] = []
    if not res.extracted:
        return out
    valence = res.extracted.get("valence", "")
    if valence in ("Fe2+", "Fe3+"):
        favoured = [n for n, p in CATALOG.items()
                    if any(SITES[k].valence == valence for k in p.mossbauer_sites if k in SITES)]
        ruled = [n for n, p in CATALOG.items()
                 if p.mossbauer_sites
                 and all(SITES[k].valence != valence for k in p.mossbauer_sites if k in SITES)]
        out.append(Evidence(
            technique="mossbauer",
            observation=f"dominant Fe valence {valence} (delta={res.extracted['delta']:.2f} mm/s)",
            weight=W["mossbauer.valence"].value, weight_key="mossbauer.valence",
            favors=tuple(favoured), rules_out=tuple(ruled)))
    if res.best is not None:
        site_key = res.best.name
        favoured = [n for n, p in CATALOG.items() if site_key in p.mossbauer_sites]
        if favoured:
            out.append(Evidence(
                technique="mossbauer", observation=f"Fe site matches {site_key}",
                weight=W["mossbauer.site_match"].value, weight_key="mossbauer.site_match",
                favors=tuple(favoured)))
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
    # Treat both SQUID literals as the same logical technique for follow-up purposes.
    seen_squid = bool(spectra_techniques & {"squid-mh", "squid-chi"})
    seen = set(spectra_techniques)
    if seen_squid:
        seen |= {"squid"}
    missing = {"raman", "xrf", "libs", "uvvis", "epr", "laicpms", "squid"} - seen
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
    calibrated: bool | str = False,
) -> DiagnosticReport:
    """Run the full diagnostic pipeline on a set of spectra and return a report.

    ``calibrated`` (opt-in) attaches an estimated P(correct) from the bundled
    Platt calibrator without altering the verdict, band, or numeric confidence;
    pass ``True``/``"platt"`` to enable it. The default path is unchanged.
    """
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
        elif s.technique in ("squid-mh", "squid-chi"):
            evidence.extend(_evidence_from_squid(s))
        elif s.technique == "pl":
            evidence.extend(_evidence_from_pl(s))
        elif s.technique == "ftir":
            evidence.extend(_evidence_from_ftir(s))
        elif s.technique == "cl":
            evidence.extend(_evidence_from_cl(s))
        elif s.technique == "mossbauer":
            evidence.extend(_evidence_from_mossbauer(s))

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
        agreement = 0
    else:
        verdict = top_name
        # Confidence = top score / (top + second), clipped to [0,1]; if no runner-up, use 1.
        denom = top_score + max(second_score, 0.0)
        confidence = float(top_score / denom) if denom > 0 else 0.0
        # Bias upward when many independent evidence pieces favour the verdict.
        agreement = sum(1 for ev in evidence if verdict in ev.favors)
        confidence = min(1.0, confidence + 0.1 * agreement)

    # Qualitative band (primary signal) + conflict detection. The numeric
    # confidence and the verdict are NOT changed by this; a contradiction only
    # demotes the band one step.
    band = scoring.confidence_band(confidence)
    conflicts = scoring.detect_conflicts(evidence, verdict, scores)
    band = scoring.band_after_conflicts(band, conflicts)
    caveats = scoring.caveats_for_techniques(spectra_techniques)

    trace = _build_trace(evidence, scores)
    # Surface each conflict as a trailing trace step (step numbering continues).
    for c in conflicts:
        trace.append(TraceStep(
            step=len(trace) + 1,
            technique=c.technique,
            finding=c.detail,
            implication="conflict — band downgraded / corroboration advised"
            if c.kind == "rules_out_verdict" else "conflict — corroboration advised",
        ))
    follow_ups = _follow_ups(spectra_techniques, top_score, second_score)
    if conflicts:
        follow_ups.append(
            f"Evidence conflict detected ({len(conflicts)}); corroborate the verdict "
            "with an additional technique before relying on it."
        )

    report = DiagnosticReport(
        verdict=verdict,
        confidence=confidence,
        candidate_scores={k: float(v) for k, v in scores.items()},
        evidence=evidence,
        reasoning_trace=trace,
        follow_up_recommendations=follow_ups,
        conflicts=conflicts,
        confidence_band=band,
        evidence_agreement=agreement,
        caveats=caveats,
    )
    if calibrated and verdict is not None:
        from checkmsg import calibrate
        method = "platt" if calibrated is True else str(calibrated)
        cal = calibrate.load_calibrator(method)
        if cal is not None:
            report.calibrated_confidence = round(cal.predict_proba(report), 4)
            report.calibration_method = method
    return report


def diagnose_calibrated(spectra, *, method: str = "platt", **kwargs) -> DiagnosticReport:
    """Convenience wrapper: ``diagnose(..., calibrated=method)``."""
    return diagnose(spectra, calibrated=method, **kwargs)


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
    if profile.squid_ordering:
        squid_meas = minerals.synthesize_squid_mh(profile, seed=seed)
        if squid_meas is not None:
            spectra.append(squid_mod.to_spectrum(squid_meas))
    if profile.pl_centers:
        sp = minerals.synthesize_pl(profile, seed=seed)
        if sp is not None:
            spectra.append(sp)
    if profile.ftir_bands or profile.diamond_type:
        sp = minerals.synthesize_ftir(profile, seed=seed)
        if sp is not None:
            spectra.append(sp)
    if profile.mossbauer_sites:
        sp = minerals.synthesize_mossbauer(profile, seed=seed)
        if sp is not None:
            spectra.append(sp)
    if profile.cl_bands:
        sp = minerals.synthesize_cl(profile, seed=seed)
        if sp is not None:
            spectra.append(sp)
    return diagnose(spectra, frequency_GHz=epr_frequency_GHz)
