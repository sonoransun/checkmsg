"""Convert a ``DiagnosticReport`` (plain dataclasses) into tier-gated JSON.

This is the single composition point for the three workstreams: the scoring
registries (`scoring`), the glossary, and the tiered renderer all feed in here.
Lower tiers omit expert-only fields (raw candidate scores, evidence weights).
"""

from __future__ import annotations

from dataclasses import asdict

from checkmsg import glossary, minerals, scoring
from checkmsg.diagnose import DiagnosticReport, Tier


def _verdict_block(report: DiagnosticReport) -> dict | None:
    if not report.verdict:
        return None
    block = {"name": report.verdict}
    try:
        p = minerals.get(report.verdict)
    except Exception:
        return block
    block.update(
        species=p.species,
        aliases=list(p.aliases),
        chemical_formula=p.chemical_formula,
        common_colors=list(p.common_colors),
    )
    return block


def _evidence_dict(ev, tier: Tier) -> dict:
    out = {
        "technique": ev.technique,
        "observation": ev.observation,
        "favors": list(ev.favors),
    }
    if tier is Tier.EXPERT:
        out["weight"] = ev.weight
        out["rules_out"] = list(ev.rules_out)
        out["weight_key"] = ev.weight_key
        spec = scoring.WEIGHTS.get(ev.weight_key)
        out["rationale"] = spec.rationale if spec else ""
    return out


def report_to_dict(report: DiagnosticReport, *, tier: Tier | str = Tier.NOVICE) -> dict:
    """Serialise a report for the API at the requested sophistication tier."""
    tier = Tier(tier)
    payload: dict = {
        "tier": tier.value,
        "verdict": _verdict_block(report),
        "confidence": {
            "value": round(float(report.confidence), 4),
            "band": report.confidence_band,
            "note": "separation ratio, not a probability of correctness",
        },
        "summary": report.render(tier),
        "follow_ups": list(report.follow_up_recommendations),
        "conflicts": [asdict(c) for c in report.conflicts],
        "glossary": [
            {"term": t.canonical, "expansion": t.expansion, "definition": t.description}
            for t in glossary.glossary_for(report)
        ],
        "disclaimer": scoring.SERVICE_DISCLAIMER,
    }
    if tier in (Tier.PRACTITIONER, Tier.EXPERT):
        payload["evidence"] = [_evidence_dict(ev, tier) for ev in report.evidence]
        payload["reasoning_trace"] = [asdict(s) for s in report.reasoning_trace]
        payload["caveats"] = [asdict(c) for c in report.caveats]
    if tier is Tier.EXPERT:
        payload["candidate_scores"] = report.candidate_scores
        payload["evidence_agreement"] = report.evidence_agreement
    return payload
