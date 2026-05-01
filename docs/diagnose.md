# The diagnose pipeline

`diagnose(spectra)` is the unified analytic entry point: feed it whatever spectra you have, and it returns a `DiagnosticReport` with a verdict, a confidence score, the per-candidate score table, every piece of evidence collected, a step-by-step reasoning trace, and follow-up technique recommendations when confidence is low.

It deliberately avoids learned classifiers. The scoring is additive and transparent so a learner can read every step and audit every weight.

## Pipeline flow

```mermaid
stateDiagram-v2
    [*] --> CollectEvidence
    CollectEvidence: per-technique analyzers run
    CollectEvidence --> ScoreCandidates: list[Evidence]
    ScoreCandidates: additive weights per profile
    ScoreCandidates --> RankCandidates
    RankCandidates --> ConfidenceCheck
    ConfidenceCheck --> EmitVerdict: top - 2nd > 0.5
    ConfidenceCheck --> RecommendFollowUp: margin too small
    RecommendFollowUp --> EmitVerdict
    EmitVerdict --> [*]
```

## Evidence collection

Every input spectrum is dispatched to the appropriate technique's analyzer. Outputs are translated into `Evidence` items with three fields:

```python
@dataclass
class Evidence:
    technique: str
    observation: str        # human-readable
    weight: float = 1.0
    favors: tuple[str, ...] = ()        # mineral names this raises score for
    rules_out: tuple[str, ...] = ()     # mineral names this lowers score for
```

For example, a Raman spectrum that matches all six catalog peaks of `corundum_Cr3plus` produces an `Evidence(technique="raman", observation="matched 6/6 catalog peaks for ruby", weight=1.0, favors=("ruby",))`.

## Scoring rules (transparent)

| Technique | Observation | Weight |
|---|---|---|
| Raman | Each catalog mineral matched (weight = matched/total) | 1.0 × match-fraction |
| Raman | Amorphous-like envelope (broad dominant + no narrow strong peaks) | 1.5 (favours amorphous, rules out crystalline) |
| UV-VIS | Each chromophore correctly identified (multi-band requires *all* bands present) | 0.6 |
| UV-VIS | No chromophores detected | 0.3 (favours minerals with empty `uvvis_bands_nm`) |
| XRF | All "major" elements of a profile present | 0.4 × N_majors |
| XRF | Each detected non-matrix element (excludes Al/Si/Ca/Mg/Na/K) | 0.3 |
| LIBS | Same major-set check | 0.3 × N_majors |
| LIBS | Each detected non-matrix element | 0.25 |
| EPR | Top center matches with cosine > 0.5 | 0.7 |
| LA-ICP-MS | Diagnostic isotope set present | 0.5 |

Light elements (H, He, Li, Be, B, C, N, O, F, Ne) are excluded from the XRF major-set check because they sit below the typical XRF detector window. A profile that lists Be as `major` (e.g. aquamarine) still passes the check on its other majors (Al + Si).

## Ranking and confidence

After all evidence is collected, the pipeline:

1. Ranks every catalog entry by total score.
2. Selects the top candidate as the verdict (or `None` if the top score is ≤ 0).
3. Computes confidence as `top / (top + second)`, then bumps it by 0.1 per favoring evidence item, clipped to 1.0.

Two candidates that tie on score produce a low-margin warning in the follow-up recommendations.

## Follow-up recommendations

The pipeline emits human-readable advice when:

- A spectrum from a particular technique was missing.
- The top two candidates differ by less than 0.5 score units.
- Confidence is below 1.0.

Examples include "to distinguish ruby from spinel, run LA-ICP-MS for Cr/Mg ratio" or "missing techniques: epr, laicpms".

## Worked example — capstone tsavorite

Synthesize Raman + UV-VIS + XRF + LIBS for a `tsavorite` profile and feed all four to `diagnose()`. The pipeline:

1. **Raman** detects six peaks at 372/549/879/1006/etc. → matches `tsavorite` and `grossular` and `andradite` (all share the garnet pattern); favors all three.
2. **UV-VIS** sees bands at 430 + 605 nm → recognises `Cr3+ d-d (emerald/alexandrite)` and `V3+ d-d (tsavorite, V-emerald)` chromophores. Tsavorite lists V³⁺ in its profile; grossular and andradite do not.
3. **XRF** detects Ca + Al + Si majors → favors all Ca-Al silicates (tsavorite, grossular, etc.).
4. **LIBS** detects Ca + Al + V → V is non-matrix, favours minerals listing V (tsavorite); the V evidence breaks the tie between tsavorite, grossular, and andradite.

Final verdict: **tsavorite** with confidence ~1.0. The reasoning trace shows every step explicitly. See `examples/19_unknown_stone_capstone.py` for the runnable version.

## Limits and caveats

The pipeline is honest about its weaknesses:

- **Quartz colour-treatment family** — rock crystal, citrine (natural and heat-treated), amethyst, and smoky quartz all share Raman, and several share EPR centres (E1', Al-hole). The pipeline frequently confuses them. The educational point is that *quartz colour-centre history is a genuinely hard case*; a single technique can't fully resolve it. Real labs use complementary FTIR + photoluminescence.
- **Synthetic data fidelity** — the catalog drives the synthesis, so a `diagnose_profile(profile)` round-trip is by construction biased toward the right answer. Real instrument data with drift, polyatomic interferences, and matrix-induced sensitivity changes will degrade scores.
- **Tied verdicts under noise** — when two minerals share most diagnostics (e.g. ruby vs sapphire_blue both being corundum + chromophore), the tie-breaker is the trace-element evidence. Rules out are *not* hard exclusions; they apply weights only.

For a top-level architectural overview see [`architecture.md`](architecture.md). For per-technique deep dives, [`techniques.md`](techniques.md). For the full curriculum walk-through, [`curriculum.md`](curriculum.md).
