# Technical accuracy & how to read a verdict

This page is the honest-limitations companion to the diagnosis pipeline. It
records, per module, how rigorous the underlying physics is, where the
reference data come from, how the additive scoring weights are justified, and —
most importantly — what the confidence number does and does not mean.

It is generated to stay in lockstep with code: the scoring weights below mirror
`checkmsg.scoring.WEIGHTS`, the per-module caveats mirror
`checkmsg.scoring.CAVEATS`, and the disclaimer is the verbatim
`checkmsg.scoring.SERVICE_DISCLAIMER`. `tests/test_accuracy_doc.py` fails if they
drift apart.

## Standing disclaimer

> Check M.S.G. is a pedagogical toolkit. The confidence band reflects how cleanly the evidence separates the top candidate from the runners-up — it is NOT a probability that the verdict is correct. Reference data come from primary literature, but several modules are physically-motivated approximations rather than measurement-calibrated instruments. Confirm any high-value stone with an accredited gemological laboratory.

## How to read confidence

`diagnose()` reports a numeric `confidence` **and** a qualitative
`confidence_band`. **The band is the signal you should read.** The float is a
*separation ratio* — `top_score / (top_score + runner_up_score)`, nudged up by
`+0.1` for each independent technique that agrees with the verdict. It measures
how cleanly the winner out-scores the field; it is **not** a calibrated
probability of being correct (there is no labelled validation set, and the
self-diagnosis round-trips are biased toward the right answer — see
[diagnose.md](diagnose.md)).

| Band | Numeric cutoff | Meaning |
|------|----------------|---------|
| `high` | ≥ 0.85 | Winner dominates and/or several techniques agree; robust to small weight changes. |
| `medium` | ≥ 0.65 | Clear winner with a non-trivial runner-up; corroborate for high-value stones. |
| `low` | ≥ 0.40 | Near-tie or a single weak technique; treat as a hypothesis. |
| `inconclusive` | < 0.40 | No verdict, or worse-than-coin-flip separation. |

A contradiction in the evidence (see *Conflict detection* below) demotes the
**band** by one step but never changes the numeric confidence or the verdict.

### What confidence does **not** mean
- It is **not** `P(verdict is correct)`.
- A `high` band on a sparse profile (Raman + XRF only) is still only as good as
  two techniques — coverage varies (see [README](../README.md) technique grid).
- Absolute numbers from formula-only modules (EPR spin counts, LA-ICP-MS
  concentrations on synthetic data) are illustrative, not measurement-grade.

## Per-module accuracy tiers

Mirrors `scoring.CAVEATS`. Tiers: **rigorous** (textbook formula + cited
constants), **approximation** (physically-motivated, demo-grade), **formula-only**
(correct equation, needs real-instrument calibration to trust the number),
**pedagogical** (deterministic teaching model).

| Module | Tier | Measurement-accurate? | In-code marker |
|--------|------|-----------------------|----------------|
| `raman` | rigorous (phase ID) | yes | diagnose.py |
| `uvvis` | approximation | yes (band centres approximate) | refdata/chromophores.py |
| `xrf` | approximation | no (relative, not concentration) | xrf.py |
| `libs` | approximation | no (relative line intensities) | libs.py |
| `epr` | formula-only | no (absolute spin count needs calibration) | epr.py |
| `laicpms` | approximation | no (sensitivity-only without IS) | laicpms.py |
| `squid` | pedagogical | no (toy mT→emu/g factor) | squid.py |
| `laser` | approximation | no (synthesis effects) | laser.py |
| `temperature` | approximation | no (heuristic phonon physics) | temperature.py |
| `muon` | pedagogical | no (textbook formulas, ±5–10%; not wired into diagnose) | muon/physics.py |
| `pl` | approximation | no (tabulated ZPLs; no phonon sidebands / temperature shifts) | pl.py |
| `ftir` | approximation | no (textbook band centres; no one-phonon deconvolution) | ftir.py |
| `mossbauer` | pedagogical | no (no thickness/texture/recoil-free-fraction correction) | mossbauer.py |
| `cl` | pedagogical | no (band centres only; growth zoning needs imaging) | cl.py |

Spot-checks against the literature passed: hc = 1239.842 eV·nm; μ_B/h ≈ 13.996
MHz/mT; λ₂₃₈ = 1.55125×10⁻¹⁰ yr⁻¹; magnetite Tc = 858 K; hematite Morin = 263 K.

## Reference-data provenance

| Dataset | Source | File |
|---------|--------|------|
| Isotope abundances | IUPAC 2021 | refdata/icpms_data.py |
| NIST SRM 610 / 612 | Pearce et al. 1997 | refdata/icpms_data.py |
| Chondrite REE | McDonough & Sun 1995 | refdata/icpms_data.py |
| U-Pb decay constants | Steiger & Jäger 1977 | refdata/icpms_data.py |
| ²³⁸U/²³⁵U ratio | Hiess et al. 2012 | refdata/icpms_data.py |
| Magnetic signatures | Dunlop & Özdemir 1997; O'Reilly 1984; Hunt et al. 1995; Morin 1950 | refdata/squid_signatures.py |
| EPR centres | Loubser & van Wyk 1978; Manenkov & Prokhorov 1956; Weil 1984; Bernstein 1979 | refdata/epr_centers.py |
| Chromophore bands | Burns 1993; Fritsch & Rossman 1987-88; Nassau 2001 | refdata/chromophores.py |
| XRF / LIBS lines | NIST X-ray transition energies / NIST ASD | refdata/nist_xray.py, refdata/nist_asd.py |
| Raman references | RRUFF (rruff.info) | refdata/rruff.py |

## Scoring-weight rationale

Every additive weight lives in `scoring.WEIGHTS` with a rationale. The numeric
values are unchanged from the original pipeline; the registry only makes the
reasoning inspectable (and tunable in one place).

| Weight key | Value | Rationale (abridged) |
|------------|-------|----------------------|
| `raman.dominant_peak` | 1.0 | Localises the phase; informational anchor, not name-specific. |
| `raman.amorphous` | 1.5 | Strong structural discriminator (rules out crystalline entries). |
| `raman.match_unit` | 1.0 | Per-mineral fingerprint credit, scaled by matched/total peaks. |
| `raman.no_peaks` | 0.5 | Absence of peaks weakly rules out crystalline entries. |
| `uvvis.chromophore` | 0.6 | Chromophore band is diagnostic of the ion, but shared across minerals. |
| `uvvis.no_chromophore` | 0.3 | Featureless spectrum weakly favours colourless profiles. |
| `xrf.detected` | 0.4 | Generic "elements detected" anchor. |
| `xrf.major_set_unit` | 0.4 | Per-major-element credit (×N) when the full major set is present. |
| `xrf.per_element` | 0.3 | Non-matrix diagnostic element present. |
| `libs.detected` | 0.3 | Generic LIBS emission anchor. |
| `libs.major_set_unit` | 0.3 | Per-major credit; LIBS sees light elements XRF cannot. |
| `libs.per_element` | 0.25 | Non-matrix LIBS line; below XRF (optical ID is noisier). |
| `epr.center_match` | 0.7 | A matched paramagnetic centre is highly specific. |
| `laicpms.isotope_set` | 0.5 | Diagnostic isotope set present; strong but mode-dependent. |
| `squid.ordering` | 0.7 | Magnetic ordering type — the single most discriminative SQUID signal. |
| `squid.tc_tn` | 0.5 | Curie/Néel temperature within 5% of a catalog value. |
| `squid.saturation` | 0.4 | Saturation moment within 20% of a catalog value. |
| `pl.line_match` | 0.7 | A sharp defect ZPL (NV/SiV/N3) is highly specific to a colour centre and synthetic-vs-natural diamond. |
| `pl.synthetic_marker` | 0.5 | The Si-V centre (~737 nm) is a near-definitive CVD/HPHT synthetic-diamond marker. |
| `ftir.diamond_type` | 0.7 | Diamond N-aggregation type (Ia/Ib/IIa/IIb) — favours the matching type, rules out others. |
| `ftir.band_match` | 0.5 | A diagnostic IR band (beryl OH, nephrite/serpentine) supports the host chemistry. |
| `ftir.polymer_flag` | 0.5 | A C-H stretch near 2870–2970 cm⁻¹ flags polymer/resin impregnation. |
| `mossbauer.valence` | 0.7 | Fe²⁺ vs Fe³⁺ from the isomer shift; rules out the wrong dominant valence. |
| `mossbauer.site_match` | 0.5 | The observed (δ, ΔEQ) doublet matches a catalog Fe site within tolerance. |
| `cl.band_match` | 0.5 | A CL emission band (band-A, Cr³⁺ red, Mn²⁺ orange, REE) corroborates the host. |

## Conflict detection

`scoring.detect_conflicts` flags transparency-relevant disagreements over the
existing `Evidence.favors`/`rules_out` sets:

- **`rules_out_verdict`** — a strong piece of evidence rules out the chosen
  verdict (e.g. a SQUID ordering mismatch contradicting a Raman-favoured name).
  This is the only conflict that demotes the confidence band.
- **`favors_non_winner`** — a strong technique points at a near-tie runner-up.
- **`technique_disagreement`** — one technique supports the verdict while another
  contradicts it.

Conflicts are surfaced in the reasoning trace and add a follow-up
recommendation, but they never silently change the verdict.

## Optional calibrated confidence & similarity search

Two opt-in tools sit *alongside* the transparent additive verdict and never change it:

- **`calibration`** — `diagnose(..., calibrated="platt")` attaches an estimated P(correct)
  from a Platt-scaling sigmoid. It is trained on noise-clean synthetic, leave-mineral-out
  cross-validated self-diagnoses: it recalibrates the pipeline's own separation ratio against
  the pipeline's own verdicts and is **NOT a real-world guarantee of correctness**.
- **`similarity`** — `checkmsg similar <technique>:<file>` retrieves the nearest synthetic
  catalog references by cosine of a normalised spectral embedding. It is a look-alike aid and
  emits no diagnostic Evidence (the additive verdict stays auditable).

## Consolidated limitations

- The pipeline deliberately avoids learned classifiers — it favours transparent,
  additive rules over maximum accuracy.
- The quartz colour-treatment family is a deliberately hard case (see
  [curriculum.md](curriculum.md)); EPR + UV-VIS alone cannot fully resolve
  treatment history.
- Synthetic round-trip diagnoses are biased toward the correct answer and are
  noise-clean by construction; real measurements carry drift, matrix effects,
  and instrument noise that these models do not simulate.
- Several entries are characterised on Raman + XRF only; a `high` band there
  rests on fewer independent lines of evidence than a multi-technique verdict.

## See also

- [architecture.md](architecture.md) — system overview and data flow.
- [techniques.md](techniques.md) — per-technique deep dives.
- [diagnose.md](diagnose.md) — pipeline state diagram and scoring rules.
