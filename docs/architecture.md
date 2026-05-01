# Architecture

`checkmsg` is a Python toolkit for gemological analysis from spectroscopic data. The architecture is intentionally flat: a single `Spectrum` data primitive, six per-technique analyzer modules, a shared mineral catalog, and a unified diagnostic pipeline that orchestrates them.

## System overview

```mermaid
flowchart TB
    subgraph Input["Input — six technique-specific spectra"]
        S1["raman cm⁻¹"]
        S2["xrf keV"]
        S3["libs nm"]
        S4["uvvis nm"]
        S5["epr mT"]
        S6["laicpms m/z"]
    end
    subgraph Analyzers["Per-technique analyzers"]
        A1["raman.analyze"]
        A2["xrf.identify_elements"]
        A3["libs.identify"]
        A4["uvvis.assign_bands"]
        A5["epr.analyze"]
        A6["laicpms.analyze"]
    end
    Catalog[("MineralProfile<br/>CATALOG (55 entries)")]
    Diagnose["diagnose.diagnose"]
    Report[["DiagnosticReport<br/>verdict + reasoning trace"]]
    S1 --> A1
    S2 --> A2
    S3 --> A3
    S4 --> A4
    S5 --> A5
    S6 --> A6
    A1 --> Diagnose
    A2 --> Diagnose
    A3 --> Diagnose
    A4 --> Diagnose
    A5 --> Diagnose
    A6 --> Diagnose
    Catalog --> Diagnose
    Diagnose --> Report
```

Every spectrum lands in `Spectrum`, an immutable dataclass that records the axis, intensity, technique tag, units, and free-form metadata. The technique-specific analyzers consume a `Spectrum` and emit either a structured result (e.g. `RamanResult`, `XrfResult`) or a list of detected features (peaks, chromophores, isotope ratios). The diagnostic pipeline turns those analyzer outputs into typed `Evidence` items, scores them against every entry in the mineral catalog, and emits a `DiagnosticReport` with a verdict, confidence, full candidate-score table, reasoning trace, and follow-up recommendations.

## Class diagram

```mermaid
classDiagram
    class Spectrum {
        +ndarray axis
        +ndarray intensity
        +str technique
        +str units
        +dict metadata
        +slice(low, high)
        +normalize(mode)
        +with_intensity(y)
    }
    class MineralProfile {
        +str name
        +str species
        +tuple raman_peaks_cm
        +tuple uvvis_bands_nm
        +dict xrf_signature
        +dict libs_signature
        +tuple epr_centers
        +tuple confusables
    }
    class IcpmsRun {
        +dict transients
        +tuple blank_window_s
        +tuple sample_window_s
        +integrate(key, window)
        +blank_subtracted_cps()
        +to_spectrum()
        +detect_segments(key)
    }
    class SpinSystem {
        +str name
        +float S
        +float|tuple g
        +float D_MHz
        +float E_MHz
        +tuple hyperfine
    }
    class LaserConfig {
        +float wavelength_nm
        +str regime
        +float photon_eV
        +shift_to_wavelength(shift)
        +resonance_enhancement(band, sigma, max)
    }
    class MicrowaveBand {
        +str name
        +float frequency_GHz
        +resonance_field_mT(g)
        +g_at(field_mT)
    }
    class DiagnosticReport {
        +str|None verdict
        +float confidence
        +dict candidate_scores
        +list evidence
        +list reasoning_trace
        +list follow_up_recommendations
        +render() str
    }
    Spectrum <.. IcpmsRun : produces (to_spectrum)
    DiagnosticReport <.. Spectrum : input
    DiagnosticReport <.. MineralProfile : scored against
```

## Repository layout

```
checkmsg/
├── src/checkmsg/
│   ├── spectrum.py            # Spectrum dataclass — every technique speaks this
│   ├── preprocess.py          # ALS / SNIP baseline, Savitzky-Golay smoothing
│   ├── peaks.py               # peak detection, Voigt fitting via lmfit
│   ├── match.py               # cosine + peak-list matching
│   ├── synthetic.py           # spectrum generators (used by examples + tests)
│   ├── raman.py               # Raman analyzer: preprocess → peaks → RRUFF match
│   ├── xrf.py                 # XRF: peak energies → NIST K/L lines → elements
│   ├── libs.py                # LIBS: emission lines → NIST ASD → elements
│   ├── uvvis.py               # UV-VIS: chromophore band assignment
│   ├── epr.py                 # EPR spin-Hamiltonian simulator + analysis
│   ├── laicpms.py             # LA-ICP-MS quant + isotopes + U-Pb + REE
│   ├── identify.py            # combined_report multi-technique fusion
│   ├── diagnose.py            # unified diagnostic pipeline + reasoning trace
│   ├── minerals.py            # MineralProfile catalog (55 entries) + helpers
│   ├── laser.py               # 8-laser catalogue (Raman excitation)
│   ├── microwave.py           # 7-band catalogue (EPR frequencies)
│   ├── temperature.py         # LN2 / room phonon physics
│   ├── cli.py                 # `checkmsg` command-line entry point
│   └── refdata/               # bundled reference data
│       ├── chromophores.py    # UV-VIS chromophore band table
│       ├── epr_centers.py     # 9 EPR centers (DPPH, P1, E1', Cr3+, ...)
│       ├── icpms_data.py      # NIST SRM 612/610, IUPAC isotopes, chondrite REE
│       ├── nist_xray.py       # K/L line table (X-ray)
│       ├── nist_asd.py        # atomic emission lines (LIBS)
│       └── rruff.py           # RRUFF Raman fetcher with on-disk cache
├── examples/                  # 19 curriculum scripts (01..19)
├── docs/                      # this directory
├── tools/                     # build_schematics.py, build_confusables_graph.py
└── tests/                     # 200+ tests covering all modules
```

## Where to look

| Question | File |
|---|---|
| "What's a Spectrum?" | `src/checkmsg/spectrum.py` |
| "How does technique X work?" | `src/checkmsg/<technique>.py` |
| "Why is this mineral confused with that one?" | `src/checkmsg/minerals.py` |
| "How does the diagnose pipeline score candidates?" | `src/checkmsg/diagnose.py` |
| "Where are reference spectra cached?" | `~/.cache/checkmsg/rruff/` (overridable via `CHECKMSG_CACHE`) |
| "How do I make a synthetic spectrum?" | `src/checkmsg/synthetic.py` |
| "Where are the worked examples?" | `examples/` (curriculum 01..19); `docs/curriculum.md` for narrated walkthroughs |

## Data-flow narrative

A typical analysis runs like this:

1. Load or measure a `Spectrum` (axis + intensity + technique tag).
2. Run `<technique>.analyze` (or call `diagnose([spec1, spec2, ...])`).
3. The analyzer pre-processes (baseline, smoothing), detects features, and matches them against bundled reference data.
4. The analyzer's output is wrapped in `Evidence` items by `diagnose`.
5. `diagnose` walks every `MineralProfile` in `CATALOG` and accumulates evidence-derived scores.
6. The top-ranked profile becomes the verdict; the runner-ups appear as confusables ruled out.
7. `DiagnosticReport.render()` produces a multi-paragraph human-readable explanation.

For a worked end-to-end walk-through, see [`curriculum.md`](curriculum.md), specifically example 19 (the capstone).
