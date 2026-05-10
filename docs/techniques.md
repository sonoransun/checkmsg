# Seven analytical techniques

Every technique answers a different question about a sample. This page summarises the physics, the bundled reference data, and the analyzer entry points for each of the seven techniques `checkmsg` supports.

| Technique | Question answered | Module | Units |
|---|---|---|---|
| Raman | What is the mineral / molecular structure? | `raman.py` | cm⁻¹ |
| XRF | What elements are in this sample? | `xrf.py` | keV |
| LIBS | What elements are at the ablated micro-spot? | `libs.py` | nm |
| UV-VIS | What gives the sample its colour? | `uvvis.py` | nm |
| EPR | Are there unpaired electrons? Where? | `epr.py` | mT |
| LA-ICP-MS | What concentrations and isotope ratios? | `laicpms.py` | m/z |
| SQUID | What is the bulk magnetic ordering and moment? | `squid.py` | mT (M-H) / K (χ-T) / Hz (χ-ac) |

The analyzers share a common shape: each accepts a `Spectrum`, calls technique-appropriate preprocessing, detects features, matches against bundled reference data, and returns a structured result.

---

## Raman spectroscopy

![Raman scattering schematic](figures/raman_schematic.png)

Raman scattering measures the inelastic frequency shift between an incident laser photon and a phonon-perturbed scattered photon. The relevant observable is the wavenumber shift ν̃_shift = 1/λ₀ − 1/λ_s (where λ₀ is the laser wavelength and λ_s the scattered wavelength), which depends only on the vibrational mode energies of the sample, not on the choice of excitation: the same material produces the same set of cm⁻¹ peaks on a 532 nm or an 830 nm laser. This invariance makes Raman the workhorse technique for *mineral identification*.

```mermaid
sequenceDiagram
    participant Laser as 532 nm laser
    participant Sample as gemstone
    participant Spectrometer as Czerny-Turner spectrometer
    participant Preprocess as preprocess_raman
    participant Match as catalog match
    Laser->>Sample: ν₀ photons
    Sample-->>Spectrometer: scattered ν₀ ± νphonon
    Spectrometer->>Preprocess: raw spectrum (with fluorescence)
    Preprocess->>Preprocess: ALS baseline + Sav-Gol smoothing
    Preprocess->>Match: cleaned spectrum
    Match-->>Match: cosine + peak-list match against RRUFF / catalog
```

**Bundled data**: RRUFF Raman reference spectra (fetched on demand to `~/.cache/checkmsg/rruff/`) plus the literature-cited Raman peak lists in every `MineralProfile`. Multi-laser corrections (1/λ⁴ scaling, Cr³⁺ resonance enhancement, fluorescence interference) live in `laser.py`; phonon population / blue-shift physics in `temperature.py`.

**Key API**:

```python
from checkmsg.raman import analyze
result = analyze(spec)            # ranks RRUFF candidates by combined cosine + peak score
result.best.mineral               # 'diamond', 'corundum', ...
result.best.cosine                # 0..1
```

**Worked example output** — `examples/01_diamond_vs_moissanite_vs_cz.py`:

![Diamond / Moissanite / CZ Raman](figures/examples/01_diamond_vs_moissanite_vs_cz.png)

Three colourless brilliants, three distinct Raman fingerprints. Diamond's razor-sharp 1332 cm⁻¹ F2g line, moissanite's 767/789 cm⁻¹ folded LO/TO doublet, and cubic zirconia's broad envelope at 269/471/641 cm⁻¹ separate the samples without ambiguity.

**Physics fundamentals.** The differential cross-section for spontaneous Raman scattering scales as ω₀⁴ |∂α/∂Q|², where ω₀ is the incident angular frequency and ∂α/∂Q is the polarizability derivative with respect to the normal coordinate Q of the vibrational mode. The ω₀⁴ dependence is what makes UV excitation ~80× more efficient than NIR for the same scatterer (cf. example 05). The Stokes/anti-Stokes intensity ratio I_AS/I_S = ((ω₀+ω)/(ω₀−ω))⁴ × exp(−ℏω/kT) provides an internal thermometer (`raman.infer_temperature`); for ℏω ≫ kT the ratio collapses toward zero, which is why anti-Stokes lines are typically weak at room temperature for vibrational modes above ~500 cm⁻¹.

---

## X-ray fluorescence (XRF)

![XRF schematic](figures/xrf_schematic.png)

XRF excites inner-shell electrons; outer-shell relaxation emits characteristic X-rays whose energies identify the element. Peaks in keV map to NIST-tabulated K/L transitions.

```mermaid
sequenceDiagram
    participant Tube as Rh X-ray tube
    participant Sample
    participant SDD as silicon-drift detector
    participant Preprocess as preprocess_xrf
    participant Match as NIST line table
    Tube->>Sample: primary X-rays
    Sample-->>SDD: K_α, K_β, L_α fluorescence
    SDD->>Preprocess: raw spectrum
    Preprocess->>Preprocess: SNIP baseline subtraction
    Preprocess->>Match: peak energies (keV)
    Match-->>Match: tolerance match against NIST K/L table
```

**Bundled data**: ~70 K/L characteristic lines for elements Z=11..82 in `refdata/data/nist_xray_lines.csv`.

**Key API**:

```python
from checkmsg.xrf import identify_elements
res = identify_elements(spec, tolerance_keV=0.05)
[e.element for e in res.elements]    # ['Al', 'Cr', 'V']
```

XRF struggles with light elements (Z < 11) — they sit below typical silicon-drift detector windows. For Be / Li / B detection, use LIBS or LA-ICP-MS.

**Physics fundamentals.** Characteristic K- and L-line energies follow Moseley's law E = R_y (Z − σ)² (1/n₁² − 1/n₂²), where R_y ≈ 13.606 eV is the Rydberg energy, σ is the screening constant (≈ 1 for K_α, ≈ 7.4 for L_α), and n₁ → n₂ is the inner-to-outer shell transition (n=2 → 1 for K_α, n=3 → 2 for L_α). The K_α₁/K_α₂ doublet from the spin-orbit-split 2p₃/₂ and 2p₁/₂ states is unresolved by typical SDD detectors (~130 eV resolution at Mn K_α), which is why the bundled NIST line table reports a single energy for K_α. The 0.05 keV tolerance used in `xrf.identify_elements` therefore covers natural detector broadening, not physical line width.

---

## LIBS (laser-induced breakdown spectroscopy)

![LIBS schematic](figures/libs_schematic.png)

LIBS uses a focused pulsed laser to ablate a small volume of material and excite it into a plasma. The plasma's atomic emission lines (in nm) reveal which elements are present at the ablated spot. LIBS detects light elements (Be, Li, B) that fall below the typical XRF detector window (Z < 11), complementing XRF for full-periodic-table coverage.

```mermaid
sequenceDiagram
    participant Laser as Nd:YAG ns pulse
    participant Sample
    participant Plasma
    participant Spectrometer as echelle
    participant Match as NIST ASD
    Laser->>Sample: 1064 nm ns pulse
    Sample-->>Plasma: ablation + breakdown
    Plasma-->>Spectrometer: atomic emission lines
    Spectrometer->>Match: gated spectrum (delay ≈ 1 µs)
    Match-->>Match: tolerance match against NIST ASD lines
```

**Bundled data**: a curated subset of NIST Atomic Spectra Database lines covering Al, Be, Cr, Cu, Fe, Ga, Mg, Mo, Ni, Pt, Si, Ti, V (gem-relevant set).

**Key API**:

![LIBS-driven sapphire origin](figures/examples/04_sapphire_origin.png)

Same gem family (corundum), four geographic origins → distinct trace-element fingerprints. Mahalanobis distance to bundled centroids classifies each sapphire.

**Physics fundamentals.** A focused ns laser pulse delivers ~10⁹–10¹⁰ W/cm² to the sample, ionising and heating it into a plasma at 8000–12000 K. Excited-state populations follow the Saha-Boltzmann distribution, with the line intensity I ∝ g_k A_ki exp(−E_k / kT) / U(T), where g_k is the upper-state degeneracy, A_ki the Einstein coefficient, E_k the upper-level energy, and U(T) the partition function. Detection is gated 0.5–2 µs after the pulse to skip the early continuum dominated by bremsstrahlung and recombination, sampling the plasma during the cooling window when atomic emission dominates. Calibration-free LIBS infers concentrations from line ratios under the local-thermodynamic-equilibrium (LTE) assumption; the toolkit's `libs.identify` takes a simpler matched-line approach, leaving quantitation to LA-ICP-MS.

---

## UV-VIS spectroscopy

![UV-VIS schematic](figures/uvvis_schematic.png)

Absorbance vs wavelength identifies the electronic transition responsible for the gem's colour: a d-d transition, an intervalence charge transfer (IVCT), or a band-edge absorption that removes a wavelength range from the transmitted light. The complementary colour to the absorbed band is what the eye sees. The bundled chromophore table (`refdata/chromophores.py`) maps band-position patterns to species like Cr³⁺ d-d, Fe²⁺/Ti⁴⁺ IVCT, and V³⁺ d-d.

```mermaid
sequenceDiagram
    participant Lamps as W + D₂ lamps
    participant Mono as monochromator
    participant Sample
    participant PMT as detector
    participant Match as chromophore table
    Lamps->>Mono: broadband UV-VIS
    Mono->>Sample: λ-selected light
    Sample-->>PMT: transmitted intensity
    PMT->>Match: -log₁₀(I/I₀) vs λ
    Match-->>Match: multi-band chromophore lookup with tolerance
```

**Key API**:

```python
from checkmsg.uvvis import assign_bands
res = assign_bands(spec)
[c.name for c in res.chromophores()]    # ['Cr3+ d-d (emerald/alexandrite)']
```

**Worked output** — `examples/03_emerald_vs_green_glass.py`:

![Emerald vs green glass](figures/examples/03_emerald_vs_green_glass.png)

Real emerald shows the Cr³⁺ d-d doublet (~430 + 605 nm in beryl host). Green glass shows a broad amorphous absorption with no chromophore signature.

**Physics fundamentals.** Absorbance follows the Beer–Lambert law A = ε(λ) c ℓ, where ε is the molar absorptivity, c the chromophore concentration, and ℓ the path length. For transition-metal d-d transitions the band centre and oscillator strength are governed by the crystal-field splitting parameter Δ_oct (or Δ_tet for tetrahedral coordination); for Cr³⁺ in octahedral coordination, the ⁴A₂g → ⁴T₂g (¹) and ⁴A₂g → ⁴T₁g (²) transitions appear near 17500 cm⁻¹ and 24000 cm⁻¹ in ruby (572 nm and 416 nm respectively). The same Cr³⁺ ion in beryl experiences a smaller Δ_oct and shifts both bands to longer wavelength (~605 + 430 nm), which is why ruby and emerald — same ion, different host — produce subtly different chromophore patterns. Spin-forbidden d-d transitions (e.g. Mn²⁺) are typically 10⁴ × weaker than spin-allowed and require higher concentrations or longer path lengths to detect.

---

## Electron paramagnetic resonance (EPR / ESR)

![EPR schematic](figures/epr_schematic.png)

EPR detects unpaired electrons by measuring microwave absorption versus magnetic field. Resonance occurs when hν = g μ_B B (the electron Zeeman condition), where g is the dimensionless Landé factor, μ_B the Bohr magneton, and B the static magnetic field. For a free electron g ≈ 2.00232; transition-metal centres deviate by 1–10 % depending on the local crystal field and spin-orbit coupling. Hyperfine coupling to nearby nuclear spins splits the resonance into multiplets whose pattern identifies the coupled nucleus (e.g. ¹⁴N gives a triplet, ⁵⁵Mn a sextet). The toolkit ships a bounded but real spin-Hamiltonian simulator.

```mermaid
sequenceDiagram
    participant Klystron as 9.5 GHz klystron
    participant Cavity as resonant cavity
    participant Sample
    participant Magnet as B₀ sweep
    participant LockIn as lock-in detector
    participant Match as bundled centers
    Klystron->>Cavity: μW
    Cavity->>Sample: tuned absorption
    Magnet->>Sample: ramp B₀
    Sample-->>LockIn: dχ"/dB
    LockIn->>Match: derivative spectrum
    Match-->>Match: simulate every center, cosine match
```

**Bundled data**: nine literature-cited paramagnetic centres (DPPH, free electron, P1 nitrogen in diamond, Ni-HPHT, E1' / Al-hole in quartz, Cr³⁺ in corundum, Fe³⁺ in corundum, Mn²⁺ in calcite).

**Key API**:

```python
from checkmsg.epr import analyze, simulate_field_sweep
result = analyze(spec, frequency_GHz=9.5)
result.best.name                    # 'diamond_P1', 'quartz_E1prime', ...
result.g_factors                    # [2.0024, 2.0026, 2.0028]
```

**Worked output** — `examples/06_epr_unpaired_electrons.py`:

![EPR unpaired electrons](figures/examples/06_epr_unpaired_electrons.png)

The P1 nitrogen triplet, smoky-quartz E1' singlet, and Mn²⁺ pearl sextet are all resolved by the same spin-Hamiltonian simulator — only the spin system parameters differ.

**Physics fundamentals.** The full spin Hamiltonian implemented in `epr.simulate_field_sweep` is

H = μ_B B · g · S + S · D · S + Σᵢ ( S · Aᵢ · Iᵢ − g_n,i μ_N B · Iᵢ )

where S is the electron spin operator, B the static magnetic field, g the (in general anisotropic) g-tensor, D the zero-field splitting tensor (with axial parameter D_MHz and rhombic E_MHz in the bundled centres), and the sum runs over hyperfine-coupled nuclear spins Iᵢ with coupling tensor Aᵢ and nuclear g-factor g_n,i. The simulator constructs H in the |M_S, M_I, …⟩ product basis using `np.kron`, batch-diagonalises across the field grid via `np.linalg.eigh`, applies the perpendicular CW transition selection rule (matrix elements of S_⊥ = S − (n·S)n where n is the field direction), and integrates over a θ/φ powder grid (auto-reduced for axial / isotropic systems). The output is the absorption derivative dχ″/dB, which is what a CW spectrometer with field-modulation lock-in detection records.

---

## LA-ICP-MS (laser ablation ICP mass spectrometry)

![LA-ICP-MS schematic](figures/laicpms_schematic.png)

The most sensitive technique in the toolkit, capable of resolving sub-ppm trace elements and individual isotopes. A 193 nm laser ablates a micro-spot; the aerosol is carried by argon into an inductively-coupled argon plasma (~8000 K) that fully ionises the sample; ions are mass-filtered by a quadrupole and counted on an electron multiplier. Outputs include ppm-level concentrations (calibrated against NIST glass standards via the Longerich internal-standard equation), Pb / Sr isotope ratios with mass-bias correction, U-Pb concordant ages, REE patterns, and time-resolved depth profiles.

```mermaid
sequenceDiagram
    participant Laser as 193 nm ArF excimer
    participant Sample
    participant ICP as Ar plasma (≈8000 K)
    participant Quad as quadrupole m/z filter
    participant Detector as electron multiplier
    participant Quant as Longerich quant
    Laser->>Sample: ablation pulse
    Sample-->>ICP: aerosol (carrier Ar)
    ICP->>Quad: ions
    Quad->>Detector: m/z-resolved counts
    Detector->>Quant: time-resolved transients
    Quant-->>Quant: ppm + isotope ratios + U-Pb age
```

**Bundled data**: 135 IUPAC isotope abundances, NIST SRM 612/610 preferred values (Pearce, Perkins, Westgate, Gorton, Jackson, Neal & Chenery 1997, *Geostandards Newsletter* 21:115), CI chondrite REE (McDonough & Sun 1995, *Chem. Geol.* 120:223), present-day terrestrial Pb composition (Stacey & Kramers 1975, *Earth Planet. Sci. Lett.* 26:207), and U-Pb decay constants (Steiger & Jäger 1977, *Earth Planet. Sci. Lett.* 36:359).

**Key API**:

```python
from checkmsg.laicpms import analyze, quantify, u_pb_age, pb_ratios, ree_pattern
result = analyze(sample_run, calibration=cal_run, internal_standard=("Ca", 400000.0))
result.concentrations["Mn"].ppm       # 25.4
result.isotope_ratios["207/206"]      # 0.836
result.u_pb_age_Ma                    # 100.0
```

**Worked output** — `examples/07_laicpms_complex_cases.py`:

![LA-ICP-MS complex cases](figures/examples/07_laicpms_complex_cases.png)

Pearl natural-vs-cultured discrimination (Mn quant + Pb isotope), HPHT-treated diamond detection (Fe+Co+Ni catalyst signature), surface-coating depth profile, and Cretaceous zircon U-Pb dating — all in one example.

**Physics fundamentals.** Quantitation follows the Longerich, Jackson & Günther 1996 internal-standard equation:

C_sample,i = (I_sample,i − I_blank,i) / (I_cal,i − I_blank,cal,i) × C_cal,i × ( C_IS,sample / C_IS,sample_apparent )

where C_sample,i is the unknown sample concentration of element i, I_sample,i and I_cal,i are the blank-subtracted count rates in the sample and calibration runs, C_cal,i is the certified concentration in the calibration glass (NIST SRM 612 or 610), and the bracketed term is the internal-standard correction: C_IS,sample is the *known* concentration of the internal-standard element (e.g. Ca = 400 000 ppm in calcite), and C_IS,sample_apparent is its *uncorrected* concentration computed from the same equation without the IS correction. The IS correction normalises out matrix-induced sensitivity changes between the calibration glass and the sample.

U-Pb concordant ages come from inverting the Steiger & Jäger 1977 decay equations: from the measured ²⁰⁶Pb*/²³⁸U ratio, t = (1/λ_238) ln(1 + ²⁰⁶Pb*/²³⁸U); from ²⁰⁷Pb*/²³⁵U, t = (1/λ_235) ln(1 + ²⁰⁷Pb*/²³⁵U). The two ages must agree (concordance) within the user-specified tolerance for the verdict to be issued; discordant data indicates Pb loss, common-Pb contamination, or inheritance and triggers a `ValueError`. Mass-bias correction for Sr isotopes uses the standard linear law normalised to the natural ⁸⁶Sr/⁸⁸Sr = 0.1194.

---

## Muon imaging (experimental)

![Muon imaging schematic](figures/muon_schematic.png)

Muon imaging — *muography* — uses a beam of negative muons to image the internal structure of large composite subjects. Three observables are extracted:

  - **Transmission**: the fraction of muons that survive a given line of sight maps the integrated mass thickness along that path.
  - **Multiple Coulomb scattering**: the RMS scattering angle after transit scales as 1/X₀, with the radiation length X₀ inversely proportional to Z². High-Z inclusions (Pt, Au, Pb, U) produce strong scattering even at low integrated mass.
  - **Muonic K_α emission**: a stopping negative muon cascades into K-shell orbits and emits muonic X-rays at energies that scale roughly as Z² × (m_µ/m_e) × 13.6 eV. These energies — 75 keV for C up to 6 MeV for Pb — identify the element non-destructively.

This module is **experimental** and explicitly assumes a "theoretical on-demand high muon source" delivering ~10⁹ µ/s collimated. Real surface-muon facilities (PSI, J-PARC) reach ~10⁸ µ/s; cosmic-ray muography uses ~1 µ cm⁻² min⁻¹.

```mermaid
sequenceDiagram
    participant Source as θ-source (10⁹ µ/s)
    participant Tracker_in as entry tracker
    participant Subject as composite subject
    participant Tracker_out as exit tracker
    participant DAQ
    Source->>Tracker_in: collimated µ⁻ beam
    Tracker_in->>Subject: tagged entry (position + direction)
    Subject-->>Tracker_out: scattered or stopped µ⁻
    Tracker_out->>DAQ: tagged exit (position + direction)
    DAQ-->>DAQ: per-ray transmission + scatter + stopping voxel
```

**Bundled data**: 18 materials (gem hosts + structural metals) with X₀ from Tsai 1974, *Rev. Mod. Phys.* 46:815; muonic K_α energies for 47 elements from Engfer, Schneuwly, Vuilleumier, Walter & Zehnder 1974, *At. Data Nucl. Data Tables* 14:509; physics constants from PDG 2024, *Phys. Rev. D* 110:030001.

**Key API**:

```python
from checkmsg.muon import VoxelGrid, MuonSource, analyze
g = VoxelGrid.filled((32, 32, 32), "corundum", spacing_mm=(2.0, 2.0, 2.0))
g.set_box((12, 12, 12), (20, 20, 20), "platinum")
src = MuonSource(mean_momentum_MeV=80, flux_per_s=1e9, polarity="negative")
img = analyze(g, src, transmission=True, scattering=True, muonic_xray=True)
img.density_map                        # 3-D ndarray, reconstructed density
img.scattering_density_map             # 3-D ndarray, Z²-weighted scatter density
img.muonic_xray_spectrum               # Spectrum (technique="muon-xray", units="keV")
```

**Physics fundamentals.** Energy loss follows Bethe-Bloch (PDG 28.31.5):

−dE/(ρ dx) = K (Z/A) (1/β²) [ ½ ln(2 m_e β²γ² T_max / I²) − β² − δ/2 ]

where K = 0.30707 MeV cm² mol⁻¹, m_e is the electron mass, I is the mean excitation potential, and δ is the density-effect correction (set to zero in our regime). The continuous-slowing-down approximation (CSDA) range is the integral of 1/(dE/dx) over kinetic energy from the entry point down to a low-energy cutoff; for surface muons (~30 MeV/c) the CSDA range in water is ≈ 1.3 mm.

Multiple Coulomb scattering follows the Highland approximation:

θ₀ = (13.6 MeV / β c p) × √(x / X₀) × (1 + 0.038 ln(x / X₀))

where p is the muon momentum, x is the path length expressed as a mass thickness (g/cm²), and X₀ is the material's radiation length. Scattering at the same x grows as Z²/A in the underlying X₀ formula, which is what makes scattering muography selectively sensitive to high-Z material.

Muonic atom K_α emission uses the hydrogenic Bohr-radius scaling: replacing the electron with a muon shrinks the Bohr radius by m_e/m_µ ≈ 1/207, so the K-shell binding scales as Z² × (m_µ/m_e) × R_y ≈ Z² × 2820 eV. Reduced-mass and finite-nuclear-size corrections lower the high-Z energies by 5–10 % from the point-Coulomb prediction; the bundled values (Engfer et al. 1974) include these corrections.

**Out of scope**: muon decay in flight, hadronic interactions in heavy materials, pion/kaon contamination of the beam, π⁻ atomic capture, full muonic-cascade transition probabilities (we model only the dominant K_α line), and muon spin rotation / relaxation (µSR is a distinct technique addressing a different physics regime).

**Worked output** — `examples/20_muon_tomography.py`:

![Muon tomography](figures/examples/20_muon_tomography.png)

Three composite subjects: a sealed reliquary (transmission), a gem geode with Pt inclusion (scattering), and a meteorite cross-section with Au inclusion (transmission + scattering + muonic K_α).

---

## SQUID magnetometry

![SQUID magnetometry schematic](figures/squid_schematic.png)

A SQUID — Superconducting Quantum Interference Device — measures magnetic flux at the **single-flux-quantum** sensitivity (Φ₀ = h / 2e ≈ 2.07 × 10⁻¹⁵ Wb), translating directly into bulk magnetic-moment sensitivity better than 10⁻⁹ emu. For mineral identification, this opens questions that no other bundled technique answers: bulk magnetic ordering type (ferri / ferro / canted-AFM / AFM / paramagnetic / diamagnetic), Curie / Néel temperatures, saturation moment, coercivity, and AC-susceptibility relaxation times. Two complementary acquisition modes are bundled, matching real instrument families:

* **dc-SQUID** — two-junction interferometer, quasi-static field sweep.
  *Voltage across the ring oscillates with applied flux at period Φ₀; flux-locked feedback turns this into a linear M readout.* Best moment sensitivity; produces M(H) hysteresis loops.
* **rf-SQUID** — single-junction tank-circuit, AC modulation.
  *The Josephson junction is read by a resonant LC tank; lock-in detection at the drive frequency yields χ′ + iχ″.* Better suited to χ(T) thermal sweeps and AC susceptibility χ′(ω) + iχ″(ω) frequency scans.

The "focused external magnetic contrast" the technique requires comes from **applied bias field sweeps**: a Helmholtz pair drives the sample to saturation in dc-mh, while an additional small AC modulation drives the lock-in for χ_ac. Sweeping H past the coercivity reveals the four classic hysteresis-loop diagnostics — saturation moment Ms, remanence Mr, coercivity Hc, and the slope dM/dH at high field — each of which discriminates a different family of magnetic minerals.

```mermaid
sequenceDiagram
    participant Sample as sample + bias coil
    participant Pickup as superconducting<br/>pickup loop
    participant SQUID as SQUID ring<br/>(dc 2-junction or rf 1-junction)
    participant FLL as flux-locked feedback<br/>or rf lock-in
    participant Analyze as squid.analyze
    Sample->>Pickup: H_bias + M(sample)
    Pickup->>SQUID: Φ ∝ M
    SQUID->>FLL: voltage / lock-in I+Q
    FLL-->>Analyze: M(H), χ(T), or χ_ac(ω)
    Analyze-->>Analyze: extract Hc, Ms, Tc, Weiss θ, loss-peak ω
    Analyze-->>Analyze: rank MagneticMineral candidates
```

**Bundled data**: 13 canonical `MagneticMineral` records in `refdata/squid_signatures.py` covering the magnetically-relevant subset of mineralogy: ferrimagnets (magnetite, pyrrhotite_4c), canted antiferromagnets (hematite — with Morin transition), pure antiferromagnets (goethite, ilmenite), the metallic FeNiCo HPHT-diamond catalyst, paramagnetic transition-metal cores (Cr³⁺, Mn²⁺, Fe²⁺, Fe³⁺), and three diamagnetic baselines (diamond, quartz, calcite). The `MineralProfile` records in `minerals.py` carry a parallel set of fields (`squid_ordering`, `squid_curie_K`, `squid_neel_K`, `squid_saturation_emu_g`, `squid_susceptibility_si`, `squid_coercivity_mT`, `squid_morin_K`) populated for every magnetically relevant entry. Sources: Dunlop & Özdemir 1997, *Rock Magnetism* (CUP); Hunt-Moskowitz-Banerjee 1995, AGU Reference Shelf 3; Morin 1950, *Phys. Rev.* 78:819.

**Key API**:

```python
from checkmsg import minerals, squid

profile = minerals.get("magnetite")
mh   = minerals.synthesize_squid_mh(profile)        # dc-SQUID hysteresis at 295 K
chi  = minerals.synthesize_squid_chi_T(profile)     # rf-SQUID χ(T) sweep
ac   = minerals.synthesize_squid_chi_ac(profile)    # rf-SQUID AC susceptibility

squid.extract_coercivity(mh)        # mT
squid.extract_saturation(mh)        # emu/g
squid.extract_curie_temperature(chi)  # K — locates the dχ/dT minimum
squid.fit_curie_weiss(chi)          # (C, θ_K)
squid.extract_loss_peak(ac)         # Hz (ωτ = 1)

result = squid.analyze(mh)          # ranks MagneticMineral candidates
result.best.name                    # 'magnetite'
```

**Physics fundamentals**

The Josephson constant K_J = 2e/h ≈ 4.836 × 10¹⁴ Hz/V quantises flux through any superconducting loop; SQUIDs exploit the resulting periodic V(Φ) curve to count flux quanta. A pickup loop wound around the sample translates ΔM into ΔΦ ≈ μ₀ ΔM A_loop, so a flux change of one quantum corresponds to ~10⁻⁷ emu — six orders of magnitude better than vibrating-sample magnetometry.

Three forward simulators capture the relevant physics:

  - `simulate_mh(ordering, ...)` — tanh saturation with explicit branch offsets at ±Hc gives the canonical hysteresis loop M(H) = Ms · tanh((H ∓ Hc) / Hk). Diamagnetic samples get a negative linear slope; paramagnets get a positive slope that scales as 1/T (Curie law).
  - `simulate_chi_T(ordering, ...)` — Curie-Weiss χ = C / (T − θ) above the ordering transition; sharp cusp at Tc / TN; optional Morin step at `morin_K` for canted-AFM systems. Sign of θ separates ferro/ferri (θ > 0, equal to Tc) from AFM (θ < 0).
  - `simulate_chi_ac(ordering, ...)` — Casimir-du Pré (Debye) relaxation: χ(ω) = χ_S + (χ_T − χ_S) / (1 + iωτ), so χ′(ω) interpolates between the high-frequency χ_S and the static χ_T while χ″(ω) peaks at ωτ = 1. The relaxation time follows Arrhenius blocking τ(T) = τ₀ · exp(blocking_K / T) when not pinned to a fixed `tau_s`.

`infer_ordering` chains these into a heuristic classifier: dc-mh uses the (Hc, Ms, Mr) triple; rf-chi-T uses the sign of θ from a Curie-Weiss fit and the steepness of dχ/dT to separate ordered from para; rf-chi-ac classifies by χ′ sign plus the presence of a χ″ loss peak.

**Diagnose-pipeline integration**: the pipeline (`diagnose.py`) collects ordering type as a single high-weight evidence row (+0.7) that favours every catalog entry whose `squid_ordering` matches and rules out every entry with a different non-empty ordering — this turns the "ferrimagnetic vs canted-AFM" SQUID readout into a 50× sharper score signal than any spectroscopic feature could deliver. Tc / TN match within 5 % adds +0.5; saturation-moment match within 20 % adds +0.4. See `docs/diagnose.md` for the full scoring rules.

**Out of scope**: instrument-side considerations (drift, calibration coil geometry, gradiometer pickup, He-3 cryogenics, SQUID flux noise floors, vortex pinning); pulsed magnetic field measurements; magneto-optical imaging; muon-spin-rotation-style local-field probes.

**Worked output** — `examples/21_squid_magnetic_minerals.py`:

![SQUID magnetic minerals](figures/examples/21_squid_magnetic_minerals.png)

Five scenarios in one script: dc-SQUID hysteresis carousel (magnetite vs hematite vs ilmenite vs diamond), rf-SQUID χ(T) Curie-Weiss fit on magnetite (Tc=858 K), AC χ_ac on a superparamagnetic FeCoNi cluster (HPHT-catalyst proxy), pearl freshwater-vs-saltwater AC screening (~22× χ′ contrast at 1 Hz), and a unified `diagnose()` integration where SQUID + Raman jointly identify magnetite.
