# Mineral catalog reference

`src/checkmsg/minerals.py` ships a catalog of 55 gemstones across 10 thematic groups. Every entry is a `MineralProfile` carrying chemical, physical, and spectroscopic fingerprints — enough to drive both the synthesis helpers and the diagnostic pipeline.

## Catalog structure

```mermaid
mindmap
  root((CATALOG<br/>55 entries))
    Diamond simulants
      diamond
      moissanite
      cubic_zirconia
      GGG
      YAG
      strontium_titanate
      white_sapphire
      white_topaz
      white_spinel
      glass_paste
    Blue stones
      sapphire_blue
      tanzanite
      iolite
      aquamarine
      blue_topaz
      blue_zircon
    Garnet group
      pyrope
      almandine
      spessartine
      grossular
      andradite
      rhodolite
      tsavorite
    Jade family
      jadeite
      nephrite
      serpentine
      aventurine_quartz
      prehnite
    Red stones
      ruby
      red_spinel
      red_beryl
      rubellite
      rhodochrosite
    Black stones
      hematite
      magnetite
      obsidian
      jet
      onyx
      schorl
    Tourmaline species
      elbaite
      dravite
      liddicoatite
    Quartz treatments
      rock_crystal
      citrine_natural
      citrine_heat_treated
      amethyst
      smoky_quartz
    Chrysoberyl
      alexandrite
      cymophane
      chrysoberyl_yellow
    Biominerals
      pearl_natural_saltwater
      pearl_akoya
      pearl_freshwater
      coral
      ivory
```

(Some minerals appear in more than one group's narrative — e.g. pyrope is both a garnet end-member and a red gem — but each has exactly one canonical entry.)

## MineralProfile schema

Every catalog entry holds:

| Field | Purpose |
|---|---|
| `name`, `species`, `aliases` | Canonical name + taxonomic species + commercial / historical names |
| `chemical_formula` | Stoichiometric formula (with dopants) |
| `crystal_system`, `mohs_hardness`, `density_g_cc`, `refractive_index` | Physical context |
| `common_colors` | Tuples for `by_color` filter |
| `raman_peaks_cm` | `((position, rel_intensity), ...)` for `synthesize_raman` and Raman matching |
| `uvvis_bands_nm` | Centers of absorption bands |
| `chromophores` | Keys into `refdata/chromophores.py` |
| `xrf_signature` / `libs_signature` | Element → `"major"`/`"minor"`/`"trace"`/`"absent"` |
| `epr_centers` | Keys into `refdata/epr_centers.py` |
| `icpms_diagnostic_isotopes` | Isotope keys diagnostic for the mineral |
| `confusables` | Names of catalog entries commonly mistaken for this one |
| `diagnostic_features` | Human-readable identifying features |
| `references` | Literature citations |

## Confusables graph (auto-generated)

Every mineral lists `confusables` — the others most easily mistaken for it. The undirected graph is built from `tools/build_confusables_graph.py` (run before each release to keep this diagram in sync with the catalog):

```mermaid
graph LR
    diamond --- moissanite
    cubic_zirconia --- diamond
    diamond --- white_sapphire
    cubic_zirconia --- moissanite
    GGG --- cubic_zirconia
    GGG --- diamond
    YAG --- diamond
    YAG --- cubic_zirconia
    GGG --- YAG
    diamond --- strontium_titanate
    moissanite --- strontium_titanate
    white_sapphire --- white_topaz
    white_sapphire --- white_spinel
    diamond --- white_topaz
    white_spinel --- white_topaz
    diamond --- white_spinel
    diamond --- glass_paste
    cubic_zirconia --- glass_paste
    sapphire_blue --- tanzanite
    iolite --- sapphire_blue
    blue_topaz --- sapphire_blue
    blue_zircon --- sapphire_blue
    iolite --- tanzanite
    aquamarine --- sapphire_blue
    aquamarine --- blue_topaz
    aquamarine --- iolite
    blue_topaz --- blue_zircon
    aquamarine --- blue_zircon
    almandine --- pyrope
    pyrope --- rhodolite
    pyrope --- ruby
    pyrope --- red_spinel
    almandine --- rhodolite
    almandine --- ruby
    almandine --- spessartine
    pyrope --- spessartine
    andradite --- grossular
    grossular --- tsavorite
    almandine --- andradite
    andradite --- tsavorite
    jadeite --- nephrite
    jadeite --- serpentine
    aventurine_quartz --- jadeite
    nephrite --- serpentine
    aventurine_quartz --- nephrite
    aventurine_quartz --- serpentine
    jadeite --- prehnite
    prehnite --- serpentine
    red_spinel --- ruby
    rubellite --- ruby
    red_spinel --- rubellite
    red_beryl --- ruby
    red_beryl --- rubellite
    almandine --- red_beryl
    rhodochrosite --- rubellite
    red_beryl --- rhodochrosite
    hematite --- magnetite
    hematite --- obsidian
    hematite --- schorl
    magnetite --- obsidian
    jet --- obsidian
    obsidian --- onyx
    hematite --- jet
    jet --- onyx
    obsidian --- schorl
    elbaite --- schorl
    dravite --- elbaite
    elbaite --- rubellite
    dravite --- schorl
    elbaite --- liddicoatite
    citrine_natural --- rock_crystal
    amethyst --- rock_crystal
    citrine_heat_treated --- citrine_natural
    amethyst --- citrine_natural
    amethyst --- citrine_heat_treated
    amethyst --- smoky_quartz
    rock_crystal --- smoky_quartz
    alexandrite --- chrysoberyl_yellow
    alexandrite --- cymophane
    chrysoberyl_yellow --- cymophane
    pearl_akoya --- pearl_natural_saltwater
    pearl_freshwater --- pearl_natural_saltwater
    coral --- pearl_natural_saltwater
    ivory --- pearl_natural_saltwater
    pearl_akoya --- pearl_freshwater
    coral --- ivory
```

The graph clusters by gem family — the catalog's confusables encode the *natural* groupings a working gemmologist would expect.

## Lookup helpers

```python
from checkmsg import minerals
minerals.get("imperial jade")        # alias-aware → returns jadeite profile
minerals.by_species("garnet")        # all 7 garnet end-members
minerals.by_color("blue")            # 6 blue gems (and any others with blue in colors)
minerals.by_confusable("ruby")       # all minerals listing ruby as a confusable
minerals.resolve_confusables("ruby") # the ruby profile's confusables resolved to entries
```

## Profile-driven synthesis

Each profile can synthesise spectra for any of the six techniques via the catalog helpers:

```python
from checkmsg.minerals import get, synthesize_raman, synthesize_uvvis, synthesize_xrf, synthesize_libs, synthesize_epr
profile = get("ruby")
raman_spec = synthesize_raman(profile, noise=0.005)
uv_spec    = synthesize_uvvis(profile)
xrf_spec   = synthesize_xrf(profile)
libs_spec  = synthesize_libs(profile)
epr_spec   = synthesize_epr(profile, frequency_GHz=9.5)
```

This is the path used by `diagnose.diagnose_profile(profile)` and by every catalog-driven curriculum example in `examples/08..19`.

## Adding a mineral

To add a new entry to `CATALOG`, populate as much of the schema as the literature supports. Required fields: `name`, `species`, `chemical_formula`, plus at least one of (`raman_peaks_cm`, `uvvis_bands_nm`, `xrf_signature`, `libs_signature`, `epr_centers`) so the diagnose pipeline has at least one diagnostic signal to score against. Recommended: `mohs_hardness`, `density_g_cc`, `refractive_index`, `common_colors`, `confusables`, `diagnostic_features`, and at least one citation in `references`.

A minimal new entry looks like:

```python
"new_mineral": MineralProfile(
    name="new_mineral",
    species="<group name>",
    aliases=("commercial name", "historical name"),
    chemical_formula="ChemicalFormula",
    crystal_system="<crystal system>",
    mohs_hardness=(low, high),
    density_g_cc=(low, high),
    refractive_index=(low, high),
    common_colors=("colour1", "colour2"),
    raman_peaks_cm=_p((peak1_cm, rel_intensity), (peak2_cm, rel_intensity), ...),
    uvvis_bands_nm=(band1_nm, band2_nm),
    chromophores=("<key from refdata.chromophores>",),
    xrf_signature={"Element1": "major", "Element2": "trace"},
    libs_signature={"Element1": "major", "Element3": "trace"},
    epr_centers=("<key from refdata.epr_centers.CENTERS>",),
    icpms_diagnostic_isotopes=("Pb208", "Sr88"),
    confusables=("similar_mineral_1", "similar_mineral_2"),
    diagnostic_features=("characteristic 1", "characteristic 2"),
    references=("Author Year, *Journal* Vol:Page",),
),
```

After adding the entry, the test suite enforces three invariants automatically:

1. `tests/test_minerals_catalog.py::test_every_profile_has_at_least_one_diagnostic_signal` — at least one of the spectroscopic-signature fields must be non-empty.
2. `tests/test_minerals_catalog.py::test_every_profile_has_at_least_one_diagnostic_feature_string` — `diagnostic_features` must be non-empty.
3. `tests/test_minerals_catalog.py::test_every_confusable_resolves` — every name listed in `confusables` must itself exist in `CATALOG` (catch typos).

If the new mineral introduces a chromophore key that isn't in `refdata/chromophores.CHROMOPHORES`, add the chromophore there first with its band positions and tolerance. Same pattern for new EPR centres.

For technique-by-technique deep dives see [`techniques.md`](techniques.md). For the diagnostic pipeline that consumes these synthesised spectra, see [`diagnose.md`](diagnose.md).
