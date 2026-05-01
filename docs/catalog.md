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

For technique-by-technique deep dives see [`techniques.md`](techniques.md). For the diagnostic pipeline that consumes these synthesised spectra, see [`diagnose.md`](diagnose.md).
