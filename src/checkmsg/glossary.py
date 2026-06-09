"""Plain-language glossary for eliminating shorthand in user-facing output.

The toolkit's reports are dense with domain shorthand — technique acronyms
(``EPR``), ion notation (``Cr3+``), physics units (``cm-1``, ``emu/g``), and
magnetic-ordering terms (``ferrimagnetic``). This module is the single place
that maps each piece of shorthand to a plain expansion + description so the CLI,
the tiered report renderer, and the web service can all explain jargon
consistently.

The glossary is assembled from a hand-authored CORE set (techniques, units,
ions, concepts, magnetic orderings) plus harvesters that pull descriptions
already present in the reference-data layer (chromophore notes, EPR centre
names, mineral species/aliases). CORE entries win on key collisions.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from checkmsg.minerals import CATALOG
from checkmsg.refdata.chromophores import CHROMOPHORES
from checkmsg.refdata.epr_centers import CENTERS

# ---------------------------------------------------------------------------
# Term model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GlossaryTerm:
    canonical: str          # the shorthand as it appears in output (e.g. "EPR")
    expansion: str          # short expansion (e.g. "electron paramagnetic resonance")
    description: str         # one-sentence plain-English explanation
    category: str           # technique | unit | ion | concept | score | magnetic | chromophore | center | mineral
    symbol: str = ""        # optional chemical formula / symbol
    citation: str = ""
    aliases: tuple[str, ...] = ()

    def inline(self) -> str:
        """First-use form: ``EPR (electron paramagnetic resonance)``."""
        return f"{self.canonical} ({self.expansion})" if self.expansion else self.canonical


# ---------------------------------------------------------------------------
# Hand-authored core terms
# ---------------------------------------------------------------------------

CORE: tuple[GlossaryTerm, ...] = (
    # --- techniques ---
    GlossaryTerm("Raman", "Raman spectroscopy", "Measures vibrational modes via inelastic light "
                 "scattering; a structural fingerprint of the crystal.", "technique"),
    GlossaryTerm("XRF", "X-ray fluorescence", "Identifies which chemical elements are present from "
                 "their characteristic X-ray emission.", "technique", aliases=("xrf",)),
    GlossaryTerm("LIBS", "laser-induced breakdown spectroscopy", "Reads element emission lines from a "
                 "laser-generated plasma; good for light and trace elements.", "technique",
                 aliases=("libs",)),
    GlossaryTerm("UV-VIS", "ultraviolet-visible spectroscopy", "Measures light absorption that gives a "
                 "gem its colour.", "technique", aliases=("uvvis", "UV-Vis")),
    GlossaryTerm("EPR", "electron paramagnetic resonance", "Detects unpaired electrons on defects or "
                 "transition-metal ions.", "technique", aliases=("epr", "ESR")),
    GlossaryTerm("LA-ICP-MS", "laser-ablation inductively-coupled-plasma mass spectrometry",
                 "Measures trace-element concentrations and isotope ratios.", "technique",
                 aliases=("laicpms", "LA-ICP-MS", "ICP-MS")),
    GlossaryTerm("SQUID", "superconducting quantum interference device", "An ultra-sensitive "
                 "magnetometer measuring a sample's magnetism.", "technique",
                 aliases=("squid", "squid-mh", "squid-chi")),
    GlossaryTerm("PL", "photoluminescence spectroscopy", "Measures sharp defect-centre emission "
                 "lines (e.g. NV/Si-V in synthetic diamond).", "technique", aliases=("pl",)),
    GlossaryTerm("FTIR", "Fourier-transform infrared spectroscopy", "Identifies chemical bonds "
                 "and defects (diamond type, water, polymer) via infrared absorption.", "technique",
                 aliases=("ftir",)),
    GlossaryTerm("Mossbauer", "Mössbauer spectroscopy", "Reads iron oxidation state (Fe²⁺ vs "
                 "Fe³⁺) and site geometry via recoilless γ-ray resonance.", "technique",
                 aliases=("mossbauer", "Mössbauer")),
    GlossaryTerm("CL", "cathodoluminescence", "Electron-beam-excited luminescence revealing "
                 "activator centres and growth zoning.", "technique", aliases=("cl",)),
    # --- units ---
    GlossaryTerm("cm-1", "wavenumbers", "Energy unit for Raman peaks (reciprocal centimetres).", "unit",
                 aliases=("cm⁻¹",)),
    GlossaryTerm("keV", "kilo-electronvolts", "Energy unit for X-ray emission lines.", "unit"),
    GlossaryTerm("nm", "nanometres", "Wavelength unit for light (UV-VIS / LIBS).", "unit"),
    GlossaryTerm("mT", "millitesla", "Unit of magnetic field strength.", "unit"),
    GlossaryTerm("emu/g", "electromagnetic units per gram", "Unit of magnetic moment per unit mass.",
                 "unit", aliases=("emu·g⁻¹",)),
    GlossaryTerm("ppm", "parts per million", "Trace-concentration unit (1 ppm = 0.0001%).", "unit"),
    GlossaryTerm("Ma", "millions of years (mega-annum)", "Geological age unit.", "unit"),
    # --- ions ---
    GlossaryTerm("Cr3+", "chromium(III) ion", "Transition-metal ion; the red colour of ruby and the "
                 "green of emerald.", "ion", symbol="Cr³⁺"),
    GlossaryTerm("Fe2+", "iron(II) ion", "Divalent iron; colours many silicates green/blue.", "ion",
                 symbol="Fe²⁺"),
    GlossaryTerm("Fe3+", "iron(III) ion", "Trivalent iron; yellow/brown colouration.", "ion",
                 symbol="Fe³⁺"),
    GlossaryTerm("Ti4+", "titanium(IV) ion", "Pairs with Fe²⁺ to give blue sapphire its colour.", "ion",
                 symbol="Ti⁴⁺"),
    GlossaryTerm("V3+", "vanadium(III) ion", "Green/colour-change chromophore (tsavorite, some "
                 "emerald).", "ion", symbol="V³⁺"),
    GlossaryTerm("Mn2+", "manganese(II) ion", "Pink/orange chromophore; a pearl/biomineral marker.",
                 "ion", symbol="Mn²⁺"),
    GlossaryTerm("Co2+", "cobalt(II) ion", "Intense blue chromophore in spinel and glass.", "ion",
                 symbol="Co²⁺"),
    # --- concepts ---
    GlossaryTerm("IVCT", "intervalence charge transfer", "Electron hopping between two metal ions "
                 "(e.g. Fe²⁺→Ti⁴⁺) that produces strong colour.", "concept"),
    GlossaryTerm("d-d", "d-orbital electronic transition", "Crystal-field transition of a "
                 "transition-metal ion that produces colour.", "concept"),
    GlossaryTerm("FWHM", "full width at half maximum", "A peak's width; broad peaks suggest "
                 "disorder/glassiness.", "concept"),
    GlossaryTerm("R-line", "ruby R-line luminescence", "Sharp red emission near 694 nm characteristic "
                 "of Cr³⁺ in corundum.", "concept"),
    GlossaryTerm("chromophore", "colour-causing centre", "The ion or defect responsible for a gem's "
                 "colour.", "concept"),
    GlossaryTerm("coercivity", "coercive field", "The reverse magnetic field needed to demagnetise a "
                 "sample.", "concept"),
    GlossaryTerm("saturation moment", "saturation magnetisation", "The maximum magnetic moment a "
                 "sample reaches in a strong field.", "concept"),
    GlossaryTerm("hysteresis", "magnetic hysteresis", "The lag between magnetisation and applied field "
                 "that traces a loop.", "concept"),
    GlossaryTerm("Curie temperature", "Curie temperature (Tc)", "Temperature above which a "
                 "ferro/ferrimagnet loses its ordering.", "concept", aliases=("Tc",)),
    GlossaryTerm("Neel temperature", "Néel temperature (TN)", "Temperature above which an "
                 "antiferromagnet loses its ordering.", "concept", aliases=("TN", "Néel temperature")),
    GlossaryTerm("Weiss constant", "Curie-Weiss intercept", "Sign/size of magnetic coupling inferred "
                 "from susceptibility vs temperature.", "concept"),
    GlossaryTerm("Morin transition", "Morin spin-flop transition", "A spin reorientation in hematite "
                 "near 263 K.", "concept"),
    GlossaryTerm("REE", "rare-earth elements", "The lanthanide series, diagnostic of a mineral's "
                 "origin.", "concept"),
    GlossaryTerm("chondrite-normalized", "chondrite-normalised", "Element abundances divided by "
                 "primitive-meteorite values to reveal patterns.", "concept"),
    GlossaryTerm("U-Pb age", "uranium-lead age", "A radiometric age from accumulated lead from "
                 "uranium decay.", "concept"),
    GlossaryTerm("isomer shift", "Mössbauer isomer shift", "Centroid of a Mössbauer doublet; "
                 "tracks iron oxidation state (Fe²⁺ high, Fe³⁺ low).", "concept"),
    GlossaryTerm("quadrupole splitting", "Mössbauer quadrupole splitting", "Separation of a "
                 "Mössbauer doublet; reflects site distortion.", "concept"),
    GlossaryTerm("zero-phonon line", "zero-phonon line (ZPL)", "The sharp purely-electronic "
                 "emission line of a defect centre.", "concept"),
    GlossaryTerm("diamond type", "diamond IR type", "Classification (Ia/Ib/IIa/IIb) by nitrogen/"
                 "boron content from FTIR — a growth-history fingerprint.", "concept"),
    # --- scoring metrics ---
    GlossaryTerm("cosine", "cosine similarity", "How closely two spectra align in shape (1.0 = "
                 "identical).", "score"),
    GlossaryTerm("peak_score", "peak-list match score", "Fraction of catalog peaks matched (0-1).",
                 "score"),
    GlossaryTerm("g_score", "g-factor match score", "How well observed EPR g-factors match a reference "
                 "centre (0-1).", "score"),
    GlossaryTerm("separation ratio", "confidence separation ratio", "How cleanly the top candidate "
                 "out-scores the runner-up; NOT a probability of being correct.", "score"),
    # --- magnetic ordering (expansion blank: the word itself is plain; the
    #     description carries the meaning in the glossary block) ---
    GlossaryTerm("diamagnetic", "", "Weakly repelled by a magnetic field; no permanent "
                 "moment (e.g. pure diamond, quartz).", "magnetic"),
    GlossaryTerm("paramagnetic", "", "Weakly attracted; isolated magnetic ions with no "
                 "long-range order.", "magnetic"),
    GlossaryTerm("ferromagnetic", "", "Spins align in parallel, giving a strong permanent "
                 "moment.", "magnetic"),
    GlossaryTerm("ferrimagnetic", "", "Opposed spins of unequal size leave a net moment "
                 "(e.g. magnetite).", "magnetic"),
    GlossaryTerm("antiferromagnetic", "", "Opposed spins cancel, leaving ~no net "
                 "moment.", "magnetic"),
    GlossaryTerm("canted-afm", "canted antiferromagnet", "Antiferromagnet whose slightly tilted spins "
                 "give a weak moment (e.g. hematite).", "magnetic", aliases=("canted-AFM",)),
)


# ---------------------------------------------------------------------------
# Harvesters (descriptions already present in the reference-data layer)
# ---------------------------------------------------------------------------


def _harvest_chromophores() -> list[GlossaryTerm]:
    out = []
    for ch in CHROMOPHORES:
        out.append(GlossaryTerm(
            canonical=ch.name, expansion="", description=ch.notes, category="chromophore"))
    return out


def _harvest_epr_centers() -> list[GlossaryTerm]:
    out = []
    for key, sys in CENTERS.items():
        out.append(GlossaryTerm(
            canonical=sys.name, expansion="",
            description=f"EPR paramagnetic centre hosted in {sys.host}.",
            category="center", aliases=(key,)))
    return out


def _harvest_minerals() -> list[GlossaryTerm]:
    out = []
    for name, p in CATALOG.items():
        colours = ", ".join(p.common_colors) if p.common_colors else "various colours"
        desc = f"A {colours} gem; species: {p.species or 'n/a'}."
        out.append(GlossaryTerm(
            canonical=name, expansion="", description=desc, category="mineral",
            symbol=p.chemical_formula, aliases=tuple(p.aliases)))
    return out


# ---------------------------------------------------------------------------
# Assembly + indexes
# ---------------------------------------------------------------------------


def _build() -> tuple[dict[str, GlossaryTerm], dict[str, str]]:
    terms: dict[str, GlossaryTerm] = {}
    # Harvested first; CORE overrides on collision.
    for gt in (*_harvest_minerals(), *_harvest_epr_centers(), *_harvest_chromophores(), *CORE):
        terms[gt.canonical] = gt
    alias_index: dict[str, str] = {}
    for gt in terms.values():
        for alias in gt.aliases:
            alias_index.setdefault(alias.lower(), gt.canonical)
    return terms, alias_index


_GLOSSARY, _ALIAS_INDEX = _build()


# ---------------------------------------------------------------------------
# Lookup API
# ---------------------------------------------------------------------------


def define(term: str) -> GlossaryTerm | None:
    """Resolve a term (canonical, alias, or case-insensitive) to its entry."""
    if term in _GLOSSARY:
        return _GLOSSARY[term]
    canon = _ALIAS_INDEX.get(term.lower())
    if canon:
        return _GLOSSARY[canon]
    for key, gt in _GLOSSARY.items():  # case-insensitive canonical fallback
        if key.lower() == term.lower():
            return gt
    return None


def has(term: str) -> bool:
    return define(term) is not None


def expand(term: str) -> str:
    """Return the inline first-use form of a term, or the term unchanged."""
    gt = define(term)
    return gt.inline() if gt else term


def all_terms() -> list[GlossaryTerm]:
    return sorted(_GLOSSARY.values(), key=lambda g: g.canonical.lower())


# ---------------------------------------------------------------------------
# Text scanning + first-use expansion
# ---------------------------------------------------------------------------

# Terms eligible for substring/word matching in free text, longest first so a
# longer term (e.g. "Curie temperature") is matched before a substring of it.
_MATCHABLE = sorted(
    (gt for gt in _GLOSSARY.values() if gt.category != "mineral"),
    key=lambda g: len(g.canonical), reverse=True,
)


def _occurs(text: str, token: str) -> bool:
    """Whether ``token`` occurs in ``text`` (boundary-aware for word tokens)."""
    if re.search(r"[^\w]", token) or any(c.isdigit() for c in token):
        return token in text  # symbol/number-bearing: exact substring, case-sensitive
    return re.search(rf"\b{re.escape(token)}\b", text, re.IGNORECASE) is not None


def find_terms(text: str) -> list[GlossaryTerm]:
    """Return the glossary terms whose shorthand appears in ``text``."""
    found: dict[str, GlossaryTerm] = {}
    for gt in _MATCHABLE:
        if gt.canonical in found:
            continue
        if _occurs(text, gt.canonical) or any(_occurs(text, a) for a in gt.aliases):
            found[gt.canonical] = gt
    return sorted(found.values(), key=lambda g: g.canonical.lower())


def expand_first_use(text: str, seen: set[str]) -> str:
    """Expand the first occurrence of each acronym-style term to its inline form.

    Only terms that carry a short ``expansion`` (technique acronyms, units,
    concepts) are inline-expanded; richer terms (chromophore/centre names) are
    surfaced through the appended glossary block instead. ``seen`` is mutated so
    a term is expanded at most once per render.
    """
    for gt in _MATCHABLE:
        if not gt.expansion or gt.canonical in seen:
            continue
        token = gt.canonical
        if re.search(r"[^\w]", token) or any(c.isdigit() for c in token):
            if token in text:
                text = text.replace(token, gt.inline(), 1)
                seen.add(token)
        else:
            pat = re.compile(rf"\b{re.escape(token)}\b")
            m = pat.search(text)
            if m:
                text = text[:m.start()] + gt.inline() + text[m.end():]
                seen.add(token)
    return text


def glossary_for(report) -> list[GlossaryTerm]:
    """Collect every glossary term referenced anywhere in a diagnostic report.

    Scans the verdict, evidence observations, and reasoning-trace text. Accepts
    any object exposing ``verdict``, ``evidence`` (with ``.observation``), and
    ``reasoning_trace`` (with ``.finding`` / ``.implication``) — typed loosely to
    avoid importing ``diagnose`` (which imports this module).
    """
    chunks: list[str] = []
    if getattr(report, "verdict", None):
        chunks.append(report.verdict)
    for ev in getattr(report, "evidence", []) or []:
        chunks.append(getattr(ev, "observation", ""))
    for st in getattr(report, "reasoning_trace", []) or []:
        chunks.append(getattr(st, "finding", ""))
        chunks.append(getattr(st, "implication", ""))
    text = "\n".join(chunks)
    return find_terms(text)
