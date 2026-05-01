"""Bundled EPR / ESR paramagnetic center library for gemological analysis.

Each center is a `SpinSystem` derived from the published spin-Hamiltonian
parameters in the cited primary source. Values are quoted in MHz (D, E,
hyperfine A) and dimensionless (g).

Sources:
  - Loubser & van Wyk 1978, *Rep. Prog. Phys.* 41:1201 (P1 / N0 substitutional nitrogen in diamond)
  - Solntsev & Yelisseyev 2000, *Mater. Sci. Forum* 354:7 (Ni-related HPHT diamond centres)
  - Weil 1984, *Phys. Chem. Minerals* 10:149 (E1' centre in α-quartz)
  - Mackey & Sander 1972, *J. Chem. Phys.* 56:5365 (Al-hole centre in smoky quartz)
  - Bernstein 1979, *Am. Mineralogist* 64:1056 (Mn²⁺ in calcite — pearl / biomineral fingerprint)
  - Manenkov & Prokhorov 1956, *Soviet Physics JETP* 1:611 (Cr³⁺ in corundum — the ruby maser system)
  - Kornienko & Prokhorov 1960, *Soviet Physics JETP* 11:1189 (Fe³⁺ in corundum)
"""

from __future__ import annotations

from checkmsg.epr import Hyperfine, SpinSystem

CENTERS: dict[str, SpinSystem] = {
    # --- Quantification & reference standards ---
    "DPPH": SpinSystem(
        name="DPPH (alpha,alpha-diphenyl-beta-picrylhydrazyl)",
        S=0.5, g=2.0036, linewidth_mT=0.15,
        host="standard",
    ),
    "free_electron": SpinSystem(
        name="free electron (g_e)",
        S=0.5, g=2.00231930, linewidth_mT=0.10,
        host="standard",
    ),

    # --- Diamond defects ---
    "diamond_P1": SpinSystem(
        name="P1 / N0 substitutional nitrogen (Ib diamond)",
        S=0.5, g=(2.0024, 2.0024, 2.0024),
        hyperfine=(Hyperfine(nucleus="14N", I=1.0,
                             A_iso_MHz=(114.03 + 81.32 + 81.32) / 3,
                             A_aniso_MHz=(81.32 - (114.03 + 81.32 + 81.32)/3,
                                          81.32 - (114.03 + 81.32 + 81.32)/3,
                                          114.03 - (114.03 + 81.32 + 81.32)/3)),),
        linewidth_mT=0.4, host="diamond",
    ),
    "diamond_Ni_HPHT": SpinSystem(
        name="Ni-related center (HPHT-grown diamond)",
        S=0.5, g=2.0319, linewidth_mT=0.6,
        host="diamond",
    ),

    # --- Quartz colour centers ---
    "quartz_E1prime": SpinSystem(
        name="E1' oxygen-vacancy center (smoky / radiation-induced quartz)",
        S=0.5, g=2.0006, linewidth_mT=0.15,
        host="quartz",
    ),
    "quartz_Al_hole": SpinSystem(
        name="[AlO4]0 hole center (smoky quartz)",
        S=0.5, g=(2.060, 2.005, 2.002), linewidth_mT=0.5,
        host="quartz",
    ),

    # --- Biominerals ---
    "calcite_Mn2plus": SpinSystem(
        name="Mn2+ in calcite / aragonite (pearl fingerprint)",
        S=2.5, g=2.0010, D_MHz=25.0, E_MHz=0.0,
        hyperfine=(Hyperfine(nucleus="55Mn", I=2.5, A_iso_MHz=245.0),),
        linewidth_mT=0.4, host="calcite",
    ),

    # --- Corundum (ruby / sapphire) ---
    "corundum_Cr3plus": SpinSystem(
        name="Cr3+ in corundum (ruby color-center / maser system)",
        S=1.5, g=(1.984, 1.984, 1.984),
        D_MHz=5715.0, E_MHz=0.0,
        linewidth_mT=0.6, host="corundum",
    ),
    "corundum_Fe3plus": SpinSystem(
        name="Fe3+ in corundum (blue-sapphire chromophore)",
        S=2.5, g=(2.003, 2.003, 2.003),
        D_MHz=12050.0, E_MHz=0.0,
        linewidth_mT=1.0, host="corundum",
    ),
}


# ---- Convenience accessors ----

def by_host(host: str) -> dict[str, SpinSystem]:
    """Filter the center registry by host material (case-insensitive)."""
    h = host.lower()
    return {k: v for k, v in CENTERS.items() if v.host.lower() == h}


def quant_standards() -> dict[str, SpinSystem]:
    """Reference centers used for spin-quantification (DPPH, free electron)."""
    return {k: v for k, v in CENTERS.items() if v.host == "standard"}


def diagnostic_for(scenario: str) -> tuple[str, ...]:
    """Names of centers that are diagnostic for a given gemological scenario."""
    table: dict[str, tuple[str, ...]] = {
        "diamond-natural-vs-synthetic": ("diamond_P1", "diamond_Ni_HPHT"),
        "smoky-vs-natural-quartz": ("quartz_E1prime", "quartz_Al_hole"),
        "pearl-fingerprint": ("calcite_Mn2plus",),
        "ruby-vs-sapphire": ("corundum_Cr3plus", "corundum_Fe3plus"),
    }
    return table.get(scenario, ())
