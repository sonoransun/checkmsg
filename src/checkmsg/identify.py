"""Multi-technique identification fusion."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

from checkmsg import epr as epr_mod
from checkmsg import laicpms as laicpms_mod
from checkmsg import libs as libs_mod
from checkmsg import raman as raman_mod
from checkmsg import uvvis as uvvis_mod
from checkmsg import xrf as xrf_mod
from checkmsg.spectrum import Spectrum


@dataclass
class IdentificationResult:
    raman: raman_mod.RamanResult | None = None
    xrf: xrf_mod.XrfResult | None = None
    libs: libs_mod.LibsResult | None = None
    uvvis: uvvis_mod.UvVisResult | None = None
    epr: epr_mod.EprResult | None = None
    laicpms: laicpms_mod.IcpmsResult | None = None
    notes: list[str] = field(default_factory=list)

    def headline(self) -> str:
        parts: list[str] = []
        if self.raman and self.raman.best:
            best = self.raman.best
            parts.append(f"Raman top match: {best.mineral} (cosine {best.cosine:.2f}, peak {best.peak_score:.2f})")
        if self.xrf and self.xrf.elements:
            top_xrf = ", ".join(f"{e.element}" for e in self.xrf.elements[:5])
            parts.append(f"XRF elements: {top_xrf}")
        if self.libs and self.libs.elements:
            top_libs = ", ".join(sorted(self.libs.elements.keys()))
            parts.append(f"LIBS lines: {top_libs}")
        if self.uvvis and self.uvvis.chromophores():
            tags = ", ".join(c.name for c in self.uvvis.chromophores())
            parts.append(f"UV-VIS chromophores: {tags}")
        if self.epr and self.epr.best:
            best = self.epr.best
            parts.append(f"EPR top center: {best.name} (cosine {best.cosine:.2f})")
        if self.laicpms:
            parts.append(f"LA-ICP-MS: {self.laicpms.headline()}")
        return "; ".join(parts) if parts else "no signal"

    def report(self) -> str:
        lines: list[str] = ["=== Check M.S.G. identification report ==="]
        lines.append(self.headline())
        if self.raman:
            lines.append("\nRaman candidates:")
            for c in self.raman.candidates:
                lines.append(
                    f"  {c.mineral:<14} cosine={c.cosine:.3f}  peak_score={c.peak_score:.3f}  matched={c.matched_peaks}")
        if self.xrf:
            lines.append("\nXRF elements (relative):")
            quant = xrf_mod.relative_quant(self.xrf)
            for e in self.xrf.elements:
                lines.append(f"  {e.element:<3} (Z={e.Z}) lines={len(e.matched_lines)} rel={quant.get(e.element, 0):.3f}")
        if self.libs:
            lines.append("\nLIBS elements:")
            for el, sig in sorted(self.libs.elements.items()):
                lines.append(f"  {el:<3} matched_lines={len(sig.matched_lines)} height={sig.summed_height:.2f}")
        if self.uvvis:
            lines.append("\nUV-VIS bands:")
            for b in self.uvvis.bands:
                lines.append(f"  {b.position:6.1f} nm  height={b.height:.3f} fwhm={b.width:.1f}")
            lines.append("  Chromophores: " + ", ".join(c.name for c in self.uvvis.chromophores()) or "  Chromophores: none")
        if self.epr:
            lines.append("\nEPR centers:")
            for c in self.epr.candidates:
                lines.append(f"  {c.name:<32} cosine={c.cosine:.3f}  g_score={c.g_score:.3f}  combined={c.combined:.3f}")
            if self.epr.g_factors:
                lines.append("  Inferred g-factors: " + ", ".join(f"{g:.4f}" for g in self.epr.g_factors[:8]))
        if self.laicpms:
            lines.append("\nLA-ICP-MS:")
            top = sorted(self.laicpms.concentrations.values(), key=lambda c: c.ppm, reverse=True)[:10]
            for c in top:
                lines.append(f"  {c.element:<3} {c.ppm:>10.3f} ppm  (LOD {c.detection_limit_ppm:.3f})")
            if self.laicpms.isotope_ratios:
                ratio_summary = ", ".join(f"{k}={v:.4f}" for k, v in self.laicpms.isotope_ratios.items())
                lines.append(f"  Isotope ratios: {ratio_summary}")
            if self.laicpms.u_pb_age_Ma is not None:
                lines.append(f"  U-Pb age: {self.laicpms.u_pb_age_Ma:.0f} Ma")
            if self.laicpms.ree_pattern:
                lines.append("  REE pattern (chondrite-normalised, top 5):")
                for el in ("La", "Ce", "Sm", "Eu", "Yb"):
                    val = self.laicpms.ree_pattern.get(el)
                    if val is not None:
                        lines.append(f"    {el}: {val:.2f}")
        if self.notes:
            lines.append("\nNotes:")
            lines.extend(f"  - {n}" for n in self.notes)
        return "\n".join(lines)


def combined_report(spectra: Iterable[Spectrum]) -> IdentificationResult:
    result = IdentificationResult()
    for s in spectra:
        if s.technique == "raman":
            result.raman = raman_mod.analyze(s)
            if raman_mod.is_amorphous(s):
                result.notes.append("Raman: no sharp crystalline peaks — sample is amorphous (e.g. glass).")
        elif s.technique == "xrf":
            result.xrf = xrf_mod.identify_elements(s)
        elif s.technique == "libs":
            result.libs = libs_mod.identify(s)
        elif s.technique == "uvvis":
            result.uvvis = uvvis_mod.assign_bands(s)
        elif s.technique == "epr":
            freq = s.metadata.get("frequency_GHz")
            if freq is None:
                result.notes.append("EPR spectrum missing frequency_GHz metadata; skipped.")
                continue
            result.epr = epr_mod.analyze(s, frequency_GHz=freq)
            if result.epr.best:
                bn = result.epr.best.name
                if bn == "diamond_P1":
                    result.notes.append("EPR: P1 nitrogen triplet — likely natural Ib or HPHT diamond.")
                elif bn == "quartz_E1prime":
                    result.notes.append("EPR: E1' oxygen-vacancy center — radiation/smoky quartz.")
                elif bn == "calcite_Mn2plus":
                    result.notes.append("EPR: Mn2+ sextet — biomineral / calcite fingerprint.")
        elif s.technique == "laicpms":
            # Combined report only handles bulk integrated spectra here. Time-resolved
            # transient runs are accessed directly via laicpms.analyze().
            result.notes.append("LA-ICP-MS bulk spectrum — call laicpms.analyze on the IcpmsRun for full results.")
    return result


def add_laicpms_result(result: IdentificationResult, lr: laicpms_mod.IcpmsResult) -> IdentificationResult:
    """Attach an already-computed `IcpmsResult` to an `IdentificationResult`."""
    result.laicpms = lr
    result.notes.extend(lr.notes)
    return result
