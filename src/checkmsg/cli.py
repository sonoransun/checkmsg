from __future__ import annotations

import argparse
import sys
from pathlib import Path

from checkmsg import diagnose as diagnose_mod
from checkmsg import epr as epr_mod
from checkmsg import glossary as glossary_mod
from checkmsg import io as io_mod
from checkmsg import laicpms as laicpms_mod
from checkmsg import squid as squid_mod
from checkmsg.identify import combined_report
from checkmsg.libs import identify as libs_identify
from checkmsg.raman import analyze as raman_analyze
from checkmsg.spectrum import Technique
from checkmsg.uvvis import assign_bands as uvvis_assign
from checkmsg.xrf import identify_elements as xrf_identify

UNITS = {
    "raman": "cm-1", "xrf": "keV", "libs": "nm", "uvvis": "nm",
    "epr": "mT", "laicpms": "m/z", "squid-mh": "mT", "squid-chi": "K",
    "pl": "nm", "ftir": "cm-1", "mossbauer": "mm/s", "cl": "nm",
}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser("checkmsg", description="Spectroscopic mineral/gem analysis.")
    sub = p.add_subparsers(dest="cmd", required=True)

    pa = sub.add_parser("analyze", help="Run a single technique on one spectrum file.")
    pa.add_argument(
        "technique",
        choices=["raman", "xrf", "libs", "uvvis", "epr", "laicpms", "squid-mh", "squid-chi",
                 "pl", "ftir", "mossbauer", "cl"],
    )
    pa.add_argument("path", type=Path)
    pa.add_argument("--frequency", type=float, default=None,
                    help="Microwave frequency in GHz (required for technique=epr).")
    pa.add_argument("--calibration", type=Path, default=None,
                    help="Calibration spectrum file (LA-ICP-MS only; NIST glass).")
    pa.add_argument("--internal-standard", type=str, default=None,
                    help="Internal standard for LA-ICP-MS, e.g. 'Al:529000' (element:ppm).")
    pa.add_argument("--temperature", type=float, default=295.0,
                    help="Sample temperature (K) for SQUID M(H) measurements.")

    pi = sub.add_parser("identify", help="Run multi-technique fusion on several files.")
    pi.add_argument("files", nargs="+", type=str,
                    help="Each entry is 'raman:path.csv' (or 'epr:path.csv:9.5' for EPR).")

    pd = sub.add_parser("diagnose", help="Run unified diagnostic pipeline on a set of spectra.")
    pd.add_argument("files", nargs="+", type=str,
                    help="Same format as 'identify': 'raman:path.csv', 'epr:path.csv:9.5', etc.")
    pd.add_argument("--tier", choices=["novice", "practitioner", "expert"], default="expert",
                    help="Report sophistication: novice (plain), practitioner (explained), "
                         "expert (full trace + caveats; default).")
    pd.add_argument("--calibrated", nargs="?", const="platt", default=None, choices=["platt", "gnb"],
                    help="Attach a calibrated P(correct) from the bundled Platt or GNB calibrator.")

    psim = sub.add_parser("similar", help="Find the nearest catalog spectra to a measurement.")
    psim.add_argument("entry", help="'<technique>:<path>', e.g. raman:spectrum.csv")
    psim.add_argument("--k", type=int, default=5, help="Number of nearest matches to show.")

    pg = sub.add_parser("glossary", help="Look up the plain-language meaning of a shorthand term.")
    pg.add_argument("term", nargs="?", default=None,
                    help="Term to define (e.g. EPR, Cr3+, ferrimagnetic). Omit to list all.")

    ps = sub.add_parser("serve", help="Launch the HTTP API (requires the 'service' extra).")
    ps.add_argument("--host", default="127.0.0.1")
    ps.add_argument("--port", type=int, default=8000)
    ps.add_argument("--reload", action="store_true", help="Auto-reload on code changes (dev).")

    args = p.parse_args(argv)

    if args.cmd == "analyze":
        spec = io_mod.read_csv(args.path, technique=args.technique, units=UNITS[args.technique])
        if args.technique == "raman":
            r = raman_analyze(spec)
            print(f"Raman peaks: {len(r.peaks)}")
            for c in r.candidates:
                print(f"  {c.mineral:<12} cos={c.cosine:.3f} peak_score={c.peak_score:.3f}")
        elif args.technique == "xrf":
            r = xrf_identify(spec)
            print(f"XRF peaks: {len(r.peaks)}")
            for e in r.elements:
                print(f"  {e.element} (Z={e.Z}) confidence={e.confidence:.2f}")
        elif args.technique == "libs":
            r = libs_identify(spec)
            for el, sig in sorted(r.elements.items()):
                print(f"  {el} matched={len(sig.matched_lines)} height={sig.summed_height:.2f}")
        elif args.technique == "uvvis":
            r = uvvis_assign(spec)
            for c in r.chromophores():
                print(f"  {c.name}")
        elif args.technique == "epr":
            if args.frequency is None:
                print("error: --frequency is required for technique=epr", file=sys.stderr)
                return 2
            spec.metadata.setdefault("frequency_GHz", args.frequency)
            r = epr_mod.analyze(spec, frequency_GHz=args.frequency)
            print(f"EPR g-factors: {[round(g, 4) for g in r.g_factors]}")
            for c in r.candidates:
                print(f"  {c.name:<32} cos={c.cosine:.3f} g_score={c.g_score:.3f}")
        elif args.technique in ("squid-mh", "squid-chi"):
            spec.metadata.setdefault("squid_mode",
                                     "dc-mh" if args.technique == "squid-mh" else "rf-chi-T")
            spec.metadata.setdefault("temperature_K", float(args.temperature))
            meas = squid_mod.from_spectrum(spec)
            r = squid_mod.analyze(meas)
            print(f"SQUID extracted: {r.headline()}")
            for c in r.candidates:
                print(f"  {c.name:<28} combined={c.combined:.2f} "
                      f"ord_match={c.ordering_match} Tc_resid={c.tc_residual_K:.0f} K")
        elif args.technique == "laicpms":
            sample_run = laicpms_mod.run_from_spectrum(spec)
            cal_run = None
            is_pair = None
            if args.calibration:
                cal_spec = io_mod.read_csv(args.calibration, technique="laicpms", units="m/z")
                cal_run = laicpms_mod.run_from_spectrum(cal_spec)
            if args.internal_standard:
                el, _, ppm_str = args.internal_standard.partition(":")
                if not el or not ppm_str:
                    print("error: --internal-standard must be 'Element:ppm'", file=sys.stderr)
                    return 2
                is_pair = (el, float(ppm_str))
            r = laicpms_mod.analyze(sample_run, calibration=cal_run, internal_standard=is_pair)
            print(r.headline())
            for note in r.notes:
                print(f"  - {note}")
            for c in sorted(r.concentrations.values(), key=lambda c: c.ppm, reverse=True)[:10]:
                print(f"  {c.element:<3} {c.ppm:>10.3f} ppm  (LOD {c.detection_limit_ppm:.3f})")
        elif args.technique == "pl":
            from checkmsg import pl as pl_mod
            r = pl_mod.analyze(spec)
            print(r.headline())
            for p, b in r.assignments:
                print(f"  {b.name:<24} at {p:.1f} nm")
        elif args.technique == "ftir":
            from checkmsg import ftir as ftir_mod
            print(ftir_mod.analyze(spec).headline())
        elif args.technique == "cl":
            from checkmsg import cl as cl_mod
            print(cl_mod.analyze(spec).headline())
        elif args.technique == "mossbauer":
            from checkmsg import mossbauer as moss_mod
            print(moss_mod.analyze(spec).headline())
        return 0

    if args.cmd == "diagnose":
        spectra = []
        for entry in args.files:
            parts = entry.split(":")
            if len(parts) < 2:
                print(f"error: '{entry}' must be in form <technique>:<path>", file=sys.stderr)
                return 2
            tech = parts[0]
            fname = parts[1]
            extra: dict = {}
            if tech == "epr" and len(parts) >= 3:
                extra["frequency_GHz"] = float(parts[2])
            if tech == "squid-mh":
                extra["squid_mode"] = "dc-mh"
                if len(parts) >= 3:
                    extra["temperature_K"] = float(parts[2])
            elif tech == "squid-chi":
                extra["squid_mode"] = "rf-chi-T"
                if len(parts) >= 3:
                    extra["applied_field_mT"] = float(parts[2])
            tech_t: Technique = tech  # type: ignore[assignment]
            s = io_mod.read_csv(fname, technique=tech_t, units=UNITS[tech])
            s.metadata.update(extra)
            spectra.append(s)
        report = diagnose_mod.diagnose(spectra, calibrated=args.calibrated or False)
        print(report.render(args.tier))
        return 0

    if args.cmd == "similar":
        from checkmsg import similarity
        parts = args.entry.split(":")
        if len(parts) < 2:
            print(f"error: '{args.entry}' must be '<technique>:<path>'", file=sys.stderr)
            return 2
        tech, fname = parts[0], parts[1]
        s = io_mod.read_csv(fname, technique=tech, units=UNITS.get(tech, ""))  # type: ignore[arg-type]
        hits = similarity.nearest(s, k=args.k)
        if not hits:
            print(f"no catalog references for technique {tech!r}")
            return 0
        print(f"Nearest {tech} catalog matches:")
        for hit in hits:
            print(f"  {hit.name:<24} similarity={hit.score:.3f}")
        return 0

    if args.cmd == "glossary":
        if args.term is None:
            for t in glossary_mod.all_terms():
                print(f"{t.canonical}: {t.description}")
            return 0
        gt = glossary_mod.define(args.term)
        if gt is None:
            print(f"error: unknown term {args.term!r}", file=sys.stderr)
            return 1
        print(gt.inline())
        print(f"  {gt.description}")
        if gt.citation:
            print(f"  source: {gt.citation}")
        return 0

    if args.cmd == "serve":
        try:
            import uvicorn
        except ImportError:
            print("error: the API server needs the 'service' extra: "
                  "pip install -e '.[service]'", file=sys.stderr)
            return 2
        uvicorn.run("checkmsg.service.app:app", host=args.host, port=args.port,
                    reload=args.reload)
        return 0

    if args.cmd == "identify":
        spectra = []
        for entry in args.files:
            parts = entry.split(":")
            if len(parts) < 2:
                print(f"error: '{entry}' must be in form <technique>:<path>", file=sys.stderr)
                return 2
            tech = parts[0]
            fname = parts[1]
            extra: dict = {}
            if tech == "epr":
                if len(parts) < 3:
                    print("error: EPR entry needs frequency: 'epr:<path>:<freq_GHz>'", file=sys.stderr)
                    return 2
                extra["frequency_GHz"] = float(parts[2])
            if tech == "squid-mh":
                extra["squid_mode"] = "dc-mh"
                if len(parts) >= 3:
                    extra["temperature_K"] = float(parts[2])
            elif tech == "squid-chi":
                extra["squid_mode"] = "rf-chi-T"
                if len(parts) >= 3:
                    extra["applied_field_mT"] = float(parts[2])
            tech_t: Technique = tech  # type: ignore[assignment]
            s = io_mod.read_csv(fname, technique=tech_t, units=UNITS[tech])
            s.metadata.update(extra)
            spectra.append(s)
        result = combined_report(spectra)
        print(result.report())
        return 0

    return 2


if __name__ == "__main__":
    sys.exit(main())
