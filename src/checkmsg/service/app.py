"""FastAPI application exposing the Check M.S.G. diagnosis pipeline.

Launch with ``checkmsg serve`` or ``uvicorn checkmsg.service.app:app``.
All endpoints are network-free (see ``config`` — offline mode is forced).
"""

from __future__ import annotations

import json
from dataclasses import asdict

from fastapi import FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from checkmsg import __version__, glossary, minerals, scoring
from checkmsg.diagnose import diagnose
from checkmsg.service.config import CONFIG
from checkmsg.service.schemas import AnalyzeRequest, DiagnoseRequest, DiagnoseResponse
from checkmsg.service.serialization import report_to_dict
from checkmsg.service.spectra import UNITS, spectrum_from_arrays, spectrum_from_upload

# Plain descriptions for the /techniques listing (logical techniques).
_TECHNIQUE_INFO = {
    "raman": "Raman", "xrf": "XRF", "libs": "LIBS", "uvvis": "UV-VIS",
    "epr": "EPR", "laicpms": "LA-ICP-MS", "squid": "SQUID",
    "pl": "PL", "ftir": "FTIR", "mossbauer": "Mössbauer", "cl": "CL",
}


def _profile_summary(p) -> dict:
    return {
        "name": p.name,
        "species": p.species,
        "aliases": list(p.aliases),
        "common_colors": list(p.common_colors),
        "confusables": list(p.confusables),
    }


def create_app() -> FastAPI:
    app = FastAPI(
        title="Check M.S.G. — mineral/gem diagnosis API",
        version=__version__,
        description=("Public-facing spectroscopic gem-diagnosis service. "
                     "Every diagnosis carries a pedagogical disclaimer; the confidence "
                     "band is a separation ratio, not a probability of correctness."),
    )
    app.add_middleware(
        CORSMiddleware,
        allow_origins=list(CONFIG.cors_origins),
        allow_methods=["GET", "POST"],
        allow_headers=["*"],
    )

    @app.exception_handler(KeyError)
    def _key_error(_request, exc):  # unknown catalog/glossary lookups
        return JSONResponse(status_code=404, content={"detail": f"not found: {exc}"})

    @app.exception_handler(ValueError)
    def _value_error(_request, exc):
        return JSONResponse(status_code=422, content={"detail": str(exc)})

    # ---- meta ----

    @app.get("/")
    def root():
        return {
            "service": "checkmsg",
            "version": __version__,
            "endpoints": ["/diagnose", "/diagnose/upload", "/analyze/{technique}",
                          "/catalog", "/catalog/{name}", "/glossary", "/glossary/{term}",
                          "/techniques", "/healthz", "/version", "/docs"],
            "disclaimer": scoring.SERVICE_DISCLAIMER,
        }

    @app.get("/healthz")
    def healthz():
        return {"status": "ok"}

    @app.get("/version")
    def version():
        return {"name": "checkmsg", "version": __version__}

    @app.get("/techniques")
    def techniques():
        out = []
        for tech, canon in _TECHNIQUE_INFO.items():
            gt = glossary.define(canon)
            cav = scoring.CAVEATS.get(tech)
            out.append({
                "technique": tech,
                "name": gt.inline() if gt else canon,
                "description": gt.description if gt else "",
                "accuracy_tier": cav.severity if cav else "unknown",
                "measurement_accurate": cav.measurement_accurate if cav else None,
                "caveat": cav.text if cav else "",
            })
        return {"count": len(out), "items": out}

    # ---- catalog ----

    @app.get("/catalog")
    def catalog(species: str | None = None, color: str | None = None,
                limit: int = Query(default=200, ge=1, le=500), offset: int = Query(default=0, ge=0)):
        if species:
            names = list(minerals.by_species(species).keys())
        elif color:
            names = list(minerals.by_color(color).keys())
        else:
            names = list(minerals.names())
        window = names[offset:offset + limit]
        items = [_profile_summary(minerals.get(n)) for n in window]
        return {"count": len(names), "items": items}

    @app.get("/catalog/{name}")
    def catalog_entry(name: str):
        p = minerals.get(name)  # KeyError -> 404
        return asdict(p)

    # ---- glossary ----

    @app.get("/glossary")
    def glossary_list(limit: int = Query(default=500, ge=1, le=1000), offset: int = Query(default=0, ge=0)):
        terms = glossary.all_terms()
        window = terms[offset:offset + limit]
        items = [{"term": t.canonical, "expansion": t.expansion,
                  "definition": t.description, "category": t.category} for t in window]
        return {"count": len(terms), "items": items}

    @app.get("/glossary/{term}")
    def glossary_define(term: str):
        gt = glossary.define(term)
        if gt is None:
            raise HTTPException(status_code=404, detail=f"unknown term: {term}")
        return {"term": gt.canonical, "expansion": gt.expansion, "expanded": gt.inline(),
                "definition": gt.description, "category": gt.category}

    # ---- diagnosis ----

    @app.post("/diagnose", response_model=DiagnoseResponse)
    def diagnose_endpoint(req: DiagnoseRequest):
        spectra = [
            spectrum_from_arrays(
                s.technique, s.axis, s.intensity,
                frequency_GHz=s.frequency_GHz, temperature_K=s.temperature_K,
                bias_mT=s.bias_mT, isotope_keys=s.isotope_keys,
            )
            for s in req.spectra
        ]
        report = diagnose(spectra)
        return report_to_dict(report, tier=req.tier)

    @app.post("/analyze/{technique}", response_model=DiagnoseResponse)
    def analyze_endpoint(technique: str, req: AnalyzeRequest):
        if technique not in UNITS:
            raise HTTPException(status_code=422, detail=f"unknown technique: {technique}")
        if technique == "epr" and req.frequency_GHz is None:
            raise HTTPException(status_code=422, detail="technique 'epr' requires frequency_GHz")
        if technique == "laicpms" and not req.isotope_keys:
            raise HTTPException(status_code=422, detail="technique 'laicpms' requires isotope_keys")
        spec = spectrum_from_arrays(
            technique, req.axis, req.intensity,
            frequency_GHz=req.frequency_GHz, temperature_K=req.temperature_K,
            bias_mT=req.bias_mT, isotope_keys=req.isotope_keys,
        )
        report = diagnose([spec])
        return report_to_dict(report, tier=req.tier)

    @app.post("/diagnose/upload", response_model=DiagnoseResponse)
    def diagnose_upload(
        files: list[UploadFile] = File(...),  # noqa: B008 (FastAPI dependency idiom)
        specs: str = Form(...),  # noqa: B008
        tier: str = Form("novice"),  # noqa: B008
    ):
        try:
            spec_meta = json.loads(specs)
        except json.JSONDecodeError as exc:
            raise HTTPException(status_code=422, detail=f"invalid specs JSON: {exc}") from exc
        if not isinstance(spec_meta, list) or len(spec_meta) != len(files):
            raise HTTPException(status_code=422,
                                detail="specs must be a JSON list parallel to the uploaded files")
        if len(files) > CONFIG.max_spectra:
            raise HTTPException(status_code=422, detail="too many spectra")
        spectra = []
        for upload, meta in zip(files, spec_meta, strict=True):
            technique = meta.get("technique")
            if technique not in UNITS:
                raise HTTPException(status_code=422, detail=f"unknown technique: {technique}")
            content = upload.file.read()
            if len(content) > CONFIG.max_upload_bytes:
                raise HTTPException(status_code=413, detail=f"{upload.filename}: file too large")
            spectra.append(spectrum_from_upload(
                technique, content,
                frequency_GHz=meta.get("frequency_GHz"), temperature_K=meta.get("temperature_K"),
                bias_mT=meta.get("bias_mT"), isotope_keys=meta.get("isotope_keys"),
            ))
        report = diagnose(spectra)
        return report_to_dict(report, tier=tier)

    return app


app = create_app()
