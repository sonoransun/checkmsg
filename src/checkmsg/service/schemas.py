"""Pydantic request/response models — the published JSON contract."""

from __future__ import annotations

import math
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from checkmsg.service.config import CONFIG

# The seven technique literals diagnose() understands (muon-xray is excluded:
# it is not wired into the diagnosis pipeline).
TechniqueLiteral = Literal[
    "raman", "xrf", "libs", "uvvis", "epr", "laicpms", "squid-mh", "squid-chi",
    "pl", "ftir", "mossbauer", "cl",
]
TierLiteral = Literal["novice", "practitioner", "expert"]


def _check_finite(values: list[float]) -> list[float]:
    if any(not math.isfinite(v) for v in values):
        raise ValueError("axis/intensity must contain only finite numbers")
    return values


class SpectrumInput(BaseModel):
    """One spectrum supplied inline as numeric arrays."""

    technique: TechniqueLiteral
    axis: list[float] = Field(min_length=2, max_length=CONFIG.max_points)
    intensity: list[float] = Field(min_length=2, max_length=CONFIG.max_points)
    frequency_GHz: float | None = Field(default=None, description="EPR microwave frequency.")
    temperature_K: float | None = Field(default=None, description="SQUID M(H) sample temperature.")
    bias_mT: float | None = Field(default=None, description="SQUID chi DC bias field.")
    isotope_keys: list[str] | None = Field(default=None, description="LA-ICP-MS measured isotopes.")

    _finite_axis = field_validator("axis")(_check_finite)
    _finite_intensity = field_validator("intensity")(_check_finite)

    @model_validator(mode="after")
    def _equal_length(self) -> SpectrumInput:
        if len(self.axis) != len(self.intensity):
            raise ValueError("axis and intensity must have equal length")
        if self.technique == "epr" and self.frequency_GHz is None:
            raise ValueError("technique 'epr' requires frequency_GHz")
        if self.technique == "laicpms" and not self.isotope_keys:
            raise ValueError("technique 'laicpms' requires isotope_keys")
        return self


class DiagnoseRequest(BaseModel):
    tier: TierLiteral = "novice"
    spectra: list[SpectrumInput] = Field(min_length=1, max_length=CONFIG.max_spectra)


class AnalyzeRequest(BaseModel):
    """Single-technique request body; the technique comes from the URL path."""

    axis: list[float] = Field(min_length=2, max_length=CONFIG.max_points)
    intensity: list[float] = Field(min_length=2, max_length=CONFIG.max_points)
    tier: TierLiteral = "novice"
    frequency_GHz: float | None = None
    temperature_K: float | None = None
    bias_mT: float | None = None
    isotope_keys: list[str] | None = None

    _finite_axis = field_validator("axis")(_check_finite)
    _finite_intensity = field_validator("intensity")(_check_finite)

    @model_validator(mode="after")
    def _equal_length(self) -> AnalyzeRequest:
        if len(self.axis) != len(self.intensity):
            raise ValueError("axis and intensity must have equal length")
        return self


# --- response models (permissive: extra fields allowed for tier variation) ---


class VerdictBlock(BaseModel):
    model_config = {"extra": "allow"}
    name: str


class ConfidenceBlock(BaseModel):
    value: float
    band: str
    note: str


class DiagnoseResponse(BaseModel):
    """Documents the response shape; extra (tier-specific) keys are allowed."""

    model_config = {"extra": "allow"}

    tier: str
    verdict: VerdictBlock | None
    confidence: ConfidenceBlock
    summary: str
    follow_ups: list[str]
    conflicts: list[dict[str, Any]]
    glossary: list[dict[str, Any]]
    disclaimer: str
