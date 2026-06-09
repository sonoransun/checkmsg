"""Service configuration and hardening limits (environment-driven).

Importing this module forces ``CHECKMSG_OFFLINE=1`` if it is unset: the public
service must never reach out to RRUFF. (``diagnose()`` with ``candidates=None``
does not use the network anyway — it calls ``raman.detect`` directly, not
``raman.analyze`` — but this is belt-and-braces.)
"""

from __future__ import annotations

import os
from dataclasses import dataclass

# Belt-and-braces: never let the service fetch reference spectra over the network.
os.environ.setdefault("CHECKMSG_OFFLINE", "1")


def _int_env(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default


@dataclass(frozen=True)
class ServiceConfig:
    max_spectra: int = 12               # per /diagnose request
    max_points: int = 200_000           # per axis/intensity array
    max_upload_bytes: int = 5 * 1024 * 1024   # per uploaded CSV file
    cors_origins: tuple[str, ...] = ()  # empty = deny cross-origin by default

    @classmethod
    def from_env(cls) -> ServiceConfig:
        origins = os.environ.get("CHECKMSG_SERVICE_CORS", "")
        return cls(
            max_spectra=_int_env("CHECKMSG_SERVICE_MAX_SPECTRA", 12),
            max_points=_int_env("CHECKMSG_SERVICE_MAX_POINTS", 200_000),
            max_upload_bytes=_int_env("CHECKMSG_SERVICE_MAX_UPLOAD_BYTES", 5 * 1024 * 1024),
            cors_origins=tuple(o.strip() for o in origins.split(",") if o.strip()),
        )


CONFIG = ServiceConfig.from_env()
