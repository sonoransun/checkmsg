"""HTTP API for the Check M.S.G. diagnosis pipeline (optional ``[service]`` extra)."""

from __future__ import annotations

from checkmsg.service.app import create_app

__all__ = ["create_app"]
