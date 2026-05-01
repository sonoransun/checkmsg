from __future__ import annotations

import os
from pathlib import Path

from platformdirs import user_cache_dir


def cache_root() -> Path:
    override = os.environ.get("CHECKMSG_CACHE")
    base = Path(override) if override else Path(user_cache_dir("checkmsg", "checkmsg"))
    base.mkdir(parents=True, exist_ok=True)
    return base


def is_offline() -> bool:
    return os.environ.get("CHECKMSG_OFFLINE", "").strip() in {"1", "true", "yes"}
