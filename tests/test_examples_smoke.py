"""Smoke-test each example script via subprocess. Verifies end-to-end behaviour."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"

SCRIPTS = [
    "01_diamond_vs_moissanite_vs_cz.py",
    "02_natural_vs_synthetic_ruby.py",
    "03_emerald_vs_green_glass.py",
    "04_sapphire_origin.py",
    "05_multi_laser_temperature.py",
    "06_epr_unpaired_electrons.py",
    "07_laicpms_complex_cases.py",
    "08_diamond_simulant_carousel.py",
    "09_blue_stone_disambiguation.py",
    "10_garnet_species_suite.py",
    "11_jade_family_confusion.py",
    "12_tourmaline_species.py",
    "13_red_gems_carousel.py",
    "14_black_opaque_stones.py",
    "15_quartz_treatment_history.py",
    "16_chrysoberyl_trio.py",
    "17_biomineral_disambiguation.py",
    "18_treatment_detection.py",
    "19_unknown_stone_capstone.py",
]


@pytest.mark.parametrize("script", SCRIPTS)
def test_example_smoke(script, tmp_path):
    env = os.environ.copy()
    env["CHECKMSG_CACHE"] = str(tmp_path / "cache")
    proc = subprocess.run(
        [sys.executable, str(EXAMPLES / script), "--smoke"],
        capture_output=True, text=True, env=env, timeout=120,
    )
    assert proc.returncode == 0, f"{script} exited {proc.returncode}\nstdout: {proc.stdout}\nstderr: {proc.stderr}"
    assert proc.stdout.strip().endswith("OK")
