"""Smoke tests for the doc-generation tools.

`build_schematics --check` verifies that the committed PNGs match what the
generator produces today (catches accidental aesthetic drift). `build_confusables_graph`
must emit syntactically-valid Mermaid that covers every `MineralProfile` listing
a confusable.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
TOOL_SCHEMATICS = REPO_ROOT / "tools" / "build_schematics.py"
TOOL_GRAPH = REPO_ROOT / "tools" / "build_confusables_graph.py"
FIGURES_DIR = REPO_ROOT / "docs" / "figures"


@pytest.mark.skipif(not FIGURES_DIR.exists(), reason="docs/figures/ not yet generated")
def test_build_schematics_check_matches_committed():
    proc = subprocess.run(
        [sys.executable, str(TOOL_SCHEMATICS), "--check", "--output", str(FIGURES_DIR)],
        capture_output=True, text=True, timeout=60,
    )
    assert proc.returncode == 0, (
        f"build_schematics --check failed; the schematics on disk no longer match "
        f"what the script produces.\nstdout:{proc.stdout}\nstderr:{proc.stderr}"
    )


def test_build_confusables_graph_is_valid_mermaid():
    proc = subprocess.run(
        [sys.executable, str(TOOL_GRAPH)],
        capture_output=True, text=True, timeout=15,
    )
    assert proc.returncode == 0
    out = proc.stdout
    assert out.startswith("graph LR\n"), "Mermaid block must open with 'graph LR'"
    edge_re = re.compile(r"^    \w+ --- \w+$")
    body_lines = out.splitlines()[1:]
    for line in body_lines:
        if not line.strip():
            continue
        assert edge_re.match(line), f"unexpected line in graph output: {line!r}"


def test_confusables_graph_covers_catalog():
    """Every catalog entry whose `confusables` lists another existing entry should
    appear in the graph."""
    from checkmsg import minerals
    proc = subprocess.run(
        [sys.executable, str(TOOL_GRAPH)],
        capture_output=True, text=True, timeout=15,
    )
    out = proc.stdout
    referenced = set()
    for line in out.splitlines()[1:]:
        if "---" in line:
            a, _, b = line.strip().partition(" --- ")
            referenced.add(a)
            referenced.add(b)
    expected = set()
    for name, profile in minerals.CATALOG.items():
        for c in profile.confusables:
            if c in minerals.CATALOG:
                expected.add(name)
                expected.add(c)
    missing = expected - referenced
    assert not missing, f"catalog entries with confusables missing from graph: {missing}"


def test_build_confusables_graph_mermaid_flag():
    """The --mermaid flag wraps output in a fenced code block."""
    proc = subprocess.run(
        [sys.executable, str(TOOL_GRAPH), "--mermaid"],
        capture_output=True, text=True, timeout=15,
    )
    assert proc.returncode == 0
    out = proc.stdout
    assert out.startswith("```mermaid\n")
    assert out.rstrip().endswith("```")


# --- Image-reference presence check ---------------------------------------------

DOC_FILES = [REPO_ROOT / "README.md", *sorted((REPO_ROOT / "docs").glob("*.md"))]
IMAGE_RE = re.compile(r"!\[[^\]]*\]\(([^)]+\.png)\)")


def test_every_doc_image_reference_resolves():
    """Every `![...](path.png)` in the README and docs/ must resolve to a file on disk.

    Catches the failure mode where a doc references a gitignored or deleted PNG —
    the case that motivated moving example plots into `docs/figures/examples/`.
    """
    missing: list[tuple[str, str]] = []
    for doc in DOC_FILES:
        if not doc.exists():
            continue
        text = doc.read_text(encoding="utf-8")
        for match in IMAGE_RE.finditer(text):
            ref = match.group(1).strip()
            # Skip absolute URLs and anchor-only fragments.
            if ref.startswith(("http://", "https://", "data:")):
                continue
            target = (doc.parent / ref).resolve()
            if not target.exists():
                missing.append((str(doc.relative_to(REPO_ROOT)), ref))
    assert not missing, (
        "broken image references in documentation:\n"
        + "\n".join(f"  {d}: {r}" for d, r in missing)
    )
