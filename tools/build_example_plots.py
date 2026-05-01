"""Regenerate the example plots committed under `docs/figures/examples/`.

Each example script (`examples/NN_*.py`) saves its plot to `examples/output/NN_*.png`
when run in non-smoke mode. This tool runs every example end-to-end and copies the
resulting PNG into the documentation asset tree at `docs/figures/examples/`, which
is the path referenced by README.md, docs/techniques.md, and docs/curriculum.md.

Total runtime is ~60–90 s for all 20 examples, so this is a one-off rebuild rather
than a CI step. The presence-check test in `tests/test_docs_tools.py` is the
guardrail that catches stale or missing image references.

Usage:
    python tools/build_example_plots.py                     # rebuild all 20
    python tools/build_example_plots.py --only 19,20        # rebuild a subset
    python tools/build_example_plots.py --output docs/figures/examples
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = REPO_ROOT / "examples"
TRANSIENT_DIR = REPO_ROOT / "examples" / "output"
DEFAULT_OUTPUT = REPO_ROOT / "docs" / "figures" / "examples"


def _resolve_only(only_arg: str | None) -> set[str] | None:
    if not only_arg:
        return None
    tokens = [t.strip() for t in only_arg.split(",") if t.strip()]
    expanded: set[str] = set()
    for tok in tokens:
        # Accept "19", "19_unknown_stone_capstone", or full filename.
        if tok.endswith(".py"):
            tok = tok[:-3]
        expanded.add(tok)
    return expanded


def _matches(stem: str, only: set[str] | None) -> bool:
    if not only:
        return True
    if stem in only:
        return True
    # Accept the leading two-digit number alone.
    prefix = stem.split("_", 1)[0]
    return prefix in only


def build_all(output_dir: Path, *, only: set[str] | None = None,
              python: str = sys.executable) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    scripts = sorted(EXAMPLES_DIR.glob("[0-9][0-9]_*.py"))
    for script in scripts:
        stem = script.stem
        if not _matches(stem, only):
            continue
        print(f"  running {script.name}", file=sys.stderr)
        result = subprocess.run([python, str(script)],
                                cwd=REPO_ROOT, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            print(result.stdout, file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            raise RuntimeError(f"{script.name} exited with {result.returncode}")
        png = TRANSIENT_DIR / f"{stem}.png"
        if not png.exists():
            print(f"  warning: {script.name} produced no PNG", file=sys.stderr)
            continue
        target = output_dir / png.name
        shutil.copy2(png, target)
        written.append(target)
        print(f"  wrote {target.relative_to(REPO_ROOT)}", file=sys.stderr)
    return written


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser("build_example_plots")
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                   help=f"target directory for committed PNGs (default: {DEFAULT_OUTPUT.relative_to(REPO_ROOT)})")
    p.add_argument("--only", type=str, default=None,
                   help="comma-separated list of example numbers or stems to rebuild")
    p.add_argument("--python", type=str, default=sys.executable,
                   help="Python executable to run the example scripts (default: current)")
    args = p.parse_args(argv)

    written = build_all(args.output, only=_resolve_only(args.only), python=args.python)
    print(f"\nrebuilt {len(written)} example plot(s) into {args.output.relative_to(REPO_ROOT)}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
