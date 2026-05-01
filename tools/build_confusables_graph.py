"""Emit a Mermaid `graph LR` block of mineral confusables.

Walks `minerals.CATALOG`, dedupes undirected pairs, and prints a Mermaid block
suitable for pasting into `docs/catalog.md`.

Usage:
    python tools/build_confusables_graph.py > /tmp/graph.mmd
"""

from __future__ import annotations

import argparse

from checkmsg import minerals


def collect_edges() -> list[tuple[str, str]]:
    seen: set[frozenset[str]] = set()
    edges: list[tuple[str, str]] = []
    for name, profile in minerals.CATALOG.items():
        for confusable in profile.confusables:
            if confusable not in minerals.CATALOG:
                continue
            pair = frozenset({name, confusable})
            if len(pair) == 2 and pair not in seen:
                seen.add(pair)
                a, b = sorted(pair)
                edges.append((a, b))
    return edges


def emit(edges: list[tuple[str, str]]) -> str:
    lines = ["graph LR"]
    # Pre-declare nodes that participate in many edges so layout is more legible.
    for a, b in edges:
        lines.append(f"    {a} --- {b}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser("build_confusables_graph")
    parser.add_argument("--mermaid", action="store_true",
                        help="wrap output in ```mermaid fenced block")
    args = parser.parse_args()

    edges = collect_edges()
    block = emit(edges)
    if args.mermaid:
        print("```mermaid")
        print(block, end="")
        print("```")
    else:
        print(block, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
