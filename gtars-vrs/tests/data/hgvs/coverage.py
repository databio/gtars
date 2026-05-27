#!/usr/bin/env python3
"""Coverage matrix for vendored HGVS test corpora.

Walks the parser-layer corpora and tallies (reference_type x edit_type)
cells, plus axis flags for forward/reverse strand (n/a here — no real
provider), intronic offset, UTR position, and gene-symbol form.

Run from the gtars repo root:

    python3 gtars-vrs/tests/data/hgvs/coverage.py

This is intentionally a one-shot audit utility. It does not import
anything from the Rust parser; it does a regex-based classification
sufficient for matrix bookkeeping. Use it to spot empty cells worth a
hand-authored gap-fill row in `gap_fill.tsv`.
"""

from __future__ import annotations

import os
import re
import sys
from collections import defaultdict
from pathlib import Path


REF_TYPES = ("g", "c", "n", "m", "r", "p")
EDIT_TYPES = ("sub", "del", "dup", "ins", "delins", "inv", "equals", "unknown")

EXPR_RE = re.compile(
    r"^(?P<acc>[A-Za-z][A-Za-z0-9_.]*)(\([^)]+\))?:(?P<rt>[gcnmrp])\.(?P<rest>.+)$"
)


def classify_edit(rest: str) -> str:
    s = rest.strip().rstrip(")")
    if "delins" in s:
        return "delins"
    if "del" in s and "ins" in s:
        return "delins"
    if s.endswith("="):
        return "equals"
    if s.endswith("?"):
        return "unknown"
    if "del" in s:
        return "del"
    if "dup" in s:
        return "dup"
    if "ins" in s:
        return "ins"
    if "inv" in s:
        return "inv"
    if ">" in s:
        return "sub"
    return "unknown"


def has_intronic(rest: str) -> bool:
    # Look for a digit followed by +N or -N (excluding lone leading minus).
    return bool(re.search(r"\d[+-]\d", rest))


def has_utr(rest: str) -> bool:
    return bool(re.match(r"^[*-]\d", rest)) or bool(re.search(r"[._][*-]\d", rest))


def is_gene_form(acc: str) -> bool:
    # Heuristic: gene symbols are short, all-caps letters/digits, no underscore
    # or dot. Conventional NM_/NC_/NR_ accessions have an underscore.
    return bool(re.match(r"^[A-Z][A-Z0-9-]{1,15}$", acc)) and "_" not in acc


def iter_lines(path: Path):
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        t = line.strip()
        if not t or t.startswith("#"):
            yield from ()
            continue
        yield lineno, t


def gather() -> list[str]:
    here = Path(__file__).parent
    sources = []
    g = here / "varfish/parser/gauntlet"
    if g.exists():
        for _, line in iter_lines(g):
            sources.append(line)
    grammar = here / "biocommons/grammar_test.tsv"
    if grammar.exists():
        for line in grammar.read_text().splitlines():
            if not line or line.startswith("#") or line.startswith("Func\t"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            func, test, valid, in_type = cols[0], cols[1], cols[2], cols[3]
            if func not in {"hgvs_variant", "c_variant", "g_variant",
                            "n_variant", "m_variant", "r_variant", "p_variant"}:
                continue
            if valid.strip() != "True":
                continue
            inputs = test.split("|") if in_type.strip() == "list" else [test]
            sources.extend(inputs)
    gap = here / "gap_fill.tsv"
    if gap.exists():
        for _, line in iter_lines(gap):
            sources.append(line.split("\t")[0])
    return sources


def main() -> int:
    matrix: dict[tuple[str, str], int] = defaultdict(int)
    flags = {"intronic": 0, "utr": 0, "gene_form": 0, "total": 0}
    for raw in gather():
        m = EXPR_RE.match(raw.strip())
        if not m:
            continue
        rt = m.group("rt")
        rest = m.group("rest")
        edit = classify_edit(rest)
        matrix[(rt, edit)] += 1
        flags["total"] += 1
        if has_intronic(rest):
            flags["intronic"] += 1
        if has_utr(rest):
            flags["utr"] += 1
        if is_gene_form(m.group("acc")):
            flags["gene_form"] += 1

    print(f"Coverage matrix (reference_type x edit_type):")
    print(f"{'':6}" + "".join(f"{e:>8}" for e in EDIT_TYPES))
    empty = []
    for rt in REF_TYPES:
        row = [f"{matrix[(rt, e)]:>8}" for e in EDIT_TYPES]
        print(f"{rt:6}" + "".join(row))
        for e in EDIT_TYPES:
            if matrix[(rt, e)] == 0:
                empty.append((rt, e))
    print()
    print(f"Total expressions classified: {flags['total']}")
    print(f"  Intronic offset present:    {flags['intronic']}")
    print(f"  UTR position (c.-N / c.*N): {flags['utr']}")
    print(f"  Gene-symbol form (no NM_):  {flags['gene_form']}")
    print()
    if empty:
        print("Empty cells (consider adding to gap_fill.tsv):")
        for rt, e in empty:
            print(f"  - {rt}.{e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
