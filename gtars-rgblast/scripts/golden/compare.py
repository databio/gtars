# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Classify the differences between two BLAST tabular outputs (ported from the PoC).

usage: compare.py <reference> <test> [--outfmt "6 std stitle"] [--json]

Lines starting with '#' (-outfmt 7 comments) are skipped; -outfmt 10 is split on commas.
Reports: per-query hit counts, exact line agreement (order-free), and for HSPs keyed by
(query, subject, qstart, qend, sstart, send) the kind of difference:
  score      evalue / bit score / raw score differ
  arrangement same endpoints and scores, other columns differ (an equal-score gap arrangement)
  format     the same numbers printed differently (e.g. 1e-05 vs 1.00e-05)
"""
import argparse
import collections
import json
import sys

STD = ["qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
       "sstart", "send", "evalue", "bitscore"]
SCORE_COLS = {"evalue", "bitscore", "score"}


def fields_of(outfmt):
    toks = outfmt.split()
    fmt, names = toks[0], toks[1:]
    out = []
    for t in names:
        out.extend(STD if t == "std" else [t])
    if not out:
        out = list(STD)
    return fmt, out


def load(path, fmt, nfields):
    sep = "," if fmt == "10" else "\t"
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split(sep)
            if len(f) >= min(nfields, 2):
                rows.append(f)
    return rows


def col(names, *cands):
    for c in cands:
        if c in names:
            return names.index(c)
    return None


def num(s):
    try:
        return float(s)
    except ValueError:
        return s


def compare(ref_path, test_path, outfmt="6"):
    fmt, names = fields_of(outfmt)
    ref, test = load(ref_path, fmt, len(names)), load(test_path, fmt, len(names))
    r = {"lines_ref": len(ref), "lines_test": len(test), "queries": {}}
    cref = collections.Counter(x[0] for x in ref)
    ctest = collections.Counter(x[0] for x in test)
    for q in sorted(set(cref) | set(ctest)):
        r["queries"][q] = (cref[q], ctest[q])
    lref = collections.Counter("\t".join(x) for x in ref)
    ltest = collections.Counter("\t".join(x) for x in test)
    r["identical_lines"] = sum((lref & ltest).values())
    iq = col(names, "qseqid", "qaccver", "qacc")
    isub = col(names, "sseqid", "saccver", "sacc")
    kcols = [col(names, "qstart"), col(names, "qend"), col(names, "sstart"), col(names, "send")]
    if iq is None or isub is None or any(k is None for k in kcols):
        r["note"] = "no HSP key columns in this -outfmt; only line agreement reported"
        return r

    def key(x):
        return (x[iq], x[isub]) + tuple(x[k] for k in kcols)

    kref = {key(x): x for x in ref}
    ktest = {key(x): x for x in test}
    both = set(kref) & set(ktest)
    r["only_ref"] = sorted(set(kref) - set(ktest))
    r["only_test"] = sorted(set(ktest) - set(kref))
    score_idx = [i for i, n in enumerate(names) if n in SCORE_COLS]
    kinds = collections.Counter()
    examples = collections.defaultdict(list)
    for k in sorted(both):
        a, b = kref[k], ktest[k]
        if a == b:
            continue
        if any(a[i] != b[i] for i in score_idx):
            kind = "format" if all(num(a[i]) == num(b[i]) for i in score_idx) else "score"
        else:
            kind = "arrangement"
        kinds[kind] += 1
        if len(examples[kind]) < 5:
            examples[kind].append(("\t".join(a), "\t".join(b)))
    r["shared"] = len(both)
    r["kinds"] = dict(kinds)
    r["examples"] = dict(examples)
    return r


def report(r):
    print(f"lines: ref={r['lines_ref']} test={r['lines_test']} identical={r.get('identical_lines')}")
    for q, (a, b) in r["queries"].items():
        flag = "" if a == b else "  <-- differs"
        print(f"  {q[:40]:40s} ref={a:6d} test={b:6d}{flag}")
    if "note" in r:
        print(r["note"])
        return
    print(f"HSP keys: shared={r['shared']} only_ref={len(r['only_ref'])} only_test={len(r['only_test'])} "
          + " ".join(f"{k}={v}" for k, v in sorted(r["kinds"].items())))
    for name in ("only_ref", "only_test"):
        for k in r[name][:5]:
            print(f"  {name}: {k}")
    for kind, ex in r["examples"].items():
        for a, b in ex[:3]:
            print(f"  {kind}:\n    ref  {a[:300]}\n    test {b[:300]}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference")
    ap.add_argument("test")
    ap.add_argument("--outfmt", default="6")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    r = compare(a.reference, a.test, a.outfmt)
    if a.json:
        json.dump(r, sys.stdout, indent=1)
        print()
    else:
        report(r)
    ok = r["lines_ref"] == r["lines_test"] == r.get("identical_lines")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
