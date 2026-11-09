# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""The rgblastn flag matrix against blastn 2.17.0 (docker) on a store selection.

usage (from the gtars workspace, after `cargo build --release -p gtars-rgblast`):
  uv run gtars-rgblast/scripts/golden/search_matrix.py --store STORE --collection hg38-primary-refgenie \
      [--tasks megablast] [--cells 40 | --smoke] [--queries 60] [--out DIR]

Reference: an `rgmakeblastdb` database of the selection (byte-identical to makeblastdb over the
plan's reference FASTA, sub-plan 1), cached under --cache/<plan hash>/; `--sequences-file`
targets get their own cached subset database. Queries come from the classes of plan 4 step 22,
read from the reference database (`blastdbcmd -range`) with a fixed seed. Every cell runs
`blastn` and `rgblastn` with identical flags; byte identity (`cmp`) is the pass rule; failures
are classified with compare.py. Results go to --out/results.tsv (one line per cell) and
--out/cells/<cell>/.
"""
import argparse
import hashlib
import itertools
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare  # noqa: E402

IMAGE = os.environ.get("BLAST_IMAGE", "ncbi/blast:2.17.0")
ALL_FIELDS = "6 std stitle sallseqid salltitles staxid staxids sallacc qlen slen sstrand btop qseq sseq"

MATRIX = {
    "task": ["megablast", "blastn", "blastn-short"],
    "evalue": ["10", "1e-5", "1000"],
    "max_target_seqs": ["1", "5", "500", "10000"],
    "strand": ["both", "plus", "minus"],
    "masking": ["default", "dust_no", "soft_false"],
    "word_size": ["default", "16", "7"],
    "outfmt": ["6", ALL_FIELDS],
    "target": ["all", "seqs1", "seqs5", "seqs50"],
}


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def docker(workdir, *args):
    return ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-v", f"{workdir}:/w", "-w", "/w", IMAGE, *args]


def plan_hash(store, selection, sequences=None):
    meta = json.loads((Path(store) / "rgstore.json").read_text()) if not store.startswith("http") else {"url": store}
    h = hashlib.sha256()
    h.update(json.dumps({k: meta.get(k) for k in ("sequences_digest", "collections_digest", "aliases_digest", "fhr_digest", "url")}, sort_keys=True).encode())
    h.update(json.dumps(selection, sort_keys=True).encode())
    if sequences:
        h.update("\n".join(sorted(sequences)).encode())
    return h.hexdigest()[:16]


def reference_db(a, selection_args, sequences=None):
    """Build or reuse the rgmakeblastdb database; returns (dir, base name 'db')."""
    key = plan_hash(a.store, selection_args, sequences)
    d = Path(a.cache) / key
    if not (d / "db.nin").exists() and not (d / "db.nal").exists():
        d.mkdir(parents=True, exist_ok=True)
        cmd = [a.rgmakeblastdb, "build", "--store", a.store, *selection_args, "--out", str(d / "db"), "--title", "golden", "--quiet", "--force"]
        if sequences:
            (d / "subset.txt").write_text("".join(f"{s}\n" for s in sequences))
            cmd += ["--sequences-file", str(d / "subset.txt")]
        sh(cmd)
    return d


def db_entries(d):
    """[(SQ. id, length)] of a database (blastdbcmd in docker)."""
    out = sh(docker(d, "blastdbcmd", "-db", "db", "-entry", "all", "-outfmt", "%a %l")).stdout
    ents = []
    for line in out.splitlines():
        i, l = line.split()
        if i.startswith("SQ."):
            ents.append((i, int(l)))
    return ents


def db_totals(d):
    out = sh(docker(d, "blastdbcmd", "-db", "db", "-info")).stdout
    m = re.search(r"([\d,]+) sequences; ([\d,]+) total bases", out)
    return int(m.group(1).replace(",", "")), int(m.group(2).replace(",", ""))


def slice_seq(d, seqid, start, end):
    """Letters start..end (1-based, inclusive) of `seqid` in the reference database `d`."""
    out = sh(docker(d, "blastdbcmd", "-db", "db", "-entry", seqid, "-range", f"{start}-{end}")).stdout
    return "".join(out.splitlines()[1:]).upper()


def rc(s):
    comp = str.maketrans("ACGTRYMKWSBDHVNacgtrymkwsbdhvn", "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn")
    return s.translate(comp)[::-1]


def mutate(rnd, s, rate):
    s = list(s)
    for i in range(len(s)):
        if rnd.random() < rate:
            s[i] = rnd.choice([c for c in "ACGT" if c != s[i]])
    return "".join(s)


def indels(rnd, s, n):
    s = list(s)
    if len(s) < 30:
        return "".join(s)
    for _ in range(n):
        p = rnd.randrange(10, len(s) - 10)
        if rnd.random() < 0.5:
            del s[p:p + rnd.randint(1, 4)]
        else:
            s[p:p] = rnd.choices("ACGT", k=rnd.randint(1, 4))
    return "".join(s)


def make_queries(d, entries, seed, n):
    """The query classes of step 22 (a fixed-seed subset), as one FASTA text."""
    rnd = random.Random(seed)
    big = [e for e in entries if e[1] >= 60000] or entries
    recs = []

    def pick(minlen):
        c = [e for e in entries if e[1] >= minlen + 10] or big
        i, L = rnd.choice(c)
        ln = max(1, min(minlen, L - 2))
        st = rnd.randrange(1, max(2, L - ln))
        return i, st, slice_seq(d, i, st, min(L, st + ln - 1))

    k = 0
    while len(recs) < n:
        cls = k % 12
        k += 1
        if cls in (0, 1):
            i, st, s = pick(rnd.randint(300, 2000))
            recs.append((f"exact_{i[3:11]}_{st}", s))
        elif cls == 2:
            i, st, s = pick(rnd.randint(50, 100))
            recs.append((f"short_{i[3:11]}_{st}", s))
        elif cls == 3:
            i, st, s = pick(rnd.randint(10000, 50000))
            recs.append((f"long_{i[3:11]}_{st}", s))
        elif cls == 4:
            i, st, s = pick(rnd.randint(300, 2000))
            r = rnd.choice([0.03, 0.08, 0.15])
            recs.append((f"sub{int(r*100)}_{i[3:11]}_{st}", mutate(rnd, s, r)))
        elif cls == 5:
            i, st, s = pick(rnd.randint(300, 2000))
            recs.append((f"indel_{i[3:11]}_{st}", indels(rnd, s, rnd.randint(1, 5))))
        elif cls == 6:
            i, st, s = pick(rnd.randint(300, 2000))
            recs.append((f"rc_{i[3:11]}_{st}", rc(s)))
        elif cls == 7:
            i, st, s = pick(rnd.randint(300, 1000))
            s = list(s)
            for p in range(5, max(5, len(s) - 5), 13):
                s[p] = "N"
            recs.append((f"withN_{i[3:11]}_{st}", "".join(s)))
        elif cls == 8:
            recs.append((f"random_{k}", "".join(rnd.choices("ACGT", k=rnd.randint(300, 800)))))
        elif cls == 9:
            i, st, s = pick(rnd.randint(300, 1000))
            recs.append((f"lower_{i[3:11]}_{st}", s.lower()))
        elif cls == 10:
            i, st, s = pick(rnd.randint(300, 1000))
            recs.append((f"rna_{i[3:11]}_{st}", s.replace("T", "U")))
        else:
            i, st, s = pick(rnd.randint(300, 1000))
            recs.append((f"odd {i} defline | with pipes", s))
    recs.append(("", slice_seq(d, entries[0][0], 1, min(200, entries[0][1]))))  # empty defline
    fa = []
    for name, s in recs:
        fa.append(f">{name}\n")
        fa.extend(s[i:i + 70] + "\n" for i in range(0, len(s), 70))
    return "".join(fa)


def cell_name(c):
    return "_".join(re.sub(r"[^A-Za-z0-9]+", "-", str(v))[:12] for v in c.values())


def flags_for(c):
    f = ["-task", c["task"], "-evalue", c["evalue"], "-max_target_seqs", c["max_target_seqs"], "-strand", c["strand"], "-outfmt", c["outfmt"]]
    if c["masking"] == "dust_no":
        f += ["-dust", "no"]
    elif c["masking"] == "soft_false":
        f += ["-soft_masking", "false"]
    if c["word_size"] != "default":
        f += ["-word_size", c["word_size"]]
    return f


def run_cell(a, c, qfa, full_db, entries, out_dir, selection_args):
    cell_dir = out_dir / "cells" / cell_name(c)
    cell_dir.mkdir(parents=True, exist_ok=True)
    (cell_dir / "q.fa").write_text(qfa)
    rnd = random.Random(f"{c['target']}:{a.seed}")
    seqs = None
    db_dir = full_db
    if c["target"] != "all":
        n = int(c["target"][4:])
        seqs = [e[0] for e in rnd.sample(entries, min(n, len(entries)))]
        db_dir = reference_db(a, selection_args, seqs)
        (cell_dir / "subset.txt").write_text("".join(f"{s}\n" for s in seqs))
    flags = flags_for(c)
    # blastn in docker: the cell dir holds the query, the db dir is mounted at /db
    bcmd = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-v", f"{cell_dir}:/w", "-v", f"{db_dir}:/db:ro", "-w", "/w", IMAGE,
            "blastn", "-db", "/db/db", "-query", "q.fa", "-num_threads", "1", *flags]
    rb = subprocess.run(bcmd, capture_output=True, text=True)
    (cell_dir / "blastn.out").write_text(rb.stdout)
    (cell_dir / "blastn.err").write_text(rb.stderr)
    rcmd = [a.rgblastn, "--store", a.store, *selection_args, "-query", str(cell_dir / "q.fa"), *flags]
    if seqs:
        rcmd += ["--sequences-file", str(cell_dir / "subset.txt")]
    rr = subprocess.run(rcmd, capture_output=True, text=True)
    (cell_dir / "rgblastn.out").write_text(rr.stdout)
    (cell_dir / "rgblastn.err").write_text(rr.stderr)
    if rb.returncode != 0 or rr.returncode != 0:
        return "error", f"blastn exit {rb.returncode}, rgblastn exit {rr.returncode}: {rr.stderr.strip()[:200]}"
    ref, mine = rb.stdout, rr.stdout
    if ref == mine:
        return "identical", f"{len(ref.splitlines())} lines"
    (cell_dir / "ref.norm").write_text(ref)
    (cell_dir / "mine.norm").write_text(mine)
    r = compare.compare(cell_dir / "ref.norm", cell_dir / "mine.norm", c["outfmt"])
    summary = f"lines {r['lines_ref']}/{r['lines_test']} identical {r.get('identical_lines')} " + " ".join(f"{k}={v}" for k, v in sorted(r.get("kinds", {}).items()))
    if r.get("only_ref") or r.get("only_test"):
        summary += f" only_ref={len(r['only_ref'])} only_test={len(r['only_test'])}"
    return "differs", summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--store", required=True)
    ap.add_argument("--collection", action="append", default=[])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--tasks", default="megablast")
    ap.add_argument("--cells", type=int, default=0, help="sample this many cells (0 = the whole matrix)")
    ap.add_argument("--smoke", action="store_true", help="8 cells, 12 queries")
    ap.add_argument("--queries", type=int, default=48)
    ap.add_argument("--seed", type=int, default=20261108)
    ap.add_argument("--cache", default="/home/nsheff/work/refget-blast/tmp/golden")
    ap.add_argument("--out", default=None)
    ap.add_argument("--rgblastn", default=None)
    ap.add_argument("--rgmakeblastdb", default=None)
    a = ap.parse_args()
    target = Path(os.environ.get("CARGO_TARGET_DIR", HERE.parents[2] / "target"))
    a.rgblastn = a.rgblastn or str(target / "release" / "rgblastn")
    a.rgmakeblastdb = a.rgmakeblastdb or str(target / "release" / "rgmakeblastdb")
    selection_args = ["--all"] if a.all else list(itertools.chain.from_iterable(["--collection", c] for c in a.collection))
    if not selection_args:
        ap.error("give --collection or --all")
    tasks = [t for t in a.tasks.split(",")]
    out_dir = Path(a.out or (Path(a.cache) / "runs" / f"matrix-{a.seed}"))
    out_dir.mkdir(parents=True, exist_ok=True)

    full_db = reference_db(a, selection_args)
    entries = db_entries(full_db)
    nseq, total = db_totals(full_db)
    print(f"reference: {full_db} ({nseq} sequences, {total} bases)", file=sys.stderr)
    if a.smoke:
        a.queries = min(a.queries, 12)
    qfa = make_queries(full_db, entries, a.seed, a.queries)
    (out_dir / "queries.fa").write_text(qfa)

    keys = list(MATRIX)
    cells = [dict(zip(keys, vals)) for vals in itertools.product(*(MATRIX[k] for k in keys))]
    cells = [c for c in cells if c["task"] in tasks]

    def runnable(c):
        return not (c["task"] == "megablast" and c["word_size"] == "7")

    rnd = random.Random(a.seed)
    if a.smoke or a.cells:
        # sample runnable cells; the skipped ones are listed, not sampled
        pool = [c for c in cells if runnable(c)]
        cells = rnd.sample(pool, min(8 if a.smoke else a.cells, len(pool)))
    results = out_dir / "results.tsv"
    with open(results, "w") as fh:
        fh.write("cell\tstatus\tsummary\t" + "\t".join(keys) + "\n")
        counts = {}
        for i, c in enumerate(cells, 1):
            if not runnable(c):
                status, summary = "skipped", "megablast needs a word size of at least 12"
            else:
                status, summary = run_cell(a, c, qfa, full_db, entries, out_dir, selection_args)
            counts[status] = counts.get(status, 0) + 1
            fh.write(f"{cell_name(c)}\t{status}\t{summary}\t" + "\t".join(str(c[k]) for k in keys) + "\n")
            fh.flush()
            print(f"[{i}/{len(cells)}] {status:9s} {cell_name(c)}  {summary}", file=sys.stderr)
    print(f"done: {counts} -> {results}", file=sys.stderr)
    sys.exit(0 if counts.get("differs", 0) + counts.get("error", 0) == 0 else 1)


if __name__ == "__main__":
    main()
