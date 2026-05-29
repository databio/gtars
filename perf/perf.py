#!/usr/bin/env python3
"""gtars-perf — occasionally-run performance regression harness (Python).

This replaces the old `gtars-perf/` Rust crate. It drives the real `gtars`
CLI as a subprocess under `/usr/bin/time -v`, so we capture true *peak RSS*
for every task it can measure (the crate returned null RSS for its in-process
extract/VRS tasks). The run-record JSON schema and the compare/seed gate are
kept identical to the crate (schema_version 1) so old tooling keeps working.

Tasks
-----
  encode : `gtars refget build <fasta>... -o DIR -j N` over a MULTI-file input
           set, swept across `--jobs` (the -j knob parallelizes across FILES,
           so >=2 files are needed for the sweep to mean anything). Real wall
           clock + real peak RSS from `/usr/bin/time -v`. Correctness is the
           store's `collections_digest` (read from rgstore.json) plus the
           bases/sequence counts parsed from the CLI's own summary line.

  extract: region/substring extraction from a refget store. *** No CLI surface
  vrs    : VCF -> VRS allele ids.                          *** exists today. ***
           The original crate measured these in-process via the library, which
           is exactly why it could not report per-task RSS. Rather than fake a
           sweep, we record an explicit "skipped" result noting the gap. Adding
           `gtars refget extract` / `gtars vrs --threads N` to the CLI would let
           this harness measure them for real (see perf/README.md).

Subcommands
-----------
  run           run measurable tasks, write a run-record JSON
  compare       gate a run JSON against targets.json (nonzero exit on regress)
  seed-targets  derive targets.json from a baseline run JSON

Exit codes (compare): 0 pass, 1 regression, 2 bad input.
Stdlib only — no pip installs.
"""

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

SCHEMA_VERSION = 1


# ─────────────────────────── machine / run metadata ───────────────────────────

def _git(repo, args):
    try:
        out = subprocess.run(["git", "-C", str(repo), *args],
                             capture_output=True, text=True)
        return out.stdout.strip() if out.returncode == 0 else None
    except OSError:
        return None


def _cpu_model():
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def _total_ram_mb():
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) // 1024
    except (OSError, ValueError):
        pass
    return 0


def _rustc_version():
    try:
        out = subprocess.run(["rustc", "--version"], capture_output=True, text=True)
        return out.stdout.strip() if out.returncode == 0 else "unknown"
    except OSError:
        return "unknown"


def collect_run_meta(repo):
    return {
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "gtars_commit": _git(repo, ["rev-parse", "--short", "HEAD"]) or "unknown",
        "gtars_commit_date": _git(repo, ["show", "-s", "--format=%cI", "HEAD"]) or "",
        "gtars_dirty": bool(_git(repo, ["status", "--porcelain"])),
        "host": platform.node() or "unknown",
        "cpu_model": _cpu_model(),
        "logical_cpus": os.cpu_count() or 0,
        "total_ram_mb": _total_ram_mb(),
        "rustc_version": _rustc_version(),
        "profile": "release",  # we always drive the release CLI binary
    }


# ──────────────────────────── /usr/bin/time -v driver ─────────────────────────

_RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")


def run_timed(argv):
    """Run `argv` under /usr/bin/time -v. Returns (seconds, peak_rss_mb, output).

    `output` is the child's combined stdout+stderr (the gtars CLI writes its
    summary lines to stderr). peak_rss_mb is None if /usr/bin/time is
    unavailable or unparseable. Raises RuntimeError on a nonzero child exit.
    """
    have_time = shutil.which("/usr/bin/time") is not None
    with tempfile.NamedTemporaryFile("r", suffix=".time", delete=False) as tf:
        time_log = tf.name
    try:
        wrapped = (["/usr/bin/time", "-v", "-o", time_log] + argv) if have_time else argv
        start = time.perf_counter()
        proc = subprocess.run(wrapped, capture_output=True, text=True)
        seconds = time.perf_counter() - start
        # /usr/bin/time -v writes ITS report to the -o file, so the child's own
        # stderr stays clean for us to mine.
        child_output = (proc.stdout or "") + (proc.stderr or "")
        if proc.returncode != 0:
            raise RuntimeError(
                f"command failed (exit {proc.returncode}): {' '.join(argv)}\n"
                f"{child_output}")
        peak_rss_mb = None
        if have_time:
            m = _RSS_RE.search(Path(time_log).read_text())
            if m:
                peak_rss_mb = int(m.group(1)) / 1024.0
        return seconds, peak_rss_mb, child_output
    finally:
        try:
            os.unlink(time_log)
        except OSError:
            pass


# ─────────────────────────────────── tasks ────────────────────────────────────

# `gtars refget build` summary lines we mine for correctness + throughput:
#   "Done: 3 sequences, 16 bases in 0.001s (...)"
_DONE_RE = re.compile(r"Done:\s*([\d,]+)\s+sequences?,\s*([\d,]+)\s+bases", re.I)


def _clear_rgsi_sidecars(fastas):
    """The importer caches `<fasta>.rgsi` sidecars; remove them so each encode
    run actually re-does the digest+encode work (matches the crate's behaviour)."""
    for f in fastas:
        f = Path(f)
        parent = f.parent
        if not parent.is_dir():
            continue
        for sib in parent.iterdir():
            if sib.name.startswith(f.name) and sib.name.endswith(".rgsi"):
                try:
                    sib.unlink()
                except OSError:
                    pass


def task_encode(gtars_bin, fastas, jobs, dataset_id):
    """Encode task: `gtars refget build <fastas> -o DIR -j jobs` under time -v."""
    _clear_rgsi_sidecars(fastas)
    out_dir = Path(tempfile.mkdtemp(prefix=f"gtars-perf-encode-j{jobs}-"))
    shutil.rmtree(out_dir, ignore_errors=True)
    argv = [str(gtars_bin), "refget", "build", *[str(f) for f in fastas],
            "-o", str(out_dir), "-j", str(jobs)]
    try:
        seconds, peak_rss_mb, output = run_timed(argv)

        # Correctness + size from the CLI's own summary and the store manifest.
        n_sequences = n_bases = 0
        m = _DONE_RE.search(output)
        if m:
            n_sequences = int(m.group(1).replace(",", ""))
            n_bases = int(m.group(2).replace(",", ""))
        collection_digest = ""
        manifest = out_dir / "rgstore.json"
        if manifest.is_file():
            try:
                collection_digest = json.loads(manifest.read_text()).get(
                    "collections_digest", "")
            except (OSError, ValueError):
                pass

        throughput = n_bases / seconds if seconds > 0 else 0.0
        return {
            "task": "encode",
            "concurrency": jobs,
            "seconds": seconds,
            "peak_rss_mb": peak_rss_mb,
            "throughput": throughput,
            "throughput_unit": "bases_per_sec",
            "extra": {
                "dataset_id": dataset_id,
                "n_bases": n_bases,
                "n_sequences": n_sequences,
                "n_files": len(fastas),
                "collection_digest": collection_digest,
            },
        }
    finally:
        shutil.rmtree(out_dir, ignore_errors=True)


def task_skipped(task, concurrency, dataset_id, reason):
    """A recognized-but-unmeasurable task. Recorded honestly, not faked.

    seconds=0 / throughput=0 / peak_rss_mb=None. `extra.skipped` is true so the
    gate and any downstream reader can distinguish it from a real measurement.
    """
    return {
        "task": task,
        "concurrency": concurrency,
        "seconds": 0.0,
        "peak_rss_mb": None,
        "throughput": 0.0,
        "throughput_unit": "n/a",
        "extra": {"dataset_id": dataset_id, "skipped": True, "reason": reason},
    }


# ──────────────────────────────────── run ─────────────────────────────────────

EXTRACT_VRS_GAP = (
    "no `gtars` CLI subcommand exposes this task today; only `refget build` "
    "(encode) is measurable via subprocess. Add a CLI extract/vrs command "
    "(vrs ideally with a --threads knob) to measure this with real peak RSS."
)


def cmd_run(a):
    fastas = [p for chunk in a.fasta for p in chunk.split(",") if p]
    if not fastas:
        sys.exit("ERROR: provide --fasta <files...> (comma- or space-separated) "
                 "for the encode jobs-sweep")
    missing = [f for f in fastas if not Path(f).is_file()]
    if missing:
        sys.exit(f"ERROR: FASTA not found: {', '.join(missing)}")
    if len(fastas) < 2:
        sys.stderr.write(
            "WARNING: encode -j only parallelizes across FILES; "
            f"{len(fastas)} file(s) supplied — the jobs sweep will be flat.\n")

    results = []
    for j in a.concurrency:
        sys.stderr.write(f"[encode] jobs={j}\n")
        results.append(task_encode(a.gtars_bin, fastas, j, a.dataset_id))

    # extract / vrs: no CLI surface — record honest skips, flag the gap.
    if a.bed:
        sys.stderr.write(f"[extract] SKIPPED — {EXTRACT_VRS_GAP}\n")
        results.append(task_skipped("extract", 1, a.dataset_id, EXTRACT_VRS_GAP))
    if a.vcf:
        for t in a.concurrency:
            sys.stderr.write(f"[vrs] threads={t} SKIPPED — {EXTRACT_VRS_GAP}\n")
            results.append(task_skipped("vrs", t, a.dataset_id, EXTRACT_VRS_GAP))

    record = {
        "schema_version": SCHEMA_VERSION,
        "run": collect_run_meta(a.repo),
        "results": results,
    }

    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = record["run"]["timestamp_utc"].replace(":", "").replace("-", "")
    json_path = out_dir / f"{record['run']['gtars_commit']}-{ts}.json"
    json_path.write_text(json.dumps(record, indent=2))
    sys.stderr.write(f"Wrote run record: {json_path}\n")

    # Summary table.
    sys.stderr.write("\ntask        conc   seconds   peak_rss_mb   throughput  unit\n")
    for r in results:
        rss = f"{r['peak_rss_mb']:.0f}" if r["peak_rss_mb"] is not None else "null"
        sys.stderr.write(
            f"{r['task']:<10}  {r['concurrency']:>4}  {r['seconds']:>8.3f}   "
            f"{rss:>11}  {r['throughput']:>11.0f}  {r['throughput_unit']}\n")
    print(json_path)


# ──────────────────────────────── compare / gate ──────────────────────────────

def _check_cell(r, c, margin):
    """Return (pass, binding-description). Mirrors the crate's check_cell."""
    msgs, ok = [], True
    if c.get("min_throughput") is not None:
        floor = c["min_throughput"] * (1.0 - margin)
        if r["throughput"] < floor:
            ok = False
            msgs.append(f"throughput {r['throughput']:.0f} < {floor:.0f} "
                        f"(min*{1.0 - margin:.2f})")
    if c.get("max_seconds") is not None:
        ceil = c["max_seconds"] * (1.0 + margin)
        if r["seconds"] > ceil:
            ok = False
            msgs.append(f"seconds {r['seconds']:.2f} > {ceil:.2f}")
    if c.get("max_peak_rss_mb") is not None and r["peak_rss_mb"] is not None:
        ceil = c["max_peak_rss_mb"] * (1.0 + margin)
        if r["peak_rss_mb"] > ceil:
            ok = False
            msgs.append(f"peak_rss {r['peak_rss_mb']:.0f} > {ceil:.0f}")
    return (ok, "all constraints OK" if ok else "; ".join(msgs))


def _load_json(path, what):
    try:
        return json.loads(Path(path).read_text())
    except OSError as e:
        sys.stderr.write(f"ERROR: cannot read {what} {path}: {e}\n")
        sys.exit(2)
    except ValueError as e:
        sys.stderr.write(f"ERROR: cannot parse {what} {path}: {e}\n")
        sys.exit(2)


def cmd_compare(a):
    run = _load_json(a.run, "run")
    if run.get("schema_version") != SCHEMA_VERSION:
        sys.stderr.write(
            f"ERROR: run schema_version {run.get('schema_version')} != harness "
            f"schema_version {SCHEMA_VERSION}; re-run / re-seed.\n")
        sys.exit(2)
    targets = _load_json(a.targets, "targets")
    if "tasks" not in targets or "margin" not in targets:
        sys.stderr.write("ERROR: targets JSON missing 'margin'/'tasks'.\n")
        sys.exit(2)

    seeded_host = targets.get("seeded_on_host", "")
    run_host = run.get("run", {}).get("host", "")
    if seeded_host and seeded_host != run_host:
        sys.stderr.write(
            f"WARNING: run host '{run_host}' differs from targets host "
            f"'{seeded_host}'. Targets are machine-specific; re-seed here.\n")

    margin = targets["margin"]
    print(f"{'task':<10} {'conc':>5} {'observed':>10} {'result':<10} binding constraint")
    any_fail = any_gated = False
    for r in run.get("results", []):
        c = targets["tasks"].get(r["task"], {}).get(str(r["concurrency"]))
        if c is None:
            print(f"{r['task']:<10} {r['concurrency']:>5} {'-':>10} "
                  f"{'WARN':<10} (no target — informational)")
            continue
        any_gated = True
        ok, binding = _check_cell(r, c, margin)
        any_fail = any_fail or not ok
        print(f"{r['task']:<10} {r['concurrency']:>5} {r['throughput']:>10.0f} "
              f"{'PASS' if ok else 'FAIL':<10} {binding}")

    if not any_gated:
        sys.stderr.write("No gated cells matched any target. "
                         "Treating as pass (informational only).\n")
    if any_fail:
        sys.stderr.write("\nGATE: FAIL — at least one cell regressed.\n")
        sys.exit(1)
    sys.stderr.write("\nGATE: PASS\n")


def cmd_seed(a):
    run = _load_json(a.run, "run")
    tasks = {}
    for r in run.get("results", []):
        # Don't seed targets from honest-skip rows (they'd gate to zero floors).
        if r.get("extra", {}).get("skipped"):
            continue
        tasks.setdefault(r["task"], {})[str(r["concurrency"])] = {
            "min_throughput": r["throughput"],
            "max_seconds": r["seconds"],
            **({"max_peak_rss_mb": r["peak_rss_mb"]}
               if r["peak_rss_mb"] is not None else {}),
        }
    targets = {
        "margin": a.margin,
        "seeded_on_host": run.get("run", {}).get("host", ""),
        "tasks": tasks,
    }
    out = Path(a.out)
    if out.parent:
        out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(targets, indent=2))
    sys.stderr.write(
        f"Seeded targets from {a.run} on host "
        f"'{targets['seeded_on_host']}' (margin {a.margin}) -> {out}\n")


# ──────────────────────────────────── cli ─────────────────────────────────────

def build_parser():
    p = argparse.ArgumentParser(
        prog="perf.py", description="gtars performance regression harness")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run tasks and emit a run-record JSON")
    r.add_argument("--fasta", action="append", default=[], metavar="F",
                   help="FASTA file(s) for the encode jobs-sweep; repeatable "
                        "and/or comma-separated (>=2 for the sweep to matter)")
    r.add_argument("--bed", help="BED file for the extract task (currently "
                                 "unmeasurable — see notes)")
    r.add_argument("--vcf", help="VCF file for the vrs task (currently "
                                 "unmeasurable — see notes)")
    r.add_argument("--concurrency", default="1,4,8",
                   type=lambda s: [int(x) for x in s.split(",") if x],
                   help="comma-separated concurrency sweep (default 1,4,8)")
    r.add_argument("--out", default="perf/runs", help="output dir for run JSON")
    r.add_argument("--gtars-bin", default="target/release/gtars",
                   help="path to the gtars CLI binary")
    r.add_argument("--repo", default=".", help="repo root for commit metadata")
    r.add_argument("--dataset-id", default="adhoc",
                   help="short label stamped into every result's extra")
    r.set_defaults(func=cmd_run)

    c = sub.add_parser("compare", help="gate a run JSON against targets")
    c.add_argument("run", help="run JSON to gate")
    c.add_argument("--targets", default="perf/targets.json", help="targets JSON")
    c.set_defaults(func=cmd_compare)

    s = sub.add_parser("seed-targets", help="derive a targets.json from a baseline")
    s.add_argument("run", help="baseline run JSON")
    s.add_argument("--out", default="perf/targets.json", help="where to write targets")
    s.add_argument("--margin", type=float, default=0.15,
                   help="tolerance applied at compare time (default 0.15)")
    s.set_defaults(func=cmd_seed)
    return p


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
