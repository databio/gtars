#!/usr/bin/env python3
"""gtars-perf — fast, refgetstore-only performance-regression suite (Python).

This drives the real `gtars` CLI (encode) and a focused in-process example
binary (`perf_task`, in gtars-refget/examples) as subprocesses under
`/usr/bin/time -v`, so we capture true *wall time* and *peak RSS* for every
task/scenario/path. The whole suite is tuned to finish in well under 30s.

It measures ONLY refgetstore — no competitor tools — across the three real
tasks of the refget benchmark:

Tasks / scenarios / paths
-------------------------
  encode  : `gtars refget build <fasta>... -o DIR -j N`. `-j` parallelizes
            across FILES, so the suite builds a SMALL-but-REAL multi-file input
            (chr20/chr21/chr22 split out of GRCh38) and sweeps jobs=1 (single
            core) and jobs=N (multi core). Records wall + peak RSS + bases/s.
            paths: "-"   scenarios: "multifile"   concurrency: 1 and N

  extract : region/substring extraction over the PREBUILT full GRCh38 store.
            Three scenarios, each measured on BOTH read paths:
              small       : ~10k narrow regions
              large_count : ~100k narrow regions
              large_width : ~2k regions of 100kb-1Mb (real, full-width slices)
            paths:
              resident : load_sequence whole touched seqs, then get_substring()
              partial  : never load_sequence; get_substring() on a stub reads
                         only the covering bytes from disk per region
              batch    : stub-only; group ranges by sequence and call
                         get_substrings() once per sequence (adaptive
                         whole-decode + unchecked-UTF-8)
            Records wall + peak RSS + bases/s.

  vrs     : the VRS/HGVS point-lookup pattern (1bp lookups), resident path, over
            a bounded real variant subset (clinvar points). Records wall +
            lookups/s.

Measurement entrypoint
----------------------
  encode  -> the `gtars` CLI (`refget build`), parsed via its summary line.
  extract / vrs -> `examples/perf_task`, which runs exactly ONE
            task+scenario+path and prints a single machine-readable line:
              RESULT task=.. scenario=.. path=.. seconds=.. items=.. \
                     bases=.. throughput=.. unit=..
            perf.py invokes it under `run_timed` (so peak RSS is captured too).

Subcommands
-----------
  run           run the whole suite, write a run-record JSON (schema_version 2)
  compare       gate a run JSON vs targets.json + HIGHLIGHT regressions in
                plain language (nonzero exit on regression beyond the margin)
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
import tarfile
import tempfile
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

SCHEMA_VERSION = 2

# ──────────────────────── bundled fixture (self-contained) ────────────────────
# The perf suite ships a SMALL, REAL subset of GRCh38 (~15 Mbp; chr20/21/22
# slices + one pure-ACGT unplaced contig) so it runs on a fresh checkout with no
# external data. The fixture is resolved at runtime in this order:
#   1. already-extracted data dir (perf/data/ with a matching VERSION marker)
#   2. local tarball (perf/perfdata-v1.tar.gz) -> extract into perf/data/
#   3. download from PERF_DATA_URL (perf/.env or process env) -> extract
# After resolution, an Encoded store is built (and cached) under perf/data/_store
# for the extract/vrs tasks. Big-store deep runs can still override with env vars
# or the matching CLI flags (PERF_STORE / PERF_BED_DIR / PERF_POINTS / ...).
HERE = Path(__file__).resolve().parent          # the perf/ directory
FIXTURE_VERSION = "perfdata-v1"
TARBALL_NAME = f"{FIXTURE_VERSION}.tar.gz"
DEF_DATA_DIR = HERE / "data"                     # extracted fixture lives here
DEF_TARBALL = HERE / TARBALL_NAME                # local cached tarball
DEF_ENV_FILE = HERE / ".env"                     # optional PERF_DATA_URL source
DEF_STORE_DIR = DEF_DATA_DIR / "_store"          # built+cached Encoded store

DEF_WORK = "/tmp/gtars-perf-work"        # cached derived BEDs/points subsets
VRS_POINTS = 200_000                      # bounded subset of the variant file


def _read_env_file(path):
    """Parse a tiny KEY=VALUE .env file (stdlib only). Returns a dict."""
    env = {}
    try:
        for line in Path(path).read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            env[k.strip()] = v.strip().strip('"').strip("'")
    except OSError:
        pass
    return env


def _fixture_version(data_dir):
    """Return the `version=...` value from a fixture's VERSION marker, or None."""
    vf = Path(data_dir) / "VERSION"
    try:
        for line in vf.read_text().splitlines():
            if line.startswith("version="):
                return line.split("=", 1)[1].strip()
    except OSError:
        pass
    return None


def _extract_tarball(tarball, data_dir):
    data_dir = Path(data_dir)
    if data_dir.exists():
        shutil.rmtree(data_dir, ignore_errors=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    sys.stderr.write(f"[data] extracting {tarball} -> {data_dir}\n")
    with tarfile.open(tarball, "r:gz") as tar:
        tar.extractall(data_dir)


def ensure_perf_data(data_dir=DEF_DATA_DIR, tarball=DEF_TARBALL,
                     env_file=DEF_ENV_FILE):
    """Resolve the bundled perf fixture, returning the extracted data dir.

    Resolution order:
      1. extracted dir already present with a matching VERSION  -> use it
      2. local tarball present                                  -> extract -> use
      3. PERF_DATA_URL set (env file or process env)            -> download -> use
      4. otherwise: a clear error explaining how to obtain it.
    """
    data_dir = Path(data_dir)
    tarball = Path(tarball)

    # 1) already extracted & not stale
    if _fixture_version(data_dir) == FIXTURE_VERSION:
        return data_dir
    if data_dir.exists() and _fixture_version(data_dir) not in (None, FIXTURE_VERSION):
        sys.stderr.write(
            f"[data] {data_dir} is a stale fixture "
            f"({_fixture_version(data_dir)} != {FIXTURE_VERSION}); re-extracting\n")

    # 2) local tarball
    if tarball.is_file():
        _extract_tarball(tarball, data_dir)
        if _fixture_version(data_dir) == FIXTURE_VERSION:
            return data_dir
        sys.exit(f"ERROR: extracted tarball {tarball} but VERSION marker did not "
                 f"match {FIXTURE_VERSION}.")

    # 3) download from PERF_DATA_URL
    url = os.environ.get("PERF_DATA_URL") or _read_env_file(env_file).get("PERF_DATA_URL")
    if url:
        sys.stderr.write(f"[data] downloading fixture from {url}\n")
        tmp = tarball.with_suffix(tarball.suffix + ".part")
        try:
            with urllib.request.urlopen(url) as resp, open(tmp, "wb") as out:
                shutil.copyfileobj(resp, out)
            tmp.replace(tarball)
        except Exception as e:  # urllib raises a zoo of exceptions
            try:
                tmp.unlink()
            except OSError:
                pass
            sys.exit(f"ERROR: failed to download fixture from {url}: {e}")
        _extract_tarball(tarball, data_dir)
        if _fixture_version(data_dir) == FIXTURE_VERSION:
            return data_dir
        sys.exit(f"ERROR: downloaded tarball but VERSION marker did not match "
                 f"{FIXTURE_VERSION}.")

    # 4) nothing worked
    sys.exit(
        "ERROR: no perf fixture available.\n"
        f"  Looked for extracted dir : {data_dir}\n"
        f"  Looked for local tarball : {tarball}\n"
        "  PERF_DATA_URL is not set (checked process env and "
        f"{env_file}).\n\n"
        "To fix, do ONE of:\n"
        f"  - place the fixture tarball at {tarball}, or\n"
        f"  - set PERF_DATA_URL=<url-to-{TARBALL_NAME}> (in {env_file} or the "
        "environment).\n"
        f"See perf/.env.example and perf/README.md for details.")


def build_fixture_store(data_dir, gtars_bin, store_dir=DEF_STORE_DIR):
    """Build (and cache) an Encoded RefgetStore from the fixture FASTAs.

    Rebuilds only if the store is missing or stale (older than the fixture's
    VERSION marker). Returns the store path."""
    data_dir = Path(data_dir)
    store_dir = Path(store_dir)
    fasta_dir = data_dir / "fasta"
    fastas = sorted(str(p) for p in fasta_dir.glob("*.fa.gz"))
    if not fastas:
        sys.exit(f"ERROR: no FASTA files found under {fasta_dir}")

    manifest = store_dir / "rgstore.json"
    version_marker = data_dir / "VERSION"
    fresh = (manifest.is_file() and version_marker.is_file()
             and manifest.stat().st_mtime >= version_marker.stat().st_mtime)
    if fresh:
        return store_dir

    sys.stderr.write(f"[data] building Encoded fixture store -> {store_dir}\n")
    shutil.rmtree(store_dir, ignore_errors=True)
    store_dir.parent.mkdir(parents=True, exist_ok=True)
    argv = [str(gtars_bin), "refget", "build", *fastas,
            "-o", str(store_dir), "-j", "0"]
    proc = subprocess.run(argv, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"ERROR: failed to build fixture store:\n{proc.stdout}{proc.stderr}")
    return store_dir


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
        "profile": "release",  # we always drive release binaries
    }


# ──────────────────────────── /usr/bin/time -v driver ─────────────────────────

_RSS_RE = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")


def run_timed(argv):
    """Run `argv` under /usr/bin/time -v. Returns (seconds, peak_rss_mb, output).

    `output` is the child's combined stdout+stderr. peak_rss_mb is None if
    /usr/bin/time is unavailable/unparseable. Raises RuntimeError on nonzero exit.
    """
    have_time = shutil.which("/usr/bin/time") is not None
    with tempfile.NamedTemporaryFile("r", suffix=".time", delete=False) as tf:
        time_log = tf.name
    try:
        wrapped = (["/usr/bin/time", "-v", "-o", time_log] + argv) if have_time else argv
        start = time.perf_counter()
        proc = subprocess.run(wrapped, capture_output=True, text=True)
        seconds = time.perf_counter() - start
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


# ───────────────────────────── cell identity helper ──────────────────────────
# Every measured result is one "cell" keyed by task/scenario/path/concurrency.
# This lets the gate and the regression highlighter address each one uniquely.

def cell_key(task, scenario, path, conc):
    return f"{task}/{scenario}/{path}/{conc}"


def cell_label(task, scenario, path, conc):
    """Human-friendly label, e.g. 'extract/large_width (partial)' or 'encode/j8'."""
    if task == "encode":
        return f"encode/j{conc}"
    p = "" if path in ("-", "", None) else f" ({path})"
    return f"{task}/{scenario}{p}"


def result_record(task, scenario, path, conc, seconds, peak_rss_mb,
                  throughput, unit, extra):
    return {
        "task": task,
        "scenario": scenario,
        "path": path,
        "concurrency": conc,
        "seconds": seconds,
        "peak_rss_mb": peak_rss_mb,
        "throughput": throughput,
        "throughput_unit": unit,
        "extra": extra,
    }


# ─────────────────────────────── data preparation ────────────────────────────

def warm_cache(*paths):
    """Read the store sequence files once so timings are page-cache-stable."""
    for p in paths:
        if not p:
            continue
        seqdir = Path(p) / "sequences"
        if not seqdir.is_dir():
            continue
        try:
            subprocess.run(
                f"find {seqdir} -type f -exec cat {{}} + > /dev/null 2>&1",
                shell=True, check=False)
        except OSError:
            pass


def fixture_fastas(data_dir):
    """The bundled, pre-split, gzipped FASTA files (>=2 so the encode `-j` sweep
    is meaningful — `-j` parallelizes across FILES)."""
    fasta_dir = Path(data_dir) / "fasta"
    files = sorted(str(p) for p in fasta_dir.glob("*.fa.gz"))
    if not files:
        sys.exit(f"ERROR: no fixture FASTAs under {fasta_dir}")
    return files


def fixture_beds(data_dir, bed_dir=None):
    """Return {scenario: bed_path} from the bundled fixture (or a `bed_dir`
    override). Files: small.bed, large_count.bed, large_width.bed."""
    bdir = Path(bed_dir) if bed_dir else Path(data_dir) / "bed"
    beds = {
        "small": bdir / "small.bed",
        "large_count": bdir / "large_count.bed",
        "large_width": bdir / "large_width.bed",
    }
    for scen, p in beds.items():
        if not p.is_file():
            sys.exit(f"ERROR: missing fixture BED for '{scen}': {p}")
    return {k: str(v) for k, v in beds.items()}


def prepare_points(src, work_dir, n):
    """First `n` lines of the variant points file, cached under work_dir.

    The cache key includes a hash of the source path so switching between the
    bundled fixture and a big real points file never reuses a stale subset."""
    import hashlib
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    tag = hashlib.sha1(str(Path(src).resolve()).encode()).hexdigest()[:8]
    dest = work / f"points_{n}_{tag}.tsv"
    if dest.is_file():
        return str(dest)
    if not Path(src).is_file():
        sys.exit(f"ERROR: variant points file not found: {src}")
    lines = []
    with open(src) as fh:
        for i, line in enumerate(fh):
            if i >= n:
                break
            lines.append(line.rstrip("\n"))
    dest.write_text("\n".join(lines) + "\n")
    return str(dest)


# ─────────────────────────────────── tasks ────────────────────────────────────

_DONE_RE = re.compile(r"Done:\s*([\d,]+)\s+sequences?,\s*([\d,]+)\s+bases", re.I)
_RESULT_RE = re.compile(
    r"RESULT\s+task=(\S+)\s+scenario=(\S+)\s+path=(\S+)\s+seconds=(\S+)\s+"
    r"items=(\S+)\s+bases=(\S+)\s+throughput=(\S+)\s+unit=(\S+)")


def _clear_rgsi_sidecars(fastas):
    """The importer caches `<fasta>.rgsi` sidecars; remove them so each encode
    run actually re-does the digest+encode work."""
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
    """encode: `gtars refget build <fastas> -o DIR -j jobs` under time -v."""
    _clear_rgsi_sidecars(fastas)
    out_dir = Path(tempfile.mkdtemp(prefix=f"gtars-perf-encode-j{jobs}-"))
    shutil.rmtree(out_dir, ignore_errors=True)
    argv = [str(gtars_bin), "refget", "build", *[str(f) for f in fastas],
            "-o", str(out_dir), "-j", str(jobs)]
    try:
        seconds, peak_rss_mb, output = run_timed(argv)
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
        return result_record(
            "encode", "multifile", "-", jobs, seconds, peak_rss_mb,
            throughput, "bases_per_sec",
            {"dataset_id": dataset_id, "n_bases": n_bases,
             "n_sequences": n_sequences, "n_files": len(fastas),
             "collection_digest": collection_digest})
    finally:
        shutil.rmtree(out_dir, ignore_errors=True)


def task_perf(perf_bin, task, scenario, path, enc, raw, bed, points, dataset_id):
    """extract/vrs: run the focused `perf_task` example, parse its RESULT line."""
    argv = [str(perf_bin), "--task", task, "--scenario", scenario,
            "--path", path, "--store", str(enc)]
    if raw:
        argv += ["--raw", str(raw)]
    if bed:
        argv += ["--bed", str(bed)]
    if points:
        argv += ["--points", str(points)]
    seconds, peak_rss_mb, output = run_timed(argv)
    m = _RESULT_RE.search(output)
    if not m:
        raise RuntimeError(f"perf_task printed no RESULT line:\n{output}")
    _t, _sc, _p, inner_secs, items, bases, throughput, unit = m.groups()
    extra = {
        "dataset_id": dataset_id,
        "items": int(items),
        "bases": int(bases),
        # `seconds` (the gated number) is the in-process sweep time, which
        # excludes process startup + any load_sequence cost; `wall_seconds` is
        # the full /usr/bin/time wall clock (kept for context / RSS provenance).
        "wall_seconds": seconds,
        "inner_seconds": float(inner_secs),
    }
    return result_record(task, scenario, path, 1, float(inner_secs),
                        peak_rss_mb, float(throughput), unit, extra)


# ──────────────────────────────────── run ─────────────────────────────────────

def cmd_run(a):
    suite_start = time.perf_counter()
    gtars_bin = a.gtars_bin
    perf_bin = a.perf_bin
    if not Path(gtars_bin).is_file():
        sys.exit(f"ERROR: gtars CLI not found: {gtars_bin} "
                 "(cargo build --release -p gtars-cli --all-features)")
    if not Path(perf_bin).is_file():
        sys.exit(f"ERROR: perf_task example not found: {perf_bin} "
                 "(cargo build --release --example perf_task -p gtars-refget)")

    # Resolve the bundled fixture (dir -> local tarball -> PERF_DATA_URL).
    data_dir = ensure_perf_data()

    # Resolve the store: env/CLI override wins (lets a big real store drive deep
    # runs); otherwise build+cache the small Encoded store from the fixture.
    enc_override = os.environ.get("PERF_STORE") or a.enc
    if enc_override:
        enc = enc_override
        if not Path(enc).is_dir():
            sys.exit(f"ERROR: PERF_STORE/--enc store not found: {enc}")
    else:
        enc = str(build_fixture_store(data_dir, gtars_bin))
    raw = os.environ.get("PERF_RAW") or a.raw  # optional; reserved/warmed

    # Inputs: fixture by default, overridable via env for deep runs.
    fastas = fixture_fastas(data_dir)
    bed_dir = os.environ.get("PERF_BED_DIR")
    beds = fixture_beds(data_dir, bed_dir)
    points_src = os.environ.get("PERF_POINTS") or a.points_src \
        or str(Path(data_dir) / "variants" / "points.tsv")

    # Warm the page cache once so numbers are stable.
    sys.stderr.write("[prep] warming page cache for stores...\n")
    warm_cache(enc, raw)

    points = prepare_points(points_src, a.work, a.vrs_points)

    results = []

    # 1) encode: single core (j1) and multi core (jN).
    ncpu = a.encode_jobs or min(os.cpu_count() or 1, 8)
    for j in (1, ncpu):
        sys.stderr.write(f"[encode] jobs={j}\n")
        results.append(task_encode(gtars_bin, fastas, j, a.dataset_id))

    # 2) extract: 3 scenarios x 3 paths.
    #    batch: stub-only, group ranges by sequence, call get_substrings() once
    #    per sequence (adaptive whole-decode + unchecked-UTF-8 levers).
    for scenario in ("small", "large_count", "large_width"):
        for path in ("resident", "partial", "batch"):
            sys.stderr.write(f"[extract] {scenario} ({path})\n")
            results.append(task_perf(
                perf_bin, "extract", scenario, path, enc, raw,
                beds[scenario], None, a.dataset_id))

    # 3) vrs: point lookups, resident path.
    sys.stderr.write("[vrs] points (resident)\n")
    results.append(task_perf(
        perf_bin, "vrs", "points", "resident", enc, None, None,
        points, a.dataset_id))

    suite_seconds = time.perf_counter() - suite_start

    record = {
        "schema_version": SCHEMA_VERSION,
        "run": collect_run_meta(a.repo),
        "suite_seconds": suite_seconds,
        "results": results,
    }

    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = record["run"]["timestamp_utc"].replace(":", "").replace("-", "")
    json_path = out_dir / f"{record['run']['gtars_commit']}-{ts}.json"
    json_path.write_text(json.dumps(record, indent=2))
    sys.stderr.write(f"\nWrote run record: {json_path}\n")

    _print_summary(results, suite_seconds)
    print(json_path)


def _fmt_throughput(r):
    unit = r["throughput_unit"]
    v = r["throughput"]
    if unit == "bases_per_sec":
        return f"{v / 1e6:>8.1f} Mbase/s"
    if unit == "lookups_per_sec":
        return f"{v / 1e6:>8.2f} Mlookup/s"
    return f"{v:>10.0f} {unit}"


def _print_summary(results, suite_seconds):
    out = sys.stderr
    out.write("\n")
    out.write(f"{'cell':<34} {'seconds':>9} {'rss_mb':>8} {'throughput':>18}\n")
    out.write("-" * 72 + "\n")
    for r in results:
        label = cell_label(r["task"], r["scenario"], r["path"], r["concurrency"])
        rss = f"{r['peak_rss_mb']:.0f}" if r["peak_rss_mb"] is not None else "null"
        out.write(f"{label:<34} {r['seconds']:>9.3f} {rss:>8} "
                  f"{_fmt_throughput(r):>18}\n")
    out.write("-" * 72 + "\n")
    out.write(f"TOTAL SUITE WALL TIME: {suite_seconds:.2f}s "
              f"({'OK <30s' if suite_seconds < 30 else 'OVER BUDGET >=30s'})\n")


# ──────────────────────────────── compare / gate ──────────────────────────────

class Color:
    def __init__(self, enabled):
        self.on = enabled
    def _w(self, code, s):
        return f"\033[{code}m{s}\033[0m" if self.on else s
    def red(self, s): return self._w("31;1", s)
    def green(self, s): return self._w("32", s)
    def yellow(self, s): return self._w("33", s)
    def cyan(self, s): return self._w("36", s)
    def bold(self, s): return self._w("1", s)


def _index_results(run):
    idx = {}
    for r in run.get("results", []):
        idx[cell_key(r["task"], r["scenario"], r["path"], r["concurrency"])] = r
    return idx


def _check_cell(r, c, margin):
    """Return (ok, [violation strings]). Same semantics as before, per-metric."""
    msgs, ok = [], True
    if c.get("min_throughput") is not None and c["min_throughput"] > 0:
        floor = c["min_throughput"] * (1.0 - margin)
        if r["throughput"] < floor:
            ok = False
            msgs.append(f"throughput {r['throughput']:.0f} < {floor:.0f} "
                        f"(min*{1.0 - margin:.2f})")
    if c.get("max_seconds") is not None:
        ceil = c["max_seconds"] * (1.0 + margin)
        if r["seconds"] > ceil:
            ok = False
            msgs.append(f"seconds {r['seconds']:.3f} > {ceil:.3f}")
    if c.get("max_peak_rss_mb") is not None and r["peak_rss_mb"] is not None:
        ceil = c["max_peak_rss_mb"] * (1.0 + margin)
        if r["peak_rss_mb"] > ceil:
            ok = False
            msgs.append(f"peak_rss {r['peak_rss_mb']:.0f} > {ceil:.0f}")
    return (ok, msgs)


def _pct(new, old):
    if old == 0:
        return 0.0
    return (new - old) / old * 100.0


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
    col = Color(sys.stdout.isatty() and not a.no_color)
    run = _load_json(a.run, "run")
    if run.get("schema_version") != SCHEMA_VERSION:
        sys.stderr.write(
            f"ERROR: run schema_version {run.get('schema_version')} != harness "
            f"schema_version {SCHEMA_VERSION}; re-run / re-seed.\n")
        sys.exit(2)
    targets = _load_json(a.targets, "targets")
    if "cells" not in targets or "margin" not in targets:
        sys.stderr.write("ERROR: targets JSON missing 'margin'/'cells'.\n")
        sys.exit(2)

    seeded_host = targets.get("seeded_on_host", "")
    run_host = run.get("run", {}).get("host", "")
    if seeded_host and seeded_host != run_host:
        sys.stderr.write(
            f"WARNING: run host '{run_host}' differs from targets host "
            f"'{seeded_host}'. Targets are machine-specific; re-seed here.\n")

    margin = targets["margin"]
    cells = targets["cells"]

    print(col.bold(f"{'cell':<34} {'seconds':>9} {'baseline':>9} "
                   f"{'Δ%':>7} {'result':>8}"))
    print("-" * 72)

    regressions = []     # (label, list of plain-language lines)
    improvements = []
    any_fail = any_gated = False

    for r in run.get("results", []):
        key = cell_key(r["task"], r["scenario"], r["path"], r["concurrency"])
        label = cell_label(r["task"], r["scenario"], r["path"], r["concurrency"])
        c = cells.get(key)
        if c is None:
            print(f"{label:<34} {r['seconds']:>9.3f} {'-':>9} {'-':>7} "
                  f"{col.yellow('WARN'):>8}  (no target)")
            continue
        any_gated = True
        ok, viol = _check_cell(r, c, margin)
        any_fail = any_fail or not ok

        base_secs = c.get("max_seconds")
        dpct = _pct(r["seconds"], base_secs) if base_secs else 0.0
        base_str = f"{base_secs:.3f}" if base_secs else "-"
        res = col.green("PASS") if ok else col.red("FAIL")
        print(f"{label:<34} {r['seconds']:>9.3f} {base_str:>9} "
              f"{dpct:>+6.1f}% {res:>8}")

        # Plain-language regression / improvement narration (per the seconds
        # metric, the primary speed signal; throughput/RSS violations append).
        if base_secs:
            if not ok and r["seconds"] > base_secs:
                # SLOWER
                slower = _pct(r["seconds"], base_secs)
                regressions.append(
                    (label,
                     f"{label}: {slower:.0f}% SLOWER "
                     f"({base_secs:.2f}s -> {r['seconds']:.2f}s)", viol))
            elif r["seconds"] < base_secs * 0.8:
                # noticeably FASTER (>20%): note as an improvement
                ratio = base_secs / r["seconds"] if r["seconds"] > 0 else 0
                improvements.append(
                    f"{label}: {ratio:.1f}x FASTER "
                    f"({base_secs:.2f}s -> {r['seconds']:.2f}s)")
        # throughput-only regressions (e.g. encode bases/s) not captured above
        if not ok and base_secs and r["seconds"] <= base_secs and viol:
            regressions.append((label, f"{label}: constraint regressed", viol))

    print("-" * 72)

    if regressions:
        print(col.red(col.bold("\n*** REGRESSIONS ***")))
        for label, line, viol in regressions:
            print("  " + col.red("REGRESSION: ") + col.red(line))
            for v in viol:
                print("      - " + v)
    if improvements:
        print(col.green(col.bold("\n*** IMPROVEMENTS ***")))
        for line in improvements:
            print("  " + col.green("FASTER: " + line))

    if not any_gated:
        sys.stderr.write("No gated cells matched any target. "
                         "Treating as pass (informational only).\n")
    if any_fail:
        print(col.red(col.bold("\nGATE: FAIL — at least one cell regressed "
                               "beyond the margin.")))
        sys.exit(1)
    print(col.green(col.bold("\nGATE: PASS")))


def cmd_seed(a):
    run = _load_json(a.run, "run")
    cells = {}
    for r in run.get("results", []):
        key = cell_key(r["task"], r["scenario"], r["path"], r["concurrency"])
        cells[key] = {
            "label": cell_label(r["task"], r["scenario"], r["path"],
                                r["concurrency"]),
            "min_throughput": r["throughput"],
            "max_seconds": r["seconds"],
            **({"max_peak_rss_mb": r["peak_rss_mb"]}
               if r["peak_rss_mb"] is not None else {}),
        }
    targets = {
        "schema_version": SCHEMA_VERSION,
        "margin": a.margin,
        "seeded_on_host": run.get("run", {}).get("host", ""),
        "cells": cells,
    }
    out = Path(a.out)
    if out.parent:
        out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(targets, indent=2))
    sys.stderr.write(
        f"Seeded {len(cells)} target cells from {a.run} on host "
        f"'{targets['seeded_on_host']}' (margin {a.margin}) -> {out}\n")


# ──────────────────────────────────── cli ─────────────────────────────────────

def build_parser():
    p = argparse.ArgumentParser(
        prog="perf.py", description="gtars refgetstore performance regression suite")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run the whole suite and emit a run-record JSON")
    r.add_argument("--enc", default=None,
                   help="ENCODED refget store dir for extract/vrs (default: build "
                        "from the bundled fixture; env PERF_STORE also honored)")
    r.add_argument("--raw", default=None,
                   help="RAW refget store dir (warmed; reserved; env PERF_RAW)")
    r.add_argument("--points-src", default=None,
                   help="variant points file (chrom<TAB>pos) for vrs (default: "
                        "fixture variants/points.tsv; env PERF_POINTS)")
    r.add_argument("--vrs-points", type=int, default=VRS_POINTS,
                   help=f"bound the vrs lookup count (default {VRS_POINTS})")
    r.add_argument("--encode-jobs", type=int, default=0,
                   help="multi-core jobs value (0 = min(cpu_count, 8))")
    r.add_argument("--work", default=DEF_WORK,
                   help="cache dir for the bounded points subset")
    r.add_argument("--out", default="perf/runs", help="output dir for run JSON")
    r.add_argument("--gtars-bin", default="target/release/gtars",
                   help="path to the gtars CLI binary")
    r.add_argument("--perf-bin",
                   default="target/release/examples/perf_task",
                   help="path to the perf_task example binary")
    r.add_argument("--repo", default=".", help="repo root for commit metadata")
    r.add_argument("--dataset-id", default=FIXTURE_VERSION,
                   help="short label stamped into every result's extra")
    r.set_defaults(func=cmd_run)

    c = sub.add_parser("compare",
                       help="gate a run JSON vs targets + highlight regressions")
    c.add_argument("run", help="run JSON to gate")
    c.add_argument("--targets", default="perf/targets.json", help="targets JSON")
    c.add_argument("--no-color", action="store_true", help="disable ANSI color")
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
