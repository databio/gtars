# gtars-perf

A fast, **refgetstore-only** performance-regression suite for `gtars`,
implemented as a single readable Python script (`perf.py`, stdlib only) plus a
focused Rust example binary (`gtars-refget/examples/perf_task.rs`).

The whole suite runs **end-to-end in well under 30 seconds** and measures
**only refgetstore** (no competitor tools) across the three real tasks of the
refget benchmark: **encode**, **extract**, and **vrs**. It captures true wall
time and **peak RSS** for every cell (via `/usr/bin/time -v`), and on `compare`
it **highlights regressions in plain language**, e.g.

```
REGRESSION: extract/large_width (partial): 22% SLOWER (1.08s -> 1.32s)
```

## Requirements

- Python 3 (stdlib only — no pip installs).
- `/usr/bin/time` (GNU time, for `-v` peak RSS). If absent, runs still work but
  `peak_rss_mb` is `null`.
- Built binaries (release):
  ```
  cargo build --release -p gtars-cli --all-features          # target/release/gtars
  cargo build --release --example perf_task -p gtars-refget  # target/release/examples/perf_task
  ```
- The bundled fixture dataset (see [Dataset](#dataset-self-contained) below). On
  a fresh checkout this is resolved automatically from a local tarball or a
  configured download URL.

## Usage

Run from the repo root. The suite is **self-contained** — it ships its own small
real dataset — so `run` is a one-liner:

```bash
# Run the whole suite (encode + extract + vrs); prints the run JSON path:
python3 perf/perf.py run --out perf/runs

# Seed a targets file from a known-good baseline run:
python3 perf/perf.py seed-targets perf/runs/<commit>-<ts>.json \
    --out perf/targets.json --margin 0.15

# Gate a new run vs targets AND highlight regressions
# (exit 0 pass / 1 regression / 2 bad input):
python3 perf/perf.py compare perf/runs/<commit>-<ts>.json --targets perf/targets.json
echo "exit: $?"
```

## Dataset (self-contained)

The suite no longer depends on a prebuilt multi-GB store or machine-specific
files. Instead it ships a **small, real subset of GRCh38** (~15 Mbp; a ~5 MB
gzipped tarball, `perf/perfdata-v1.tar.gz`) and resolves it at runtime.

**Tarball contents** (source data only — the store is built+cached at runtime,
not shipped):

| path                    | what                                                        |
|-------------------------|-------------------------------------------------------------|
| `fasta/part1.fa.gz`     | chr20 slice (first 5 Mbp; has N runs → `dna3bit`)           |
| `fasta/part2.fa.gz`     | chr21 slice (5 Mbp; has N runs → `dna3bit`)                 |
| `fasta/part3.fa.gz`     | chr22 slice (5 Mbp) + `GL000226.1` (pure ACGT → `dna2bit`)  |
| `bed/small.bed`         | ~10k narrow regions                                         |
| `bed/large_count.bed`   | ~100k narrow regions                                        |
| `bed/large_width.bed`   | ~2k regions 100kb–1Mb (sized to fit the small genome)       |
| `variants/points.tsv`   | ~200k `chrom<TAB>pos` (clinvar chr20/21/22 + synthetic fill)|
| `VERSION`               | fixture version marker (`version=perfdata-v1`)              |

Three separate FASTA files keep the encode `-j` single-vs-multi-core comparison
meaningful (`-j` parallelizes across **files**). Both store alphabets
(`dna2bit`, `dna3bit`) are present so encode/extract exercise both decoders.

### Resolution order

`run` calls `ensure_perf_data()`, which resolves the fixture in this order:

1. **Extracted dir** — `perf/data/` already present with a matching `VERSION`
   marker → use it.
2. **Local tarball** — `perf/perfdata-v1.tar.gz` present → extract into
   `perf/data/` → use it.
3. **Download** — `PERF_DATA_URL` set (read from `perf/.env`, else the process
   environment) → download (stdlib `urllib`) to `perf/perfdata-v1.tar.gz` →
   extract → use it.
4. Otherwise: a clear error explaining how to obtain the data.

After resolution, an **Encoded** store is built once from the bundled FASTAs and
cached under `perf/data/_store/` (rebuilt only if missing/stale) for the
extract/vrs tasks.

`perf/data/`, `perf/*.tar.gz`, and `perf/.env` are all gitignored.

### Configuring the download URL

The tarball is too large to commit, so upload it to cloud storage and point the
suite at it:

```bash
cp perf/.env.example perf/.env
# edit perf/.env and set PERF_DATA_URL=https://<bucket>/perfdata-v1.tar.gz
```

On a machine that already has the tarball or an extracted `perf/data/`, no URL
is needed.

### Overrides (deep runs against the big real store)

The bundled fixture is the default, but you can point any task at full-size real
data via env vars (or the matching `run` flags):

| env var         | overrides                                              |
|-----------------|--------------------------------------------------------|
| `PERF_STORE`    | the Encoded store dir for extract/vrs (`--enc`)        |
| `PERF_RAW`      | the Raw store dir (warmed; reserved) (`--raw`)         |
| `PERF_BED_DIR`  | a dir holding `small.bed` / `large_count.bed` / `large_width.bed` |
| `PERF_POINTS`   | a `chrom<TAB>pos` variant file (`--points-src`)        |
| `PERF_DATA_URL` | where to download the fixture tarball                  |

The bounded points subset is cached under `--work` (`/tmp/gtars-perf-work` by
default; cache key includes the source path so switching sources never reuses a
stale subset).

## Tasks, scenarios, and paths

Each measured **cell** is keyed by `task / scenario / path / concurrency`.

| task    | scenario(s)                          | path(s)              | knob | what's measured            |
|---------|--------------------------------------|----------------------|------|----------------------------|
| encode  | `multifile`                          | `-`                  | `-j` | wall, peak RSS, bases/s    |
| extract | `small`, `large_count`, `large_width`| `resident`, `partial`, `batch`| —    | wall, peak RSS, bases/s    |
| vrs     | `points`                             | `resident`           | —    | wall, lookups/s            |

> **Caveats on the `batch` path and fixture size.** `batch` (`get_substrings`,
> one call per sequence) exists mainly to amortize **FFI** overhead for callers
> like the Python binding on high-region-count workloads; in pure-Rust it is at
> best a wash and can be *slower* on `large_width` because it materializes a
> whole sequence's output at once. The benchmark **contestant uses the per-region
> partial path**, not `batch`. Also note the bundled fixture's sequences are only
> ~5 Mbp, so `large_width` here understates real chromosome-scale costs — for
> decode/whole-vs-partial decisions, spot-check against a full-size store (set
> `PERF_STORE`) rather than trusting the fixture alone.

### encode

`gtars refget build <fasta>... -o DIR -j N`. The `-j` knob parallelizes across
**FILES**, so the suite builds the three bundled FASTA files (`fasta/part*.fa.gz`
— real chr20/chr21/chr22 slices plus a pure-ACGT contig) and runs the build at
**jobs=1 (single core)** and **jobs=N (multi core, `min(cpu_count, 8)`)**. The
multi-file fixture keeps the 1-core-vs-N-core comparison meaningful while staying
fast. Correctness is the store's `collections_digest` plus the parsed
`n_bases`/`n_sequences`.

### extract

Region/substring extraction over the **cached Encoded fixture store**
(`perf/data/_store/`, built once from the bundled FASTAs). Three scenarios mirror
the real benchmark, all against the small fixture genome:

- **small**       — ~10k narrow regions
- **large_count** — ~100k narrow regions
- **large_width** — ~2k regions of 100kb–1Mb (full-width slices sized to fit the
  ~5 Mbp fixture chromosomes)

Each scenario is measured on **both read paths**, because they have different
perf profiles and a regression in either must be caught:

- **resident** — `load_sequence` the whole touched sequences, then
  `get_substring()` slices/decodes each span.
- **partial**  — never `load_sequence`; `get_substring()` on a stub reads only
  the bytes covering each queried region straight from the `.seq` file.

### vrs

The VRS/HGVS point-lookup pattern: many **1bp lookups** on the resident path,
over a **bounded** variant subset (~200k of the bundled
`variants/points.tsv`, `chrom<TAB>pos` — real clinvar positions on chr20/21/22
plus synthetic fill to reach the count). Records lookups/s.

## Measurement entrypoint

- **encode** is driven by the real `gtars` CLI (`refget build`) and parsed from
  its `Done:` summary line.
- **extract** and **vrs** are driven by `examples/perf_task`, which runs exactly
  **one** task+scenario+path and prints a single machine-readable line:

  ```
  RESULT task=extract scenario=large_width path=partial seconds=1.316 \
         items=2000 bases=1112931389 throughput=845400000.0 unit=bases_per_sec
  ```

  `perf.py` invokes it under `/usr/bin/time -v` (so peak RSS is captured too)
  and parses the `RESULT` line. The gated `seconds` is the in-process sweep
  time (excludes process startup and, for resident, the one-time
  `load_sequence`); the full wall clock is kept in `extra.wall_seconds`.

Before measuring, the suite **warms the page cache** once (reads the encoded and
raw store sequence files to `/dev/null`) so numbers are stable.

## Run-record schema (`schema_version: 2`)

```jsonc
{
  "schema_version": 2,
  "run": { "timestamp_utc": "...", "gtars_commit": "...", "host": "...",
           "cpu_model": "...", "logical_cpus": 16, "rustc_version": "...",
           "profile": "release", ... },
  "suite_seconds": 6.56,            // total wall time of the whole suite
  "results": [
    {
      "task": "extract", "scenario": "large_width", "path": "partial",
      "concurrency": 1,
      "seconds": 1.316,             // gated number (in-process sweep)
      "peak_rss_mb": 5.0,           // from /usr/bin/time -v; null if unmeasured
      "throughput": 845400000.0,
      "throughput_unit": "bases_per_sec",
      "extra": { "items": 2000, "bases": 1112931389, "wall_seconds": 1.33, ... }
    }
    // encode rows additionally carry n_bases/n_sequences/collection_digest
  ]
}
```

## Gate / targets.json (`schema_version: 2`)

`seed-targets` writes one entry per cell keyed by `task/scenario/path/conc`:

```jsonc
{
  "schema_version": 2,
  "margin": 0.15,                    // tolerance applied to every constraint
  "seeded_on_host": "machine",       // compare warns if run host differs
  "cells": {
    "extract/large_width/partial/1": {
      "label": "extract/large_width (partial)",
      "min_throughput": 845400000.0, "max_seconds": 1.316, "max_peak_rss_mb": 5.0
    }
    // ...one per measured cell
  }
}
```

`compare` applies `margin` symmetrically (fails if `seconds > max_seconds*(1+m)`,
`throughput < min_throughput*(1-m)`, or `peak_rss_mb > max_peak_rss_mb*(1+m)`),
prints a per-cell table with the per-cell **Δ%** vs baseline, and then prints a
plain-language summary:

- **REGRESSIONS** — e.g. `REGRESSION: extract/large_width (partial): 22% SLOWER
  (1.08s -> 1.32s)`, with the specific violated constraints listed. ANSI red on
  a tty (use `--no-color` to disable).
- **IMPROVEMENTS** — e.g. `FASTER: extract/large_width (resident): 3.2x FASTER`.

Cells with no matching target are informational (`WARN`). Exit codes: **0**
pass, **1** regression, **2** bad input. Targets are machine-specific — commit
`targets.json`, but re-seed it on the machine you gate on.

## Adding a task/scenario/path

1. For a CLI-measurable task, add a `task_*` that builds an `argv` and calls
   `run_timed(...)`. For an in-process refgetstore path, add a branch to
   `examples/perf_task.rs` that emits a `RESULT` line and call it via
   `task_perf(...)`.
2. Wire it into `cmd_run` (mind the < 30s budget — shrink region/lookup counts
   if needed; draw inputs from the bundled fixture under `perf/data/`).
3. Re-seed `targets.json` so the new cell gets a floor/ceiling.
