# gtars-perf

An occasionally-run performance regression harness for `gtars`, implemented as
a single readable Python script (`perf.py`, stdlib only). It replaces the old
`gtars-perf/` Rust crate while keeping the **same run-record JSON schema
(`schema_version: 1`)** and the **same compare/seed gate design**.

The script drives the real `gtars` CLI as a subprocess under `/usr/bin/time -v`,
so it captures the child's **true peak RSS** for every task it measures. (The
old crate ran extract/VRS in-process and therefore reported `null` RSS for them.)

## Requirements

- Python 3 (stdlib only — no pip installs).
- `/usr/bin/time` (GNU time, for `-v` peak RSS). If absent, runs still work but
  `peak_rss_mb` is `null`.
- A built `gtars` CLI:
  ```
  cargo build --release -p gtars-cli --all-features
  # binary: target/release/gtars
  ```

## Usage

Run from the repo root.

```bash
# Run the encode jobs-sweep over a multi-file FASTA set:
python3 perf/perf.py run \
    --fasta tests/data/fasta/a.fa --fasta tests/data/fasta/b.fa \
    --concurrency 1,4,8 \
    --gtars-bin target/release/gtars \
    --dataset-id mydataset \
    --out perf/runs
# (--fasta is repeatable and/or comma-separated.) Prints the run JSON path on stdout.

# Seed a targets file from a known-good baseline run:
python3 perf/perf.py seed-targets perf/runs/<commit>-<ts>.json \
    --out perf/targets.json --margin 0.15

# Gate a new run against targets (exit 0 pass / 1 regression / 2 bad input):
python3 perf/perf.py compare perf/runs/<commit>-<ts>.json --targets perf/targets.json
echo "exit: $?"
```

## Tasks

| task    | how measured | concurrency knob | real peak RSS? |
|---------|--------------|------------------|----------------|
| encode  | `gtars refget build <fastas> -o DIR -j N` | `-j` (FASTA files in parallel) | yes |
| extract | *not measurable* | — | n/a (skipped) |
| vrs     | *not measurable* | — | n/a (skipped) |

The `-j` knob parallelizes across **files**, so use **>=2 FASTA files** for the
encode sweep to mean anything (a single file gives a flat sweep; the harness
warns).

### CLI gap: extract and vrs

As of this writing the `gtars` CLI exposes **only** `refget build`. There is no
subcommand for region/substring **extract** from a refget store, nor for
**VCF -> VRS** id computation. The old crate measured those by calling the
library functions in-process (`SubstringsFromRegions`,
`compute_vrs_ids_parallel_encoded`), which is exactly why it could not report
per-task peak RSS for them.

To keep this harness honest it does **not** fake those measurements. If you pass
`--bed`/`--vcf`, the run records an explicit skipped row
(`extra.skipped == true`, with a `reason`) instead of a fabricated number.
`seed-targets` ignores skipped rows, and `compare` treats them as informational
unless a target exists for them.

To measure extract/VRS with real peak RSS, add CLI subcommands that wrap the
existing library calls, e.g.:

- `gtars refget extract --store DIR --bed FILE` (substrings from regions)
- `gtars vrs --store DIR --vcf FILE --threads N` (VRS allele ids; the thread
  knob is what lets us do a 1/4/8 concurrency sweep)

Then add a task function in `perf.py` (see "Adding a task") that shells out to it
under `run_timed(...)`.

## Run-record schema (`schema_version: 1`)

```jsonc
{
  "schema_version": 1,
  "run": {
    "timestamp_utc":   "2026-05-29T18:58:01Z",
    "gtars_commit":    "1a1a694",
    "gtars_commit_date": "2026-05-29T...",
    "gtars_dirty":     false,
    "host":            "machine",
    "cpu_model":       "...",
    "logical_cpus":    16,
    "total_ram_mb":    64000,
    "rustc_version":   "rustc 1.x.y ...",
    "profile":         "release"
  },
  "results": [
    {
      "task":            "encode",
      "concurrency":     4,          // -j for encode; threads for vrs; 1 for extract
      "seconds":         1.234,      // wall clock measured around the child
      "peak_rss_mb":     321.0,      // from /usr/bin/time -v; null if unmeasured
      "throughput":      1234567.0,
      "throughput_unit": "bases_per_sec",
      "extra": {                     // task-specific; includes correctness checks
        "dataset_id": "mydataset",
        "n_bases": 16, "n_sequences": 3, "n_files": 2,
        "collection_digest": "3e0213b5..."   // store collections_digest (regression-visible)
      }
    }
  ]
}
```

Correctness signals live in `extra`: encode records the store
`collection_digest` plus `n_bases`/`n_sequences`, so a change in *output*
(not just speed) is visible across runs.

## Gate / targets.json

```jsonc
{
  "margin": 0.15,                    // tolerance applied to every constraint
  "seeded_on_host": "machine",       // compare warns if run host differs
  "tasks": {
    "encode": {
      "1": { "min_throughput": 1.0e6, "max_seconds": 2.0, "max_peak_rss_mb": 400.0 },
      "4": { "min_throughput": 3.0e6, "max_seconds": 1.0, "max_peak_rss_mb": 450.0 }
    }
  }
}
```

`compare` applies `margin` symmetrically:
- fails if `throughput < min_throughput * (1 - margin)`
- fails if `seconds     > max_seconds     * (1 + margin)`
- fails if `peak_rss_mb > max_peak_rss_mb * (1 + margin)` (only when both present)

Cells with no matching target are informational (`WARN`, not a failure). Exit
codes: **0** pass, **1** regression, **2** bad input (unreadable/unparseable
JSON, schema mismatch, malformed targets).

Targets are machine-specific — commit `targets.json`, but re-seed it on the
machine you gate on.

## Adding a task

1. Write a `task_<name>(...)` function that builds an `argv` for the `gtars`
   CLI and calls `run_timed(argv)` -> `(seconds, peak_rss_mb, stdout)`.
2. Parse correctness/size signals from `stdout` (or a store manifest) into
   `extra`, compute `throughput`, and return a result dict with the standard
   keys (`task, concurrency, seconds, peak_rss_mb, throughput, throughput_unit,
   extra`).
3. Call it from `cmd_run` (sweep over `a.concurrency` if it has a thread knob).
4. Re-seed `targets.json` so the new cells get a floor/ceiling.
