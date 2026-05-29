# gtars-perf

An occasionally-run **performance regression harness** for gtars. It exercises
three tasks at several concurrency levels, emits a stable machine-readable run
record (wall time + peak RSS + throughput per cell), and ships a compare/gate
that exits nonzero when a run misses per-cell targets.

This is the lightweight gtars-internal regression gate. It is **not** run in CI
on every push; run it occasionally on a reference machine.
`refgetstore-benchmark` remains the heavyweight cross-format comparison.

## Tasks

| task    | what it does                                              | concurrency knob | throughput unit  |
|---------|-----------------------------------------------------------|------------------|------------------|
| encode  | `gtars refget build <fasta>... -o <dir> -j <N>`           | `-j/--jobs` = FASTA files imported **concurrently** | `bases_per_sec` |
| extract | iterate a BED over a store (`SubstringsFromRegions`)      | fixed at 1 (single-threaded in lib today) | `regions_per_sec` (`bytes_per_sec` in extra) |
| vrs     | `compute_vrs_ids_parallel_encoded` over a VCF             | worker threads   | `variants_per_sec` |

> **Encode `--jobs` reconciliation:** the FASTA-import parallelism is a single
> `--jobs` knob = number of FASTA **files** imported concurrently (the old
> per-sequence `--threads` knob was removed). The sweep therefore only matters
> over a **multi-file** input set; a single FASTA gives a flat sweep (the harness
> warns). The encode task is driven through the **real CLI binary** as a
> subprocess wrapped in `/usr/bin/time -v` for true peak RSS.

## Measuring peak RSS

Each (task, concurrency) cell runs in its **own process** so peak RSS is
per-task:
- library tasks (extract, vrs) re-invoke this binary's hidden `run-one` mode and
  read the child's `VmHWM` from `/proc/self/status`;
- encode reads the CLI child's `Maximum resident set size` from `/usr/bin/time -v`.

Pass `--no-subprocess` to run library cells in-process (peak RSS becomes
process-cumulative and is recorded as `null`; the gate skips that constraint).

## Usage

Build the CLI it drives, then build the harness:

```bash
cargo build --release -p gtars-cli --all-features   # provides target/release/gtars
cargo build --release -p gtars-perf
```

Run all tasks and write a run record:

```bash
cargo run --release -p gtars-perf -- run \
  --fasta a.fa.gz,b.fa.gz,c.fa.gz \
  --bed regions.bed \
  --vcf clinvar.vcf.gz \
  --concurrency 1,4,8 \
  --out perf/runs --csv \
  --gtars-bin target/release/gtars \
  --dataset-id grch38
```

(If `--store` is omitted, a setup store is built once from `--fasta` at `-j 1`
and used for extract + vrs, then deleted.)

Seed targets from a baseline run (do this once per machine):

```bash
cargo run --release -p gtars-perf -- seed-targets perf/runs/<baseline>.json \
  --out perf/targets.json --margin 0.15
```

Gate a later run (exit `0` pass, `1` regression, `2` bad input/schema):

```bash
cargo run --release -p gtars-perf -- compare perf/runs/<latest>.json \
  --targets perf/targets.json
```

## Output schema (`schema_version: 1`)

One JSON file per run at `perf/runs/<commit>-<timestamp>.json`:

```json
{
  "schema_version": 1,
  "run": {
    "timestamp_utc": "2026-05-29T18:22:04Z",
    "gtars_commit": "6f48d9b",
    "gtars_commit_date": "2026-05-20T11:03:00+00:00",
    "gtars_dirty": false,
    "host": "hostname",
    "cpu_model": "AMD Ryzen ...",
    "logical_cpus": 16,
    "total_ram_mb": 64000,
    "rustc_version": "rustc 1.86.0",
    "profile": "release"
  },
  "results": [
    {
      "task": "encode",
      "concurrency": 4,
      "seconds": 41.2,
      "peak_rss_mb": 1834.0,
      "throughput": 75100000.0,
      "throughput_unit": "bases_per_sec",
      "extra": { "n_bases": 3100000000, "n_sequences": 24, "n_files": 3, "dataset_id": "grch38" }
    }
  ]
}
```

- `results` is a flat array, one object per `(task, concurrency)` cell.
- `peak_rss_mb` is `null` when measured in-process (`--no-subprocess`).
- `extra` holds task-specific correctness fields (digest, region count,
  VRS-id checksum). The gate ignores fields it does not know.
- `--csv` writes a flattened one-row-per-cell mirror beside the JSON.

## Targets + gate

`perf/targets.json`:

```json
{
  "margin": 0.15,
  "seeded_on_host": "reference-host",
  "tasks": {
    "encode":  { "1": {"min_throughput": 18000000, "max_seconds": 200, "max_peak_rss_mb": 2500} },
    "extract": { "1": {"min_throughput": 20000, "max_seconds": 30, "max_peak_rss_mb": 4000} },
    "vrs":     { "8": {"min_throughput": 48000, "max_seconds": 25, "max_peak_rss_mb": 2000} }
  }
}
```

A cell PASSES when every present constraint holds, relaxed by `margin`:
`throughput >= min_throughput*(1-margin)`, `seconds <= max_seconds*(1+margin)`,
`peak_rss_mb <= max_peak_rss_mb*(1+margin)`. Cells with no matching target are
informational (warn, never fail) so new tasks/concurrencies don't break the gate
before targets are seeded.

**Targets are machine-specific.** The gate warns loudly (does not fail) when the
run host differs from the host the targets were seeded on. Re-seed per machine.

## Adding a task

Implement the cell logic in `src/tasks.rs` and wire it into `run_lib_task`
(library tasks) or add a CLI-driven branch like `run_encode_cli`, then add it to
the loop in `cmd_run`. The output schema and gate need no changes — a new
`(task, concurrency)` row is informational until you seed a target for it.

## Deviations from the original plan

- **Encode via the CLI, with `--jobs` not `--threads`.** The plan predated the
  removal of the per-sequence `--threads` knob and the addition of the
  `gtars refget build -j` CLI. Encode is now driven through the real binary over
  a multi-file FASTA set; the sweep is over `--jobs` (concurrent files).
- **Targets/config are JSON, not YAML/TOML.** Keeps the crate dependency-light
  (no yaml/toml parser); the schema is otherwise identical to the plan.
- **No auto-download datasets module.** Datasets are passed by path
  (`--fasta/--bed/--vcf`) with a `--dataset-id` label stamped into every result.
  Wire `refgetstore-benchmark` data dirs in via those flags.
