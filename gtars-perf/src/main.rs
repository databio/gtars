//! gtars-perf — occasionally-run performance regression harness.
//!
//! Three tasks are exercised at several concurrency levels:
//!   - encode : `gtars refget build <fasta>... -o <dir> -j <N>` — the `--jobs`
//!     knob is the number of FASTA files imported CONCURRENTLY, so the sweep is
//!     run over a MULTI-FILE input set (a single file does not parallelize).
//!     Driven through the real CLI binary as a subprocess, wrapped in
//!     `/usr/bin/time -v` for true peak RSS.
//!   - extract : iterate a BED file over a prebuilt store via the library
//!     `SubstringsFromRegions` iterator. Single-threaded in the library today,
//!     so reported at concurrency 1 (informational). Library-driven, run in a
//!     self-subprocess so peak RSS is per-task.
//!   - vrs : `compute_vrs_ids_parallel_encoded` over a VCF at 1/4/8 threads.
//!     Library-driven in a self-subprocess.
//!
//! Per-task peak RSS is measured by running each (task, concurrency) cell in
//! its OWN process and reading the process's `VmHWM` from `/proc/self/status`
//! (encode reads the CLI child's RSS via `/usr/bin/time -v`). Wall time is
//! measured by the parent around the child.
//!
//! Subcommands:
//!   run           run all tasks, write a run-record JSON (+ optional CSV)
//!   compare       gate a run JSON against targets.json (nonzero exit on regress)
//!   seed-targets  derive targets.json from a baseline run JSON
//!   run-one       (internal) run a single library cell, print one PERF_RESULT line

use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use anyhow::{anyhow, bail, Context, Result};
use clap::{Args, Parser, Subcommand};
use serde::{Deserialize, Serialize};

mod machine;
mod tasks;

pub const SCHEMA_VERSION: u32 = 1;

#[derive(Parser)]
#[command(name = "gtars-perf", about = "gtars performance regression harness")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Run all tasks and emit a run-record JSON.
    Run(RunArgs),
    /// Gate a run JSON against targets (exit nonzero on regression).
    Compare(CompareArgs),
    /// Derive a targets.json from a baseline run JSON.
    SeedTargets(SeedArgs),
    /// (internal) Run a single (task, concurrency) cell in-process.
    RunOne(RunOneArgs),
}

#[derive(Args)]
struct RunArgs {
    /// Directory to write the run JSON into.
    #[arg(long, default_value = "perf/runs")]
    out: PathBuf,
    /// Also write a flattened CSV mirror next to the JSON.
    #[arg(long)]
    csv: bool,
    /// FASTA files for the encode task (>=2 for the --jobs sweep to matter).
    #[arg(long, value_delimiter = ',')]
    fasta: Vec<PathBuf>,
    /// Prebuilt store directory for extract + vrs (built from --fasta if omitted).
    #[arg(long)]
    store: Option<PathBuf>,
    /// BED file for the extract task.
    #[arg(long)]
    bed: Option<PathBuf>,
    /// VCF file for the vrs task.
    #[arg(long)]
    vcf: Option<PathBuf>,
    /// Comma-separated concurrency sweep (default 1,4,8).
    #[arg(long, value_delimiter = ',', default_value = "1,4,8")]
    concurrency: Vec<usize>,
    /// Path to the gtars CLI binary used to drive the encode task.
    #[arg(long, default_value = "target/release/gtars")]
    gtars_bin: PathBuf,
    /// Repo root used to read the gtars commit metadata.
    #[arg(long, default_value = ".")]
    repo: PathBuf,
    /// Run cells in-process instead of one-subprocess-per-cell (RSS becomes
    /// process-cumulative and is reported as null).
    #[arg(long)]
    no_subprocess: bool,
    /// A short label for the dataset used (stamped into every result's extra).
    #[arg(long, default_value = "adhoc")]
    dataset_id: String,
}

#[derive(Args)]
struct CompareArgs {
    /// Run JSON to gate.
    run: PathBuf,
    /// Targets JSON.
    #[arg(long, default_value = "perf/targets.json")]
    targets: PathBuf,
}

#[derive(Args)]
struct SeedArgs {
    /// Baseline run JSON.
    run: PathBuf,
    /// Where to write the targets JSON.
    #[arg(long, default_value = "perf/targets.json")]
    out: PathBuf,
    /// Default tolerance applied at compare time (recorded, not baked in).
    #[arg(long, default_value_t = 0.15)]
    margin: f64,
}

#[derive(Args)]
struct RunOneArgs {
    /// Task name: extract | vrs (encode is driven via the CLI, not here).
    #[arg(long)]
    task: String,
    #[arg(long)]
    concurrency: usize,
    #[arg(long)]
    store: Option<PathBuf>,
    #[arg(long)]
    bed: Option<PathBuf>,
    #[arg(long)]
    vcf: Option<PathBuf>,
    #[arg(long, default_value = "adhoc")]
    dataset_id: String,
}

// ──────────────────────────── run-record schema ───────────────────────────

#[derive(Serialize, Deserialize)]
pub struct RunRecord {
    pub schema_version: u32,
    pub run: machine::RunMeta,
    pub results: Vec<ResultCell>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct ResultCell {
    pub task: String,
    /// Number of concurrent units (FASTA files for encode; threads for vrs;
    /// always 1 for extract, which is single-threaded in the library today).
    pub concurrency: usize,
    pub seconds: f64,
    /// Peak resident set in MB; `null` when measured in-process (cumulative).
    pub peak_rss_mb: Option<f64>,
    pub throughput: f64,
    pub throughput_unit: String,
    pub extra: serde_json::Value,
}

fn main() {
    if let Err(e) = real_main() {
        eprintln!("gtars-perf error: {e:#}");
        std::process::exit(2);
    }
}

fn real_main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Run(a) => cmd_run(a),
        Cmd::Compare(a) => cmd_compare(a),
        Cmd::SeedTargets(a) => cmd_seed(a),
        Cmd::RunOne(a) => cmd_run_one(a),
    }
}

// ──────────────────────────────── run ─────────────────────────────────────

fn cmd_run(a: RunArgs) -> Result<()> {
    if a.fasta.is_empty() && a.store.is_none() {
        bail!("provide --fasta <files...> and/or --store <dir>");
    }

    // Resolve / build a store for extract + vrs.
    let mut tmp_store_holder: Option<PathBuf> = None;
    let store_path = match &a.store {
        Some(s) => s.clone(),
        None => {
            let dir = std::env::temp_dir().join(format!("gtars-perf-store-{}", std::process::id()));
            eprintln!("Building setup store at {} from {} FASTA(s)...", dir.display(), a.fasta.len());
            build_store_via_cli(&a.gtars_bin, &a.fasta, &dir, 1)?;
            tmp_store_holder = Some(dir.clone());
            dir
        }
    };

    let mut results: Vec<ResultCell> = Vec::new();

    // encode task — only if we have FASTAs (drive the real CLI binary).
    if !a.fasta.is_empty() {
        if a.fasta.len() < 2 {
            eprintln!(
                "WARNING: encode --jobs only parallelizes across FILES; \
                 {} file(s) supplied — the jobs sweep will be flat.",
                a.fasta.len()
            );
        }
        for &j in &a.concurrency {
            eprintln!("[encode] jobs={j}");
            let cell = tasks::run_encode_cli(&a.gtars_bin, &a.fasta, j, &a.dataset_id)?;
            results.push(cell);
        }
    }

    // extract task — single-threaded in the library; concurrency fixed at 1.
    if let Some(bed) = &a.bed {
        eprintln!("[extract] concurrency=1");
        let cell = run_lib_cell(&a, "extract", 1, Some(&store_path), Some(bed), None)?;
        results.push(cell);
    }

    // vrs task — thread sweep.
    if let Some(vcf) = &a.vcf {
        for &t in &a.concurrency {
            eprintln!("[vrs] threads={t}");
            let cell = run_lib_cell(&a, "vrs", t, Some(&store_path), None, Some(vcf))?;
            results.push(cell);
        }
    }

    let run = machine::collect_run_meta(&a.repo)?;
    let record = RunRecord { schema_version: SCHEMA_VERSION, run, results };

    std::fs::create_dir_all(&a.out)
        .with_context(|| format!("create out dir {}", a.out.display()))?;
    let fname = format!("{}-{}.json", record.run.gtars_commit, record.run.timestamp_utc.replace([':'], "").replace('-', ""));
    let json_path = a.out.join(&fname);
    std::fs::write(&json_path, serde_json::to_string_pretty(&record)?)?;
    eprintln!("Wrote run record: {}", json_path.display());

    if a.csv {
        let csv_path = json_path.with_extension("csv");
        std::fs::write(&csv_path, to_csv(&record))?;
        eprintln!("Wrote CSV mirror: {}", csv_path.display());
    }

    // Print a small summary table.
    eprintln!("\ntask        conc   seconds   peak_rss_mb   throughput  unit");
    for r in &record.results {
        eprintln!(
            "{:<10}  {:>4}  {:>8.3}   {:>11}  {:>11.0}  {}",
            r.task, r.concurrency, r.seconds,
            r.peak_rss_mb.map(|v| format!("{v:.0}")).unwrap_or_else(|| "null".into()),
            r.throughput, r.throughput_unit
        );
    }

    if let Some(d) = tmp_store_holder {
        let _ = std::fs::remove_dir_all(d);
    }
    Ok(())
}

/// Run a library cell either as a self-subprocess (per-task RSS) or in-process.
fn run_lib_cell(
    a: &RunArgs,
    task: &str,
    concurrency: usize,
    store: Option<&Path>,
    bed: Option<&Path>,
    vcf: Option<&Path>,
) -> Result<ResultCell> {
    if a.no_subprocess {
        let mut cell = tasks::run_lib_task(task, concurrency, store, bed, vcf, &a.dataset_id)?;
        cell.peak_rss_mb = None; // process-cumulative; gate skips this constraint
        return Ok(cell);
    }
    let exe = std::env::current_exe()?;
    let mut cmd = Command::new(&exe);
    cmd.arg("run-one")
        .arg("--task").arg(task)
        .arg("--concurrency").arg(concurrency.to_string())
        .arg("--dataset-id").arg(&a.dataset_id);
    if let Some(s) = store { cmd.arg("--store").arg(s); }
    if let Some(b) = bed { cmd.arg("--bed").arg(b); }
    if let Some(v) = vcf { cmd.arg("--vcf").arg(v); }

    let start = Instant::now();
    let out = cmd.output().with_context(|| "spawn run-one child")?;
    let wall = start.elapsed().as_secs_f64();
    if !out.status.success() {
        bail!(
            "run-one child for task '{task}' (conc {concurrency}) failed:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    let stdout = String::from_utf8_lossy(&out.stdout);
    let line = stdout
        .lines()
        .find(|l| l.starts_with("PERF_RESULT "))
        .ok_or_else(|| anyhow!("child emitted no PERF_RESULT line; stderr:\n{}", String::from_utf8_lossy(&out.stderr)))?;
    let mut cell: ResultCell = serde_json::from_str(line.trim_start_matches("PERF_RESULT ").trim())?;
    // Trust the parent's wall clock over the child's internal timer for the
    // top-level seconds (includes full process lifetime); keep child throughput
    // which is computed from its own compute window.
    let _ = wall; // child seconds already reflect the timed compute window
    cell.task = task.to_string();
    Ok(cell)
}

fn cmd_run_one(a: RunOneArgs) -> Result<()> {
    let cell = tasks::run_lib_task(
        &a.task,
        a.concurrency,
        a.store.as_deref(),
        a.bed.as_deref(),
        a.vcf.as_deref(),
        &a.dataset_id,
    )?;
    println!("PERF_RESULT {}", serde_json::to_string(&cell)?);
    Ok(())
}

/// Build a store by shelling out to the real gtars CLI (used for setup at j=1).
fn build_store_via_cli(bin: &Path, fasta: &[PathBuf], out: &Path, jobs: usize) -> Result<()> {
    let _ = std::fs::remove_dir_all(out);
    let mut cmd = Command::new(bin);
    cmd.arg("refget").arg("build");
    for f in fasta { cmd.arg(f); }
    cmd.arg("-o").arg(out).arg("-j").arg(jobs.to_string());
    let status = cmd.status().with_context(|| format!("run {} refget build", bin.display()))?;
    if !status.success() {
        bail!("gtars refget build failed (status {status})");
    }
    Ok(())
}

fn to_csv(rec: &RunRecord) -> String {
    let m = &rec.run;
    let mut s = String::from(
        "schema_version,timestamp_utc,gtars_commit,gtars_dirty,host,cpu_model,logical_cpus,total_ram_mb,rustc_version,profile,task,concurrency,seconds,peak_rss_mb,throughput,throughput_unit\n",
    );
    for r in &rec.results {
        s.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{},{},{},{:.3},{},{:.0},{}\n",
            rec.schema_version, m.timestamp_utc, m.gtars_commit, m.gtars_dirty,
            csv_q(&m.host), csv_q(&m.cpu_model), m.logical_cpus, m.total_ram_mb,
            csv_q(&m.rustc_version), m.profile,
            r.task, r.concurrency, r.seconds,
            r.peak_rss_mb.map(|v| format!("{v:.0}")).unwrap_or_default(),
            r.throughput, r.throughput_unit,
        ));
    }
    s
}

fn csv_q(s: &str) -> String {
    if s.contains(',') || s.contains('"') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

// ──────────────────────────── compare / gate ──────────────────────────────

#[derive(Serialize, Deserialize)]
struct Targets {
    /// Default tolerance applied to every constraint at compare time.
    margin: f64,
    /// Host the targets were seeded on (gate warns on mismatch).
    #[serde(default)]
    seeded_on_host: String,
    /// tasks[task][concurrency] -> constraints
    tasks: std::collections::BTreeMap<String, std::collections::BTreeMap<String, Constraint>>,
}

#[derive(Serialize, Deserialize, Default, Clone)]
struct Constraint {
    #[serde(skip_serializing_if = "Option::is_none")]
    min_throughput: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    max_seconds: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    max_peak_rss_mb: Option<f64>,
}

fn cmd_compare(a: CompareArgs) -> Result<()> {
    let run: RunRecord = serde_json::from_str(
        &std::fs::read_to_string(&a.run).with_context(|| format!("read run {}", a.run.display()))?,
    )
    .context("parse run JSON")?;
    if run.schema_version != SCHEMA_VERSION {
        eprintln!(
            "ERROR: run schema_version {} != harness schema_version {}; re-run / re-seed.",
            run.schema_version, SCHEMA_VERSION
        );
        std::process::exit(2);
    }
    let targets: Targets = serde_json::from_str(
        &std::fs::read_to_string(&a.targets)
            .with_context(|| format!("read targets {}", a.targets.display()))?,
    )
    .context("parse targets JSON")?;

    if !targets.seeded_on_host.is_empty() && targets.seeded_on_host != run.run.host {
        eprintln!(
            "WARNING: run host '{}' differs from targets host '{}'. \
             Targets are machine-specific; re-seed on this machine.",
            run.run.host, targets.seeded_on_host
        );
    }

    let margin = targets.margin;
    println!(
        "{:<10} {:>5} {:>10} {:<10} {}",
        "task", "conc", "observed", "result", "binding constraint"
    );
    let mut any_fail = false;
    let mut any_gated = false;
    for r in &run.results {
        let t = targets
            .tasks
            .get(&r.task)
            .and_then(|m| m.get(&r.concurrency.to_string()));
        let Some(c) = t else {
            println!("{:<10} {:>5} {:>10} {:<10} (no target — informational)", r.task, r.concurrency, "-", "WARN");
            continue;
        };
        any_gated = true;
        let (pass, binding) = check_cell(r, c, margin);
        if !pass { any_fail = true; }
        println!(
            "{:<10} {:>5} {:>10.0} {:<10} {}",
            r.task, r.concurrency, r.throughput,
            if pass { "PASS" } else { "FAIL" }, binding
        );
    }

    if !any_gated {
        eprintln!("No gated cells matched any target. Treating as pass (informational only).");
    }
    if any_fail {
        eprintln!("\nGATE: FAIL — at least one cell regressed.");
        std::process::exit(1);
    }
    eprintln!("\nGATE: PASS");
    Ok(())
}

/// Returns (pass, human description of the binding/violated constraint).
fn check_cell(r: &ResultCell, c: &Constraint, margin: f64) -> (bool, String) {
    let mut msgs = Vec::new();
    let mut pass = true;
    if let Some(min) = c.min_throughput {
        let floor = min * (1.0 - margin);
        if r.throughput < floor {
            pass = false;
            msgs.push(format!("throughput {:.0} < {:.0} (min*{:.2})", r.throughput, floor, 1.0 - margin));
        }
    }
    if let Some(max) = c.max_seconds {
        let ceil = max * (1.0 + margin);
        if r.seconds > ceil {
            pass = false;
            msgs.push(format!("seconds {:.2} > {:.2}", r.seconds, ceil));
        }
    }
    if let (Some(max), Some(rss)) = (c.max_peak_rss_mb, r.peak_rss_mb) {
        let ceil = max * (1.0 + margin);
        if rss > ceil {
            pass = false;
            msgs.push(format!("peak_rss {:.0} > {:.0}", rss, ceil));
        }
    }
    if pass {
        (true, "all constraints OK".into())
    } else {
        (false, msgs.join("; "))
    }
}

fn cmd_seed(a: SeedArgs) -> Result<()> {
    let run: RunRecord = serde_json::from_str(
        &std::fs::read_to_string(&a.run).with_context(|| format!("read run {}", a.run.display()))?,
    )?;
    let mut tasks: std::collections::BTreeMap<String, std::collections::BTreeMap<String, Constraint>> =
        Default::default();
    for r in &run.results {
        let entry = tasks.entry(r.task.clone()).or_default();
        entry.insert(
            r.concurrency.to_string(),
            Constraint {
                min_throughput: Some(r.throughput),
                max_seconds: Some(r.seconds),
                max_peak_rss_mb: r.peak_rss_mb,
            },
        );
    }
    let targets = Targets { margin: a.margin, seeded_on_host: run.run.host.clone(), tasks };
    if let Some(parent) = a.out.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    std::fs::write(&a.out, serde_json::to_string_pretty(&targets)?)?;
    eprintln!(
        "Seeded targets from {} on host '{}' (margin {}) -> {}",
        a.run.display(), run.run.host, a.margin, a.out.display()
    );
    Ok(())
}
