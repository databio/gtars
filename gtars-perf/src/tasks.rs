//! Task implementations. Each produces a [`ResultCell`]. `peak_rss_mb` for the
//! library tasks is the process `VmHWM` (authoritative when each cell runs in
//! its own subprocess); for the encode task it comes from `/usr/bin/time -v`.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use anyhow::{anyhow, bail, Context, Result};
use serde_json::json;

use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore};
use gtars_vrs::vcf::compute_vrs_ids_parallel_encoded;

use crate::ResultCell;

/// Peak resident set of the current process, in MB, via `/proc/self/status`.
pub fn peak_rss_mb() -> f64 {
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            let kb: f64 = rest.split_whitespace().next().and_then(|s| s.parse().ok()).unwrap_or(0.0);
            return kb / 1024.0;
        }
    }
    0.0
}

// ──────────────────────────── encode (CLI) ────────────────────────────────

/// Drive `gtars refget build <fasta>... -o <tmp> -j <jobs>` through the real CLI
/// binary, wrapped in `/usr/bin/time -v` so we capture the child's true peak RSS.
/// The `.rgsi` sidecars next to each FASTA are cleared first so the
/// digest+encode stage is actually exercised, and a fresh output dir is used.
pub fn run_encode_cli(
    bin: &Path,
    fasta: &[PathBuf],
    jobs: usize,
    dataset_id: &str,
) -> Result<ResultCell> {
    for f in fasta {
        clear_rgsi_sidecar(f);
    }
    let out_dir = std::env::temp_dir().join(format!(
        "gtars-perf-encode-j{}-{}",
        jobs,
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&out_dir);

    // /usr/bin/time -v writes its report to stderr; -o keeps it out of the
    // child's own stderr so parsing is robust.
    let time_log = std::env::temp_dir().join(format!("gtars-perf-time-j{}-{}.log", jobs, std::process::id()));
    let mut cmd = Command::new("/usr/bin/time");
    cmd.arg("-v").arg("-o").arg(&time_log).arg(bin);
    cmd.arg("refget").arg("build");
    for f in fasta {
        cmd.arg(f);
    }
    cmd.arg("-o").arg(&out_dir).arg("-j").arg(jobs.to_string());

    let start = Instant::now();
    let status = cmd.status().with_context(|| format!("run {} via /usr/bin/time", bin.display()))?;
    let seconds = start.elapsed().as_secs_f64();
    if !status.success() {
        bail!("encode (jobs={jobs}) failed with status {status}");
    }

    let peak_rss_mb = parse_time_v_max_rss_kb(&time_log).map(|kb| kb / 1024.0);

    // Read store to get total bases, sequence count, and a stable digest for the
    // correctness cross-check across jobs.
    let mut store = RefgetStore::open_local(&out_dir).context("open built store")?;
    store.load_all_collections().context("load_all_collections")?;
    let n_bases: u64 = store.list_sequences().iter().map(|m| m.length as u64).sum();
    let n_sequences = store.list_sequences().len();
    let coll_digest = store
        .list_collections(0, 1, &[])
        .ok()
        .and_then(|p| p.results.first().map(|c| c.digest.clone()))
        .unwrap_or_default();

    let throughput = if seconds > 0.0 { n_bases as f64 / seconds } else { 0.0 };

    let _ = std::fs::remove_dir_all(&out_dir);
    let _ = std::fs::remove_file(&time_log);

    Ok(ResultCell {
        task: "encode".into(),
        concurrency: jobs,
        seconds,
        peak_rss_mb,
        throughput,
        throughput_unit: "bases_per_sec".into(),
        extra: json!({
            "dataset_id": dataset_id,
            "n_bases": n_bases,
            "n_sequences": n_sequences,
            "n_files": fasta.len(),
            "collection_digest": coll_digest,
        }),
    })
}

fn clear_rgsi_sidecar(fasta: &Path) {
    // The import caches a `<fasta>.rgsi` (and variants) sidecar; remove any
    // sibling whose name starts with the fasta filename and ends in `.rgsi`.
    if let Some(dir) = fasta.parent() {
        if let Some(stem) = fasta.file_name().and_then(|s| s.to_str()) {
            if let Ok(rd) = std::fs::read_dir(dir) {
                for e in rd.flatten() {
                    let name = e.file_name();
                    let name = name.to_string_lossy();
                    if name.starts_with(stem) && name.ends_with(".rgsi") {
                        let _ = std::fs::remove_file(e.path());
                    }
                }
            }
        }
    }
}

fn parse_time_v_max_rss_kb(log: &Path) -> Option<f64> {
    let text = std::fs::read_to_string(log).ok()?;
    for line in text.lines() {
        let l = line.trim();
        if let Some(rest) = l.strip_prefix("Maximum resident set size (kbytes):") {
            return rest.trim().parse::<f64>().ok();
        }
    }
    None
}

// ─────────────────── library tasks (extract / vrs) ────────────────────────

/// Dispatch a single library cell. Called both in-process and from the
/// `run-one` subprocess. Reports `peak_rss_mb` from this process's `VmHWM`.
pub fn run_lib_task(
    task: &str,
    concurrency: usize,
    store: Option<&Path>,
    bed: Option<&Path>,
    vcf: Option<&Path>,
    dataset_id: &str,
) -> Result<ResultCell> {
    match task {
        "extract" => run_extract(
            store.ok_or_else(|| anyhow!("extract needs --store"))?,
            bed.ok_or_else(|| anyhow!("extract needs --bed"))?,
            dataset_id,
        ),
        "vrs" => run_vrs(
            store.ok_or_else(|| anyhow!("vrs needs --store"))?,
            vcf.ok_or_else(|| anyhow!("vrs needs --vcf"))?,
            concurrency,
            dataset_id,
        ),
        other => bail!("unknown library task '{other}' (encode is driven via the CLI)"),
    }
}

/// Iterate a BED over a store via the `SubstringsFromRegions` iterator.
/// Single-threaded in the library today, so concurrency is fixed at 1 and
/// reported informationally. Throughput is `regions_per_sec`; extracted bytes
/// and an order-independent checksum land in `extra` for correctness checks.
fn run_extract(store_path: &Path, bed: &Path, dataset_id: &str) -> Result<ResultCell> {
    let mut store = RefgetStore::open_local(store_path).context("open store for extract")?;
    store.load_all_collections().context("load_all_collections")?;
    let coll_digest = store
        .list_collections(0, 1, &[])?
        .results
        .first()
        .ok_or_else(|| anyhow!("no collections in store"))?
        .digest
        .clone();
    // Make all sequences resident so extraction does not pay per-region disk IO.
    let digests: Vec<String> = store.list_sequences().iter().map(|m| m.sha512t24u.clone()).collect();
    for d in &digests {
        store.load_sequence(d).ok();
    }
    let store: ReadonlyRefgetStore = store.into_readonly();

    let bed_str = bed.to_string_lossy().to_string();
    let start = Instant::now();
    let mut n_regions: u64 = 0;
    let mut n_bytes: u64 = 0;
    let mut checksum: u64 = 0;
    for item in store.substrings_from_regions(coll_digest.as_bytes(), &bed_str)? {
        let seq = item?;
        let bytes = seq.sequence.as_bytes();
        n_bytes += bytes.len() as u64;
        // order-independent per-region hash, summed
        let mut h: u64 = 1469598103934665603; // FNV offset
        for &b in bytes {
            h ^= b as u64;
            h = h.wrapping_mul(1099511628211);
        }
        checksum = checksum.wrapping_add(h);
        n_regions += 1;
    }
    let seconds = start.elapsed().as_secs_f64();
    let throughput = if seconds > 0.0 { n_regions as f64 / seconds } else { 0.0 };
    let bytes_per_sec = if seconds > 0.0 { n_bytes as f64 / seconds } else { 0.0 };

    Ok(ResultCell {
        task: "extract".into(),
        concurrency: 1,
        seconds,
        peak_rss_mb: Some(peak_rss_mb()),
        throughput,
        throughput_unit: "regions_per_sec".into(),
        extra: json!({
            "dataset_id": dataset_id,
            "n_regions": n_regions,
            "n_bytes": n_bytes,
            "bytes_per_sec": bytes_per_sec,
            "checksum": checksum,
            "note": "single-threaded in library; concurrency fixed at 1",
        }),
    })
}

/// Compute VRS Allele ids over a VCF at `threads` workers via the library's
/// parallel encoded path. Throughput is `variants_per_sec`; a VRS-id checksum
/// is recorded in `extra` so a regression that also changes output is caught.
fn run_vrs(store_path: &Path, vcf: &Path, threads: usize, dataset_id: &str) -> Result<ResultCell> {
    let mut store = RefgetStore::open_local(store_path).context("open store for vrs")?;
    store.load_all_collections().context("load_all_collections")?;
    let coll_digest = store
        .list_collections(0, 1, &[])?
        .results
        .first()
        .ok_or_else(|| anyhow!("no collections in store"))?
        .digest
        .clone();
    let collection = store.get_collection(&coll_digest).context("get_collection")?;
    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    for rec in &collection.sequences {
        let m = rec.metadata();
        name_to_digest.insert(m.name.clone(), m.sha512t24u.clone());
    }
    // Make referenced sequences resident, then freeze to readonly.
    for d in name_to_digest.values() {
        store.load_sequence(d).ok();
    }
    let store: ReadonlyRefgetStore = store.into_readonly();

    let vcf_str = vcf.to_string_lossy().to_string();
    let mut checksum: u64 = 0;
    let start = Instant::now();
    let count = compute_vrs_ids_parallel_encoded(
        &store,
        &name_to_digest,
        &vcf_str,
        threads,
        |res| {
            for &b in res.vrs_id.as_bytes() {
                checksum = checksum.wrapping_add(b as u64);
            }
        },
    )
    .context("compute_vrs_ids_parallel_encoded")?;
    let seconds = start.elapsed().as_secs_f64();
    let throughput = if seconds > 0.0 { count as f64 / seconds } else { 0.0 };

    Ok(ResultCell {
        task: "vrs".into(),
        concurrency: threads,
        seconds,
        peak_rss_mb: Some(peak_rss_mb()),
        throughput,
        throughput_unit: "variants_per_sec".into(),
        extra: json!({
            "dataset_id": dataset_id,
            "n_variants": count,
            "checksum": checksum,
        }),
    })
}
