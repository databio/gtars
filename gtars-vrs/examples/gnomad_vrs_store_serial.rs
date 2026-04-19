//! Single-threaded pure-Rust VRS-ID benchmark — the serial counterpart to
//! `gnomad_vrs_store.rs`. Preloads exactly the same way (open store →
//! `load_all_sequences` → decode every sequence in the first collection)
//! so the open+preload phase is directly comparable, then runs the
//! zero-copy, single-threaded `compute_vrs_ids_streaming_readonly` with a
//! callback that writes each variant's TSV row straight to a
//! `BufWriter<File>`. No worker pool, no channels, no `Vec<VrsResult>`
//! accumulation.
//!
//! This exists so we can answer the "does CPU parallelism actually help
//! a preloaded-store VCF → VRS workload?" question with a pure-Rust
//! serial-vs-parallel comparison that isolates the parallelism decision
//! from any pyo3 or Python-binding overhead.
//!
//! Usage:
//!   cargo run --release --example gnomad_vrs_store_serial -- \
//!       <store_dir> <vcf_bgz_path> <results_tsv> <variants_tsv> [num_workers]
//!
//! Arguments:
//!   store_dir     — RefgetStore directory (as written by RefgetStore::save_local)
//!   vcf_bgz_path  — BGZF-compressed VCF (plain gzip also accepted by open_vcf)
//!   results_tsv   — Output path for the single-row phase-timing TSV
//!   variants_tsv  — Output path for per-variant VRS IDs (chrom/pos/ref/alt/vrs_id)
//!   num_workers   — Optional: accepted for CLI parity with gnomad_vrs_store, ignored
//!                   (this binary is single-threaded by construction).

use std::io::{BufWriter, Write};
use std::time::Instant;

use gtars_refget::store::RefgetStore;
use gtars_vrs::vcf::{build_name_to_digest_readonly, compute_vrs_ids_streaming_readonly};

/// Peak resident set size of the current process, in MB.
///
/// Reads /proc/self/status VmHWM (Linux-only). Matches the Python
/// contestants' `max_rss_mb()` helper so values are comparable across
/// the tourney.
fn max_rss_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|n| n.parse::<u64>().ok())
                .map(|kb| (kb as f64 / 1024.0 * 10.0).round() / 10.0)
        })
        .unwrap_or(0.0)
}

fn usage_and_exit(argv0: &str) -> ! {
    eprintln!(
        "Usage: {} <store_dir> <vcf_bgz_path> <results_tsv> <variants_tsv> [num_workers]",
        argv0
    );
    std::process::exit(2);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if !(5..=6).contains(&args.len()) {
        usage_and_exit(&args[0]);
    }
    let store_dir = &args[1];
    let vcf_path = &args[2];
    let results_tsv = &args[3];
    let variants_tsv = &args[4];
    // num_workers is accepted for CLI parity but intentionally ignored.
    let _num_workers_arg: usize = if args.len() == 6 {
        args[5].parse().unwrap_or(0)
    } else {
        0
    };

    // --- Phase 1: open store ---
    let t0 = Instant::now();
    let mut store = RefgetStore::open_local(store_dir).expect("Failed to open RefgetStore");
    let t_open = t0.elapsed();

    // Identify the collection (first one, matching Python contestants).
    let paged = store
        .list_collections(0, 10, &[])
        .expect("Failed to list collections");
    let collections = &paged.results;
    assert!(!collections.is_empty(), "No collections in store");
    let collection_digest = collections[0].digest.clone();

    // --- Phase 2: preload + decode every sequence ---
    // Mirrors gnomad_vrs_store.rs exactly so open+preload time is comparable.
    let t1 = Instant::now();
    store
        .load_all_sequences()
        .expect("Failed to load_all_sequences");
    let collection = store
        .get_collection(&collection_digest)
        .expect("Failed to get collection")
        .clone();
    for seq_record in &collection.sequences {
        let digest = seq_record.metadata().sha512t24u.clone();
        store
            .ensure_decoded(digest.as_str())
            .unwrap_or_else(|e| panic!("Failed to decode {}: {}", digest, e));
    }
    let t_preload = t1.elapsed();

    // --- Phase 3: convert to readonly and compute serially via streaming callback ---
    let readonly = store.into_readonly();
    let name_to_digest = build_name_to_digest_readonly(&readonly, &collection_digest)
        .expect("Failed to build name_to_digest");

    // Open the per-variant TSV up front; the streaming callback writes
    // each row directly — no Vec<VrsResult> accumulation.
    let vf = std::fs::File::create(variants_tsv).expect("Failed to open variants TSV");
    let mut vw = BufWriter::with_capacity(1 << 20, vf);
    writeln!(vw, "chrom\tpos\tref\talt\tvrs_id").unwrap();

    let t2 = Instant::now();
    let n = compute_vrs_ids_streaming_readonly(
        &readonly,
        &name_to_digest,
        vcf_path,
        |r| {
            writeln!(
                vw,
                "{}\t{}\t{}\t{}\t{}",
                r.chrom, r.pos, r.ref_allele, r.alt_allele, r.vrs_id
            )
            .expect("Failed to write variants TSV row");
        },
    )
    .expect("compute_vrs_ids_streaming_readonly failed");
    let t_compute = t2.elapsed();

    vw.flush().expect("Failed to flush variants TSV");

    let open_s = t_open.as_secs_f64();
    let preload_s = t_preload.as_secs_f64();
    let compute_s = t_compute.as_secs_f64();
    let total_s = open_s + preload_s + compute_s;
    let vps = if compute_s > 0.0 {
        n as f64 / compute_s
    } else {
        0.0
    };

    let vcf_basename = std::path::Path::new(vcf_path)
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| vcf_path.to_string());

    // Match the TSV schema used by gnomad_vrs_store / bgzf_tsv.py exactly:
    //   contestant, vcf, n_variants, num_workers, open_s, preload_s,
    //   compute_s, total_s, variants_per_sec, max_rss_mb
    // num_workers is 1 — this path is single-threaded.
    let rss = max_rss_mb();
    let rf = std::fs::File::create(results_tsv).expect("Failed to open results TSV");
    let mut rw = BufWriter::new(rf);
    writeln!(
        rw,
        "contestant\tvcf\tn_variants\tnum_workers\topen_s\tpreload_s\tcompute_s\ttotal_s\tvariants_per_sec\tmax_rss_mb"
    )
    .unwrap();
    writeln!(
        rw,
        "rust_serial\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.0}\t{:.1}",
        vcf_basename, n, 1, open_s, preload_s, compute_s, total_s, vps, rss
    )
    .unwrap();
    rw.flush().unwrap();

    eprintln!(
        "[rust_serial] n={} workers=1 open={:.2}s preload={:.2}s compute={:.2}s total={:.2}s  {:.0} v/s  rss={:.1} MB  (wrote {})",
        n, open_s, preload_s, compute_s, total_s, vps, rss, variants_tsv
    );
}
