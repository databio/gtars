//! Pure-Rust VRS-ID benchmark matching the `bgzf_tsv` Python contestant.
//!
//! Opens a pre-built RefgetStore directory (produced by the benchmark's
//! `build_store` task), preloads all sequences + decodes them, then runs
//! either the parallel or serial VRS computation. Emits the same TSV output
//! format as the Python contestants expect: a single header row plus a single
//! data row with phase-split timings.
//!
//! Usage:
//!   cargo run --release --example gnomad_vrs_store -- \
//!       <store_dir> <vcf_bgz_path> <results_tsv> <variants_tsv> [--serial] [num_workers]
//!
//! Arguments:
//!   store_dir     — RefgetStore directory (as written by RefgetStore::save_local)
//!   vcf_bgz_path  — BGZF-compressed VCF (plain gzip is rejected by parallel path)
//!   results_tsv   — Output path for the single-row phase-timing TSV
//!   variants_tsv  — Output path for per-variant VRS IDs (chrom/pos/ref/alt/vrs_id)
//!   --serial      — Optional flag: run single-threaded streaming path instead of parallel
//!   num_workers   — Optional: number of VRS worker threads (parallel path only;
//!                   default: available_parallelism - 2)

#[cfg(target_os = "linux")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use std::io::{BufWriter, Write};
use std::time::Instant;

use gtars_refget::store::RefgetStore;
use gtars_vrs::vcf::{
    build_name_to_digest_readonly, compute_vrs_ids_parallel_bgzf_with_sink,
    compute_vrs_ids_streaming_readonly, decode_vcf_chroms,
};

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
        "Usage: {} <store_dir> <vcf_bgz_path> <results_tsv> <variants_tsv> [--serial] [num_workers]",
        argv0
    );
    std::process::exit(2);
}

fn main() {
    let raw_args: Vec<String> = std::env::args().collect();
    let serial = raw_args.iter().any(|a| a == "--serial");
    let args: Vec<&String> = raw_args.iter().filter(|a| *a != "--serial").collect();

    if !(5..=6).contains(&args.len()) {
        usage_and_exit(&raw_args[0]);
    }
    let store_dir = args[1];
    let vcf_path = args[2];
    let results_tsv = args[3];
    let variants_tsv = args[4];
    let num_workers: usize = if args.len() == 6 {
        args[5].parse().unwrap_or(0)
    } else {
        0
    };
    let num_workers = if num_workers == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get().saturating_sub(2).max(1))
            .unwrap_or(1)
    } else {
        num_workers
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

    // --- Phase 2: decode only the sequences the VCF references ---
    // Lazy per-chromosome decode — avoids eager `load_all_sequences` +
    // full-collection `ensure_decoded`, which for a single-chromosome VCF
    // against GRCh38 would waste ~3-9 GB of RSS and ~15 s of wall-clock.
    let t1 = Instant::now();
    decode_vcf_chroms(&mut store, &collection_digest, vcf_path)
        .expect("Failed to decode VCF-referenced sequences");
    let t_preload = t1.elapsed();

    // --- Phase 3: compute VRS IDs (parallel or serial) ---
    let readonly = store.into_readonly();
    let name_to_digest = build_name_to_digest_readonly(&readonly, &collection_digest)
        .expect("Failed to build name_to_digest");

    // Open the per-variant TSV up front; the sink/callback writes each
    // row directly — no Vec<VrsResult> accumulation.
    let vf = std::fs::File::create(variants_tsv).expect("Failed to open variants TSV");
    let mut vw = BufWriter::with_capacity(1 << 20, vf);
    writeln!(vw, "chrom\tpos\tref\talt\tvrs_id").unwrap();

    let t2 = Instant::now();
    let (n, contestant, workers_used) = if serial {
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
        (n, "rust_serial", 1usize)
    } else {
        let n = compute_vrs_ids_parallel_bgzf_with_sink(
            &readonly,
            &name_to_digest,
            vcf_path,
            num_workers,
            |r| {
                writeln!(
                    vw,
                    "{}\t{}\t{}\t{}\t{}",
                    r.chrom, r.pos, r.ref_allele, r.alt_allele, r.vrs_id
                )
                .expect("Failed to write variants TSV row");
            },
        )
        .expect("compute_vrs_ids_parallel_bgzf_with_sink failed");
        (n, "rust_native", num_workers)
    };
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

    // Match the same header as bgzf_tsv.py, with num_workers and max_rss_mb included.
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
        "{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.0}\t{:.1}",
        contestant, vcf_basename, n, workers_used, open_s, preload_s, compute_s, total_s, vps, rss
    )
    .unwrap();
    rw.flush().unwrap();

    eprintln!(
        "[{}] n={} workers={} open={:.2}s preload={:.2}s compute={:.2}s total={:.2}s  {:.0} v/s  rss={:.1} MB  (wrote {})",
        contestant, n, workers_used, open_s, preload_s, compute_s, total_s, vps, rss, variants_tsv
    );
}
