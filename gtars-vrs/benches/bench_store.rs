//! Benchmark: Load a RefgetStore from disk and compute VRS IDs from a VCF.
//! Reports throughput and peak RSS for mmap vs in-memory comparison.
//!
//! Usage:
//!   cargo bench --bench bench_store --features filesystem -- <store_path> <vcf_path>

use std::time::Instant;

use gtars_refget::store::RefgetStore;
use gtars_vrs::vcf::compute_vrs_ids_streaming;

fn peak_rss_mb() -> f64 {
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    for line in status.lines() {
        if line.starts_with("VmHWM:") {
            let kb: f64 = line
                .split_whitespace()
                .nth(1)
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0);
            return kb / 1024.0;
        }
    }
    0.0
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <store_path> <vcf_path>", args[0]);
        std::process::exit(1);
    }
    let store_path = &args[1];
    let vcf_path = &args[2];

    // Step 1: Open existing RefgetStore from disk
    eprintln!("Opening store: {}", store_path);
    let t0 = Instant::now();
    let mut store = RefgetStore::open_local(store_path).expect("Failed to open store");
    let load_time = t0.elapsed();
    eprintln!("Opened in {:.3}s", load_time.as_secs_f64());
    eprintln!("Peak RSS after open: {:.0} MB", peak_rss_mb());

    // Get collection digest
    let paged = store.list_collections(0, 10, &[]).expect("Failed to list collections");
    let collections = &paged.results;
    assert!(!collections.is_empty(), "No collections found in store");
    let collection_digest = &collections[0].digest;
    eprintln!(
        "Collection: {} ({} sequences)",
        collection_digest, collections[0].n_sequences
    );

    // Preload all sequences
    eprintln!("Preloading sequences...");
    let preload_start = Instant::now();
    store.load_all_sequences().expect("Failed to preload sequences");
    eprintln!("Preloaded in {:.3}s", preload_start.elapsed().as_secs_f64());
    eprintln!("Peak RSS after preload: {:.0} MB", peak_rss_mb());

    // Step 2: Stream VRS IDs from VCF (suppress output for speed)
    eprintln!("\nProcessing VCF: {}", vcf_path);
    let t1 = Instant::now();
    let count = compute_vrs_ids_streaming(&mut store, collection_digest, vcf_path, |_r| {
        // suppress output for benchmark
    })
    .expect("Failed to compute VRS IDs");
    let vrs_time = t1.elapsed();

    let throughput = count as f64 / vrs_time.as_secs_f64();
    let peak_mb = peak_rss_mb();

    eprintln!("\n=== BENCHMARK RESULTS ===");
    eprintln!("Variants:    {}", count);
    eprintln!("VRS time:    {:.3}s", vrs_time.as_secs_f64());
    eprintln!("Throughput:  {:.0} variants/sec", throughput);
    eprintln!("Peak RSS:    {:.0} MB", peak_mb);
    eprintln!("Total time:  {:.1}s", t0.elapsed().as_secs_f64());
}
