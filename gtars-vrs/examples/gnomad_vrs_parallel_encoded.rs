//! Compute VRS Allele IDs from a (gnomAD) VCF using the immutable 2-bit-encoded
//! refget store with in-process parallel workers (decode-on-the-fly).
//!
//! Mirrors the store setup of `bench_encoded.rs`: open a local store, load all
//! collections, build the name->digest map, make the needed sequences resident
//! with `load_sequence`, freeze with `into_readonly`, then call the parallel
//! encoded path. Prints throughput (variants/sec) and the thread count.
//!
//! Usage:
//!   cargo run --release --example gnomad_vrs_parallel_encoded -- <store_path> <vcf_path> [threads] [--print] [--bgzf]
//!     threads: worker count (default: available_parallelism)
//!     --print: also print each result row (chrom pos ref alt vrs_id) to stdout
//!     --bgzf:  use the BGZF-block-parallel path (raw-block reader + workers own
//!              decompression). Requires BGZF (`.bgz`) input. Scales past the
//!              single-reader default path. Without it, the single-reader
//!              `compute_vrs_ids_parallel_encoded` path is used.

use std::collections::HashMap;
use std::time::Instant;

use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore};
use gtars_vrs::vcf::{
    compute_vrs_ids_parallel_bgzf_encoded_with_sink, compute_vrs_ids_parallel_encoded,
    parse_vcf_record,
};

/// Peak resident set size (RSS) in kilobytes, read from /proc/self/status
/// (VmHWM = "high water mark"). Returns 0 if unavailable.
fn peak_rss_kb() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|v| v.parse().ok())
        })
        .unwrap_or(0)
}

fn open_vcf_reader(path: &str) -> Box<dyn std::io::BufRead> {
    use flate2::read::MultiGzDecoder;
    use std::fs::File;
    use std::io::BufReader;
    let file = File::open(path).unwrap_or_else(|e| panic!("open {path}: {e}"));
    let cap = 1 << 20;
    if path.ends_with(".gz") || path.ends_with(".bgz") {
        Box::new(BufReader::with_capacity(cap, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(cap, file))
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <store_path> <vcf_path> [threads] [--print]",
            args[0]
        );
        std::process::exit(1);
    }
    let store_path = &args[1];
    let vcf_path = &args[2];
    let default_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8);
    let threads: usize = args
        .get(3)
        .and_then(|s| if s.starts_with("--") { None } else { s.parse().ok() })
        .unwrap_or(default_threads);
    let do_print = args.iter().any(|a| a == "--print");
    let use_bgzf = args.iter().any(|a| a == "--bgzf");

    // ── Setup (untimed) ───────────────────────────────────────────────────
    eprintln!("Opening store: {store_path}");
    let t0 = Instant::now();
    let mut store = RefgetStore::open_local(store_path).expect("open_local");
    store.load_all_collections().expect("load_all_collections");
    let paged = store
        .list_collections(0, 10, &[])
        .expect("list_collections");
    let coll_digest = paged.results.first().expect("no collections").digest.clone();
    let collection = store.get_collection(&coll_digest).expect("get_collection");
    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    for rec in &collection.sequences {
        let m = rec.metadata();
        name_to_digest.insert(m.name.clone(), m.sha512t24u.clone());
    }
    eprintln!(
        "Store opened in {:.1}s ({} sequences in collection {})",
        t0.elapsed().as_secs_f64(),
        collection.sequences.len(),
        coll_digest
    );

    // Determine which sequences the VCF actually references, then make only
    // those resident (encoded bytes) before freezing to readonly.
    eprintln!("Scanning VCF for referenced chromosomes: {vcf_path}");
    let mut needed: HashMap<String, String> = HashMap::new();
    {
        use std::io::BufRead;
        let mut reader = open_vcf_reader(vcf_path);
        let mut line = String::new();
        loop {
            line.clear();
            match reader.read_line(&mut line) {
                Ok(0) => break,
                Ok(_) => {}
                Err(e) if e.kind() == std::io::ErrorKind::InvalidInput => break,
                Err(e) => panic!("read vcf line: {e}"),
            }
            let Some(rec) = parse_vcf_record(&line) else {
                continue;
            };
            if let Some(d) = name_to_digest.get(rec.chrom) {
                needed
                    .entry(rec.chrom.to_string())
                    .or_insert_with(|| d.clone());
            }
        }
    }
    eprintln!("VCF references {} known chromosome(s)", needed.len());

    for d in needed.values() {
        store
            .load_sequence(d)
            .unwrap_or_else(|e| panic!("load_sequence {d}: {e}"));
    }
    let store: ReadonlyRefgetStore = store.into_readonly();

    // ── Timed parallel pass ───────────────────────────────────────────────
    let path_name = if use_bgzf {
        "BGZF-block-parallel (workers own decompression)"
    } else {
        "single-reader parallel"
    };
    eprintln!("\nComputing VRS IDs with {threads} thread(s) via {path_name}...");
    if do_print {
        println!("chrom\tpos\tref\talt\tvrs_id");
    }
    let on_result = |r: gtars_vrs::vcf::VrsResult| {
        if do_print {
            println!(
                "{}\t{}\t{}\t{}\t{}",
                r.chrom, r.pos, r.ref_allele, r.alt_allele, r.vrs_id
            );
        }
    };
    let t1 = Instant::now();
    let count = if use_bgzf {
        compute_vrs_ids_parallel_bgzf_encoded_with_sink(&store, &needed, vcf_path, threads, on_result)
            .expect("compute_vrs_ids_parallel_bgzf_encoded_with_sink")
    } else {
        compute_vrs_ids_parallel_encoded(&store, &needed, vcf_path, threads, on_result)
            .expect("compute_vrs_ids_parallel_encoded")
    };
    let secs = t1.elapsed().as_secs_f64();

    eprintln!(
        "path={}\tthreads={threads}\tn_variants={count}\tcompute_s={secs:.3}\tvariants_per_sec={:.0}",
        if use_bgzf { "bgzf-block" } else { "single-reader" },
        count as f64 / secs
    );
    eprintln!(
        "peak_rss={:.1} MB\tTotal time: {:.1}s",
        peak_rss_kb() as f64 / 1024.0,
        t0.elapsed().as_secs_f64()
    );
}
