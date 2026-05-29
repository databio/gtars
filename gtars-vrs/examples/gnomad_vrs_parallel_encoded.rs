//! Compute VRS Allele IDs from a (gnomAD) VCF using the immutable 2-bit-encoded
//! refget store with in-process parallel workers (decode-on-the-fly).
//!
//! Mirrors the store setup of `bench_encoded.rs`: open a local store, load all
//! collections, build the name->digest map, make the needed sequences resident
//! with `load_sequence`, freeze with `into_readonly`, then call the parallel
//! encoded path. Prints throughput (variants/sec) and the thread count.
//!
//! Usage:
//!   cargo run --release --example gnomad_vrs_parallel_encoded -- <store_path> <vcf_path> [threads] [--print]
//!     threads: worker count (default: available_parallelism)
//!     --print: also print each result row (chrom pos ref alt vrs_id) to stdout

use std::collections::HashMap;
use std::time::Instant;

use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore};
use gtars_vrs::vcf::compute_vrs_ids_parallel_encoded;

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
            let l = line.trim_end_matches(['\n', '\r']);
            if l.is_empty() || l.starts_with('#') {
                continue;
            }
            let chrom = l.split('\t').next().unwrap_or("");
            if let Some(d) = name_to_digest.get(chrom) {
                needed.entry(chrom.to_string()).or_insert_with(|| d.clone());
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
    eprintln!("\nComputing VRS IDs with {threads} thread(s)...");
    if do_print {
        println!("chrom\tpos\tref\talt\tvrs_id");
    }
    let t1 = Instant::now();
    let count = compute_vrs_ids_parallel_encoded(
        &store,
        &needed,
        vcf_path,
        threads,
        |r| {
            if do_print {
                println!(
                    "{}\t{}\t{}\t{}\t{}",
                    r.chrom, r.pos, r.ref_allele, r.alt_allele, r.vrs_id
                );
            }
        },
    )
    .expect("compute_vrs_ids_parallel_encoded");
    let secs = t1.elapsed().as_secs_f64();

    eprintln!(
        "threads={threads}\tn_variants={count}\tcompute_s={secs:.3}\tvariants_per_sec={:.0}",
        count as f64 / secs
    );
    eprintln!("Total time: {:.1}s", t0.elapsed().as_secs_f64());
}
