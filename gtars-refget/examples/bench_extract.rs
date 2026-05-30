//! Benchmark: extraction read-path strategies on a RefgetStore.
//!
//! Isolates the cost of the on-the-fly decode in the extract path versus
//! decoding each touched sequence once (resident) and slicing, versus a
//! Raw-mode store, versus the new partial-read path that fetches only the
//! bytes covering each queried region straight from the `.seq` file. Also
//! exercises the VRS/HGVS-style access pattern (many tiny point lookups) from
//! an optional variants file. No Python, no FASTA writing -- pure in-process
//! extraction over real regions.
//!
//! Usage:
//!   cargo run --release --example bench_extract -- \
//!       <encoded_store> <bed> [raw_store] [variants_vcf_or_bed]
//!
//! Strategies (all run after collection metadata is loaded):
//!   A  on-the-fly    : load whole touched seqs, get_substring() decodes each span
//!   B  decode-once   : load whole touched seqs, decode each fully once, then slice
//!   C  raw-resident  : load whole touched seqs from a Raw store, direct slice
//!   D  partial-read  : NEVER load_sequence; get_substring() on a stub reads only
//!                      the covering bytes from disk per query
//!   V  vrs-pattern   : many 1bp point lookups (the VRS/HGVS access shape), run
//!                      under both the resident path and the partial-read path

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;

use gtars_refget::store::RefgetStore;

fn read_bed(path: &str) -> anyhow::Result<Vec<(String, usize, usize)>> {
    let f = File::open(path)?;
    let mut out = Vec::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split('\t');
        let (Some(c), Some(s), Some(e)) = (it.next(), it.next(), it.next()) else {
            continue;
        };
        let (Ok(s), Ok(e)) = (s.parse::<usize>(), e.parse::<usize>()) else {
            continue;
        };
        if e > s {
            out.push((c.to_string(), s, e));
        }
    }
    Ok(out)
}

/// Parse a variants source as 1bp point lookups (chrom, pos0, pos0+1).
/// Accepts VCF (chrom \t pos1 ...) or BED (chrom \t start \t end); for VCF we
/// use the 1-based POS converted to a 0-based single-base window.
fn read_points(path: &str) -> anyhow::Result<Vec<(String, usize, usize)>> {
    let f = File::open(path)?;
    let mut out = Vec::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split('\t');
        let (Some(c), Some(p)) = (it.next(), it.next()) else {
            continue;
        };
        let third = it.next();
        // VCF: second column is 1-based POS. BED: second is 0-based start, third is end.
        if let (Ok(p), Some(Ok(e))) = (p.parse::<usize>(), third.map(|t| t.parse::<usize>())) {
            // looks like BED (has a valid integer third column): take a 1bp probe at start
            let _ = e;
            out.push((c.to_string(), p, p + 1));
        } else if let Ok(p1) = p.parse::<usize>() {
            // VCF POS (1-based) -> 0-based single base
            if p1 >= 1 {
                out.push((c.to_string(), p1 - 1, p1));
            }
        }
    }
    Ok(out)
}

/// Resolve chrom -> seq digest for every chrom referenced, optionally calling
/// load_sequence (resident) or not (stub / partial-read).
fn resolve(
    store: &mut RefgetStore,
    coll_digest: &str,
    unique: &[String],
    do_load: bool,
) -> anyhow::Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    for chrom in unique {
        let digest = match store.get_sequence_by_name(coll_digest, chrom) {
            Ok(rec) => rec.metadata().sha512t24u.clone(),
            Err(_) => continue,
        };
        if do_load {
            store.load_sequence(&digest)?;
        }
        map.insert(chrom.clone(), digest);
    }
    Ok(map)
}

fn unique_chroms(regions: &[(String, usize, usize)]) -> Vec<String> {
    let mut v: Vec<String> = regions.iter().map(|(c, _, _)| c.clone()).collect();
    v.sort();
    v.dedup();
    v
}

/// Time a get_substring sweep over regions, returning (secs, sink).
fn sweep(store: &RefgetStore, regions: &[(String, usize, usize)], cd: &HashMap<String, String>) -> (f64, u64) {
    let t = Instant::now();
    let mut sink = 0u64;
    for (chrom, s, e) in regions {
        if let Some(d) = cd.get(chrom) {
            let sub = store.get_substring(d, *s, *e).unwrap();
            sink = sink
                .wrapping_add(sub.len() as u64)
                .wrapping_add(sub.as_bytes()[0] as u64);
        }
    }
    (t.elapsed().as_secs_f64(), sink)
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: bench_extract <encoded_store> <bed> [raw_store] [variants]");
        std::process::exit(1);
    }
    let enc_path = &args[1];
    let bed_path = &args[2];
    let raw_path = args.get(3);
    let variants_path = args.get(4);

    let regions = read_bed(bed_path)?;
    let total_width: usize = regions.iter().map(|(_, s, e)| e - s).sum();
    let unique = unique_chroms(&regions);
    println!(
        "BED regions: {}   total width: {} Mbp   unique chroms: {}",
        regions.len(),
        total_width / 1_000_000,
        unique.len()
    );

    // ---- ENCODED store: resident strategies A and B ----
    {
        let mut store = RefgetStore::open_local(enc_path)?;
        store.set_quiet(true);
        store.load_all_collections()?;
        let coll_digest = store.list_collections(0, 10_000, &[])?.results[0].digest.clone();

        let t = Instant::now();
        let cd = resolve(&mut store, &coll_digest, &unique, true)?;
        let load = t.elapsed().as_secs_f64();
        let regions: Vec<_> = regions.iter().filter(|(c, _, _)| cd.contains_key(c)).cloned().collect();
        println!("\nload_sequence (whole, {} seqs): {:.3}s", cd.len(), load);

        // A: on-the-fly
        let (a, sink) = sweep(&store, &regions, &cd);
        println!("[A] on-the-fly  (resident): {:.3}s  ({:.2} us/region)  sink={}", a, a * 1e6 / regions.len() as f64, sink);

        // B: decode whole once, slice
        let t = Instant::now();
        let mut decoded: HashMap<String, String> = HashMap::new();
        let mut bases = 0usize;
        for (chrom, d) in &cd {
            let len = store.get_sequence_metadata(d).unwrap().length;
            let full = store.get_substring(d, 0, len)?;
            bases += full.len();
            decoded.insert(chrom.clone(), full);
        }
        let b_dec = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let mut sink_b = 0u64;
        for (chrom, s, e) in &regions {
            let sl = &decoded[chrom].as_bytes()[*s..*e];
            sink_b = sink_b.wrapping_add(sl.len() as u64).wrapping_add(sl[0] as u64);
        }
        let b_sl = t.elapsed().as_secs_f64();
        println!("[B] decode-once (resident): {:.3}s  (decode-whole {:.3}s of {} Mbp + slice {:.3}s)  sink={}", b_dec + b_sl, b_dec, bases / 1_000_000, b_sl, sink_b);
    }

    // ---- ENCODED store: partial-read strategy D (no load_sequence) ----
    {
        let mut store = RefgetStore::open_local(enc_path)?;
        store.set_quiet(true);
        store.load_all_collections()?;
        let coll_digest = store.list_collections(0, 10_000, &[])?.results[0].digest.clone();
        let cd = resolve(&mut store, &coll_digest, &unique, false)?; // stubs only
        let regions: Vec<_> = regions.iter().filter(|(c, _, _)| cd.contains_key(c)).cloned().collect();
        let (d, sink) = sweep(&store, &regions, &cd);
        println!("[D] partial-read (no load): {:.3}s  ({:.2} us/region)  sink={}", d, d * 1e6 / regions.len() as f64, sink);
    }

    // ---- RAW store: strategy C ----
    if let Some(raw_path) = raw_path {
        let mut store = RefgetStore::open_local(raw_path)?;
        store.set_quiet(true);
        store.load_all_collections()?;
        let coll_digest = store.list_collections(0, 10_000, &[])?.results[0].digest.clone();
        let t = Instant::now();
        let cd = resolve(&mut store, &coll_digest, &unique, true)?;
        let load = t.elapsed().as_secs_f64();
        let regions: Vec<_> = regions.iter().filter(|(c, _, _)| cd.contains_key(c)).cloned().collect();
        let (c, sink) = sweep(&store, &regions, &cd);
        println!("\n[raw] load_sequence (whole): {:.3}s", load);
        println!("[C] raw-resident          : {:.3}s  ({:.2} us/region)  sink={}", c, c * 1e6 / regions.len() as f64, sink);
    }

    // ---- VRS/HGVS pattern: many 1bp point lookups ----
    if let Some(vp) = variants_path {
        let points = read_points(vp)?;
        let upoints = unique_chroms(&points);
        let total_pts = points.len();
        println!("\n=== VRS/HGVS pattern: {} point lookups (1bp), {} chroms ===", total_pts, upoints.len());

        // resident
        {
            let mut store = RefgetStore::open_local(enc_path)?;
            store.set_quiet(true);
            store.load_all_collections()?;
            let coll_digest = store.list_collections(0, 10_000, &[])?.results[0].digest.clone();
            let t = Instant::now();
            let cd = resolve(&mut store, &coll_digest, &upoints, true)?;
            let load = t.elapsed().as_secs_f64();
            let points: Vec<_> = points.iter().filter(|(c, _, _)| cd.contains_key(c)).cloned().collect();
            let (v, sink) = sweep(&store, &points, &cd);
            println!("[V-resident] load {:.3}s + lookups {:.3}s  ({:.2} us/lookup)  sink={}", load, v, v * 1e6 / points.len() as f64, sink);
        }
        // partial-read
        {
            let mut store = RefgetStore::open_local(enc_path)?;
            store.set_quiet(true);
            store.load_all_collections()?;
            let coll_digest = store.list_collections(0, 10_000, &[])?.results[0].digest.clone();
            let cd = resolve(&mut store, &coll_digest, &upoints, false)?;
            let points: Vec<_> = points.iter().filter(|(c, _, _)| cd.contains_key(c)).cloned().collect();
            let (v, sink) = sweep(&store, &points, &cd);
            println!("[V-partial ] no load + lookups {:.3}s  ({:.2} us/lookup)  sink={}", v, v * 1e6 / points.len() as f64, sink);
        }
    }

    Ok(())
}
