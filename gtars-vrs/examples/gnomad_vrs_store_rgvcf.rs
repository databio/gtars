//! Pure-Rust rgvcf contestant for the VRS-ID bench. Drop-in for
//! `gnomad_vrs_store` but takes an rgvcf file instead of a BGZF VCF.
//!
//! Arguments mirror `gnomad_vrs_store`'s positional interface:
//!
//!   gnomad_vrs_store_rgvcf <store_dir> <rgvcf_path> <results_tsv>
//!       <variants_tsv> [num_workers]
//!
//! Emits the same single-row `results_tsv` so the tourney wrapper doesn't
//! care which contestant produced it; contestant string is `rust_rgvcf`.

#[cfg(target_os = "linux")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use gtars_refget::store::RefgetStore;
use gtars_rgvcf::RgvcfReader;
use gtars_vrs::rgvcf_input::compute_vrs_ids_from_rgvcf_with_sink;

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

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if !(5..=6).contains(&args.len()) {
        eprintln!(
            "Usage: {} <store_dir> <rgvcf_path> <results_tsv> <variants_tsv> [num_workers]",
            args[0]
        );
        std::process::exit(2);
    }
    let store_dir = &args[1];
    let rgvcf_path = &args[2];
    let results_tsv = &args[3];
    let variants_tsv = &args[4];
    let num_workers_arg: usize = if args.len() == 6 { args[5].parse().unwrap_or(0) } else { 0 };
    let effective_workers = if num_workers_arg == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get().saturating_sub(2).max(1))
            .unwrap_or(1)
    } else {
        num_workers_arg
    };

    let t0 = Instant::now();
    let mut store = RefgetStore::open_local(store_dir).context("open store")?;
    let t_open = t0.elapsed();

    let paged = store.list_collections(0, 10, &[]).context("list collections")?;
    let collections = &paged.results;
    if collections.is_empty() {
        bail!("no collections");
    }
    let collection_digest = collections[0].digest.clone();
    store
        .load_collection(collection_digest.as_str())
        .context("load_collection")?;

    // Read rgvcf header (cheap) to learn which chromosomes to decode.
    let zero_digest = "\0".repeat(32);
    let t1 = Instant::now();
    {
        let reader = RgvcfReader::open(rgvcf_path)?;
        for c in reader.chromosomes() {
            if c.seq_digest == zero_digest {
                eprintln!("[rust_rgvcf] WARN chrom {} zero digest; skip", c.name);
                continue;
            }
            store.load_sequence(c.seq_digest.as_str())?;
            store.ensure_decoded(c.seq_digest.as_str())?;
        }
    }
    let t_preload = t1.elapsed();

    let readonly = store.into_readonly();
    // Build name_to_digest from the readonly store for the entry point's
    // validation. Names come from the collection; we use the encoded
    // rgvcf's declared chrom names to map to raw digests.
    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    let reader2 = RgvcfReader::open(rgvcf_path)?;
    for c in reader2.chromosomes() {
        if c.seq_digest == zero_digest {
            continue;
        }
        name_to_digest.insert(c.name.clone(), c.seq_digest.clone());
    }

    let vf = File::create(variants_tsv)?;
    let mut vw = BufWriter::with_capacity(1 << 20, vf);
    writeln!(vw, "chrom\tpos\tref\talt\tvrs_id")?;

    let t2 = Instant::now();
    let n = compute_vrs_ids_from_rgvcf_with_sink(
        &readonly,
        &name_to_digest,
        rgvcf_path,
        effective_workers,
        |r| {
            let _ = writeln!(
                vw,
                "{}\t{}\t{}\t{}\t{}",
                r.chrom, r.pos, r.ref_allele, r.alt_allele, r.vrs_id
            );
        },
    )?;
    let t_compute = t2.elapsed();
    vw.flush()?;

    let open_s = t_open.as_secs_f64();
    let preload_s = t_preload.as_secs_f64();
    let compute_s = t_compute.as_secs_f64();
    let total_s = open_s + preload_s + compute_s;
    let vps = if compute_s > 0.0 { n as f64 / compute_s } else { 0.0 };
    let rss = max_rss_mb();

    let vcf_basename = Path::new(rgvcf_path)
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| rgvcf_path.to_string());

    let rf = File::create(results_tsv)?;
    let mut rw = BufWriter::new(rf);
    writeln!(
        rw,
        "contestant\tvcf\tn_variants\tnum_workers\topen_s\tpreload_s\tcompute_s\ttotal_s\tvariants_per_sec\tmax_rss_mb"
    )?;
    writeln!(
        rw,
        "rust_rgvcf\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.0}\t{:.1}",
        vcf_basename, n, effective_workers, open_s, preload_s, compute_s, total_s, vps, rss
    )?;
    rw.flush()?;

    eprintln!(
        "[rust_rgvcf] n={} workers={} open={:.2}s preload={:.2}s compute={:.2}s total={:.2}s  {:.0} v/s  rss={:.1} MB  (wrote {})",
        n, effective_workers, open_s, preload_s, compute_s, total_s, vps, rss, variants_tsv
    );
    Ok(())
}
