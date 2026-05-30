//! Focused, single-shot performance task for the perf regression suite.
//!
//! Unlike `bench_extract` (which runs ALL read strategies at once and is too
//! heavy for a tight regression suite), this binary runs exactly ONE
//! task+scenario+path and prints a SINGLE machine-readable RESULT line that
//! `perf/perf.py` parses. It is invoked under `/usr/bin/time -v` so the
//! harness can also capture peak RSS.
//!
//! Usage:
//!   perf_task --task extract --path partial --bed FILE --store ENC [--raw RAW]
//!   perf_task --task vrs     --path resident --points FILE --store ENC
//!
//! Tasks / paths:
//!   extract, path=resident : load_sequence whole touched seqs, get_substring()
//!   extract, path=partial  : never load_sequence; get_substring() reads only
//!                            the covering bytes from disk per region
//!   vrs,     path=resident : many 1bp point lookups over resident sequences
//!   vrs,     path=partial  : many 1bp point lookups, partial-read path
//!
//! Output (one line, on stdout):
//!   RESULT task=extract scenario=- path=partial seconds=0.612 items=10000 \
//!          bases=5507919716 throughput=8999869480.0 unit=bases_per_sec

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

/// Parse `chrom \t pos` (1-based POS, the clinvar_points.tsv shape) into 1bp
/// 0-based probes. Tolerates a BED-like third column too.
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
        if let Ok(p1) = p.parse::<usize>() {
            if p1 >= 1 {
                out.push((c.to_string(), p1 - 1, p1));
            }
        }
    }
    Ok(out)
}

fn unique_chroms(regions: &[(String, usize, usize)]) -> Vec<String> {
    let mut v: Vec<String> = regions.iter().map(|(c, _, _)| c.clone()).collect();
    v.sort();
    v.dedup();
    v
}

/// Resolve chrom -> seq digest, optionally calling load_sequence (resident).
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

/// Time a get_substring sweep over regions, returning (secs, total_bases).
fn sweep(
    store: &RefgetStore,
    regions: &[(String, usize, usize)],
    cd: &HashMap<String, String>,
) -> (f64, u64) {
    let t = Instant::now();
    let mut sink = 0u64;
    let mut bases = 0u64;
    for (chrom, s, e) in regions {
        if let Some(d) = cd.get(chrom) {
            let sub = store.get_substring(d, *s, *e).unwrap();
            bases += sub.len() as u64;
            sink = sink
                .wrapping_add(sub.len() as u64)
                .wrapping_add(sub.as_bytes()[0] as u64);
        }
    }
    let secs = t.elapsed().as_secs_f64();
    // keep the optimizer honest
    if secs < 0.0 {
        eprintln!("{}", sink);
    }
    (secs, bases)
}

struct Args {
    task: String,
    path: String,
    scenario: String,
    store: String,
    raw: Option<String>,
    bed: Option<String>,
    points: Option<String>,
}

fn parse_args() -> Args {
    let mut a = Args {
        task: String::new(),
        path: "resident".into(),
        scenario: "-".into(),
        store: String::new(),
        raw: None,
        bed: None,
        points: None,
    };
    let argv: Vec<String> = std::env::args().collect();
    let mut i = 1;
    while i < argv.len() {
        let k = argv[i].as_str();
        let mut next = || {
            i += 1;
            argv.get(i).cloned().unwrap_or_default()
        };
        match k {
            "--task" => a.task = next(),
            "--path" => a.path = next(),
            "--scenario" => a.scenario = next(),
            "--store" => a.store = next(),
            "--raw" => a.raw = Some(next()),
            "--bed" => a.bed = Some(next()),
            "--points" => a.points = Some(next()),
            other => {
                eprintln!("unknown arg: {other}");
                std::process::exit(2);
            }
        }
        i += 1;
    }
    if a.task.is_empty() || a.store.is_empty() {
        eprintln!("usage: perf_task --task extract|vrs --path resident|partial --store ENC [--raw RAW] [--bed F | --points F] [--scenario NAME]");
        std::process::exit(2);
    }
    a
}

fn open_with_collection(path: &str) -> anyhow::Result<(RefgetStore, String)> {
    let mut store = RefgetStore::open_local(path)?;
    store.set_quiet(true);
    store.load_all_collections()?;
    let coll_digest = store.list_collections(0, 10_000, &[])?.results[0]
        .digest
        .clone();
    Ok((store, coll_digest))
}

fn main() -> anyhow::Result<()> {
    let a = parse_args();

    let (items, bases, seconds, unit): (u64, u64, f64, &str) = match a.task.as_str() {
        "extract" => {
            let bed = a.bed.as_ref().expect("extract task requires --bed");
            let regions = read_bed(bed)?;
            let unique = unique_chroms(&regions);
            let resident = a.path == "resident";
            // strategy C uses the raw store; otherwise the encoded store.
            let store_path = if a.path == "raw" {
                a.raw.as_ref().expect("raw path requires --raw")
            } else {
                &a.store
            };
            let (mut store, coll) = open_with_collection(store_path)?;
            let do_load = resident || a.path == "raw";
            let cd = resolve(&mut store, &coll, &unique, do_load)?;
            let regions: Vec<_> = regions
                .iter()
                .filter(|(c, _, _)| cd.contains_key(c))
                .cloned()
                .collect();
            let (secs, bases) = sweep(&store, &regions, &cd);
            (regions.len() as u64, bases, secs, "bases_per_sec")
        }
        "vrs" => {
            let pts = a.points.as_ref().expect("vrs task requires --points");
            let points = read_points(pts)?;
            let unique = unique_chroms(&points);
            let (mut store, coll) = open_with_collection(&a.store)?;
            let do_load = a.path != "partial";
            let cd = resolve(&mut store, &coll, &unique, do_load)?;
            let points: Vec<_> = points
                .iter()
                .filter(|(c, _, _)| cd.contains_key(c))
                .cloned()
                .collect();
            let (secs, _bases) = sweep(&store, &points, &cd);
            (points.len() as u64, points.len() as u64, secs, "lookups_per_sec")
        }
        other => {
            eprintln!("unknown task: {other}");
            std::process::exit(2);
        }
    };

    let throughput = if seconds > 0.0 {
        match unit {
            "lookups_per_sec" => items as f64 / seconds,
            _ => bases as f64 / seconds,
        }
    } else {
        0.0
    };

    println!(
        "RESULT task={} scenario={} path={} seconds={:.6} items={} bases={} throughput={:.1} unit={}",
        a.task, a.scenario, a.path, seconds, items, bases, throughput, unit
    );
    Ok(())
}
