//! Catalog memory and lookup-latency harness for many-sequence stores.
//!
//! Builds (once) a synthetic packed store with N unique sequences spread over
//! K collections, then opens it the way the seqcol server does
//! (`open_local` + `load_all_collections` + `into_readonly`) and reports:
//!
//! - RSS delta of the open (from `/proc/self/statm`)
//! - open time
//! - 100k `get_sequence_by_name`, 100k sha lookups, 100k md5 lookups
//! - one `get_collection_level2` of the largest collection
//!
//! Each collection holds a shared pool of sequences (in every collection) plus
//! its own unique sequences, with `share` setting the shared fraction of each
//! collection. Defaults: 1M sequences, 10 collections, 50% shared, 20-char names.
//!
//! Usage:
//!   catalog_memory --dir /path/to/store [--n 1000000] [--k 10] [--share 0.5]
//!
//! The store is built only if `--dir` does not exist yet. Build and measure in
//! separate runs if you want a clean RSS number: `--build-only`, then again.

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use gtars_refget::digest::types::digest_sequence;
use gtars_refget::store::{FastaImportOptions, RefgetStore};

fn rss_bytes() -> u64 {
    let s = std::fs::read_to_string("/proc/self/statm").unwrap_or_default();
    let pages: u64 = s.split_whitespace().nth(1).and_then(|x| x.parse().ok()).unwrap_or(0);
    pages * 4096
}

/// Deterministic, unique 64-bp sequence for index `i`.
fn seq_for(i: usize) -> Vec<u8> {
    let mut x = (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xD1B5_4A32_D192_ED03;
    let mut out = Vec::with_capacity(64);
    // First 32 bases encode i exactly (2 bits each), so sequences are unique.
    let mut v = i as u64;
    for _ in 0..32 {
        out.push(b"ACGT"[(v & 3) as usize]);
        v >>= 2;
    }
    for _ in 0..32 {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        out.push(b"ACGT"[(x & 3) as usize]);
    }
    out
}

fn name_for(i: usize) -> String {
    format!("scaffold_{:011}", i) // 20 chars
}

struct Layout {
    shared: usize,
    unique_per: usize,
}

fn layout(n: usize, k: usize, share: f64) -> Layout {
    // n = shared + k * unique_per; shared / (shared + unique_per) = share
    let ratio = share / (1.0 - share); // shared = ratio * unique_per
    let unique_per = (n as f64 / (ratio + k as f64)).floor() as usize;
    let shared = n - k * unique_per;
    Layout { shared, unique_per }
}

fn build(dir: &Path, n: usize, k: usize, share: f64) -> anyhow::Result<()> {
    let l = layout(n, k, share);
    eprintln!(
        "building store: {} unique seqs, {} collections of {} members ({} shared)",
        n,
        k,
        l.shared + l.unique_per,
        l.shared
    );
    let fa_dir = dir.with_extension("fastas");
    std::fs::create_dir_all(&fa_dir)?;
    let mut fastas = Vec::new();
    for c in 0..k {
        let p = fa_dir.join(format!("c{c}.fa"));
        let mut w = std::io::BufWriter::new(std::fs::File::create(&p)?);
        let unique_start = l.shared + c * l.unique_per;
        for i in (0..l.shared).chain(unique_start..unique_start + l.unique_per) {
            writeln!(w, ">{}", name_for(i))?;
            w.write_all(&seq_for(i))?;
            w.write_all(b"\n")?;
        }
        w.flush()?;
        fastas.push(p);
    }
    let t = Instant::now();
    let mut store = RefgetStore::on_disk_packed(dir)?;
    store.set_quiet(true);
    for p in &fastas {
        store.add_sequence_collection_from_fasta(p, FastaImportOptions::new())?;
    }
    store.write()?;
    eprintln!("built in {:.1}s", t.elapsed().as_secs_f64());
    std::fs::remove_dir_all(&fa_dir)?;
    Ok(())
}

fn main() -> anyhow::Result<()> {
    let mut args: HashMap<String, String> = HashMap::new();
    let mut flags = Vec::new();
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let a = argv[i].trim_start_matches("--").to_string();
        if a == "build-only" {
            flags.push(a);
            i += 1;
        } else {
            args.insert(a, argv.get(i + 1).cloned().unwrap_or_default());
            i += 2;
        }
    }
    let dir = PathBuf::from(args.get("dir").expect("--dir required"));
    let n: usize = args.get("n").map(|s| s.parse().unwrap()).unwrap_or(1_000_000);
    let k: usize = args.get("k").map(|s| s.parse().unwrap()).unwrap_or(10);
    let share: f64 = args.get("share").map(|s| s.parse().unwrap()).unwrap_or(0.5);
    let iters: usize = args.get("iters").map(|s| s.parse().unwrap()).unwrap_or(100_000);

    if !dir.exists() {
        build(&dir, n, k, share)?;
    }
    if flags.iter().any(|f| f == "build-only") {
        return Ok(());
    }
    let l = layout(n, k, share);

    // Compute probe keys BEFORE measuring RSS so they don't count.
    // sha/md5 probes span all sequences; name probes use the shared pool, which
    // every collection holds under the same names.
    let probe: Vec<usize> = (0..iters).map(|j| (j * 7919) % n).collect();
    let probe_digests: Vec<(String, String)> = probe
        .iter()
        .map(|&i| {
            let m = digest_sequence("x", &seq_for(i));
            let m = m.metadata();
            (m.sha512t24u.clone(), m.md5.clone())
        })
        .collect();
    let probe_names: Vec<String> =
        (0..iters).map(|j| name_for((j * 7919) % l.shared.max(1))).collect();

    let rss0 = rss_bytes();
    let t = Instant::now();
    let mut store = RefgetStore::open_local(&dir)?;
    store.set_quiet(true);
    let open_local_s = t.elapsed().as_secs_f64();
    store.load_all_collections()?;
    let ro = store.into_readonly();
    let open_s = t.elapsed().as_secs_f64();
    let rss1 = rss_bytes();
    let delta = rss1.saturating_sub(rss0);

    let colls = ro.list_collections(0, usize::MAX, &[])?;
    let first = colls.results[0].digest.clone();

    let t = Instant::now();
    let mut hits = 0usize;
    for name in &probe_names {
        if ro.get_sequence_by_name(&first, name).is_ok() {
            hits += 1;
        }
    }
    let by_name_s = t.elapsed().as_secs_f64();
    assert_eq!(hits, iters, "name lookups missed");

    let t = Instant::now();
    let mut hits = 0usize;
    for (sha, _) in &probe_digests {
        if ro.get_sequence_metadata(sha).is_some() {
            hits += 1;
        }
    }
    let by_sha_s = t.elapsed().as_secs_f64();
    assert_eq!(hits, iters, "sha lookups missed");

    let t = Instant::now();
    let mut hits = 0usize;
    for (_, md5) in &probe_digests {
        if ro.get_sequence(md5).is_ok() {
            hits += 1;
        }
    }
    let by_md5_s = t.elapsed().as_secs_f64();
    assert_eq!(hits, iters, "md5 lookups missed");

    // Borrowed-view lookups (no owned metadata built per call).
    let t = Instant::now();
    let mut hits = 0usize;
    for (sha, _) in &probe_digests {
        if ro.sequence_view(sha).is_some() {
            hits += 1;
        }
    }
    let view_sha_s = t.elapsed().as_secs_f64();
    assert!(hits > 0);
    let t = Instant::now();
    let mut hits = 0usize;
    for (_, md5) in &probe_digests {
        if ro.sequence_view(md5).is_some() {
            hits += 1;
        }
    }
    let view_md5_s = t.elapsed().as_secs_f64();
    assert!(hits > 0);
    let t = Instant::now();
    let mut hits = 0usize;
    for name in &probe_names {
        if let Ok(m) = ro.member_by_name(&first, name) {
            hits += m.seq().length();
        }
    }
    let view_name_s = t.elapsed().as_secs_f64();
    assert!(hits > 0);

    let t = Instant::now();
    let lvl2 = ro.get_collection_level2(&first)?;
    let lvl2_s = t.elapsed().as_secs_f64();

    // First packed read loads the local pack location table; measure its RSS.
    let rss2 = rss_bytes();
    let _ = ro.get_substring(&probe_digests[0].0, 0, 10)?;
    let pack_delta = rss_bytes().saturating_sub(rss2);

    let total_members = k * (l.shared + l.unique_per);
    println!(
        "RESULT n={} k={} members={} rss_delta_mb={:.1} bytes_per_seq={:.0} open_s={:.2} \
         open_local_s={:.2} view_name_100k_ms={:.1} view_sha_100k_ms={:.1} view_md5_100k_ms={:.1} \
         by_name_100k_ms={:.1} by_sha_100k_ms={:.1} by_md5_100k_ms={:.1} level2_ms={:.1} level2_n={} pack_index_mb={:.1}",
        n,
        k,
        total_members,
        delta as f64 / 1e6,
        delta as f64 / n as f64,
        open_s,
        open_local_s,
        view_name_s * 1e3 * 100_000.0 / iters as f64,
        view_sha_s * 1e3 * 100_000.0 / iters as f64,
        view_md5_s * 1e3 * 100_000.0 / iters as f64,
        by_name_s * 1e3 * 100_000.0 / iters as f64,
        by_sha_s * 1e3 * 100_000.0 / iters as f64,
        by_md5_s * 1e3 * 100_000.0 / iters as f64,
        lvl2_s * 1e3,
        lvl2.names.len(),
        pack_delta as f64 / 1e6,
    );
    Ok(())
}
