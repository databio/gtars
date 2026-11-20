//! Open a store the way the seqcol server does (`open_local` +
//! `load_all_collections` + `into_readonly`) and report the RSS it costs.
//!
//! Usage: open_rss <store_dir>
use std::time::Instant;

use gtars_refget::store::RefgetStore;

fn rss_bytes() -> u64 {
    let s = std::fs::read_to_string("/proc/self/statm").unwrap_or_default();
    let pages: u64 = s.split_whitespace().nth(1).and_then(|x| x.parse().ok()).unwrap_or(0);
    pages * 4096
}

fn main() -> anyhow::Result<()> {
    let dir = std::env::args().nth(1).expect("usage: open_rss <store_dir>");
    let rss0 = rss_bytes();
    let t = Instant::now();
    let mut store = RefgetStore::open_local(&dir)?;
    store.set_quiet(true);
    store.load_all_collections()?;
    let ro = store.into_readonly();
    let open_s = t.elapsed().as_secs_f64();
    let delta = rss_bytes().saturating_sub(rss0);
    let stats = ro.stats();
    println!(
        "RESULT n_sequences={} n_collections={} rss_delta_mb={:.1} bytes_per_seq={:.0} open_s={:.2}",
        stats.n_sequences,
        stats.n_collections,
        delta as f64 / 1e6,
        delta as f64 / stats.n_sequences.max(1) as f64,
        open_s
    );
    Ok(())
}
