//! A coalescing block cache over pack files.
//!
//! Pack files are read through a bounded LRU of fixed-size, block-aligned chunks
//! shared across every pack of a store. This cache — not the pack layout alone —
//! is what collapses a 10K-region batch extraction from ~10K positioned reads to
//! a few dozen: overlapping and adjacent regions land in the same 1 MiB block,
//! and a run of missing blocks is fetched with ONE call to the backend (a single
//! `pread` locally, a single HTTP `Range` remotely).

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use anyhow::{anyhow, Result};

/// Default aligned block size: 1 MiB.
pub(crate) const PACK_BLOCK_SIZE: u64 = 1 << 20;

/// Default cache capacity in blocks (128 MiB at the default block size).
pub(crate) const DEFAULT_PACK_CACHE_BLOCKS: usize = 128;

/// A bounded LRU of block-aligned pack chunks.
///
/// Keyed by `(pack name, block index)`. Blocks are immutable (`Arc<Vec<u8>>`),
/// so a cached block never goes stale and concurrent inserters of the same block
/// are harmless. Recency is a monotonically increasing stamp; eviction scans for
/// the minimum stamp (cheap at the default cap, like the fd cache's LRU).
/// Key: `(pack name, block index)`. Value: `(block bytes, last-use stamp)`.
type BlockKey = (Arc<str>, u64);
type BlockSlot = (Arc<Vec<u8>>, u64);

#[derive(Debug)]
struct Lru {
    cap: usize,
    tick: u64,
    entries: HashMap<BlockKey, BlockSlot>,
}

impl Lru {
    fn new(cap: usize) -> Self {
        Lru {
            cap: cap.max(1),
            tick: 0,
            entries: HashMap::new(),
        }
    }

    fn get(&mut self, key: &(Arc<str>, u64)) -> Option<Arc<Vec<u8>>> {
        self.tick += 1;
        let tick = self.tick;
        if let Some(slot) = self.entries.get_mut(key) {
            slot.1 = tick;
            Some(Arc::clone(&slot.0))
        } else {
            None
        }
    }

    fn insert(&mut self, key: (Arc<str>, u64), block: Arc<Vec<u8>>) {
        self.tick += 1;
        let tick = self.tick;
        if !self.entries.contains_key(&key)
            && self.entries.len() >= self.cap
            && let Some(lru_key) = self
                .entries
                .iter()
                .min_by_key(|(_, (_, stamp))| *stamp)
                .map(|(k, _)| k.clone())
        {
            self.entries.remove(&lru_key);
        }
        self.entries.insert(key, (block, tick));
    }
}

/// The coalescing block cache. Behind a `Mutex` because reads take `&self`.
#[derive(Debug)]
pub(crate) struct PackBlockCache {
    lru: Mutex<Lru>,
    block_size: u64,
}

impl PackBlockCache {
    pub(crate) fn new(cap_blocks: usize, block_size: u64) -> Self {
        PackBlockCache {
            lru: Mutex::new(Lru::new(cap_blocks)),
            block_size: block_size.max(1),
        }
    }

    /// Test-only: number of distinct (pack, block) entries currently cached.
    #[cfg(all(test, feature = "filesystem"))]
    pub(crate) fn cached_block_count(&self) -> usize {
        self.lru.lock().unwrap().entries.len()
    }

    /// Read bytes `[start, end)` of `pack` (whose total size is `pack_len`).
    ///
    /// Missing blocks are identified under the lock; each CONTIGUOUS RUN of
    /// missing blocks is fetched with ONE call to `fetch(byte_off, byte_len)`.
    /// The fetched bytes are split into blocks and inserted, then the requested
    /// window is assembled from the covered slices.
    ///
    /// A request spanning more bytes than the cache can hold BYPASSES the cache
    /// (one direct fetch), so a single large read cannot evict everything.
    pub(crate) fn read_range(
        &self,
        pack: &Arc<str>,
        pack_len: u64,
        start: u64,
        end: u64,
        fetch: &dyn Fn(u64, u64) -> Result<Vec<u8>>,
    ) -> Result<Vec<u8>> {
        if start > end || end > pack_len {
            return Err(anyhow!(
                "pack read range [{}, {}) out of bounds for pack len {}",
                start,
                end,
                pack_len
            ));
        }
        if start == end {
            return Ok(Vec::new());
        }

        let bs = self.block_size;
        let first_block = start / bs;
        let last_block = (end - 1) / bs;
        let n_blocks = last_block - first_block + 1;

        // Bypass: a span wider than the cache capacity would thrash it.
        let cap = { self.lru.lock().unwrap().cap as u64 };
        if n_blocks > cap {
            return fetch(start, end - start);
        }

        // 1. Under the lock: collect hits, note which blocks are missing.
        let mut blocks: HashMap<u64, Arc<Vec<u8>>> = HashMap::new();
        let mut missing: Vec<u64> = Vec::new();
        {
            let mut lru = self.lru.lock().unwrap();
            for b in first_block..=last_block {
                let key = (Arc::clone(pack), b);
                if let Some(block) = lru.get(&key) {
                    blocks.insert(b, block);
                } else {
                    missing.push(b);
                }
            }
        }

        // 2. Outside the lock: fetch each contiguous run of missing blocks once.
        let runs = coalesce_runs(&missing);
        let mut fetched: Vec<(u64, Arc<Vec<u8>>)> = Vec::new();
        for (run_start_block, run_end_block) in runs {
            let byte_off = run_start_block * bs;
            let byte_end = ((run_end_block + 1) * bs).min(pack_len);
            let data = fetch(byte_off, byte_end - byte_off)?;
            // Split into blocks.
            for b in run_start_block..=run_end_block {
                let block_start = (b * bs - byte_off) as usize;
                if block_start >= data.len() {
                    break;
                }
                let block_end = (block_start + bs as usize).min(data.len());
                let block = Arc::new(data[block_start..block_end].to_vec());
                blocks.insert(b, Arc::clone(&block));
                fetched.push((b, block));
            }
        }

        // 3. Under the lock: insert freshly fetched blocks.
        if !fetched.is_empty() {
            let mut lru = self.lru.lock().unwrap();
            for (b, block) in fetched {
                lru.insert((Arc::clone(pack), b), block);
            }
        }

        // 4. Assemble [start, end) from the covered slices.
        let mut out = Vec::with_capacity((end - start) as usize);
        for b in first_block..=last_block {
            let block = blocks
                .get(&b)
                .ok_or_else(|| anyhow!("internal: pack block {} missing after fetch", b))?;
            let block_base = b * bs;
            let s = start.max(block_base);
            let e = end.min(block_base + bs);
            let off_in_block = (s - block_base) as usize;
            let len = (e - s) as usize;
            if off_in_block + len > block.len() {
                return Err(anyhow!(
                    "pack block {} too short: need [{}, {}), have {} bytes",
                    b,
                    off_in_block,
                    off_in_block + len,
                    block.len()
                ));
            }
            out.extend_from_slice(&block[off_in_block..off_in_block + len]);
        }
        Ok(out)
    }
}

/// Collapse a sorted list of block indices into contiguous `[start, end]` runs
/// (inclusive block indices).
fn coalesce_runs(missing: &[u64]) -> Vec<(u64, u64)> {
    let mut runs = Vec::new();
    let mut iter = missing.iter().copied();
    let Some(mut run_start) = iter.next() else {
        return runs;
    };
    let mut prev = run_start;
    for b in iter {
        if b == prev + 1 {
            prev = b;
        } else {
            runs.push((run_start, prev));
            run_start = b;
            prev = b;
        }
    }
    runs.push((run_start, prev));
    runs
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[test]
    fn coalesce_runs_splits_on_gaps() {
        assert_eq!(coalesce_runs(&[]), vec![]);
        assert_eq!(coalesce_runs(&[3]), vec![(3, 3)]);
        assert_eq!(coalesce_runs(&[0, 1, 2]), vec![(0, 2)]);
        assert_eq!(coalesce_runs(&[0, 1, 3, 4, 7]), vec![(0, 1), (3, 4), (7, 7)]);
    }

    fn backing(len: usize) -> Vec<u8> {
        (0..len).map(|i| (i % 251) as u8).collect()
    }

    #[test]
    fn read_range_matches_backing_and_caches() {
        let data = backing(10_000);
        let pack_len = data.len() as u64;
        let cache = PackBlockCache::new(64, 1024); // 1 KiB blocks
        let pack: Arc<str> = Arc::from("deadbeef");
        let fetch_count = AtomicUsize::new(0);
        let fetch = |off: u64, len: u64| -> Result<Vec<u8>> {
            fetch_count.fetch_add(1, Ordering::SeqCst);
            Ok(data[off as usize..(off + len) as usize].to_vec())
        };

        // First read of [2000, 3500): spans blocks 1,2,3 -> one coalesced fetch.
        let got = cache.read_range(&pack, pack_len, 2000, 3500, &fetch).unwrap();
        assert_eq!(got, &data[2000..3500]);
        assert_eq!(fetch_count.load(Ordering::SeqCst), 1);

        // Overlapping read [2500, 3000): fully cached -> zero new fetches.
        let got = cache.read_range(&pack, pack_len, 2500, 3000, &fetch).unwrap();
        assert_eq!(got, &data[2500..3000]);
        assert_eq!(fetch_count.load(Ordering::SeqCst), 1);

        // A read touching a fresh block -> exactly one more fetch.
        let got = cache.read_range(&pack, pack_len, 5000, 5100, &fetch).unwrap();
        assert_eq!(got, &data[5000..5100]);
        assert_eq!(fetch_count.load(Ordering::SeqCst), 2);
    }

    #[test]
    fn read_range_clamps_last_block_to_pack_len() {
        let data = backing(1500);
        let pack_len = data.len() as u64;
        let cache = PackBlockCache::new(64, 1024);
        let pack: Arc<str> = Arc::from("cafef00d");
        let fetch = |off: u64, len: u64| -> Result<Vec<u8>> {
            // Never asked to read past EOF.
            assert!(off + len <= pack_len, "fetch past EOF: {}+{}", off, len);
            Ok(data[off as usize..(off + len) as usize].to_vec())
        };
        let got = cache.read_range(&pack, pack_len, 1400, 1500, &fetch).unwrap();
        assert_eq!(got, &data[1400..1500]);
    }

    #[test]
    fn oversize_request_bypasses_cache() {
        let data = backing(10_000);
        let pack_len = data.len() as u64;
        let cache = PackBlockCache::new(2, 1024); // cap 2 blocks
        let pack: Arc<str> = Arc::from("beadfeed");
        let fetch = |off: u64, len: u64| -> Result<Vec<u8>> {
            Ok(data[off as usize..(off + len) as usize].to_vec())
        };
        // 4000 bytes over 1KiB blocks == ~4 blocks > cap 2 -> bypass, nothing cached.
        let got = cache.read_range(&pack, pack_len, 1000, 5000, &fetch).unwrap();
        assert_eq!(got, &data[1000..5000]);
        assert_eq!(cache.cached_block_count(), 0);
    }
}
