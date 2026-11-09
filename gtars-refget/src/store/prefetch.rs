//! Bulk sequence prefetch: plan where each sequence comes from, fetch them in
//! parallel within a memory budget, check them, optionally cache them, and
//! hand them to the caller as they arrive.
//!
//! Planning reads store metadata on the caller's thread. The dispatcher and
//! workers own everything they touch (an `Arc<Fetcher>`, paths, sizes), so a
//! [`Prefetch`] is `'static`: the caller may change the store while it
//! iterates.
//!
//! Peak memory is about [`DEFAULT_MEMORY_BUDGET`] plus `jobs x` the ranged
//! piece size (8 MiB): each range body is read into a per-request buffer, then
//! copied into its sequence's buffer.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use anyhow::{Result, anyhow};
use crossbeam_channel::{Receiver, RecvTimeoutError, Sender, unbounded};

use super::StorageMode;
use super::cache::SeqCache;
#[cfg(feature = "http")]
use super::fetch::Fetcher;
use super::fetch::{AtomicCounters, FetchCounters, Field};
use super::readonly::ReadonlyRefgetStore;
use crate::digest::{SequenceMetadata, SequenceRecord, lookup_alphabet};
use crate::hashkeyable::{DigestKey, key_to_digest_string};

/// The most downloaded bytes a prefetch holds in memory at once: 1 GiB.
pub(crate) const DEFAULT_MEMORY_BUDGET: u64 = 1 << 30;

/// How a local store reads its files at most in parallel.
const LOCAL_THREADS: usize = 8;

// ============================================================================
// Public types
// ============================================================================

/// The order in which a [`Prefetch`] yields items.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PrefetchOrder {
    /// Each unique sequence once, in completion order. Largest start first.
    #[default]
    Any,
    /// One item per requested position, in request order. A repeated digest
    /// is read (or downloaded) once per position.
    Input,
}

/// How a prefetch runs.
#[derive(Debug, Clone, Default)]
pub struct PrefetchOptions {
    pub order: PrefetchOrder,
    /// Publish sequences downloaded from a remote store to the shared cache.
    pub write_cache: bool,
    /// Print a progress line to stderr, at most one every 5 s.
    pub progress: bool,
}

/// Where an item's bytes came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ItemSource {
    /// Already in memory (a loaded record or the resident overlay).
    Memory,
    /// A local store's own file.
    Local,
    /// The shared cache.
    Cache,
    /// Downloaded now.
    Remote,
}

/// One sequence's stored bytes (packed, raw or zstd, as on disk).
pub struct PrefetchItem {
    /// Position in the request (`Input`), or the first position of this
    /// sequence (`Any`).
    pub index: usize,
    /// sha512t24u digest, without the `SQ.` prefix.
    pub digest: String,
    pub source: ItemSource,
    buf: Buffer,
}

impl PrefetchItem {
    pub fn bytes(&self) -> &[u8] {
        match &self.buf.data {
            BufData::Owned(v) => v,
            BufData::Shared(a) => a,
        }
    }

    /// The bytes as an owned vector (copied only when shared with the store).
    pub(crate) fn into_vec(self) -> Vec<u8> {
        match self.buf.data {
            BufData::Owned(v) => v,
            BufData::Shared(a) => Arc::try_unwrap(a).unwrap_or_else(|a| a.as_ref().clone()),
        }
    }

    /// The bytes as a shared buffer (for the resident overlay).
    pub(crate) fn into_shared(self) -> Arc<Vec<u8>> {
        match self.buf.data {
            BufData::Owned(v) => Arc::new(v),
            BufData::Shared(a) => a,
        }
    }
}

impl std::fmt::Debug for PrefetchItem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PrefetchItem")
            .field("index", &self.index)
            .field("digest", &self.digest)
            .field("source", &self.source)
            .field("len", &self.bytes().len())
            .finish()
    }
}

/// A sequence's bytes, holding its share of the memory budget until dropped.
struct Buffer {
    data: BufData,
    _token: Option<BudgetToken>,
}

enum BufData {
    Owned(Vec<u8>),
    Shared(Arc<Vec<u8>>),
}

// ============================================================================
// Memory budget
// ============================================================================

struct Budget {
    limit: u64,
    state: Mutex<BudgetState>,
    cv: Condvar,
}

#[derive(Default)]
struct BudgetState {
    reserved: u64,
    max_reserved: u64,
}

struct BudgetToken {
    budget: Arc<Budget>,
    bytes: u64,
}

impl Drop for BudgetToken {
    fn drop(&mut self) {
        let mut st = self.budget.state.lock().unwrap();
        st.reserved -= self.bytes;
        drop(st);
        self.budget.cv.notify_all();
    }
}

impl Budget {
    /// Reserve `bytes`, waiting until it fits, nothing is reserved (a
    /// sequence bigger than the budget runs alone), or `bypass()` says this is
    /// the next item to deliver. `None` when cancelled.
    fn acquire(
        self: &Arc<Self>,
        bytes: u64,
        bypass: &dyn Fn() -> bool,
        cancel: &AtomicBool,
    ) -> Option<BudgetToken> {
        let mut st = self.state.lock().unwrap();
        loop {
            if cancel.load(Ordering::Relaxed) {
                return None;
            }
            if st.reserved == 0 || st.reserved.saturating_add(bytes) <= self.limit || bypass() {
                st.reserved += bytes;
                st.max_reserved = st.max_reserved.max(st.reserved);
                return Some(BudgetToken {
                    budget: Arc::clone(self),
                    bytes,
                });
            }
            st = self
                .cv
                .wait_timeout(st, Duration::from_millis(50))
                .unwrap()
                .0;
        }
    }
}

// ============================================================================
// Plan
// ============================================================================

/// Everything the dispatcher needs about one sequence, owned.
struct Entry {
    digest: String,
    meta: SequenceMetadata,
    /// Packed size; `None` for Zstd.
    size: Option<u64>,
    /// Size estimate for budgets and progress.
    est: u64,
    /// First position of this sequence in the request.
    first_index: usize,
    memory: Option<Arc<Vec<u8>>>,
    /// Local store file, or the cache file of a remote store.
    local: Option<PathBuf>,
    /// Path relative to the remote store root.
    rel: Option<String>,
}

/// Owned context shared by the dispatcher and workers.
struct Ctx {
    mode: StorageMode,
    #[cfg(feature = "http")]
    fetcher: Option<Arc<Fetcher>>,
    cache: Option<SeqCache>,
    write_cache: bool,
    verify: bool,
    input: bool,
    stats: Arc<AtomicCounters>,
    cancel: Arc<AtomicBool>,
    budget: Arc<Budget>,
    tx: Sender<Msg>,
    next_deliver: Arc<AtomicUsize>,
}

impl Ctx {
    fn is_remote(&self) -> bool {
        #[cfg(feature = "http")]
        {
            self.fetcher.is_some()
        }
        #[cfg(not(feature = "http"))]
        {
            false
        }
    }

    /// Count into this run's counters and the store's cumulative ones.
    fn count(&self, field: Field, n: u64) {
        self.stats.add(field, n);
        #[cfg(feature = "http")]
        if let Some(f) = &self.fetcher {
            f.stats().add(field, n);
        }
    }

    fn fail(&self, e: anyhow::Error) {
        if !self.cancel.swap(true, Ordering::Relaxed) {
            let _ = self.tx.send(Msg::Err(e));
        }
        self.budget.cv.notify_all();
    }
}

/// One fetch of one sequence: position `id` of the dispatch order.
struct Unit {
    id: usize,
    entry: Arc<Entry>,
    source: ItemSource,
    /// Range downloads assemble here.
    assembly: Mutex<Vec<u8>>,
    /// Ranges still in flight.
    remaining: AtomicUsize,
    token: Mutex<Option<BudgetToken>>,
}

enum Job {
    Read {
        unit: Arc<Unit>,
        path: PathBuf,
    },
    #[cfg(feature = "http")]
    Whole {
        unit: Arc<Unit>,
    },
    #[cfg(feature = "http")]
    Range {
        unit: Arc<Unit>,
        start: u64,
        end: u64,
    },
}

enum Msg {
    Done {
        unit_id: usize,
        first_index: usize,
        digest: String,
        source: ItemSource,
        est: u64,
        buf: Buffer,
    },
    Err(anyhow::Error),
}

// ============================================================================
// Prefetch
// ============================================================================

/// A running prefetch. Iterate it for the items; dropping it cancels the
/// downloads and joins its threads.
pub struct Prefetch {
    rx: Option<Receiver<Msg>>,
    input: bool,
    total: usize,
    delivered: usize,
    bytes_total: u64,
    bytes_done: u64,
    /// `Input`: finished units waiting for their turn, by position.
    done_units: HashMap<usize, (Buffer, ItemSource, String, u64)>,
    next_deliver: Arc<AtomicUsize>,
    stats: Arc<AtomicCounters>,
    started: Instant,
    cancel: Arc<AtomicBool>,
    budget: Arc<Budget>,
    progress: bool,
    last_print: Option<Instant>,
    threads: Vec<JoinHandle<()>>,
    failed: bool,
}

impl Prefetch {
    /// What this run has done so far.
    pub fn counters(&self) -> FetchCounters {
        self.stats.snapshot(self.started.elapsed())
    }

    #[cfg(test)]
    fn max_reserved(&self) -> u64 {
        self.budget.state.lock().unwrap().max_reserved
    }

    /// The progress line, e.g. `fetch: 812/6359 sequences, 3.21/15.90 GB,
    /// 27.4 MB/s, 1 retry, 120 cache hits`: first after 5 s, then every 5 s,
    /// and once at the end.
    fn report(&mut self) {
        if !self.progress {
            return;
        }
        let finished = self.delivered == self.total;
        match self.last_print {
            None if !finished => {
                self.last_print = Some(Instant::now());
                return;
            }
            Some(t) if !finished && t.elapsed() < Duration::from_secs(5) => return,
            _ => {}
        }
        self.last_print = Some(Instant::now());
        let c = self.counters();
        eprintln!(
            "fetch: {}/{} sequences, {:.2}/{:.2} GB, {:.1} MB/s, {} {}, {} cache hits",
            self.delivered,
            self.total,
            self.bytes_done as f64 / 1e9,
            self.bytes_total as f64 / 1e9,
            c.rate_mb_s(),
            c.retries,
            if c.retries == 1 { "retry" } else { "retries" },
            c.cache_hits
        );
    }

    fn stop(&mut self) {
        self.cancel.store(true, Ordering::Relaxed);
        self.budget.cv.notify_all();
        self.rx = None;
        for t in self.threads.drain(..) {
            let _ = t.join();
        }
    }

    fn deliver(
        &mut self,
        index: usize,
        buf: Buffer,
        source: ItemSource,
        digest: String,
        est: u64,
    ) -> PrefetchItem {
        self.delivered += 1;
        self.bytes_done += est;
        self.next_deliver.store(self.delivered, Ordering::Relaxed);
        self.budget.cv.notify_all();
        self.report();
        PrefetchItem {
            index,
            digest,
            source,
            buf,
        }
    }
}

impl Iterator for Prefetch {
    type Item = Result<PrefetchItem>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed || self.delivered >= self.total {
            if !self.threads.is_empty() {
                self.stop();
            }
            return None;
        }
        loop {
            if self.input {
                let pos = self.delivered;
                if let Some((buf, source, digest, est)) = self.done_units.remove(&pos) {
                    return Some(Ok(self.deliver(pos, buf, source, digest, est)));
                }
            }
            let msg = match self.rx.as_ref()?.recv_timeout(Duration::from_secs(1)) {
                Ok(m) => m,
                Err(RecvTimeoutError::Timeout) => {
                    self.report();
                    continue;
                }
                Err(RecvTimeoutError::Disconnected) => {
                    self.failed = true;
                    self.stop();
                    return Some(Err(anyhow!(
                        "prefetch stopped before delivering every sequence"
                    )));
                }
            };
            match msg {
                Msg::Err(e) => {
                    self.failed = true;
                    self.stop();
                    return Some(Err(e));
                }
                Msg::Done {
                    unit_id,
                    first_index,
                    digest,
                    source,
                    est,
                    buf,
                } => {
                    if self.input {
                        self.done_units.insert(unit_id, (buf, source, digest, est));
                    } else {
                        return Some(Ok(self.deliver(first_index, buf, source, digest, est)));
                    }
                }
            }
        }
    }
}

impl Drop for Prefetch {
    fn drop(&mut self) {
        self.stop();
    }
}

// ============================================================================
// Store entry points
// ============================================================================

fn packed_size(mode: StorageMode, meta: &SequenceMetadata) -> Option<u64> {
    match mode {
        StorageMode::Encoded => {
            let bits = lookup_alphabet(&meta.alphabet).bits_per_symbol as u64;
            Some((meta.length as u64 * bits).div_ceil(8))
        }
        StorageMode::Raw => Some(meta.length as u64),
        StorageMode::Zstd => None,
    }
}

impl ReadonlyRefgetStore {
    /// Fetch the stored bytes of many sequences in parallel.
    ///
    /// Digests may be sha512t24u or md5, with or without `SQ.`; an unknown
    /// digest fails before any download. Sources, in order: memory (a loaded
    /// record or the resident overlay), the local store's file, the shared
    /// cache (remote stores), then the remote. Whole sequences downloaded
    /// over the network are length- and digest-checked before delivery
    /// (`FetchConfig::verify`), and published to the cache when
    /// `opts.write_cache` is on. A remote store needs its sequence index loaded
    /// ([`RefgetStore::prefetch`](super::RefgetStore::prefetch) loads it).
    pub fn prefetch(&self, digests: &[String], opts: PrefetchOptions) -> Result<Prefetch> {
        let mut keys = Vec::with_capacity(digests.len());
        for d in digests {
            keys.push(
                self.resolve_sequence_key(d)
                    .ok_or_else(|| anyhow!("Sequence not found: {}", d))?,
            );
        }
        self.prefetch_keys(&keys, opts, DEFAULT_MEMORY_BUDGET)
    }

    /// Prefetch with cache writes on, dropping the bytes. Sequences already
    /// in the cache are not downloaded again. Returns the counters.
    pub fn cache_sequences(
        &self,
        digests: &[String],
        opts: PrefetchOptions,
    ) -> Result<FetchCounters> {
        let opts = PrefetchOptions {
            write_cache: true,
            order: PrefetchOrder::Any,
            ..opts
        };
        let mut run = self.prefetch(digests, opts)?;
        for item in run.by_ref() {
            item?;
        }
        Ok(run.counters())
    }

    /// Prefetch resolved keys, holding at most `budget` downloaded bytes.
    pub(crate) fn prefetch_keys(
        &self,
        keys: &[DigestKey],
        opts: PrefetchOptions,
        budget: u64,
    ) -> Result<Prefetch> {
        #[cfg(feature = "http")]
        let fetcher = self.fetcher.clone();
        #[cfg(feature = "http")]
        let (jobs, chunk, verify) = match &fetcher {
            Some(f) => (f.config().jobs, f.chunk(), f.config().verify),
            None => (LOCAL_THREADS, 0, false),
        };
        #[cfg(not(feature = "http"))]
        let (jobs, chunk, verify) = (LOCAL_THREADS, 0u64, false);
        let remote = self.remote_source.is_some();

        // Plan on the caller's thread: owned entries, one per unique sequence.
        let resident = self
            .resident
            .read()
            .map_err(|_| anyhow!("resident overlay lock poisoned"))?;
        let mut entries: HashMap<DigestKey, Arc<Entry>> = HashMap::new();
        let mut uniq: Vec<DigestKey> = Vec::new();
        for (i, key) in keys.iter().enumerate() {
            if entries.contains_key(key) {
                continue;
            }
            let record = self
                .sequence_store
                .get(key)
                .ok_or_else(|| anyhow!("Sequence not found: {}", key_to_digest_string(key)))?;
            let meta = record.metadata().clone();
            let memory = match record {
                SequenceRecord::Full { sequence, .. } => Some(Arc::clone(sequence)),
                SequenceRecord::Stub(_) => resident.get(key).cloned(),
            };
            let local = if memory.is_some() {
                None
            } else {
                self.local_seq_path(&meta)
            };
            let rel = if remote && memory.is_none() {
                Some(
                    self.sequence_relative_path(&meta.sha512t24u)?
                        .to_string_lossy()
                        .into_owned(),
                )
            } else {
                None
            };
            if memory.is_none() && local.is_none() && rel.is_none() {
                return Err(anyhow!(
                    "Sequence data for {} not found locally and no remote source configured",
                    meta.sha512t24u
                ));
            }
            let size = packed_size(self.mode, &meta);
            let est = size.unwrap_or(meta.length as u64);
            entries.insert(
                *key,
                Arc::new(Entry {
                    digest: meta.sha512t24u.clone(),
                    meta,
                    size,
                    est,
                    first_index: i,
                    memory,
                    local,
                    rel,
                }),
            );
            uniq.push(*key);
        }
        drop(resident);

        let input = opts.order == PrefetchOrder::Input;
        let positions: Vec<DigestKey> = if input {
            keys.to_vec()
        } else {
            // Largest first, so the long poles start first.
            let mut u = uniq.clone();
            u.sort_by(|a, b| {
                entries[b]
                    .est
                    .cmp(&entries[a].est)
                    .then(entries[a].first_index.cmp(&entries[b].first_index))
            });
            u
        };
        let total = positions.len();
        let bytes_total = positions.iter().map(|k| entries[k].est).sum();

        // Thread counts sized to the work, so a one-sequence load stays cheap.
        let need: Vec<&Arc<Entry>> = uniq
            .iter()
            .map(|k| &entries[k])
            .filter(|e| e.memory.is_none())
            .collect();
        let local_threads = LOCAL_THREADS
            .min(jobs)
            .min(need.len())
            .max(usize::from(!need.is_empty()));
        let remote_jobs: u64 = if remote {
            need.iter()
                .map(|e| match e.size {
                    Some(s) if chunk > 0 && s > chunk => s.div_ceil(chunk),
                    _ => 1,
                })
                .sum()
        } else {
            0
        };
        let remote_threads = (jobs as u64).min(remote_jobs) as usize;

        let (tx, rx) = unbounded::<Msg>();
        let stats = Arc::new(AtomicCounters::default());
        let cancel = Arc::new(AtomicBool::new(false));
        let budget = Arc::new(Budget {
            limit: budget.max(1),
            state: Mutex::new(BudgetState::default()),
            cv: Condvar::new(),
        });
        let next_deliver = Arc::new(AtomicUsize::new(0));
        let ctx = Arc::new(Ctx {
            mode: self.mode,
            #[cfg(feature = "http")]
            fetcher,
            cache: self.seq_cache.clone(),
            write_cache: opts.write_cache,
            verify,
            input,
            stats: Arc::clone(&stats),
            cancel: Arc::clone(&cancel),
            budget: Arc::clone(&budget),
            tx,
            next_deliver: Arc::clone(&next_deliver),
        });

        let mut threads = Vec::new();
        let (local_tx, local_rx) = unbounded::<Job>();
        for _ in 0..local_threads {
            let rx = local_rx.clone();
            let ctx = Arc::clone(&ctx);
            threads.push(std::thread::spawn(move || worker(ctx, rx)));
        }
        let (remote_tx, remote_rx) = unbounded::<Job>();
        for _ in 0..remote_threads {
            let rx = remote_rx.clone();
            let ctx = Arc::clone(&ctx);
            threads.push(std::thread::spawn(move || worker(ctx, rx)));
        }
        drop((local_rx, remote_rx));
        {
            let ctx = Arc::clone(&ctx);
            threads.push(std::thread::spawn(move || {
                dispatch(ctx, entries, positions, local_tx, remote_tx)
            }));
        }
        drop(ctx);

        Ok(Prefetch {
            rx: Some(rx),
            input,
            total,
            delivered: 0,
            bytes_total,
            bytes_done: 0,
            done_units: HashMap::new(),
            next_deliver,
            stats,
            started: Instant::now(),
            cancel,
            budget,
            progress: opts.progress,
            last_print: None,
            threads,
            failed: false,
        })
    }
}

// ============================================================================
// Dispatcher
// ============================================================================

fn dispatch(
    ctx: Arc<Ctx>,
    entries: HashMap<DigestKey, Arc<Entry>>,
    positions: Vec<DigestKey>,
    local_tx: Sender<Job>,
    remote_tx: Sender<Job>,
) {
    for (pos, key) in positions.iter().enumerate() {
        if ctx.cancel.load(Ordering::Relaxed) {
            return;
        }
        let entry = Arc::clone(&entries[key]);
        if let Err(e) = start_unit(&ctx, pos, entry, &local_tx, &remote_tx) {
            ctx.fail(e);
            return;
        }
    }
}

fn send_done(ctx: &Ctx, unit_id: usize, entry: &Entry, source: ItemSource, buf: Buffer) {
    let _ = ctx.tx.send(Msg::Done {
        unit_id,
        first_index: entry.first_index,
        digest: entry.digest.clone(),
        source,
        est: entry.est,
        buf,
    });
}

fn start_unit(
    ctx: &Arc<Ctx>,
    id: usize,
    entry: Arc<Entry>,
    local_tx: &Sender<Job>,
    remote_tx: &Sender<Job>,
) -> Result<()> {
    // 1. In memory: no worker, no budget.
    if let Some(mem) = &entry.memory {
        let buf = Buffer {
            data: BufData::Shared(Arc::clone(mem)),
            _token: None,
        };
        send_done(ctx, id, &entry, ItemSource::Memory, buf);
        return Ok(());
    }
    // An empty sequence has no bytes to read.
    if entry.size == Some(0) {
        let buf = Buffer {
            data: BufData::Owned(Vec::new()),
            _token: None,
        };
        send_done(ctx, id, &entry, ItemSource::Memory, buf);
        return Ok(());
    }

    let next_deliver = Arc::clone(&ctx.next_deliver);
    let input = ctx.input;
    let bypass = move || input && next_deliver.load(Ordering::Relaxed) == id;
    let Some(token) = ctx.budget.acquire(entry.est, &bypass, &ctx.cancel) else {
        return Ok(());
    };
    let unit = |source: ItemSource| {
        Arc::new(Unit {
            id,
            entry: Arc::clone(&entry),
            source,
            assembly: Mutex::new(Vec::new()),
            remaining: AtomicUsize::new(0),
            token: Mutex::new(Some(token)),
        })
    };

    // 2. Local store: read its file.
    if !ctx.is_remote() {
        let path = entry.local.clone().ok_or_else(|| {
            anyhow!(
                "Sequence data for {} not found locally and no remote source configured",
                entry.digest
            )
        })?;
        let _ = local_tx.send(Job::Read {
            unit: unit(ItemSource::Local),
            path,
        });
        return Ok(());
    }

    // 3. Remote store, cache file present (and the right length).
    if let Some(path) = &entry.local
        && let Ok(md) = std::fs::metadata(path)
    {
        if entry.size.is_none_or(|s| s == md.len()) {
            let _ = local_tx.send(Job::Read {
                unit: unit(ItemSource::Cache),
                path: path.clone(),
            });
            return Ok(());
        }
        let _ = std::fs::remove_file(path);
    }

    // 4. Remote: one GET, or ranged GETs for a big file.
    #[cfg(feature = "http")]
    {
        let chunk = ctx
            .fetcher
            .as_ref()
            .expect("remote store has a fetcher")
            .chunk();
        let unit = unit(ItemSource::Remote);
        match entry.size {
            Some(size) if size > chunk => {
                let n = size.div_ceil(chunk);
                *unit.assembly.lock().unwrap() = vec![0u8; size as usize];
                unit.remaining.store(n as usize, Ordering::Relaxed);
                for i in 0..n {
                    let _ = remote_tx.send(Job::Range {
                        unit: Arc::clone(&unit),
                        start: i * chunk,
                        end: ((i + 1) * chunk).min(size),
                    });
                }
            }
            _ => {
                let _ = remote_tx.send(Job::Whole { unit });
            }
        }
        Ok(())
    }
    #[cfg(not(feature = "http"))]
    {
        let _ = (remote_tx, unit);
        Err(anyhow!("Remote fetching requires the 'http' feature"))
    }
}

// ============================================================================
// Workers
// ============================================================================

fn worker(ctx: Arc<Ctx>, rx: Receiver<Job>) {
    while let Ok(job) = rx.recv() {
        if ctx.cancel.load(Ordering::Relaxed) {
            continue; // drain without work
        }
        if let Err(e) = run_job(&ctx, job) {
            ctx.fail(e);
        }
    }
}

fn run_job(ctx: &Arc<Ctx>, job: Job) -> Result<()> {
    match job {
        Job::Read { unit, path } => {
            let bytes = std::fs::read(&path)
                .map_err(|e| anyhow!("Failed to read sequence file {}: {}", path.display(), e))?;
            check_len(&unit, &bytes, &path.display().to_string())?;
            let n = bytes.len() as u64;
            if unit.source == ItemSource::Cache {
                ctx.count(|s| &s.cache_hits, 1);
                ctx.count(|s| &s.cache_hit_bytes, n);
            } else {
                ctx.count(|s| &s.local_reads, 1);
                ctx.count(|s| &s.local_read_bytes, n);
            }
            deliver(ctx, &unit, bytes);
            Ok(())
        }
        #[cfg(feature = "http")]
        Job::Whole { unit } => {
            let fetcher = ctx.fetcher.as_ref().expect("remote store has a fetcher");
            let rel = unit.entry.rel.as_deref().expect("remote entry has a path");
            let bytes = fetcher.get_file(rel, unit.entry.size, Some(&ctx.stats))?;
            finish_download(ctx, &unit, bytes)
        }
        #[cfg(feature = "http")]
        Job::Range { unit, start, end } => {
            let fetcher = ctx.fetcher.as_ref().expect("remote store has a fetcher");
            let rel = unit.entry.rel.as_deref().expect("remote entry has a path");
            let body = fetcher.get_range_with(rel, start, end, Some(&ctx.stats))?;
            unit.assembly.lock().unwrap()[start as usize..end as usize].copy_from_slice(&body);
            drop(body);
            if unit.remaining.fetch_sub(1, Ordering::AcqRel) == 1 {
                let bytes = std::mem::take(&mut *unit.assembly.lock().unwrap());
                finish_download(ctx, &unit, bytes)?;
            }
            Ok(())
        }
    }
}

fn check_len(unit: &Unit, bytes: &[u8], what: &str) -> Result<()> {
    if let Some(s) = unit.entry.size
        && bytes.len() as u64 != s
    {
        return Err(anyhow!(
            "{}: got {} bytes, expected {} for sequence {}",
            what,
            bytes.len(),
            s,
            unit.entry.digest
        ));
    }
    Ok(())
}

fn deliver(ctx: &Ctx, unit: &Unit, bytes: Vec<u8>) {
    let buf = Buffer {
        data: BufData::Owned(bytes),
        _token: unit.token.lock().unwrap().take(),
    };
    send_done(ctx, unit.id, &unit.entry, unit.source, buf);
}

/// A whole sequence came over the network: check it, cache it, deliver it.
#[cfg(feature = "http")]
fn finish_download(ctx: &Ctx, unit: &Unit, bytes: Vec<u8>) -> Result<()> {
    let fetcher = ctx.fetcher.as_ref().expect("remote store has a fetcher");
    let rel = unit.entry.rel.as_deref().expect("remote entry has a path");
    let url = fetcher.url(rel)?;
    check_len(unit, &bytes, &url)?;
    if ctx.verify {
        let got = super::verify::digest_packed(&bytes, &unit.entry.meta, ctx.mode)?;
        if got != unit.entry.digest {
            return Err(anyhow!(
                "{}: downloaded bytes hash to {}, expected {}; nothing was cached",
                url,
                got,
                unit.entry.digest
            ));
        }
    }
    if let (Some(cache), Some(dest), true) = (&ctx.cache, &unit.entry.local, ctx.write_cache) {
        cache.publish(dest, &bytes)?;
    }
    ctx.count(|s| &s.sequences_fetched, 1);
    deliver(ctx, unit, bytes);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hashkeyable::HashKeyable;
    use crate::store::{FastaImportOptions, RefgetStore};
    use tempfile::TempDir;

    fn fasta(sizes: &[usize]) -> String {
        sizes
            .iter()
            .enumerate()
            .map(|(i, n)| {
                let bases: String = (0..*n)
                    .map(|j| b"ACGT"[(i + j * 3 + j / 5) % 4] as char)
                    .collect();
                format!(">s{}\n{}\n", i, bases)
            })
            .collect()
    }

    /// A local store holding `text`, plus its sequence digests in FASTA order.
    fn build(text: &str, mode: StorageMode) -> (TempDir, std::path::PathBuf, Vec<String>) {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("s.fa");
        std::fs::write(&fa, text).unwrap();
        let path = dir.path().join("store");
        let mut w = RefgetStore::on_disk(&path).unwrap();
        w.set_encoding_mode(mode).unwrap();
        let (meta, _) = w
            .add_sequence_collection_from_fasta(&fa, FastaImportOptions::new())
            .unwrap();
        w.write().unwrap();
        drop(w);
        let mut s = RefgetStore::open_local(&path).unwrap();
        let coll = s.get_collection(&meta.digest).unwrap();
        let digests = coll
            .sequences
            .iter()
            .map(|r| r.metadata().sha512t24u.clone())
            .collect();
        (dir, path, digests)
    }

    /// A local Raw store whose sequences are `sizes` bytes each, plus their digests.
    fn store(sizes: &[usize]) -> (TempDir, RefgetStore, Vec<String>) {
        let (dir, path, digests) = build(&fasta(sizes), StorageMode::Raw);
        (dir, RefgetStore::open_local(&path).unwrap(), digests)
    }

    fn run(s: &RefgetStore, ds: &[String], order: PrefetchOrder, budget: u64) -> Prefetch {
        let keys: Vec<DigestKey> = ds.iter().map(|d| d.to_key()).collect();
        let opts = PrefetchOptions {
            order,
            ..PrefetchOptions::default()
        };
        s.inner.prefetch_keys(&keys, opts, budget).unwrap()
    }

    #[test]
    fn budget_caps_reserved_bytes_in_both_orders() {
        let (_d, s, ds) = store(&[10, 20, 30, 40]);
        for order in [PrefetchOrder::Any, PrefetchOrder::Input] {
            let mut run = run(&s, &ds, order, 50);
            let mut n = 0;
            for item in run.by_ref() {
                assert_eq!(item.unwrap().source, ItemSource::Local);
                n += 1;
            }
            assert_eq!(n, 4);
            assert!(
                run.max_reserved() <= 50,
                "{order:?}: {}",
                run.max_reserved()
            );
        }
    }

    #[test]
    fn budget_smaller_than_every_sequence_runs_one_at_a_time() {
        let (_d, s, ds) = store(&[10, 20, 30, 40]);
        let mut run = run(&s, &ds, PrefetchOrder::Any, 5);
        assert_eq!(run.by_ref().filter(|i| i.is_ok()).count(), 4);
        assert!(run.max_reserved() <= 40, "{}", run.max_reserved());
    }

    #[test]
    fn input_yields_in_request_order_with_repeats() {
        let (_d, s, ds) = store(&[10, 20, 30, 40]);
        let req = vec![
            ds[3].clone(),
            ds[0].clone(),
            ds[0].clone(),
            ds[2].clone(),
            ds[3].clone(),
        ];
        let got: Vec<(usize, String, usize)> = run(&s, &req, PrefetchOrder::Input, 1 << 20)
            .map(|i| {
                let i = i.unwrap();
                (i.index, i.digest.clone(), i.bytes().len())
            })
            .collect();
        let want: Vec<(usize, String, usize)> = req
            .iter()
            .enumerate()
            .map(|(i, d)| (i, d.clone(), [40, 10, 10, 30, 40][i]))
            .collect();
        assert_eq!(got, want);
    }

    #[test]
    fn held_items_block_new_dispatch() {
        let (_d, s, ds) = store(&[40, 40]);
        let mut run = run(&s, &ds, PrefetchOrder::Any, 50);
        let first = run.next().unwrap().unwrap();
        let (tx, rx) = std::sync::mpsc::channel();
        let handle = std::thread::spawn(move || {
            let second = run.next().map(|r| r.map(|i| i.digest));
            tx.send(second.is_some()).unwrap();
            run
        });
        assert!(
            rx.recv_timeout(Duration::from_millis(300)).is_err(),
            "dispatch must wait while the item is held"
        );
        drop(first);
        assert!(rx.recv_timeout(Duration::from_secs(5)).unwrap());
        let run = handle.join().unwrap();
        assert!(run.max_reserved() <= 50);
    }

    #[test]
    fn unknown_digest_fails_before_any_read() {
        let (_d, mut s, _ds) = store(&[10]);
        let err = s
            .prefetch(&["nope".to_string()], PrefetchOptions::default())
            .err()
            .unwrap();
        assert!(err.to_string().contains("Sequence not found: nope"));
    }

    // ------------------------------------------------------------------------
    // Ranged downloads, with a 7-byte piece so tiny files split
    // ------------------------------------------------------------------------

    #[cfg(feature = "http")]
    mod remote {
        use super::*;
        use crate::store::FetchConfig;
        use crate::store::http_fixture::{FaultKind, Fixture};

        const CHUNK: u64 = 7;

        /// Open `fx` as a remote store whose fetcher splits files every 7 bytes.
        fn open(fx: &Fixture, cache: &std::path::Path) -> RefgetStore {
            let cfg = FetchConfig {
                jobs: 4,
                ..FetchConfig::default()
            };
            let mut s = RefgetStore::open_remote(cache, fx.url("s"), &cfg).unwrap();
            s.set_quiet(true);
            s.inner.fetcher = Some(Arc::new(
                Fetcher::new(&fx.url("s"), cfg).unwrap().with_chunk(CHUNK),
            ));
            s.load_sequence_index().unwrap();
            s
        }

        /// Stored sizes chunk-1, chunk, chunk+1 and 3*chunk+2 (in Raw mode; 2-bit
        /// Encoded packs 4 bases per byte).
        const SIZES: [usize; 4] = [6, 7, 8, 23];

        #[test]
        fn ranged_downloads_match_local_bytes_and_run_concurrently() {
            for (mode, per_byte) in [(StorageMode::Raw, 1), (StorageMode::Encoded, 4)] {
                let bases: Vec<usize> = SIZES.iter().map(|n| n * per_byte).collect();
                let (_d, path, ds) = build(&fasta(&bases), mode);
                let local: HashMap<String, Vec<u8>> = run(
                    &RefgetStore::open_local(&path).unwrap(),
                    &ds,
                    PrefetchOrder::Any,
                    DEFAULT_MEMORY_BUDGET,
                )
                .map(|i| {
                    let i = i.unwrap();
                    (i.digest.clone(), i.bytes().to_vec())
                })
                .collect();
                let planned: usize = local
                    .values()
                    .map(|v| (v.len() as u64).div_ceil(CHUNK) as usize)
                    .sum();
                let fx = Fixture::serve(&path, &["s"]);
                fx.add_fault(".seq", FaultKind::Delay(Duration::from_millis(30)), 1000);
                let cache = tempfile::tempdir().unwrap();
                for order in [PrefetchOrder::Any, PrefetchOrder::Input] {
                    let s = open(&fx, cache.path());
                    fx.clear_log();
                    for item in run(&s, &ds, order, DEFAULT_MEMORY_BUDGET) {
                        let item = item.unwrap();
                        assert_eq!(item.source, ItemSource::Remote);
                        assert_eq!(item.bytes(), local[&item.digest].as_slice(), "{mode:?}");
                    }
                    let reqs = fx.requests_matching(".seq");
                    assert_eq!(reqs.len(), planned, "{mode:?} {order:?}");
                    let ranged = reqs.iter().filter(|r| r.range.is_some()).count();
                    assert!(ranged > 0 && ranged < planned, "{ranged} of {planned}");
                    assert!(fx.peak_in_flight() > 1 && fx.peak_in_flight() <= 4);
                }
            }
        }

        #[test]
        fn big_files_need_a_server_that_honors_range() {
            let (_d, path, ds) = build(&fasta(&SIZES), StorageMode::Raw);
            let fx = Fixture::with_options(&path, &["s"], false, false);
            let cache = tempfile::tempdir().unwrap();
            let s = open(&fx, cache.path());
            // At most one piece: one whole-file GET, which needs no Range.
            assert_eq!(run(&s, &ds[..2], PrefetchOrder::Any, 1 << 20).count(), 2);
            let err = run(&s, &ds[3..], PrefetchOrder::Any, 1 << 20)
                .next()
                .unwrap()
                .unwrap_err();
            assert!(
                err.to_string().contains("ignored the Range header"),
                "{err}"
            );
        }
    }
}
