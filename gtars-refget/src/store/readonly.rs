//! ReadonlyRefgetStore struct definition and core methods.

use super::*;
use super::alias::AliasManager;
use super::lock::{LockOptions, StoreLock};
#[cfg(feature = "filesystem")]
use super::prefetch::{PrefetchOptions, PrefetchOrder, DEFAULT_MEMORY_BUDGET};

use std::collections::{HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::path::{Path, PathBuf};

use super::catalog::{Catalog, CollRow, MemberView, NewMember, SeqCore, SeqRow, SeqView};

use anyhow::{anyhow, Context, Result};

use crate::collection::SequenceRecordExt;
use crate::digest::lookup_alphabet;
use crate::digest::{SequenceCollectionMetadata, SequenceMetadata, SequenceRecord};
use crate::digest::{decode_string_from_bytes, decode_substring_from_bytes, encode_sequence, StreamingDecoder};
use crate::digest::{byte_range_for_bases, decode_substring_from_bytes_at_offset};
use crate::hashkeyable::{DigestKey, HashKeyable, key_to_digest_string};
use crate::seqcol::metadata_matches_attribute;

use std::fs::{self, create_dir_all, File};
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::sync::{Arc, Mutex, RwLock};

/// Fallback capacity for the per-store open-file-descriptor cache, used when
/// the process's soft `RLIMIT_NOFILE` cannot be read at all (non-unix
/// targets, or the `getrlimit` syscall failing).
const SEQ_FD_CACHE_CAP_FLOOR: usize = 256;

/// Ceiling for the adaptive fd-cache capacity. Even under a very high or
/// unlimited `RLIMIT_NOFILE`, the cache stays bounded so its O(n) LRU eviction
/// scan (see `FdCache::insert`) stays cheap and a single store instance does
/// not gratuitously hold an unbounded number of fds open (a process may have
/// several stores open at once, plus its own stdio/sockets/other files).
const SEQ_FD_CACHE_CAP_CEILING: usize = 8192;

/// Headroom subtracted from the soft `RLIMIT_NOFILE` when sizing the fd
/// cache, reserved for the rest of the process (stdio, sockets, other open
/// files, other store instances).
const SEQ_FD_CACHE_HEADROOM: usize = 64;

/// Capacity of the per-store open-file-descriptor cache used by the
/// partial-read path, adapted to the process's soft `RLIMIT_NOFILE`.
///
/// A fixed cap of 256 (seqrepo's `fd_cache_size`) is too small: a
/// genome-scale refget store has tens of thousands of segment `.seq` files
/// (e.g. ~52k for a human genome at typical segment sizes; a transcriptome
/// store can have ~122k sequences). A scattered/genome-wide batch of queries
/// touches far more than 256 distinct files, so a fixed cap thrashes:
/// nearly every region evicts and reopens its file, costing an extra
/// `openat`+`close` pair per region. Sizing the cache from the process's own
/// fd limit lets it hold far more open handles when the ulimit allows it,
/// while still leaving headroom and an upper bound so it can't exhaust fds or
/// let LRU eviction get expensive.
fn default_fd_cache_cap() -> usize {
    #[cfg(unix)]
    {
        let mut rlim = libc::rlimit {
            rlim_cur: 0,
            rlim_max: 0,
        };
        if unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &mut rlim) } == 0 {
            let soft = rlim.rlim_cur;
            if soft == libc::RLIM_INFINITY {
                // No real ceiling to respect; use our own bound.
                return SEQ_FD_CACHE_CAP_CEILING;
            }
            let usable = (soft as usize).saturating_sub(SEQ_FD_CACHE_HEADROOM).max(1);
            // Never size the cache above what the rlimit (minus headroom)
            // actually allows, even when that is below the floor -- an
            // unusually small ulimit must be respected, not overridden by a
            // "typical minimum". Above the floor, cap at the ceiling.
            return usable.min(SEQ_FD_CACHE_CAP_CEILING);
        }
    }
    SEQ_FD_CACHE_CAP_FLOOR
}

/// True when `err` (as returned by [`ReadonlyRefgetStore::get_substring_from_disk`]
/// / [`ReadonlyRefgetStore::get_cached_seq_file`]) is ultimately caused by the
/// `.seq` file not existing, as opposed to some other failure (permission
/// denied, a bad range, a genuine I/O error). Used to fall through to the
/// remote fallback exactly on a missing file, replacing an upfront
/// `path.exists()` stat with inspecting the result of an attempted open --
/// the common case (file present) then costs no extra syscall at all, since
/// the fd cache would have to open it anyway.
fn is_missing_seq_file(err: &anyhow::Error) -> bool {
    err.chain()
        .any(|cause| matches!(cause.downcast_ref::<std::io::Error>(), Some(io_err) if io_err.kind() == std::io::ErrorKind::NotFound))
}

/// Error unless `[start, end)` is a non-empty range inside a sequence of
/// `length` bases.
fn check_range(start: usize, end: usize, length: usize) -> Result<()> {
    if start >= length || end > length || start >= end {
        return Err(anyhow!(
            "Invalid substring range: start={}, end={}, sequence length={}",
            start,
            end,
            length
        ));
    }
    Ok(())
}

/// Batch size at which `get_substrings` against a *remote-only* sequence stops
/// issuing one HTTP byte-range request per range and instead downloads the whole
/// stored record once (cached locally as a `.seq`, like `load_sequence`), then
/// serves every range from the local file. Each remote range is a separate
/// round-trip, so for bulk extraction (e.g. BED region sets) the per-request
/// latency dominates and a single cached download wins; below the threshold,
/// fetching only the touched bytes is cheaper than pulling a whole chromosome.
/// Heuristic, count-based: small/ad-hoc reads stay lean, bulk reads promote.
const REMOTE_BULK_FETCH_THRESHOLD: usize = 16;

/// Convert decoded sequence bytes into a `String` without re-validating UTF-8.
///
/// # Why this exists
/// gtars sequence bytes are *always* ASCII: the Encoded decoder emits only
/// bytes drawn from `alphabet.decoding_array` (A/C/G/T/N/IUPAC/protein), and Raw
/// mode stores one ASCII base per byte per the refget data model. So
/// `String::from_utf8`'s validation scan can never fail on this data — it is
/// pure overhead.
///
/// # Why it justifies a (single, audited) `unsafe`
/// Benchmarked (AMD Ryzen 7 3800X, release build): on the whole-sequence Encoded
/// decode path the `from_utf8` validation scan costs a steady ~6–10% of total
/// decode-and-convert time — roughly 10–20 ms per 250 Mbp chromosome. That is a
/// real, repeatable cost on the whole-genome read hot path, so we skip it here.
///
/// This is the crate's only `from_utf8_unchecked`. Remote data is validated in
/// `get_substring_from_remote` before it reaches this conversion.
///
/// The `debug_assert` re-validates the ASCII invariant in debug and test builds,
/// so any future break (a non-ASCII alphabet, an offset bug) fails loudly under
/// test while release builds stay fast.
#[inline]
fn decoded_sequence_to_string(bytes: Vec<u8>) -> String {
    debug_assert!(
        bytes.is_ascii(),
        "decoded sequence bytes must be ASCII (decoder invariant violated)"
    );
    // SAFETY: gtars sequence bytes are guaranteed ASCII (see doc above), hence
    // valid UTF-8, so skipping validation is sound.
    unsafe { String::from_utf8_unchecked(bytes) }
}

/// A tiny bounded LRU cache mapping a sequence digest to an open, shared
/// `.seq` file handle. `.seq` content for a digest is immutable, so a cached
/// handle never goes stale and needs no invalidation. When the cache is full,
/// inserting a new handle evicts (closes) the least-recently-used one.
///
/// Recency is tracked with a monotonically increasing counter stamped on each
/// access; eviction scans for the minimum stamp. With a small cap (default 256)
/// this O(n) eviction is negligible compared to the avoided `open`/`seek`
/// syscalls on the hot path.
#[derive(Debug)]
struct FdCache {
    cap: usize,
    tick: u64,
    /// relative-or-absolute file path -> (shared file handle, last-use stamp).
    /// Keyed by PATH (not digest) so pack files and `.seq` files share one
    /// bounded LRU of open handles.
    entries: HashMap<PathBuf, (Arc<File>, u64)>,
}

impl FdCache {
    fn new(cap: usize) -> Self {
        FdCache {
            cap: cap.max(1),
            tick: 0,
            entries: HashMap::with_capacity(cap.max(1)),
        }
    }

    /// Return a cloned handle on hit, bumping recency. `None` on miss.
    fn get(&mut self, key: &Path) -> Option<Arc<File>> {
        self.tick += 1;
        let tick = self.tick;
        if let Some(slot) = self.entries.get_mut(key) {
            slot.1 = tick;
            Some(Arc::clone(&slot.0))
        } else {
            None
        }
    }

    /// Insert a freshly opened handle, evicting (closing) the LRU entry if full.
    fn insert(&mut self, key: PathBuf, file: Arc<File>) {
        self.tick += 1;
        let tick = self.tick;
        if !self.entries.contains_key(&key) && self.entries.len() >= self.cap {
            // Evict the least-recently-used entry. Dropping the Arc closes the
            // fd once no in-flight read still holds a clone.
            if let Some(lru_key) = self
                .entries
                .iter()
                .min_by_key(|(_, (_, stamp))| *stamp)
                .map(|(k, _)| k.clone())
            {
                self.entries.remove(&lru_key);
            }
        }
        self.entries.insert(key, (file, tick));
    }
}


/// A `Read` adapter that yields a byte range from a shared `Arc<Vec<u8>>`
/// buffer without cloning the backing data.
///
/// Used by `stream_sequence` when the target sequence is already resident
/// in memory as `SequenceRecord::Full`. The alternative of slicing+cloning
/// would allocate a second copy of the (possibly chromosome-sized) buffer
/// and defeat the O(1) memory guarantee of streaming.
pub(crate) struct ArcSliceReader {
    buf: Arc<Vec<u8>>,
    pos: usize,
    end: usize,
}

impl ArcSliceReader {
    pub(crate) fn new(buf: Arc<Vec<u8>>, start: usize, end: usize) -> Self {
        debug_assert!(start <= end);
        debug_assert!(end <= buf.len());
        Self { buf, pos: start, end }
    }
}

impl Read for ArcSliceReader {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        let remaining = self.end - self.pos;
        if remaining == 0 {
            return Ok(0);
        }
        let n = remaining.min(out.len());
        out[..n].copy_from_slice(&self.buf[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

/// Stand-in for the fetcher when the `http` feature is off; never constructed.
#[cfg(not(feature = "http"))]
struct NoFetcher;

#[cfg(not(feature = "http"))]
impl NoFetcher {
    fn get_range(&self, _rel: &str, _start: u64, _end: u64) -> Result<Vec<u8>> {
        Err(anyhow!("Remote access requires the 'http' feature to be enabled"))
    }

    fn open_range(&self, _rel: &str, _start: u64, _end: u64) -> Result<Box<dyn Read + Send>> {
        Err(anyhow!("Remote access requires the 'http' feature to be enabled"))
    }
}

/// Core refget store with `&self` read methods, suitable for `Arc` sharing in servers.
///
/// Mutating methods are used during the setup/loading phase; once wrapped in `Arc`,
/// only `&self` reads are accessible, making concurrent access thread-safe.
///
/// All sequence and collection metadata lives in one columnar [`Catalog`]:
/// sequences deduplicated store-wide, each collection's members (with their
/// per-collection names) as a contiguous range, and digest/name lookups as
/// tables of row numbers. Sequence bytes are never part of the catalog; they
/// live in [`Self::loaded`] (owned bytes of in-memory stores) or the
/// evictable `resident` overlay, both keyed by sequence row.
#[derive(Debug)]
pub struct ReadonlyRefgetStore {
    /// Every sequence, membership and collection fact, each stored once.
    pub(crate) catalog: Catalog,
    /// Bytes owned by the store: sequences of an in-memory store, or ones
    /// loaded with `load_sequence`. Keyed by catalog sequence row.
    pub(crate) loaded: HashMap<SeqRow, Arc<Vec<u8>>>,
    /// Storage strategy for sequences
    pub(crate) mode: StorageMode,
    /// Where the store lives on disk (local store or cache directory)
    pub(crate) local_path: Option<PathBuf>,
    /// Where to pull sequences from (if remote-backed)
    pub(crate) remote_source: Option<String>,
    /// Template for sequence file paths (e.g., "sequences/%s2/%s.seq")
    pub(crate) seqdata_path_template: Option<String>,
    /// Whether to persist sequences to disk (write-through caching)
    pub(crate) persist_to_disk: bool,
    /// Whether to suppress progress output
    pub(crate) quiet: bool,
    /// Whether to compute ancillary digests (nlp, snlp, sorted_sequences).
    /// Default: true for new stores.
    pub(crate) ancillary_digests: bool,
    /// Whether on-disk attribute reverse index is enabled.
    /// Default: false.
    pub(crate) attribute_index: bool,
    /// Human-readable aliases for sequences and collections.
    pub(crate) aliases: AliasManager,
    /// FHR metadata for collections, keyed by collection digest.
    pub(crate) fhr_metadata: HashMap<DigestKey, super::fhr_metadata::FhrMetadata>,
    /// Available sequence alias namespaces (from manifest, for remote discovery).
    pub(crate) available_sequence_alias_namespaces: Vec<String>,
    /// Available collection alias namespaces (from manifest, for remote discovery).
    pub(crate) available_collection_alias_namespaces: Vec<String>,
    /// Bounded LRU cache of open `.seq` file handles for the partial-read path,
    /// behind a `Mutex` because `get_substring` takes `&self`. Avoids re-opening
    /// the same chromosome file on every per-region `get_substring` call.
    seq_fd_cache: Mutex<FdCache>,
    /// `&self` resident overlay: whole, on-disk-packed `.seq` byte arrays keyed by
    /// SHA512t24u digest. Populated via [`Self::load_resident`] to promote a
    /// sequence that is not loaded to RAM without taking `&mut self` (the store is
    /// `Arc`-shared behind a readonly handle). Consulted in `get_substring` /
    /// `get_substrings` *before* the on-disk partial-read path: on a hit the
    /// covering span is decoded from these bytes exactly like loaded bytes, so
    /// output is byte-identical to the disk path but every region query becomes an
    /// in-memory decode instead of a positioned read. The bytes are the raw file
    /// contents (2-bit packed in `Encoded` mode, one byte/base in `Raw`), so a hit
    /// is mode-agnostic and only the sliced span is ever decoded. `RwLock` so many
    /// readers proceed concurrently; the `Arc` is cloned out under a brief read
    /// lock and decoded after the lock is dropped, so decodes never serialize.
    pub(crate) resident: RwLock<HashMap<SeqRow, Arc<Vec<u8>>>>,
    /// Whether the sequence index (sequences.rgsi) has been loaded.
    /// For remote stores, this starts as `false` and is lazily loaded on first
    /// sequence access, avoiding the costly download when only browsing collections.
    pub(crate) sequence_index_loaded: bool,
    /// The relative path to the sequence index file (from rgstore.json),
    /// stored for deferred loading in remote stores.
    pub(crate) sequence_index_path: Option<String>,
    /// What this handle has actually changed since its last successful commit.
    /// See [`PendingChanges`]. Behind a `Mutex` because the commit path is
    /// `&self`.
    pub(crate) pending: Mutex<PendingChanges>,
    /// A store-directory write lock held across several mutations, set by
    /// [`Self::lock_for_batch`]. When present, individual commits skip acquiring
    /// their own lock — otherwise a batch guard would deadlock against the
    /// `write_index_files` each mutation triggers.
    pub(crate) commit_lock: Option<Arc<StoreLock>>,
    /// Timeout/staleness settings used when acquiring the store write lock.
    pub(crate) lock_options: LockOptions,
    /// Re-point an existing alias (in memory, or already published on disk) to
    /// a different target, instead of erroring on the conflict. Off by default:
    /// binding one human-readable name to different content is a semantic
    /// conflict, not a merge detail. Never overrides a clash within one load
    /// (one TSV, or one FASTA's header tags).
    pub(crate) force_alias: bool,

    // --- Packed layout -------------------------------------------------------
    /// Physical layout of THIS store's sequence bytes. Invariant:
    /// `layout == Packed` implies `remote_source == None` (the local cache of a
    /// remote-backed store is always standalone).
    pub(crate) layout: Layout,
    /// Layout of the REMOTE this store caches from (if any). A remote packed
    /// store is read through its sidecar while the local cache stays standalone.
    pub(crate) remote_layout: Layout,
    /// Per-pack size cap in bytes for new packs (default 256 MiB).
    pub(crate) pack_cap_bytes: u64,
    /// Optional block cache in front of local pack reads, shared across all
    /// packs of this store. Off (`None`) by default; see
    /// [`Self::set_pack_cache_mb`].
    pub(crate) pack_cache: Option<super::packcache::PackBlockCache>,
    /// Pack location of each sequence row on the LOCAL side, built lazily from
    /// the sidecar on the first packed read. Invalidated (set to `None`)
    /// whenever a commit changes the sidecar or the catalog is renumbered.
    pub(crate) pack_cols: Mutex<Option<Arc<super::catalog::PackColumns>>>,
    /// The REMOTE side's `sequences.pack.idx` bytes, fetched whole on the
    /// first remote packed read and searched by digest. Not row-keyed like
    /// `pack_cols`, because a remote store's catalog grows after open as
    /// collections and the sequence index load.
    pub(crate) remote_pack_index: Mutex<Option<Arc<Vec<u8>>>>,
    /// Cache of pack file byte lengths (needed to clamp the final block of a
    /// local block-cache read), keyed by content-hash pack name.
    pub(crate) pack_lengths: Mutex<HashMap<Arc<str>, u64>>,
    /// The open, not-yet-sealed pack for an incremental (`add_sequence_record`)
    /// build session on a packed store. Sealed at cap and always at commit.
    pub(crate) open_session_pack: Mutex<Option<SessionPack>>,
    /// Downloads files from the remote store (remote stores only).
    #[cfg(feature = "http")]
    pub(crate) fetcher: Option<Arc<super::fetch::Fetcher>>,
    /// The shared download cache that holds this remote store's `.seq` files
    /// (remote stores only). Metadata lives in `local_path`, the remote's
    /// origin folder inside the same cache.
    pub(crate) seq_cache: Option<super::cache::SeqCache>,
}

/// The open, not-yet-sealed pack file a packed store writes new sequences
/// into, from FASTA import, `add_sequence_record` and `compact`. Bytes are
/// appended in call order, the content hash is fed the same bytes in the same
/// order, and the pack is sealed (fsync + rename to its content-hash name) at
/// the size cap, at the end of an import, and at commit.
#[derive(Debug)]
pub(crate) struct SessionPack {
    /// Directory `packs/` under the store.
    pub(crate) packs_dir: PathBuf,
    /// Names this session's temp packs, so no other writer opens the same file.
    pub(crate) token: String,
    /// Ordinal of the current open sub-pack (for temp-file naming).
    pub(crate) ordinal: usize,
    /// Temp file path `packs/.rgstore.tmp.pack.session.<token>.<ordinal>`.
    pub(crate) tmp_path: PathBuf,
    /// Open handle to the temp pack.
    pub(crate) file: File,
    /// Running byte length == next offset.
    pub(crate) len: u64,
    /// Running content hash of the bytes appended so far.
    pub(crate) hasher: sha2::Sha256,
    /// `(digest, offset, len)` rows for sequences appended to THIS open pack.
    pub(crate) rows: Vec<(String, u64, u64)>,
    /// Sidecar rows of the packs this session already sealed.
    pub(crate) sealed: Vec<(String, super::packidx::PackRow)>,
    /// Names of packs this session renamed into place and has not handed to
    /// a commit yet. A pack that already existed under the same name is not
    /// listed: another writer owns it.
    pub(crate) created: Vec<String>,
}

impl SessionPack {
    /// Create the temp file for sub-pack `ordinal`. `create_new`, so a name
    /// another writer already holds is an error instead of a shared file.
    fn create_temp(packs_dir: &Path, token: &str, ordinal: usize) -> Result<(PathBuf, File)> {
        let tmp_path = packs_dir.join(format!(".rgstore.tmp.pack.session.{}.{}", token, ordinal));
        let file = std::fs::OpenOptions::new()
            .create_new(true)
            .read(true)
            .write(true)
            .open(&tmp_path)
            .with_context(|| format!("creating temp pack {}", tmp_path.display()))?;
        Ok((tmp_path, file))
    }

    /// Seal the open pack to its content-hash name and open a fresh empty
    /// temp for continued writing. No-op on an empty pack.
    ///
    /// On error the session is unchanged (same open temp, rows and hash), so
    /// the seal can be retried. The next temp is created before the rename for
    /// that reason: after the rename nothing can fail.
    pub(crate) fn seal(&mut self) -> Result<()> {
        use sha2::Digest as _;
        if self.len == 0 {
            return Ok(());
        }
        // Cut off bytes past the tracked length (left by a failed append), so
        // the sealed file is exactly the bytes that were hashed.
        self.file
            .set_len(self.len)
            .context("truncating session pack to its written length")?;
        let _ = self.file.sync_all();
        let digest = self.hasher.clone().finalize();
        let name: String = digest[..16].iter().map(|b| format!("{:02x}", b)).collect();
        let final_path = self.packs_dir.join(format!("{}.pack", name));
        let (next_tmp, next_file) =
            Self::create_temp(&self.packs_dir, &self.token, self.ordinal + 1)?;
        if final_path.exists() {
            let _ = std::fs::remove_file(&self.tmp_path);
        } else if let Err(e) = std::fs::rename(&self.tmp_path, &final_path) {
            let _ = std::fs::remove_file(&next_tmp);
            return Err(e).with_context(|| format!("sealing session pack {}", final_path.display()));
        } else {
            self.created.push(name.clone());
        }
        if let Ok(dir) = File::open(&self.packs_dir) {
            let _ = dir.sync_all();
        }
        let pack: Arc<str> = Arc::from(name.as_str());
        for (digest, offset, len) in self.rows.drain(..) {
            let row = super::packidx::PackRow { pack: Arc::clone(&pack), offset, len };
            self.sealed.push((digest, row));
        }

        // Continue in the fresh temp.
        self.ordinal += 1;
        self.file = next_file;
        self.tmp_path = next_tmp;
        self.len = 0;
        self.hasher = sha2::Sha256::new();
        Ok(())
    }
}

/// A session pack that is dropped without being sealed (the store handle went
/// away before a commit) leaves no temp file behind. After a seal, `tmp_path`
/// names the fresh empty temp, or a path already renamed away. The packs it
/// sealed but never handed to a commit (`created`) are deleted too, since no
/// row will ever point at them.
impl Drop for SessionPack {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.tmp_path);
        remove_uncommitted_packs(&self.packs_dir, &self.created, &HashSet::new());
    }
}

/// Delete each pack in `created` (packs a writer renamed into place whose
/// rows never reached a commit) that neither `staged` nor the on-disk
/// sidecar references. Deletes nothing if the sidecar cannot be read.
///
/// Runs without the store lock. Another writer may have reused one of these
/// packs by name (same content hash). If it has committed, the sidecar lists
/// the pack and it is kept. If it has not, its commit finds the pack gone and
/// refuses to publish (see `commit_pack_idx`), a clear error rather than a
/// bad read. A commit landing between the sidecar read here and the unlink is
/// not caught: its rows then name a missing pack and their reads fail.
pub(crate) fn remove_uncommitted_packs(
    packs_dir: &Path,
    created: &[String],
    staged: &HashSet<String>,
) {
    if created.is_empty() {
        return;
    }
    let Some(local_path) = packs_dir.parent() else { return };
    let Ok(listed) = ReadonlyRefgetStore::sidecar_pack_names(local_path) else { return };
    for name in created {
        if !staged.contains(name) && !listed.contains(name) {
            let _ = std::fs::remove_file(packs_dir.join(format!("{}.pack", name)));
        }
    }
}

/// Which side of a (possibly remote-backed) store a resolution targets. The
/// local cache of a remote packed store is standalone — sealed packs cannot be
/// appended to — so the two sides can have different layouts.
#[derive(Clone, Copy, Debug)]
pub(crate) enum StoreSide {
    Local,
    Remote,
}

/// Where one sequence's stored bytes sit on this machine (see
/// [`ReadonlyRefgetStore::stored_span`]): a whole `.seq` file, or a span of a
/// pack file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StoredSpan {
    /// The `.seq` file or the pack file.
    pub path: PathBuf,
    /// Byte offset of the record in `path` (0 for a `.seq` file).
    pub offset: u64,
    /// Exact record length for a pack span; `None` = the whole file.
    pub len: Option<u64>,
}

/// Path of a sealed pack file, relative to the store root.
pub(crate) fn pack_relpath(pack: &str) -> PathBuf {
    PathBuf::from(format!("packs/{}.pack", pack))
}

/// Result of [`ReadonlyRefgetStore::compact`].
#[derive(Debug, Clone)]
pub struct CompactReport {
    /// Pack files on disk before compaction.
    pub packs_before: usize,
    /// Pack files after compaction.
    pub packs_after: usize,
    /// Dead bytes reclaimed.
    pub bytes_reclaimed: u64,
    /// Unreferenced pack files deleted by
    /// [`CompactOptions::prune_unreferenced`]; 0 without it.
    pub packs_pruned: usize,
}

/// Options for [`ReadonlyRefgetStore::compact_with`].
#[derive(Debug, Clone, Copy, Default)]
pub struct CompactOptions {
    /// Also delete every pack file that the new index does not reference,
    /// including packs no index ever listed (left by a crashed writer, for
    /// example). Such a pack can belong to another writer that has sealed it
    /// but not committed yet, so set this only when no other process is
    /// writing to the store.
    pub prune_unreferenced: bool,
}

/// What this handle has actually created, modified, or deleted since its last
/// successful commit.
///
/// A store handle's in-memory maps hold two different things mixed together:
/// rows it produced, and rows it merely LOADED so they could be read
/// (`open_local` loads a stub for every collection in `collections.rgci`, and
/// every row of `sequences.rgsi`). A commit that cannot tell those apart has to
/// write the whole pile back — and the pile is an open-time snapshot, which on a
/// long import is hours stale.
///
/// Recording changes explicitly is what makes a commit a DELTA: take the lock,
/// re-read the index from disk, apply exactly these sets, publish. A row this
/// handle never touched is never rewritten, so a concurrent writer's additions
/// cannot be clobbered and its deletions cannot be resurrected.
///
/// In-memory only, never written to disk: this describes an intent that is fully
/// discharged by the commit that applies it, and is cleared as soon as that
/// commit succeeds. A handle dropped without committing abandons its pending
/// changes — see [`ReadonlyRefgetStore::has_uncommitted_changes`].
#[derive(Debug, Default, Clone)]
pub struct PendingChanges {
    /// Collection digests to publish into `collections.rgci`.
    pub(crate) collections: HashSet<DigestKey>,
    /// Sequence digests to publish into `sequences.rgsi`.
    pub(crate) sequences: HashSet<DigestKey>,
    /// `(kind, namespace, alias)` triples to publish into the alias TSVs.
    pub(crate) aliases: HashSet<(AliasKind, String, String)>,
    /// Collection digests to drop from `collections.rgci`.
    pub(crate) removed_collections: HashSet<DigestKey>,
    /// Sequence digests to drop from `sequences.rgsi`.
    pub(crate) removed_sequences: HashSet<DigestKey>,
    /// `(kind, namespace, alias)` triples to drop from the alias TSVs.
    pub(crate) removed_aliases: HashSet<(AliasKind, String, String)>,
    /// `(kind, namespace)` pairs this handle explicitly emptied. Recorded so the
    /// commit visits the namespace at all; whether the TSV is actually deleted
    /// depends on the MERGED result being empty, not on this handle's intent.
    pub(crate) emptied_alias_namespaces: HashSet<(AliasKind, String)>,
    /// Subset of `aliases` the caller asked to force: on commit these overwrite
    /// a conflicting binding another writer already published, as if
    /// `force_alias` were set for that one entry.
    pub(crate) forced_aliases: HashSet<(AliasKind, String, String)>,
    /// Sidecar rows to publish into `sequences.pack.idx` (packed layout only):
    /// digest key -> its `(pack, offset, len)` span. Merged and republished
    /// under the store lock in `write_index_files`.
    pub(crate) pack_rows: HashMap<DigestKey, super::packidx::PackRow>,
}

impl PendingChanges {
    pub(crate) fn is_empty(&self) -> bool {
        self.collections.is_empty()
            && self.sequences.is_empty()
            && self.aliases.is_empty()
            && self.removed_collections.is_empty()
            && self.removed_sequences.is_empty()
            && self.removed_aliases.is_empty()
            && self.emptied_alias_namespaces.is_empty()
            && self.pack_rows.is_empty()
    }

    /// Record a collection this handle added. An add supersedes a pending
    /// removal of the same digest (and vice versa), so the two sets stay
    /// disjoint and the commit never has to guess which wins.
    pub(crate) fn add_collection(&mut self, key: DigestKey) {
        self.removed_collections.remove(&key);
        self.collections.insert(key);
    }

    pub(crate) fn remove_collection(&mut self, key: DigestKey) {
        self.collections.remove(&key);
        self.removed_collections.insert(key);
    }

    pub(crate) fn add_sequence(&mut self, key: DigestKey) {
        self.removed_sequences.remove(&key);
        self.sequences.insert(key);
    }

    pub(crate) fn remove_sequence(&mut self, key: DigestKey) {
        self.sequences.remove(&key);
        self.removed_sequences.insert(key);
    }

    /// Drop sequences from the pending set WITHOUT recording a removal.
    ///
    /// For rows that were staged but never published — the failed-import
    /// cleanup path. Recording a removal instead would be wrong: it would delete
    /// an identical row another writer had legitimately published, since
    /// sequence rows are content-addressed and therefore shared.
    pub(crate) fn forget_sequences(&mut self, keys: &HashSet<DigestKey>) {
        self.sequences.retain(|k| !keys.contains(k));
        self.pack_rows.retain(|k, _| !keys.contains(k));
    }

    pub(crate) fn add_alias(&mut self, kind: AliasKind, namespace: &str, alias: &str) {
        let entry = (kind, namespace.to_string(), alias.to_string());
        self.removed_aliases.remove(&entry);
        self.aliases.insert(entry);
    }

    pub(crate) fn add_forced_alias(&mut self, kind: AliasKind, namespace: &str, alias: &str) {
        self.add_alias(kind, namespace, alias);
        self.forced_aliases
            .insert((kind, namespace.to_string(), alias.to_string()));
    }

    pub(crate) fn is_forced(&self, kind: AliasKind, namespace: &str, alias: &str) -> bool {
        self.forced_aliases
            .contains(&(kind, namespace.to_string(), alias.to_string()))
    }

    pub(crate) fn remove_alias(&mut self, kind: AliasKind, namespace: &str, alias: &str) {
        let entry = (kind, namespace.to_string(), alias.to_string());
        self.aliases.remove(&entry);
        self.forced_aliases.remove(&entry);
        self.removed_aliases.insert(entry);
    }

    /// Namespaces of one kind this handle has touched — the only ones a commit
    /// needs to rewrite. Namespaces nobody changed are left alone on disk.
    pub(crate) fn touched_alias_namespaces(&self, kind: AliasKind) -> HashSet<String> {
        let mut out = HashSet::new();
        for (k, ns, _) in self.aliases.iter().chain(self.removed_aliases.iter()) {
            if *k == kind {
                out.insert(ns.clone());
            }
        }
        for (k, ns) in &self.emptied_alias_namespaces {
            if *k == kind {
                out.insert(ns.clone());
            }
        }
        out
    }

    /// Drop exactly the changes a completed commit applied.
    ///
    /// Set difference, not `clear()`: another thread may have recorded a change
    /// after the commit snapshotted, and that one still has to be applied.
    pub(crate) fn clear_matching(&mut self, applied: &PendingChanges) {
        self.collections.retain(|k| !applied.collections.contains(k));
        self.sequences.retain(|k| !applied.sequences.contains(k));
        self.aliases.retain(|k| !applied.aliases.contains(k));
        self.forced_aliases.retain(|k| !applied.aliases.contains(k));
        self.removed_collections
            .retain(|k| !applied.removed_collections.contains(k));
        self.removed_sequences
            .retain(|k| !applied.removed_sequences.contains(k));
        self.removed_aliases
            .retain(|k| !applied.removed_aliases.contains(k));
        self.emptied_alias_namespaces
            .retain(|k| !applied.emptied_alias_namespaces.contains(k));
        self.pack_rows
            .retain(|k, _| !applied.pack_rows.contains_key(k));
    }

    /// Drop the alias changes for one namespace, after the cheap
    /// single-namespace commit path published it.
    pub(crate) fn clear_alias_namespace(&mut self, kind: AliasKind, namespace: &str) {
        self.aliases
            .retain(|(k, ns, _)| !(*k == kind && ns == namespace));
        self.forced_aliases
            .retain(|(k, ns, _)| !(*k == kind && ns == namespace));
        self.removed_aliases
            .retain(|(k, ns, _)| !(*k == kind && ns == namespace));
        self.emptied_alias_namespaces
            .retain(|(k, ns)| !(*k == kind && ns == namespace));
    }
}

impl ReadonlyRefgetStore {
    /// Generic constructor. Creates a new, empty `ReadonlyRefgetStore`.
    /// Internal only - users should go through RefgetStore.
    pub(crate) fn new(mode: StorageMode) -> Self {
        ReadonlyRefgetStore {
            catalog: Catalog::new(),
            loaded: HashMap::new(),
            mode,
            local_path: None,
            remote_source: None,
            seqdata_path_template: None,
            persist_to_disk: false,
            quiet: false,
            ancillary_digests: true,
            attribute_index: false,
            aliases: AliasManager::default(),
            fhr_metadata: HashMap::new(),
            available_sequence_alias_namespaces: Vec::new(),
            available_collection_alias_namespaces: Vec::new(),
            seq_fd_cache: Mutex::new(FdCache::new(default_fd_cache_cap())),
            resident: RwLock::new(HashMap::new()),
            sequence_index_loaded: true,
            sequence_index_path: None,
            pending: Mutex::new(PendingChanges::default()),
            commit_lock: None,
            lock_options: LockOptions::default(),
            force_alias: false,
            layout: Layout::Standalone,
            remote_layout: Layout::Standalone,
            pack_cap_bytes: 256 * 1024 * 1024,
            pack_cache: None,
            pack_cols: Mutex::new(None),
            remote_pack_index: Mutex::new(None),
            pack_lengths: Mutex::new(HashMap::new()),
            open_session_pack: Mutex::new(None),
            #[cfg(feature = "http")]
            fetcher: None,
            seq_cache: None,
        }
    }

    /// Network and cache counters of this store since it opened. Zeros for a
    /// local store.
    pub fn fetch_stats(&self) -> super::fetch::FetchCounters {
        #[cfg(feature = "http")]
        if let Some(f) = &self.fetcher {
            return f.counters();
        }
        super::fetch::FetchCounters::default()
    }

    /// Resolve a sha512t24u or md5 digest, with or without `SQ.`.
    #[cfg_attr(not(feature = "filesystem"), allow(dead_code))]
    pub(crate) fn resolve_sequence_key(&self, digest: &str) -> Option<DigestKey> {
        self.seq_row(digest)
            .map(|r| self.catalog.seq(r).sha512t24u().to_key())
    }

    /// Where a sequence's `.seq` file lives on this machine: the local store's
    /// own file, or (remote stores) its file in the shared cache. `None` when
    /// the store has no such file (no path, or a packed local store).
    #[cfg_attr(not(feature = "filesystem"), allow(dead_code))]
    pub(crate) fn local_seq_path(&self, meta: &SequenceMetadata) -> Option<PathBuf> {
        self.local_seq_path_parts(&meta.sha512t24u, meta.alphabet)
    }

    /// [`Self::local_seq_path`] from a digest and alphabet.
    pub(crate) fn local_seq_path_parts(
        &self,
        sha512t24u: &str,
        alphabet: crate::digest::AlphabetType,
    ) -> Option<PathBuf> {
        if let Some(cache) = &self.seq_cache {
            return Some(cache.seq_path(self.mode, alphabet, sha512t24u));
        }
        if self.layout == Layout::Packed {
            return None;
        }
        let local_path = self.local_path.as_ref()?;
        let template = self.seqdata_path_template.as_ref()?;
        Some(local_path.join(Self::expand_template(sha512t24u, template)))
    }

    // =========================================================================
    // Catalog access
    // =========================================================================

    /// Resolve a sequence digest (sha512t24u, `SQ.`-prefixed, or MD5) to its
    /// catalog row.
    pub(crate) fn seq_row<K: AsRef<[u8]>>(&self, digest: K) -> Option<SeqRow> {
        self.catalog.resolve(digest)
    }

    /// A borrowed view of one sequence's metadata, by any digest form.
    pub fn sequence_view<K: AsRef<[u8]>>(&self, digest: K) -> Option<SeqView<'_>> {
        self.seq_row(digest).map(|r| self.catalog.seq(r))
    }

    /// Views of every sequence in the store, in insertion order.
    pub fn sequence_views(&self) -> impl Iterator<Item = SeqView<'_>> + '_ {
        self.catalog.seq_rows().map(|r| self.catalog.seq(r))
    }

    /// Row of a collection whose member list is loaded; an error if the
    /// collection is unknown or not loaded.
    pub(crate) fn loaded_coll_row<K: AsRef<[u8]>>(&self, digest: K) -> Result<CollRow> {
        let d = digest.as_ref();
        match self.catalog.find_coll(d) {
            Some(c) if self.catalog.is_loaded(c) => Ok(c),
            Some(_) => Err(anyhow!(
                "Collection {} not loaded. Call load_collection() or load_all_collections() first.",
                String::from_utf8_lossy(d)
            )),
            None => Err(anyhow!("Collection not found: {}", String::from_utf8_lossy(d))),
        }
    }

    /// Per-collection members (name, description, sequence) in FASTA order.
    /// This is the source of truth for names/headers; a sequence's store-wide
    /// label only names the shared bytes.
    pub fn collection_members<K: AsRef<[u8]>>(
        &self,
        digest: K,
    ) -> Result<impl ExactSizeIterator<Item = MemberView<'_>> + '_> {
        let c = self.loaded_coll_row(digest)?;
        Ok(self.catalog.members(c))
    }

    /// Borrow the store-owned (loaded) bytes of a sequence, in the store's
    /// native encoding. `None` if the sequence is unknown or not loaded.
    pub fn loaded_bytes<K: AsRef<[u8]>>(&self, digest: K) -> Option<&[u8]> {
        let row = self.seq_row(digest)?;
        self.loaded.get(&row).map(|b| b.as_slice())
    }

    /// Store-owned bytes of a sequence, if it is loaded.
    pub(crate) fn loaded_arc(&self, row: SeqRow) -> Option<Arc<Vec<u8>>> {
        self.loaded.get(&row).cloned()
    }

    /// Owned record for a sequence row, with its store-wide label. Bytes are an
    /// `Arc` bump, not a copy.
    pub(crate) fn seq_record(&self, row: SeqRow) -> SequenceRecord {
        let metadata = self.catalog.seq(row).to_metadata();
        match self.loaded_arc(row) {
            Some(sequence) => SequenceRecord::Full { metadata, sequence },
            None => SequenceRecord::Stub(metadata),
        }
    }

    /// Owned record for a collection member: THIS collection's name and
    /// description, the shared bytes.
    pub(crate) fn member_record(&self, m: MemberView<'_>) -> SequenceRecord {
        let metadata = m.to_metadata();
        match self.loaded_arc(m.seq_row()) {
            Some(sequence) => SequenceRecord::Full { metadata, sequence },
            None => SequenceRecord::Stub(metadata),
        }
    }

    /// Insert (or, with `force`, relabel) a sequence and optionally its bytes.
    /// Returns the row and whether the store changed (new row or forced).
    pub(crate) fn insert_sequence(
        &mut self,
        meta: &SequenceMetadata,
        bytes: Option<Arc<Vec<u8>>>,
        force: bool,
    ) -> Result<(SeqRow, bool)> {
        let (row, new) = self.catalog.add_sequence(meta)?;
        if !new && !force {
            return Ok((row, false));
        }
        if !new {
            self.catalog.set_sequence_label(row, meta)?;
        }
        match bytes {
            Some(b) => {
                self.loaded.insert(row, b);
            }
            None => {
                self.loaded.remove(&row);
            }
        }
        Ok((row, true))
    }

    /// Install a collection with its members (sequences added as needed,
    /// without bytes). Returns the collection row and the sequence rows that
    /// were new to the store.
    pub(crate) fn insert_collection(
        &mut self,
        meta: SequenceCollectionMetadata,
        members: &[SequenceMetadata],
    ) -> Result<(CollRow, Vec<SeqRow>)> {
        let mut rows = Vec::with_capacity(members.len());
        let mut new_rows = Vec::new();
        for m in members {
            let (row, new) = self.catalog.add_sequence(m)?;
            if new {
                new_rows.push(row);
            }
            rows.push(row);
        }
        let c = self.catalog.upsert_collection(meta)?;
        self.catalog.set_members(
            c,
            members.iter().zip(&rows).map(|(m, &seq)| NewMember {
                seq,
                name: &m.name,
                description: m.description.as_deref(),
                fai: m.fai.as_ref(),
            }),
        )?;
        Ok((c, new_rows))
    }

    /// Tombstone sequences and drop their bytes.
    pub(crate) fn remove_sequence_rows(&mut self, rows: &[SeqRow]) {
        let mut resident = self.resident.write().ok();
        for &r in rows {
            self.catalog.remove_sequence(r);
            self.loaded.remove(&r);
            if let Some(res) = resident.as_mut() {
                res.remove(&r);
            }
        }
    }

    /// Renumber catalog rows once enough are dead, rekeying the byte maps.
    pub(crate) fn maybe_compact_catalog(&mut self) {
        if !self.catalog.needs_compaction() {
            return;
        }
        let remap = self.catalog.compact();
        let rekey = |m: HashMap<SeqRow, Arc<Vec<u8>>>| -> HashMap<SeqRow, Arc<Vec<u8>>> {
            m.into_iter()
                .filter_map(|(r, b)| remap.get(r as usize).copied().flatten().map(|n| (n, b)))
                .collect()
        };
        self.loaded = rekey(std::mem::take(&mut self.loaded));
        if let Ok(mut res) = self.resident.write() {
            *res = rekey(std::mem::take(&mut *res));
        }
        self.invalidate_pack_index();
    }

    // =========================================================================
    // Write locking
    // =========================================================================

    /// Hold the store's exclusive writer lock across several mutations.
    ///
    /// Individual mutations already commit under the lock, so this is only needed
    /// when a sequence of them must be atomic with respect to other writers — the
    /// registry's alias step, for instance, which adds many aliases and expects
    /// them to appear together. The guard is dropped by
    /// [`Self::release_batch_lock`] or when the store is dropped.
    ///
    /// Nested commits are re-entrant: they see the held guard and skip acquiring.
    pub fn lock_for_batch(&mut self, operation: &str) -> Result<()> {
        if self.commit_lock.is_some() {
            return Ok(());
        }
        let local_path = self
            .local_path
            .as_ref()
            .context("lock_for_batch requires a disk-backed store")?;
        let lock = StoreLock::acquire(local_path, operation, self.lock_options.clone())?;
        self.commit_lock = Some(Arc::new(lock));
        Ok(())
    }

    /// Release a lock taken by [`Self::lock_for_batch`].
    pub fn release_batch_lock(&mut self) {
        self.commit_lock = None;
    }

    /// Whether this store currently holds a batch write lock.
    pub fn holds_batch_lock(&self) -> bool {
        self.commit_lock.is_some()
    }

    /// Whether this store has additions or removals that have not yet been
    /// committed to the shared index files.
    ///
    /// Pending changes live only in memory, so dropping a store without
    /// committing silently abandons them: added rows never reach
    /// `sequences.rgsi`/`collections.rgci`, and removed rows are still on disk.
    pub fn has_uncommitted_changes(&self) -> bool {
        !self.pending.lock().unwrap().is_empty()
    }

    /// Record a change this handle made, so the next commit applies it to the
    /// index it reads fresh from disk.
    pub(crate) fn record<F: FnOnce(&mut PendingChanges)>(&self, f: F) {
        f(&mut self.pending.lock().unwrap());
    }

    /// Override the timeout/staleness settings used when acquiring the write lock.
    pub fn set_lock_options(&mut self, options: LockOptions) {
        self.lock_options = options;
    }

    /// Allow re-pointing an existing alias (in memory, or one another writer
    /// already published) to a different digest, instead of erroring on the
    /// conflict. A load that maps one alias to two digests still errors.
    pub fn set_force_alias(&mut self, force: bool) {
        self.force_alias = force;
    }

    /// Test-only: shrink the partial-read fd-cache to a tiny capacity so the
    /// eviction/reopen path is exercised deterministically. Only used by the
    /// filesystem-gated partial-read tests.
    #[cfg(all(test, feature = "filesystem"))]
    pub(crate) fn set_seq_fd_cache_cap_for_test(&self, cap: usize) {
        let mut cache = self.seq_fd_cache.lock().unwrap();
        *cache = FdCache::new(cap);
    }

    /// Set whether to suppress progress output.
    pub fn set_quiet(&mut self, quiet: bool) {
        self.quiet = quiet;
    }

    /// Returns whether the store is in quiet mode.
    pub fn is_quiet(&self) -> bool {
        self.quiet
    }

    /// Check whether a valid RefgetStore exists at the given path.
    pub fn store_exists<P: AsRef<Path>>(path: P) -> bool {
        path.as_ref().join("rgstore.json").exists()
    }

    /// Change the storage mode, re-encoding/decoding existing sequences as needed.
    ///
    /// Switching to [`StorageMode::Encoded`] fails if a sequence holds a byte
    /// that its stored alphabet cannot represent. That means the store's
    /// alphabet metadata is wrong (written by an older gtars) and the store
    /// must be re-imported. On error the store is left unchanged.
    ///
    /// For a disk-backed PACKED store that already holds sequences, changing the
    /// mode is refused with an error: the packed bytes on disk are in the
    /// current mode, and re-encoding only the resident records would desync
    /// them from the packs. Mode is fixed at creation for a populated packed
    /// store. (An empty packed store, before any import, may still switch,
    /// since there are no packed bytes to desync.)
    pub fn set_encoding_mode(&mut self, new_mode: StorageMode) -> Result<()> {
        if self.mode == new_mode {
            return Ok(());
        }
        if self.persist_to_disk
            && self.layout == Layout::Packed
            && self.catalog.n_sequences() > 0
        {
            return Err(anyhow!(
                "refusing to change storage mode on a populated packed store \
                 (would desync the packs); mode stays {:?}",
                self.mode
            ));
        }
        let old_mode = self.mode;

        let encode = |metadata: &SequenceMetadata, ascii: &[u8]| -> Result<Vec<u8>> {
            let alphabet = lookup_alphabet(&metadata.alphabet);
            encode_sequence(ascii, alphabet).map_err(|e| {
                anyhow::anyhow!(
                    "cannot encode sequence {} ({}): {e}. The stored alphabet \
                     metadata is wrong; re-import this store.",
                    metadata.sha512t24u,
                    metadata.name
                )
            })
        };

        // Convert everything first, then swap in, so an error leaves the store
        // untouched rather than half converted.
        let mut converted: Vec<(SeqRow, Vec<u8>)> = Vec::new();
        for (key, sequence) in self.loaded.iter() {
            let owned_meta = self.catalog.seq(*key).to_metadata();
            let metadata = &owned_meta;
            let new_bytes = match (old_mode, new_mode) {
                (StorageMode::Raw, StorageMode::Encoded) => encode(metadata, &sequence[..])?,
                (StorageMode::Encoded, StorageMode::Raw) => {
                    let alphabet = lookup_alphabet(&metadata.alphabet);
                    decode_string_from_bytes(&sequence[..], metadata.length, alphabet)
                }
                // Zstd conversions: normalize through ASCII. The stored
                // body for Zstd is the whole-record zstd frame.
                (StorageMode::Raw, StorageMode::Zstd) => {
                    zstd::encode_all(&sequence[..], super::ZSTD_STORAGE_LEVEL)
                        .expect("zstd encode during mode switch")
                }
                (StorageMode::Zstd, StorageMode::Raw) => {
                    zstd::decode_all(&sequence[..]).expect("zstd decode during mode switch")
                }
                (StorageMode::Encoded, StorageMode::Zstd) => {
                    let alphabet = lookup_alphabet(&metadata.alphabet);
                    let ascii =
                        decode_string_from_bytes(&sequence[..], metadata.length, alphabet);
                    zstd::encode_all(&ascii[..], super::ZSTD_STORAGE_LEVEL)
                        .expect("zstd encode during mode switch")
                }
                (StorageMode::Zstd, StorageMode::Encoded) => {
                    let ascii = zstd::decode_all(&sequence[..])
                        .expect("zstd decode during mode switch");
                    encode(metadata, &ascii[..])?
                }
                _ => continue,
            };
            converted.push((*key, new_bytes));
        }

        for (key, new_bytes) in converted {
            self.loaded.insert(key, std::sync::Arc::new(new_bytes));
        }

        self.mode = new_mode;
        Ok(())
    }

    /// Enable 2-bit encoding for space efficiency.
    ///
    /// Fails, leaving the store unchanged, if a sequence cannot be encoded
    /// with its stored alphabet (see [`set_encoding_mode`](Self::set_encoding_mode)).
    pub fn enable_encoding(&mut self) -> Result<()> {
        self.set_encoding_mode(StorageMode::Encoded)
    }

    /// Disable encoding, use raw byte storage.
    pub fn disable_encoding(&mut self) {
        self.set_encoding_mode(StorageMode::Raw)
            .expect("switching to Raw never encodes, so it cannot fail");
    }

    /// Enable disk persistence for this store.
    pub fn enable_persistence<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let path = path.as_ref();

        self.local_path = Some(path.to_path_buf());
        self.persist_to_disk = true;
        self.seqdata_path_template
            .get_or_insert_with(|| DEFAULT_SEQDATA_PATH_TEMPLATE.to_string());

        create_dir_all(path.join("sequences"))?;
        create_dir_all(path.join("collections"))?;

        // Materializing an in-memory store onto disk: EVERYTHING resident is a
        // change, because none of it came from this directory. This is the one
        // place where "write back the whole in-memory pile" is the correct
        // commit, and it is correct precisely because nothing here was loaded.
        {
            let mut pending = self.pending.lock().unwrap();
            for v in self.sequence_views() {
                pending.add_sequence(v.sha512t24u().to_key());
            }
            for c in self.catalog.coll_rows() {
                pending.add_collection(self.catalog.coll_meta(c).digest.to_key());
            }
            for kind in [AliasKind::Sequence, AliasKind::Collection] {
                for ns in self.aliases.namespaces_for(kind) {
                    for alias in self.aliases.namespace_map(kind, &ns).keys() {
                        pending.add_alias(kind, &ns, alias);
                    }
                }
            }
        }

        let rows: Vec<SeqRow> = self.loaded.keys().copied().collect();
        for row in rows {
            if let Some(sequence) = self.loaded.remove(&row) {
                let metadata = self.catalog.seq(row).to_metadata();
                self.write_sequence_to_disk_single(&metadata, &sequence)?;
            }
        }

        for c in self.catalog.coll_rows() {
            self.write_collection_to_disk_single(
                self.catalog.coll_meta(c),
                self.catalog.members(c).map(|m| m.to_metadata()),
            )?;
        }

        self.write_index_files()?;

        Ok(())
    }

    /// Disable disk persistence for this store.
    pub fn disable_persistence(&mut self) {
        self.persist_to_disk = false;
    }

    /// Check if persistence to disk is enabled.
    pub fn is_persisting(&self) -> bool {
        self.persist_to_disk
    }

    /// Adds a collection, and all sequences in it, to the store.
    pub fn add_sequence_collection(
        &mut self,
        collection: crate::digest::SequenceCollection,
    ) -> Result<()> {
        self.add_sequence_collection_internal(collection, false)
    }

    /// Adds a collection, overwriting existing data.
    pub fn add_sequence_collection_force(
        &mut self,
        collection: crate::digest::SequenceCollection,
    ) -> Result<()> {
        self.add_sequence_collection_internal(collection, true)
    }

    /// Internal implementation for adding a sequence collection.
    pub(crate) fn add_sequence_collection_internal(
        &mut self,
        collection: crate::digest::SequenceCollection,
        force: bool,
    ) -> Result<()> {
        let coll_digest = collection.metadata.digest.to_key();

        if !force && self.catalog.find_coll(&collection.metadata.digest).is_some() {
            return Ok(());
        }

        let crate::digest::SequenceCollection { metadata, sequences } = collection;
        let members: Vec<SequenceMetadata> =
            sequences.iter().map(|s| s.metadata().clone()).collect();

        if self.persist_to_disk && self.local_path.is_some() {
            self.write_collection_to_disk_single(&metadata, &members)?;
        }

        for sequence_record in sequences {
            self.add_sequence_record(sequence_record, force)?;
        }
        self.insert_collection(metadata, &members)?;
        self.record(|p| p.add_collection(coll_digest));

        if self.persist_to_disk && self.local_path.is_some() {
            self.write_index_files()?;
        }

        Ok(())
    }

    /// Pack a raw ASCII record when the store uses encoded storage.
    ///
    /// Raw and packed forms are distinguished by their expected lengths, so a
    /// single-base sequence is ambiguous. Use [`ingest_sequence`](Self::ingest_sequence)
    /// when the input is known to be raw; it handles that case explicitly.
    fn pack_raw_ascii_for_encoded_mode(&self, sr: SequenceRecord) -> Result<SequenceRecord> {
        if self.mode != StorageMode::Encoded {
            return Ok(sr);
        }
        match sr {
            SequenceRecord::Full { metadata, sequence } => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                let bps = alphabet.bits_per_symbol;
                let expected_encoded = metadata.length.saturating_mul(bps).div_ceil(8);
                let looks_raw = sequence.len() == metadata.length
                    && sequence.len() != expected_encoded
                    && sequence.is_ascii();
                if looks_raw {
                    let packed = encode_sequence(&sequence[..], alphabet).with_context(|| {
                        format!(
                            "encoding sequence {} ({}) with alphabet {}",
                            metadata.sha512t24u, metadata.name, metadata.alphabet
                        )
                    })?;
                    Ok(SequenceRecord::Full {
                        metadata,
                        sequence: std::sync::Arc::new(packed),
                    })
                } else {
                    Ok(SequenceRecord::Full { metadata, sequence })
                }
            }
            stub @ SequenceRecord::Stub(_) => Ok(stub),
        }
    }

    /// Digest and insert raw bytes using the store's native encoding.
    ///
    /// Because the input is known to be raw, this method also handles
    /// single-base sequences, whose raw and packed lengths are identical.
    /// Input is normalized in the same way as `digest_sequence`.
    pub fn ingest_sequence(&mut self, name: &str, bytes: &[u8]) -> Result<String> {
        let record = crate::digest::digest_sequence(name, bytes);
        let digest = record.metadata().sha512t24u.clone();

        let record = match record {
            SequenceRecord::Full { metadata, sequence } if self.mode == StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                let packed = encode_sequence(&sequence[..], alphabet).with_context(|| {
                    format!("encoding sequence {} ({})", metadata.sha512t24u, metadata.name)
                })?;
                SequenceRecord::Full {
                    metadata,
                    sequence: std::sync::Arc::new(packed),
                }
            }
            other => other,
        };

        self.add_sequence_record(record, false)?;
        Ok(digest)
    }

    /// Adds a SequenceRecord directly to the store without collection association.
    ///
    /// Raw ASCII records are packed when the store uses encoded storage.
    /// Already-packed records are stored unchanged.
    pub fn add_sequence_record(&mut self, sr: SequenceRecord, force: bool) -> Result<()> {
        let sr = self.pack_raw_ascii_for_encoded_mode(sr)?;

        let metadata = sr.metadata();
        let key = metadata.sha512t24u.to_key();

        // Records must now match the store's native representation: ASCII in
        // Raw mode or the alphabet's packed size in Encoded mode.
        debug_assert!(
            match &sr {
                SequenceRecord::Full { metadata, sequence } => match self.mode {
                    StorageMode::Raw => sequence.is_ascii(),
                    StorageMode::Encoded => {
                        let bps = lookup_alphabet(&metadata.alphabet).bits_per_symbol;
                        sequence.len() == metadata.length.saturating_mul(bps).div_ceil(8)
                    }
                    StorageMode::Zstd => true, // zstd bytes have no simple invariant to check
                },
                SequenceRecord::Stub(_) => true,
            },
            "add_sequence_record: sequence bytes must be stored in the store's native encoding (ASCII in Raw mode, packed to the alphabet's encoded size in Encoded mode)"
        );

        // A dedup hit records NOTHING: the row is already in memory, which for a
        // disk-backed store means it is already in `sequences.rgsi`. Re-asserting
        // it would put a row this handle did not produce back into the delta.
        if !force && self.catalog.find_sha(&metadata.sha512t24u).is_some() {
            return Ok(());
        }

        // In Zstd mode the caller passes the ASCII record body (per contract,
        // "do not 2-bit pre-pack, let zstd own the entropy"); the store owns the
        // whole-record zstd compression here so the on-disk / in-memory body is
        // the compressed frame. Reads decompress the whole frame and slice.
        let sr = self.compress_record_for_storage(sr)?;
        let (metadata, bytes) = match sr {
            SequenceRecord::Full { metadata, sequence } => (metadata, Some(sequence)),
            SequenceRecord::Stub(metadata) => (metadata, None),
        };

        let writes_bytes = self.persist_to_disk && self.local_path.is_some() && bytes.is_some();
        // What the insert below changes, so a failed byte write can put it
        // back: the prior label and bytes of a sequence already in the catalog
        // (only a forced re-add gets here with one), and whether this handle
        // had already staged the sequence or a removal of it.
        let undo = writes_bytes.then(|| {
            let prior = self
                .catalog
                .find_sha(&metadata.sha512t24u)
                .map(|r| (self.catalog.seq(r).to_metadata(), self.loaded.get(&r).cloned()));
            let p = self.pending.lock().unwrap();
            (prior, p.sequences.contains(&key), p.removed_sequences.contains(&key))
        });

        let (row, _) = self.insert_sequence(&metadata, bytes.clone(), true)?;
        self.record(|p| p.add_sequence(key));

        if let (Some(sequence), Some((prior, was_staged, was_removed))) = (&bytes, undo) {
            let packed = self.layout == Layout::Packed;
            let written = if packed {
                // Packed incremental write (panget's path): append to the open
                // session pack. The bytes stay loaded so mid-session reads are
                // served from RAM; on reopen the sequence is resolved through
                // the sidecar.
                self.session_pack_append(&metadata.sha512t24u, sequence)
            } else {
                self.write_sequence_to_disk_single(&metadata, sequence)
            };
            if let Err(e) = written {
                // Undo the insert, so a later commit publishes no row whose
                // bytes are not on disk. A sequence this call added is
                // removed; one already in the catalog keeps its row and gets
                // its prior label and bytes back.
                match prior {
                    None => self.remove_sequence_rows(&[row]),
                    Some((meta, prior_bytes)) => {
                        let _ = self.catalog.set_sequence_label(row, &meta);
                        match prior_bytes {
                            Some(b) => self.loaded.insert(row, b),
                            None => self.loaded.remove(&row),
                        };
                    }
                }
                let mut p = self.pending.lock().unwrap();
                if !was_staged {
                    p.sequences.remove(&key);
                }
                if was_removed {
                    p.removed_sequences.insert(key);
                }
                return Err(e);
            }
            if !packed {
                self.loaded.remove(&row);
            }
        }
        Ok(())
    }

    /// For a `Zstd`-mode store, compress a `Full` record's ASCII body into a zstd
    /// frame (the stored representation). A no-op for every other mode and for
    /// `Stub` records. `metadata.length` stays the logical (ASCII) base count.
    fn compress_record_for_storage(&self, sr: SequenceRecord) -> Result<SequenceRecord> {
        if self.mode != StorageMode::Zstd {
            return Ok(sr);
        }
        match sr {
            SequenceRecord::Full { metadata, sequence } => {
                let compressed = zstd::encode_all(&sequence[..], super::ZSTD_STORAGE_LEVEL)
                    .context("zstd compression of sequence record failed")?;
                Ok(SequenceRecord::Full {
                    metadata,
                    sequence: std::sync::Arc::new(compressed),
                })
            }
            other => Ok(other),
        }
    }

    /// Inserter-side half of [`add_sequence_record`](Self::add_sequence_record)
    /// for the parallel streaming import: do the FAST in-memory dedup decision
    /// and catalog insert here (single inserter, so no race on the
    /// dedup truth), but DO NOT write the bytes to disk. Instead return the
    /// expanded full path + bytes, so the caller can dispatch the (independent,
    /// per-digest) `.seq` write to a bounded writer pool, or append the bytes
    /// to the session pack of a packed store.
    ///
    /// Returns `Some((full_path, bytes))` when a write must be performed
    /// (disk-backed store, `Full` record, and either a first occurrence or a
    /// forced overwrite), or `None` when nothing needs to be written to disk
    /// (in-memory mode, a `Stub`, or a dedup hit). The in-memory state is fully
    /// updated before returning, so dedup stays the single source of truth on
    /// the inserter thread.
    #[cfg_attr(not(feature = "filesystem"), allow(dead_code))]
    pub(crate) fn add_sequence_record_deferred_write(
        &mut self,
        sr: SequenceRecord,
        force: bool,
    ) -> Result<Option<(PathBuf, Vec<u8>)>> {
        let metadata = sr.metadata();
        let key = metadata.sha512t24u.to_key();

        // Dedup hit: nothing inserted, so nothing to record. See
        // [`Self::add_sequence_record`].
        if !force && self.catalog.find_sha(&metadata.sha512t24u).is_some() {
            return Ok(None);
        }

        if self.persist_to_disk && self.local_path.is_some() {
            match sr {
                SequenceRecord::Full { metadata, sequence } => {
                    // Insert the metadata NOW (reserve), so the dedup
                    // decision is committed on the inserter thread. The caller
                    // does the actual byte write.
                    let full_path = self
                        .sequence_file_path(&metadata.sha512t24u)
                        .context("seqdata_path_template/local_path not set")?;
                    self.insert_sequence(&metadata, None, true)?;
                    self.record(|p| p.add_sequence(key));
                    // INVARIANT: strong_count == 1 here — the Arc was just
                    // constructed in the import pipeline (Arc::new in import.rs)
                    // and no other owner exists at the point of this call, so
                    // try_unwrap always succeeds and reclaims the buffer with zero
                    // copy. The clone fallback is dead code / future-proofing only.
                    debug_assert_eq!(
                        std::sync::Arc::strong_count(&sequence),
                        1,
                        "sequence Arc must be uniquely owned for zero-copy reclaim"
                    );
                    let bytes = std::sync::Arc::try_unwrap(sequence)
                        .unwrap_or_else(|arc| (*arc).clone());
                    // Zstd mode: the caller handed us the ASCII body; compress the
                    // whole record so the zstd frame is what gets written.
                    let bytes = if self.mode == StorageMode::Zstd {
                        zstd::encode_all(&bytes[..], super::ZSTD_STORAGE_LEVEL)
                            .context("zstd compression of sequence record failed")?
                    } else {
                        bytes
                    };
                    return Ok(Some((full_path, bytes)));
                }
                SequenceRecord::Stub(s) => {
                    self.insert_sequence(&s, None, true)?;
                    self.record(|p| p.add_sequence(key));
                    return Ok(None);
                }
            }
        }

        // In-memory (non-disk-backed): keep the resident body in the store's
        // native form, so a Zstd store holds the compressed frame just like the
        // disk-backed path returns compressed bytes above.
        let (metadata, bytes) = match self.compress_record_for_storage(sr)? {
            SequenceRecord::Full { metadata, sequence } => (metadata, Some(sequence)),
            SequenceRecord::Stub(metadata) => (metadata, None),
        };
        self.insert_sequence(&metadata, bytes, true)?;
        self.record(|p| p.add_sequence(key));
        Ok(None)
    }

    // =========================================================================
    // Sequence query methods
    // =========================================================================

    /// Returns an iterator over all sequence digests in the store
    pub fn sequence_digests(&self) -> impl Iterator<Item = DigestKey> + '_ {
        self.sequence_views().map(|v| v.sha512t24u().to_key())
    }

    /// Returns an iterator over (owned) sequence metadata for all sequences in
    /// the store. Use [`Self::sequence_views`] to avoid the per-item copy.
    pub fn sequence_metadata(&self) -> impl Iterator<Item = SequenceMetadata> + '_ {
        self.sequence_views().map(|v| v.to_metadata())
    }

    /// Logical size in bytes of all sequence *payloads*, computed from metadata
    /// (length x storage mode). Cheap: O(n_sequences) arithmetic, no I/O, no walk,
    /// works on stubs. Does NOT include index files, sidecars, or rgstore.json.
    /// This is the value cached in the manifest at write time (see write path); for
    /// the exact full on-disk footprint use `actual_disk_usage()` (walks the dir).
    pub fn logical_sequence_bytes(&self) -> usize {
        self.sequence_views().map(|v| v.disk_size(&self.mode)).sum()
    }

    /// Exact on-disk footprint of the store: recursively walks the store directory
    /// and sums every file's byte length (`.seq` data, `.rgsi`/`.rgci` indexes, FHR
    /// sidecars, alias data, `rgstore.json`). Accurate and complete, but does I/O on
    /// every call. Use this when an exact footprint is needed; for the cheap cached
    /// sequence-content size use `logical_sequence_bytes()`.
    pub fn actual_disk_usage(&self) -> usize {
        let Some(path) = &self.local_path else {
            return 0;
        };

        fn dir_size(path: &std::path::Path) -> usize {
            let mut total = 0;
            if let Ok(entries) = std::fs::read_dir(path) {
                for entry in entries.flatten() {
                    let path = entry.path();
                    if path.is_file() {
                        total += entry.metadata().map(|m| m.len() as usize).unwrap_or(0);
                    } else if path.is_dir() {
                        total += dir_size(&path);
                    }
                }
            }
            total
        }

        dir_size(path)
    }

    // =========================================================================
    // Collection API
    // =========================================================================

    /// List collections with pagination and optional attribute filtering.
    pub fn list_collections(
        &self,
        page: usize,
        page_size: usize,
        filters: &[(&str, &str)],
    ) -> Result<PagedResult<SequenceCollectionMetadata>> {
        let mut filtered: Vec<SequenceCollectionMetadata> = Vec::new();
        for c in self.catalog.coll_rows() {
            let meta = self.catalog.coll_meta(c);
            let mut passes = true;
            for &(attr_name, attr_digest) in filters {
                if !metadata_matches_attribute(meta, attr_name, attr_digest)? {
                    passes = false;
                    break;
                }
            }
            if passes {
                filtered.push(meta.clone());
            }
        }

        filtered.sort_by(|a, b| a.digest.cmp(&b.digest));

        let total = filtered.len();
        let start = page * page_size;
        let results = if start < total {
            filtered.into_iter().skip(start).take(page_size).collect()
        } else {
            Vec::new()
        };

        Ok(PagedResult {
            results,
            pagination: Pagination {
                page,
                page_size,
                total,
            },
        })
    }

    /// Get metadata for a single collection by digest (no sequence data).
    pub fn get_collection_metadata<K: AsRef<[u8]>>(
        &self,
        collection_digest: K,
    ) -> Option<&SequenceCollectionMetadata> {
        self.catalog
            .find_coll(collection_digest)
            .map(|c| self.catalog.coll_meta(c))
    }

    /// Get a collection with all its sequences loaded.
    pub fn get_collection(&self, collection_digest: &str) -> Result<crate::digest::SequenceCollection> {
        let c = self.loaded_coll_row(collection_digest)?;
        // Each result carries the per-collection name AND description, not the
        // first-imported store-wide label of a sequence shared across
        // collections. Only the bytes (and `fai`, if the member lacks it) come
        // from the sequence row.
        let sequences: Vec<SequenceRecord> =
            self.catalog.members(c).map(|m| self.member_record(m)).collect();
        Ok(crate::digest::SequenceCollection {
            metadata: self.catalog.coll_meta(c).clone(),
            sequences,
        })
    }

    /// Dry-run of [`Self::remove_collection`]'s orphan cleanup: the sequence
    /// digests that WOULD be deleted, without touching anything.
    ///
    /// Uses the same disk-derived live set as the real thing, so it fails with the
    /// same error if the store is not in a state where orphan GC is safe.
    ///
    /// # Advisory, by design
    ///
    /// This takes NO lock. It exists for confirmation prompts ("this will free
    /// 21 sequences, proceed?"), where blocking every concurrent writer for the
    /// length of a full-store scan to answer a question the user may say no to
    /// would be absurd. The authoritative scan runs again inside
    /// `remove_collection`, under the lock, and the two can legitimately differ:
    /// a collection committed in between makes some planned orphan live again,
    /// and the real removal will correctly spare it. Callers comparing the two
    /// should treat a shortfall as normal concurrency, not corruption.
    pub fn plan_orphan_removal(&self, digest: &str) -> Result<Vec<String>> {
        let key = digest.to_key();
        let candidates = self.orphan_candidates(digest)?;
        let live = self.live_sequence_set(&key)?;
        Ok(candidates
            .into_iter()
            .filter(|k| !live.contains(k))
            .map(|k| key_to_digest_string(&k))
            .collect())
    }

    /// Sequence digests belonging to the collection being removed.
    fn orphan_candidates(&self, digest: &str) -> Result<Vec<DigestKey>> {
        if self.local_path.is_some() {
            self.collection_sequence_digests(digest)
        } else {
            Ok(self.loaded_member_seq_keys(digest))
        }
    }

    /// Sequence digest keys of a loaded collection's members (empty if the
    /// collection is unknown or not loaded).
    pub(crate) fn loaded_member_seq_keys(&self, digest: &str) -> Vec<DigestKey> {
        match self.catalog.find_coll(digest) {
            Some(c) => self
                .catalog
                .members(c)
                .map(|m| m.seq().sha512t24u().to_key())
                .collect(),
            None => Vec::new(),
        }
    }

    /// Sequence digests still referenced by some OTHER collection.
    ///
    /// Disk-backed stores derive this from disk. Only a pure in-memory store
    /// (`local_path == None`) falls back to the loaded collection members:
    /// there, the in-memory catalog is all there is, so it is authoritative.
    fn live_sequence_set(&self, exclude: &DigestKey) -> Result<HashSet<DigestKey>> {
        if self.local_path.is_some() {
            return self.live_sequence_digests_from_disk(exclude);
        }
        let mut live = HashSet::new();
        for c in self.catalog.coll_rows() {
            if self.catalog.coll_meta(c).digest.to_key() == *exclude {
                continue;
            }
            live.extend(self.catalog.members(c).map(|m| m.seq().sha512t24u().to_key()));
        }
        Ok(live)
    }

    /// Remove a collection from the store.
    ///
    /// When `remove_orphan_sequences` is true, the orphan-cleanup path is
    /// O(sequences + collections): each orphan is one catalog lookup, never a
    /// scan. Do not reintroduce a per-orphan scan here; on large stores that is
    /// quadratic and effectively unbounded.
    ///
    /// The live set is derived from DISK, not from the loaded members; see
    /// [`Self::live_sequence_digests_from_disk`] for why that distinction is the
    /// difference between a GC and a data-loss bug.
    ///
    /// # The lock is held across the orphan scan, deliberately
    ///
    /// The scan is the expensive part — on the plantref store (147 collections,
    /// ~1.5M sequence digests) it reads every `collections/*.rgsi` and takes
    /// roughly a minute. Everywhere else in this module the rule is "do the
    /// expensive work outside the lock"; here that rule is wrong. Computing the
    /// live set before taking the lock leaves a window in which another writer
    /// commits a collection referencing one of the digests the scan just
    /// classified as an orphan, and the unlink that follows deletes live data.
    /// Removals are rare and hand-run, so blocking concurrent builds for a
    /// minute is the cheap side of that trade.
    pub fn remove_collection(
        &mut self,
        digest: &str,
        remove_orphan_sequences: bool,
    ) -> Result<bool> {
        let key = digest.to_key();

        if self.catalog.find_coll(digest).is_none() {
            return Ok(false);
        }

        // Take the lock BEFORE the orphan scan and hold it through the scan, the
        // index commit, and the unlinks. All three have to see the same store: a
        // collection another writer adds in between would make an "orphan" live
        // again, and we would unlink it anyway. The nested `write_index_files`
        // sees the guard and does not re-acquire.
        let acquired = self.acquire_commit_lock("remove_collection")?;
        let restore = match acquired {
            Some(lock) => {
                let previous = self.commit_lock.replace(Arc::new(lock));
                Some(previous)
            }
            // Already inside a caller's `lock_for_batch`; leave it alone.
            None => None,
        };

        let result = self.remove_collection_locked(digest, key, remove_orphan_sequences);

        if let Some(previous) = restore {
            self.commit_lock = previous;
        }
        result
    }

    fn remove_collection_locked(
        &mut self,
        digest: &str,
        key: DigestKey,
        remove_orphan_sequences: bool,
    ) -> Result<bool> {
        // Compute both sets BEFORE mutating anything: the candidates come from
        // the collection's own `.rgsi`, which the cleanup below unlinks.
        let orphan_candidates: Vec<DigestKey> = if remove_orphan_sequences {
            self.orphan_candidates(digest)?
        } else {
            Vec::new()
        };
        let still_referenced: HashSet<DigestKey> = if remove_orphan_sequences {
            self.live_sequence_set(&key)?
        } else {
            HashSet::new()
        };

        if let Some(c) = self.catalog.find_coll(digest) {
            self.catalog.remove_collection(c);
        }
        self.fhr_metadata.remove(&key);
        self.record(|p| p.remove_collection(key));

        // Remove collection aliases pointing to this digest. The in-memory
        // alias table is an open-time snapshot; an alias another process
        // published for this collection since then is only on disk, and
        // leaving it behind makes it resolve to nothing. Under the lock, also
        // scan the published TSVs for every entry targeting this digest.
        let mut alias_pairs = self.aliases.reverse_lookup_collection(digest);
        if self.persist_to_disk {
            if let Some(local_path) = &self.local_path {
                let dir = local_path.join("aliases").join(AliasKind::Collection.subdir());
                if let Ok(entries) = fs::read_dir(&dir) {
                    for entry in entries.filter_map(|e| e.ok()) {
                        let path = entry.path();
                        if path.extension().map_or(true, |e| e != "tsv") {
                            continue;
                        }
                        let Some(ns) = path.file_stem().and_then(|s| s.to_str()) else {
                            continue;
                        };
                        for (alias, target) in super::persistence::read_alias_tsv(&path)? {
                            if target == key {
                                alias_pairs.push((ns.to_string(), alias));
                            }
                        }
                    }
                }
            }
        }
        alias_pairs.sort();
        alias_pairs.dedup();
        let affected_namespaces: HashSet<String> = alias_pairs
            .iter()
            .map(|(ns, _)| ns.clone())
            .collect();
        for (ns, alias) in &alias_pairs {
            self.aliases.remove_collection(ns, alias);
            self.record_alias_removal(AliasKind::Collection, ns, alias);
        }
        for ns in &affected_namespaces {
            self.persist_alias_namespace(AliasKind::Collection, ns)?;
        }

        // Which sequences to unlink, decided under the lock from the live set
        // computed under the same lock. A sequence that a loaded collection
        // still uses is kept even when the on-disk index does not list that
        // collection: the removed collection's members are already dropped
        // above, so any remaining reference comes from another collection.
        let orphans: Vec<DigestKey> = if remove_orphan_sequences {
            orphan_candidates
                .into_iter()
                .filter(|k| !still_referenced.contains(k))
                .filter(|k| {
                    self.catalog
                        .find_sha(key_to_digest_string(k))
                        .is_none_or(|r| self.catalog.seq_refcount(r) == 0)
                })
                .collect()
        } else {
            Vec::new()
        };

        if !orphans.is_empty() {
            {
                let mut pending = self.pending.lock().unwrap();
                for orphan_key in &orphans {
                    pending.remove_sequence(*orphan_key);
                }
            }

            let rows: Vec<SeqRow> = orphans
                .iter()
                .filter_map(|k| self.catalog.find_sha(key_to_digest_string(k)))
                .collect();
            self.remove_sequence_rows(&rows);
        }
        self.maybe_compact_catalog();

        // INDEX FIRST, FILES SECOND. Both orderings can be interrupted by a
        // crash; only one of the two resulting states is recoverable. Rows gone,
        // files present is garbage on disk that a later GC can reclaim. Rows
        // present, files gone is a store that lies: `get_substring` raises,
        // `get_collection` cannot open its sequences, and an alias resolves to a
        // collection that no longer exists.
        if self.persist_to_disk {
            // The packs the sidecar lists before this commit, read under the
            // lock this removal holds: the only packs it may unlink.
            let packs_before = match (&self.local_path, self.layout) {
                (Some(local_path), Layout::Packed) => Self::sidecar_pack_names(local_path)?,
                _ => HashSet::new(),
            };
            self.write_index_files()?;

            if self.layout == Layout::Packed {
                // Packed: no per-`.seq` unlink. Orphan rows were already dropped
                // from the sidecar by the commit above; here we eagerly unlink
                // any listed pack that no live sidecar row references any more.
                if let Some(local_path) = self.local_path.clone() {
                    self.remove_fully_dead_packs(&local_path, &packs_before)?;
                }
            } else {
                self.remove_seq_files(&orphans);
            }

            if let Some(local_path) = &self.local_path {
                let rgsi_path = local_path.join(format!("collections/{}.rgsi", digest));
                let _ = fs::remove_file(&rgsi_path);
                let fhr_path = local_path.join(format!("fhr/{}.fhr.json", digest));
                let _ = fs::remove_file(&fhr_path);
            }
        }

        Ok(true)
    }

    /// Remove every sequence this handle staged (and its `.seq` file on disk)
    /// that no collection in the catalog references.
    ///
    /// Used as a best-effort cleanup after a failed import:
    /// [`add_sequence_collections_from_fastas`](Self::add_sequence_collections_from_fastas)
    /// may have inserted sequences for a partially-built collection that never
    /// received an `End` message (and was therefore never registered). Those
    /// sequences are orphans, with no owning collection, and this method
    /// removes them so the store stays internally consistent. The operation is
    /// O(sequences + collections): each orphan is one catalog lookup.
    ///
    /// Sequences shared with a successfully-finalized collection are preserved.
    ///
    /// Only sequences THIS import staged (still pending, never committed) are
    /// candidates. Anything already in the on-disk index is out of bounds:
    /// `open_local` loads every sequence row but no collection members, so
    /// "a sequence no loaded collection references" describes the entire
    /// pre-existing store, not an orphan. Deriving candidates from that would
    /// unlink every published `.seq` file after one failed import into a
    /// reopened store.
    #[cfg_attr(not(feature = "filesystem"), allow(dead_code))]
    pub(crate) fn remove_orphan_seq_files(&mut self) {
        // Build the set of all sequence digest keys that are referenced by at
        // least one surviving (finalized) collection.
        let mut referenced: std::collections::HashSet<DigestKey> =
            std::collections::HashSet::new();
        for c in self.catalog.coll_rows() {
            referenced.extend(self.catalog.members(c).map(|m| m.seq().sha512t24u().to_key()));
        }

        // Candidates: staged by this handle and not yet published.
        let staged: Vec<DigestKey> = self
            .pending
            .lock()
            .unwrap()
            .sequences
            .iter()
            .copied()
            .collect();

        // Collect the orphan keys (staged but not referenced).
        let orphans: Vec<DigestKey> = staged
            .into_iter()
            .filter(|k| !referenced.contains(k))
            .collect();

        if orphans.is_empty() {
            return;
        }

        // Remove from in-memory structures.
        let rows: Vec<SeqRow> = orphans
            .iter()
            .filter_map(|k| self.catalog.find_sha(key_to_digest_string(k)))
            .collect();
        self.remove_sequence_rows(&rows);
        let orphan_set: std::collections::HashSet<DigestKey> = orphans.iter().copied().collect();

        // FORGET these, do not tombstone them. They were staged by an import
        // that failed before its single commit, so no row of ours was ever
        // published; recording a removal would delete the identical row a
        // different writer legitimately published, since sequence rows are
        // content-addressed and therefore shared between collections.
        self.record(|p| p.forget_sequences(&orphan_set));
        self.forget_staged_aliases(AliasKind::Sequence, &orphan_set);

        // Best-effort remove on-disk `.seq` files.
        if self.persist_to_disk {
            self.remove_seq_files(&orphans);
        }
    }

    /// Best-effort unlink of the standalone `.seq` files of `keys`, then of
    /// each shard directory they leave empty.
    ///
    /// Collect parent shard dirs while unlinking, then rmdir them ONCE at the
    /// end. Do NOT move the rmdir back into the loop: `.seq` files are sharded
    /// over ~4096 two-char prefix dirs, so a per-file rmdir issues N syscalls
    /// (nearly all failing with ENOTEMPTY) to remove at most ~4096 dirs,
    /// roughly doubling the syscall count in a loop that is filesystem-metadata
    /// bound. The rmdir pass must run AFTER all unlinks, or dirs that only
    /// become empty later would be skipped.
    fn remove_seq_files(&self, keys: &[DigestKey]) {
        let (Some(local_path), Some(template)) = (&self.local_path, &self.seqdata_path_template)
        else {
            return;
        };
        let mut parent_dirs: HashSet<PathBuf> = HashSet::new();
        for key in keys {
            let full = local_path.join(Self::expand_template(&key_to_digest_string(key), template));
            let _ = fs::remove_file(&full);
            if let Some(parent) = full.parent() {
                parent_dirs.insert(parent.to_path_buf());
            }
        }
        for parent in &parent_dirs {
            let _ = fs::remove_dir(parent); // ignore if non-empty
        }
    }

    /// Drop every collection this handle staged that has a member in `seqs`,
    /// with the collection aliases staged for it, without recording a removal.
    ///
    /// Used after a failed import for sequences whose bytes never reached disk
    /// (no sealed pack, or no `.seq` file): a collection that names them cannot
    /// be committed.
    /// [`Self::remove_orphan_seq_files`] then drops the sequences themselves.
    #[cfg_attr(not(feature = "filesystem"), allow(dead_code))]
    pub(crate) fn forget_collections_using(&mut self, seqs: &HashSet<DigestKey>) {
        let staged = self.pending.lock().unwrap().collections.clone();
        let doomed: Vec<CollRow> = self
            .catalog
            .coll_rows()
            .filter(|&c| staged.contains(&self.catalog.coll_meta(c).digest.to_key()))
            .filter(|&c| {
                self.catalog
                    .members(c)
                    .any(|m| seqs.contains(&m.seq().sha512t24u().to_key()))
            })
            .collect();
        let keys: HashSet<DigestKey> = doomed
            .iter()
            .map(|&c| self.catalog.coll_meta(c).digest.to_key())
            .collect();
        for c in doomed {
            self.catalog.remove_collection(c);
        }
        self.pending
            .lock()
            .unwrap()
            .collections
            .retain(|k| !keys.contains(k));
        self.forget_staged_aliases(AliasKind::Collection, &keys);

        // Each collection's `.rgsi` was written when it finalized. Remove it
        // unless the published index lists the collection (another writer
        // committed the same one): a store whose index lists no collections
        // loads them from `collections/` on open.
        if let (true, Some(local_path)) = (self.persist_to_disk, &self.local_path) {
            let index = local_path.join("collections.rgci");
            if let Ok(published) = super::persistence::read_collections_rgci_rows(&index) {
                for key in keys.iter().filter(|k| !published.contains_key(*k)) {
                    let rgsi = format!("collections/{}.rgsi", key_to_digest_string(key));
                    let _ = fs::remove_file(local_path.join(rgsi));
                }
            }
        }
    }

    /// Unbind the aliases of `kind` this handle staged but never published whose
    /// target is in `targets`, without recording a removal.
    fn forget_staged_aliases(&mut self, kind: AliasKind, targets: &HashSet<DigestKey>) {
        let mut pending = self.pending.lock().unwrap();
        let doomed: Vec<(AliasKind, String, String)> = pending
            .aliases
            .iter()
            .filter(|(k, ns, alias)| {
                *k == kind
                    && match kind {
                        AliasKind::Sequence => self.aliases.resolve_sequence(ns, alias),
                        AliasKind::Collection => self.aliases.resolve_collection(ns, alias),
                    }
                    .is_some_and(|t| targets.contains(&t))
            })
            .cloned()
            .collect();
        for entry in doomed {
            match kind {
                AliasKind::Sequence => self.aliases.remove_sequence(&entry.1, &entry.2),
                AliasKind::Collection => self.aliases.remove_collection(&entry.1, &entry.2),
            };
            pending.forced_aliases.remove(&entry);
            pending.aliases.remove(&entry);
        }
    }

    // =========================================================================
    // Import from another store
    // =========================================================================

    /// Import a single collection (with all its sequences, aliases, and FHR
    /// metadata) from another store into this store.
    ///
    /// Both stores must be disk-backed with matching storage modes.
    /// The source store must have the collection loaded (call
    /// `load_collection()` or `load_all_collections()` first).
    pub fn import_collection(&mut self, source: &ReadonlyRefgetStore, digest: &str) -> Result<()> {
        // Store-to-store copy does not handle packed stores.
        if self.layout == Layout::Packed || source.layout == Layout::Packed {
            return Err(anyhow!(
                "import_collection does not support packed stores: copying a \
                 collection into or out of a packed store is not implemented. \
                 Export it to FASTA and import that instead."
            ));
        }
        // Both stores must be disk-backed with same mode
        if source.local_path.is_none() || self.local_path.is_none() || !self.persist_to_disk {
            return Err(anyhow!("import_collection requires both stores to be disk-backed"));
        }
        if source.mode != self.mode {
            return Err(anyhow!(
                "import_collection requires matching storage modes (source={:?}, dest={:?})",
                source.mode,
                self.mode,
            ));
        }

        self.import_collection_file_copy(source, digest)?;

        // Copy sequence aliases for every sequence in the imported collection
        for seq_key in source.loaded_member_seq_keys(digest) {
            let seq_digest = key_to_digest_string(&seq_key);
            for (ns, alias) in source.aliases.reverse_lookup_sequence(&seq_digest) {
                self.add_sequence_alias(&ns, &alias, &seq_digest)?;
            }
        }

        // Copy collection aliases
        for (ns, alias) in source.aliases.reverse_lookup_collection(digest) {
            self.add_collection_alias(&ns, &alias, digest)?;
        }

        // Copy FHR metadata
        if let Some(fhr) = source.get_fhr_metadata(digest) {
            self.set_fhr_metadata(digest, fhr.clone())?;
        }

        Ok(())
    }

    /// File-copy based import: copies RGSI and .seq files directly from
    /// source to destination, then registers the collection in memory.
    fn import_collection_file_copy(
        &mut self,
        source: &ReadonlyRefgetStore,
        digest: &str,
    ) -> Result<()> {
        use crate::collection::read_rgsi_file;

        // 1. Read the source RGSI file to get collection metadata
        let src_rgsi_path = source
            .collection_file_path(digest)
            .ok_or_else(|| anyhow!("Source store has no local path for collection {}", digest))?;
        let collection = read_rgsi_file(&src_rgsi_path)
            .with_context(|| format!("Failed to read source RGSI file: {}", src_rgsi_path.display()))?;

        let mut metadata = collection.metadata.clone();

        // 2. Handle ancillary digests and copy/write the RGSI file
        let dst_rgsi_path = self
            .collection_file_path(digest)
            .ok_or_else(|| anyhow!("Dest store has no local path for collection {}", digest))?;
        if let Some(parent) = dst_rgsi_path.parent() {
            create_dir_all(parent)?;
        }

        let needs_ancillary_rewrite = self.ancillary_digests
            && metadata.name_length_pairs_digest.is_none();

        if needs_ancillary_rewrite {
            // Source lacks ancillary digests but destination wants them --
            // compute them and write a new RGSI file.
            metadata.compute_ancillary_digests(&collection.sequences);
            crate::collection::write_collection_rgsi_rows(
                &dst_rgsi_path,
                &metadata,
                collection.sequences.iter().map(|s| s.metadata()),
            )?;
        } else {
            // Byte-for-byte copy of the RGSI file
            fs::copy(&src_rgsi_path, &dst_rgsi_path)
                .with_context(|| format!(
                    "Failed to copy RGSI {} -> {}",
                    src_rgsi_path.display(),
                    dst_rgsi_path.display(),
                ))?;
        }

        // 3. Copy sequence data files
        for seq_record in &collection.sequences {
            let seq_meta = seq_record.metadata();
            let seq_digest = &seq_meta.sha512t24u;

            let src_seq_path = source
                .sequence_file_path(seq_digest)
                .ok_or_else(|| anyhow!("Source has no path for sequence {}", seq_digest))?;
            let dst_seq_path = self
                .sequence_file_path(seq_digest)
                .ok_or_else(|| anyhow!("Dest has no path for sequence {}", seq_digest))?;

            // Skip if destination already has this sequence (dedup across collections)
            if dst_seq_path.exists() {
                // Still need to register in memory below
            } else {
                if let Some(parent) = dst_seq_path.parent() {
                    create_dir_all(parent)?;
                }
                // Cache and store files are never modified in place, so a hard
                // link is safe and keeps a pull from doubling disk use. Fall
                // back to a copy (for example across devices).
                if fs::hard_link(&src_seq_path, &dst_seq_path).is_err() {
                    fs::copy(&src_seq_path, &dst_seq_path).with_context(|| {
                        format!(
                            "Failed to copy sequence {} -> {}",
                            src_seq_path.display(),
                            dst_seq_path.display(),
                        )
                    })?;
                }
            }
        }

        // 4. Register in memory (sequences deduped against the catalog).
        let coll_key = digest.to_key();
        let members: Vec<SequenceMetadata> =
            collection.sequences.iter().map(|s| s.metadata().clone()).collect();
        let (_, new_rows) = self.insert_collection(metadata.clone(), &members)?;
        self.record(|p| p.add_collection(coll_key));
        for row in new_rows {
            let seq_key = self.catalog.seq(row).sha512t24u().to_key();
            self.record(|p| p.add_sequence(seq_key));
        }

        // 5. Update index files
        self.write_index_files()?;

        Ok(())
    }

    // =========================================================================
    // Sequence API
    // =========================================================================

    /// List all sequences in the store (metadata only, no sequence data).
    pub fn list_sequences(&self) -> Vec<SequenceMetadata> {
        let mut result: Vec<_> = self.sequence_metadata().collect();
        result.sort_by(|a, b| a.sha512t24u.cmp(&b.sha512t24u));
        result
    }

    /// Get metadata for a single sequence by digest (no sequence data).
    /// Accepts sha512t24u (bare or `SQ.`-prefixed) or MD5.
    pub fn get_sequence_metadata<K: AsRef<[u8]>>(
        &self,
        seq_digest: K,
    ) -> Option<SequenceMetadata> {
        self.sequence_view(seq_digest).map(|v| v.to_metadata())
    }

    /// Get a sequence by its SHA512t24u (or MD5) digest. Bytes, if loaded, are
    /// shared (`Arc`), not copied.
    pub fn get_sequence<K: AsRef<[u8]>>(&self, seq_digest: K) -> Result<SequenceRecord> {
        let row = self.seq_row(seq_digest.as_ref()).ok_or_else(|| {
            anyhow!(
                "Sequence not found: {}",
                String::from_utf8_lossy(seq_digest.as_ref())
            )
        })?;
        Ok(self.seq_record(row))
    }

    /// Clear sequence data from the store to free memory.
    ///
    /// Drops every sequence's loaded bytes. Metadata is preserved: sequences,
    /// collections, name lookups, MD5 lookups, aliases, and FHR metadata
    /// remain intact.
    ///
    /// For on-disk stores, data on disk is unaffected.
    pub fn clear(&mut self) {
        self.loaded.clear();
        if let Ok(mut cache) = self.seq_fd_cache.lock() {
            cache.entries.clear();
        }
    }

    /// Check whether a sequence's bytes are loaded (owned by the store) rather
    /// than metadata only.
    pub fn is_sequence_loaded<K: AsRef<[u8]>>(&self, seq_digest: K) -> bool {
        self.seq_row(seq_digest)
            .is_some_and(|r| self.loaded.contains_key(&r))
    }

    /// Get a sequence by collection digest and name.
    ///
    /// The returned record carries THIS collection's name and description for
    /// the sequence, not the first-imported label of a sequence shared across
    /// collections. Bytes (and `fai`, if the member lacks it) come from the
    /// sequence row; the bytes are an `Arc` bump, not a copy.
    pub fn get_sequence_by_name<K: AsRef<[u8]>>(
        &self,
        collection_digest: K,
        sequence_name: &str,
    ) -> Result<SequenceRecord> {
        let m = self.member_by_name(collection_digest, sequence_name)?;
        Ok(self.member_record(m))
    }

    /// The member of a loaded collection named `sequence_name`.
    pub fn member_by_name<K: AsRef<[u8]>>(
        &self,
        collection_digest: K,
        sequence_name: &str,
    ) -> Result<MemberView<'_>> {
        let c = self.loaded_coll_row(collection_digest)?;
        self.catalog
            .member_by_name(c, sequence_name)
            .ok_or_else(|| anyhow!("Sequence '{}' not found in collection", sequence_name))
    }

    // =========================================================================
    // Loading methods
    // =========================================================================

    /// Ensure the sequence index (sequences.rgsi) is loaded.
    /// For remote stores this is deferred until first sequence access.
    /// Returns Ok(()) immediately if already loaded.
    pub(crate) fn ensure_sequence_index_loaded(&mut self) -> Result<()> {
        if self.sequence_index_loaded {
            return Ok(());
        }

        let seq_index_path = self.sequence_index_path.clone()
            .ok_or_else(|| anyhow!("No sequence_index_path set for deferred loading"))?;

        let cached = self.local_path.as_ref().is_some_and(|p| p.join(&seq_index_path).exists());
        if !self.quiet && !cached {
            eprintln!("Downloading sequence index {}...", seq_index_path);
        }
        let expected = self.manifest_digest(|m| m.sequences_digest.clone());
        let data = self
            .fetch_meta_checked(&seq_index_path, expected.as_deref())?
            .ok_or_else(|| anyhow!("Sequence index {} not found", seq_index_path))?;
        let data_str = String::from_utf8(data)
            .context("sequence index contains invalid UTF-8")?;

        Self::load_sequences_from_reader(self, data_str.as_bytes())?;
        self.sequence_index_loaded = true;
        Ok(())
    }

    /// Explicitly load the sequence index. For servers that want to preload
    /// all data during startup before converting to `ReadonlyRefgetStore`.
    pub fn load_sequence_index(&mut self) -> Result<()> {
        self.ensure_sequence_index_loaded()
    }

    /// Load the member lists of all collections. A remote store downloads
    /// the missing collection files in parallel.
    pub fn load_all_collections(&mut self) -> Result<()> {
        let digests: Vec<String> = self
            .catalog
            .coll_rows()
            .map(|c| self.catalog.coll_meta(c).digest.clone())
            .collect();
        self.ensure_collections_loaded(&digests)?;
        self.catalog.shrink_to_fit();
        Ok(())
    }

    /// Load several collections by digest. A remote store downloads the
    /// missing collection files in parallel.
    pub fn load_collections(&mut self, digests: &[String]) -> Result<()> {
        self.ensure_collections_loaded(digests)
    }

    /// Eagerly load the bytes of every sequence. A remote store downloads
    /// them in parallel and writes them to the shared cache.
    #[cfg(not(feature = "filesystem"))]
    pub fn load_all_sequences(&mut self) -> Result<()> {
        let rows: Vec<SeqRow> = self.catalog.seq_rows().collect();
        for row in rows {
            self.ensure_sequence_loaded(row)?;
        }
        Ok(())
    }

    /// Eagerly load all Stub sequences to Full. A remote store downloads them
    /// in parallel and writes them to the shared cache.
    #[cfg(feature = "filesystem")]
    pub fn load_all_sequences(&mut self) -> Result<()> {
        let keys: Vec<DigestKey> = self
            .catalog
            .seq_rows()
            .filter(|r| !self.loaded.contains_key(r))
            .map(|r| self.catalog.seq(r).sha512t24u().to_key())
            .collect();
        if keys.is_empty() {
            return Ok(());
        }
        let opts = PrefetchOptions {
            order: PrefetchOrder::Any,
            write_cache: true,
            progress: false,
        };
        // Everything stays in RAM anyway: no budget.
        let mut run = self.prefetch_keys(&keys, opts, u64::MAX)?;
        let mut got: Vec<(SeqRow, Arc<Vec<u8>>)> = Vec::with_capacity(keys.len());
        for item in run.by_ref() {
            let item = item?;
            if let Some(row) = self.catalog.find_sha(&item.digest) {
                got.push((row, item.into_shared()));
            }
        }
        drop(run);
        self.loaded.extend(got);
        Ok(())
    }

    /// Load a single collection by digest.
    pub fn load_collection(&mut self, digest: &str) -> Result<()> {
        self.ensure_collection_loaded(digest)
    }

    /// Load a single sequence by digest.
    pub fn load_sequence(&mut self, digest: &str) -> Result<()> {
        let row = self
            .seq_row(digest)
            .ok_or_else(|| anyhow!("Sequence not found in store"))?;
        self.ensure_sequence_loaded(row)
    }

    // =========================================================================
    // Resident overlay (&self): promote whole `.seq` files to RAM without &mut
    // =========================================================================

    /// Promote a set of sequences to the `&self` resident overlay, reading each
    /// whole `.seq` file into RAM as its raw on-disk (packed) bytes.
    ///
    /// Unlike [`load_sequence`](Self::load_sequence) / [`load_all_sequences`](
    /// Self::load_all_sequences) — which take `&mut self` and are therefore
    /// unusable once the store is wrapped in an `Arc` — this takes `&self` and
    /// stores into an `RwLock`-guarded overlay consulted by `get_substring` /
    /// `get_substrings` before the on-disk partial-read path.
    ///
    /// The loads are made cheap and predictable: digests are deduped (md5 aliases
    /// resolved), any already loaded or already resident are skipped, and the rest
    /// are **sorted by `.seq` path** so the reads walk the sequence directory in
    /// order rather than in `HashMap` hash order (the pathology of
    /// [`load_all_sequences`](Self::load_all_sequences)). Each file is one
    /// `File::open` + `read_to_end` (a few large sequential reads), not a scatter
    /// of positioned reads.
    ///
    /// Returns the number of bytes newly loaded into the overlay by this call
    /// (already resident or loaded skips contribute 0), for RAM-budget accounting.
    pub fn load_resident(&self, digests: &[String]) -> Result<u64> {
        // Plan: resolve to rows, skip loaded / already-resident, dedupe.
        let mut plan: Vec<SeqRow> = Vec::new();
        let mut seen: HashSet<SeqRow> = HashSet::new();
        {
            let resident = self
                .resident
                .read()
                .map_err(|_| anyhow!("resident overlay lock poisoned"))?;
            for d in digests {
                match self.seq_row(d) {
                    Some(key) if !seen.insert(key) => continue, // duplicate in this request
                    Some(key) if resident.contains_key(&key) => continue, // already resident
                    Some(key) if self.loaded.contains_key(&key) => continue, // already loaded
                    Some(key) => plan.push(key),
                    None => {
                        return Err(anyhow!(
                            "Cannot load resident: sequence not found: {}",
                            d
                        ));
                    }
                }
            }
        }
        if plan.is_empty() {
            return Ok(0);
        }

        // One prefetch: parallel reads (local files, packs, cache) or
        // downloads. Downloads go into the cache, as `load_sequence` does.
        #[cfg(feature = "filesystem")]
        let fetched = {
            let keys: Vec<DigestKey> =
                plan.iter().map(|&r| self.catalog.seq(r).sha512t24u().to_key()).collect();
            let opts = PrefetchOptions {
                order: PrefetchOrder::Any,
                write_cache: true,
                progress: false,
            };
            self.prefetch_keys(&keys, opts, u64::MAX)?.map(|item| {
                item.and_then(|item| {
                    let row = self
                        .catalog
                        .find_sha(&item.digest)
                        .ok_or_else(|| anyhow!("Sequence not found: {}", item.digest))?;
                    Ok((row, item.into_shared()))
                })
            })
        };
        #[cfg(not(feature = "filesystem"))]
        let fetched = plan.iter().map(|&row| self.read_stored_bytes(row, true).map(|b| (row, b)));
        let mut loaded: Vec<(SeqRow, Arc<Vec<u8>>)> = Vec::with_capacity(plan.len());
        let mut total: u64 = 0;
        for item in fetched {
            let (key, data) = item?;
            total += data.len() as u64;
            loaded.push((key, data));
        }

        // Publish into the overlay under a single write lock.
        {
            let mut resident = self
                .resident
                .write()
                .map_err(|_| anyhow!("resident overlay lock poisoned"))?;
            for (key, arc) in loaded {
                resident.insert(key, arc);
            }
        }

        Ok(total)
    }

    /// Drop the given digests from the resident overlay. Subsequent queries for
    /// them fall back to the on-disk partial-read path. Unknown / non-resident
    /// digests are ignored.
    pub fn evict_resident(&self, digests: &[String]) {
        if let Ok(mut resident) = self.resident.write() {
            for d in digests {
                if let Some(key) = self.seq_row(d) {
                    resident.remove(&key);
                }
            }
        }
    }

    /// Drop the entire resident overlay.
    pub fn evict_all_resident(&self) {
        if let Ok(mut resident) = self.resident.write() {
            resident.clear();
        }
    }

    /// Total bytes currently held in the resident overlay (sum of the packed
    /// `.seq` byte arrays). Used for RAM-budget accounting.
    pub fn resident_bytes(&self) -> u64 {
        self.resident
            .read()
            .map(|r| r.values().map(|b| b.len() as u64).sum())
            .unwrap_or(0)
    }

    /// Number of sequences currently resident.
    pub fn resident_len(&self) -> usize {
        self.resident.read().map(|r| r.len()).unwrap_or(0)
    }

    /// Estimate, from in-memory metadata only (no I/O), the additional bytes
    /// [`load_resident`](Self::load_resident) would pull into RAM for `digests`.
    /// Already loaded or already resident digests, and duplicates, contribute 0.
    /// Used for RAM-budget planning before committing to a load.
    pub fn resident_estimate_bytes(&self, digests: &[String]) -> u64 {
        let resident = self.resident.read().ok();
        let mut seen: HashSet<SeqRow> = HashSet::new();
        let mut total: u64 = 0;
        for d in digests {
            let Some(key) = self.seq_row(d) else { continue };
            if !seen.insert(key) {
                continue;
            }
            if let Some(r) = &resident {
                if r.contains_key(&key) {
                    continue;
                }
            }
            if self.loaded.contains_key(&key) {
                continue;
            }
            // Zstd: compressed size varies; disk_size uses uncompressed as an
            // upper bound.
            total += self.catalog.seq(key).disk_size(&self.mode) as u64;
        }
        total
    }

    /// Clone the resident packed bytes for `key`, if present. The `Arc` is cloned
    /// out under a brief read lock so the lock is released before the caller
    /// decodes -- decodes therefore never serialize against each other.
    fn resident_lookup(&self, key: SeqRow) -> Option<Arc<Vec<u8>>> {
        self.resident.read().ok().and_then(|r| r.get(&key).cloned())
    }

    /// Decode `[start, end)` from a whole stored sequence buffer (loaded or
    /// resident bytes), after the bounds check.
    fn decode_full_span(
        &self,
        metadata: &SeqCore,
        sequence: &[u8],
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>> {
        check_range(start, end, metadata.length)?;
        match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                Ok(decode_substring_from_bytes(sequence, start, end, alphabet))
            }
            StorageMode::Raw => Ok(sequence[start..end].to_vec()),
            StorageMode::Zstd => {
                let ascii = zstd::decode_all(sequence)
                    .context("zstd decompression of sequence record failed")?;
                Ok(ascii[start..end].to_vec())
            }
        }
    }

    /// Iterate over all collections with their sequences loaded.
    pub fn iter_collections(&self) -> impl Iterator<Item = crate::digest::SequenceCollection> + '_ {
        let mut digests: Vec<String> = self
            .catalog
            .coll_rows()
            .map(|c| self.catalog.coll_meta(c).digest.clone())
            .collect();
        digests.sort();

        digests.into_iter().filter_map(move |digest| {
            self.get_collection(&digest).ok()
        })
    }

    /// Iterate over all sequences with their data loaded.
    pub fn iter_sequences(&self) -> impl Iterator<Item = SequenceRecord> + '_ {
        let mut rows: Vec<SeqRow> = self.catalog.seq_rows().collect();
        rows.sort_by_key(|&r| self.catalog.seq(r).sha512t24u());
        rows.into_iter().map(move |r| self.seq_record(r))
    }

    /// Check if a collection is fully loaded.
    pub fn is_collection_loaded<K: AsRef<[u8]>>(&self, collection_digest: K) -> bool {
        self.catalog
            .find_coll(collection_digest)
            .is_some_and(|c| self.catalog.is_loaded(c))
    }

    /// Returns the local path where the store is located (if any)
    pub fn local_path(&self) -> Option<&PathBuf> {
        self.local_path.as_ref()
    }

    /// Returns the remote source URL (if any)
    pub fn remote_source(&self) -> Option<&str> {
        self.remote_source.as_deref()
    }

    /// Returns the storage mode used by this store
    pub fn storage_mode(&self) -> StorageMode {
        self.mode
    }

    /// Returns the physical layout of this store's sequence bytes, fixed when
    /// the store was created.
    pub fn layout(&self) -> Layout {
        self.layout
    }

    // =========================================================================
    // Substring retrieval
    // =========================================================================

    /// Retrieves a substring from an encoded sequence by its SHA512t24u digest.
    ///
    /// Thin `String` wrapper over [`get_substring_bytes`](Self::get_substring_bytes);
    /// see it for the actual retrieval logic.
    pub fn get_substring<K: AsRef<[u8]>>(
        &self,
        sha512_digest: K,
        start: usize,
        end: usize,
    ) -> Result<String> {
        self.get_substring_bytes(sha512_digest, start, end)
            .map(decoded_sequence_to_string)
    }

    /// Retrieves a substring from an encoded sequence by its SHA512t24u digest,
    /// as raw bytes rather than a `String`.
    pub fn get_substring_bytes<K: AsRef<[u8]>>(
        &self,
        sha512_digest: K,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>> {
        let actual_key = self.seq_row(sha512_digest.as_ref()).ok_or_else(|| {
            anyhow!(
                "Sequence not found: {}",
                String::from_utf8_lossy(sha512_digest.as_ref())
            )
        })?;

        // Zero-length range: a request for `[start, start)` is a well-formed
        // empty interval, so return an empty result instead of erroring. This
        // is the natural substring semantics, matches the wasm/JS
        // `RemoteRefgetStore.getSubstring` contract, and is what bulk/region
        // callers expect for an empty interval. Handling it here, before record
        // resolution, means every downstream path (loaded or resident bytes,
        // on-disk, remote) agrees without each needing its own guard.
        // Inverted ranges are rejected by the bounds checks below.
        if start == end {
            return Ok(Vec::new());
        }

        let core = self.catalog.seq(actual_key).core();
        // Loaded bytes, or the resident overlay (whole `.seq` promoted to RAM
        // via `load_resident`): decode the covering span from memory, with no
        // disk read.
        if let Some(bytes) = self.loaded_arc(actual_key).or_else(|| self.resident_lookup(actual_key)) {
            return self.decode_full_span(&core, &bytes, start, end);
        }
        // Not in memory: for a disk-backed store, read only the bytes covering
        // [start, end) instead of loading the whole sequence into RAM. This is
        // the random-access extract path: sparse region extraction touches a
        // tiny fraction of each sequence, so a per-query partial read avoids
        // the cost of `load_sequence` pulling entire chromosomes into memory.
        //
        // Local `.seq` (if present) -> remote byte-range (if configured).
        // Attempt the disk read directly rather than `path.exists()`-ing
        // first: the fd cache already amortizes opens across queries, so a
        // pre-check would be a pure-overhead `statx` on the overwhelmingly
        // common case (the file exists). A missing local file still falls
        // through to the remote fallback rather than erroring, detected from
        // the open failing with `NotFound`.
        if self.local_path.is_some() && self.seqdata_path_template.is_some() {
            match self.get_substring_from_disk(&core, start, end) {
                Ok(s) => return Ok(s),
                Err(e) if is_missing_seq_file(&e) => {} // fall through to remote below
                Err(e) => return Err(e),
            }
        }
        if self.remote_source.is_some() {
            return self.get_substring_from_remote(&core, start, end);
        }
        Err(anyhow!("Sequence data not loaded (stub only)"))
    }

    /// Batch substring retrieval for many ranges of ONE sequence.
    ///
    /// Resolves the record once and validates all ranges up front, then:
    /// - **Loaded or resident**: decode/slice each range from the bytes in RAM.
    /// - **Disk-backed, not loaded**: partial-read each range independently via the
    ///   fd-cached [`get_substring_from_disk`](Self::get_substring_from_disk)
    ///   path (reads only the bytes covering each range; reuses the open `.seq`
    ///   handle across ranges).
    ///
    /// The main benefit over a loop of [`get_substring`](Self::get_substring) is
    /// for FFI callers (e.g. the Python binding): resolving the record once and
    /// crossing the boundary once per sequence instead of once per range.
    ///
    /// A whole-sequence "decode once, slice all" strategy was tried for
    /// wide/overlapping range sets but measured *slower* on real chromosome-scale
    /// sequences (decoding a 250 Mbp chromosome into one buffer, plus the large
    /// allocation, outweighs the avoided re-decode of scattered spans), so the
    /// partial path is used uniformly. Callers wanting resident behavior should
    /// load the sequence (`load_sequence`) first; the in-memory arm then serves it.
    ///
    /// Output strings are byte-identical to a loop of
    /// [`get_substring`](Self::get_substring).
    ///
    /// Thin `String` wrapper over
    /// [`get_substrings_bytes`](Self::get_substrings_bytes); see it for the
    /// actual retrieval logic.
    pub fn get_substrings<K: AsRef<[u8]>>(
        &self,
        sha512_digest: K,
        ranges: &[(usize, usize)],
    ) -> Result<Vec<String>> {
        self.get_substrings_bytes(sha512_digest, ranges).map(|v| {
            v.into_iter().map(decoded_sequence_to_string).collect()
        })
    }

    /// Batch substring retrieval for many ranges of ONE sequence, as raw bytes
    /// rather than `String`s. See [`get_substrings`](Self::get_substrings) for
    /// the full behavior description.
    pub fn get_substrings_bytes<K: AsRef<[u8]>>(
        &self,
        sha512_digest: K,
        ranges: &[(usize, usize)],
    ) -> Result<Vec<Vec<u8>>> {
        let actual_key = self.seq_row(sha512_digest.as_ref()).ok_or_else(|| {
            anyhow!(
                "Sequence not found: {}",
                String::from_utf8_lossy(sha512_digest.as_ref())
            )
        })?;

        let meta = &self.catalog.seq(actual_key).core();
        let seq_length = meta.length;

        // Validate every range up front (same rule as get_substring). A
        // zero-length range (`start == end`) is a valid empty interval that
        // yields "" (see get_substring); only an inverted range (`start > end`)
        // or an out-of-bounds end is rejected. `start == seq_length` is allowed
        // only when it is also an empty range (start == end == seq_length).
        for &(start, end) in ranges {
            if start != end || end > seq_length {
                check_range(start, end, seq_length)?;
            }
        }

        match self.loaded_arc(actual_key).or_else(|| self.resident_lookup(actual_key)) {
            // Loaded bytes, or the resident overlay (whole `.seq` promoted to
            // RAM via `load_resident`): decode every range from memory.
            Some(sequence) => {
                let seq: &[u8] = sequence.as_slice();
                let mut out = Vec::with_capacity(ranges.len());
                match self.mode {
                    StorageMode::Encoded => {
                        let alphabet = lookup_alphabet(&meta.alphabet);
                        for &(start, end) in ranges {
                            out.push(decode_substring_from_bytes(seq, start, end, alphabet));
                        }
                    }
                    StorageMode::Raw => {
                        for &(start, end) in ranges {
                            out.push(seq[start..end].to_vec());
                        }
                    }
                    StorageMode::Zstd => {
                        // Decompress the whole record once, slice every range.
                        let ascii = zstd::decode_all(seq)
                            .context("zstd decompression of sequence record failed")?;
                        for &(start, end) in ranges {
                            out.push(ascii[start..end].to_vec());
                        }
                    }
                }
                Ok(out)
            }
            None => {
                // Resolve once for the whole batch: local `.seq` (if present) ->
                // remote byte-range (if configured). Whether the local file
                // exists is determined lazily below from an actual open/read
                // attempt rather than a `path.exists()` stat here -- the fd
                // cache already amortizes opens, so pre-checking existence is
                // pure overhead on the common (file present) case. A missing
                // local file still falls through to the remote fallback rather
                // than erroring.
                let has_local_path = self.local_path.is_some() && self.seqdata_path_template.is_some();
                let has_remote = self.remote_source.is_some();
                if !has_local_path && !has_remote {
                    return Err(anyhow!("Sequence data not loaded (stub only)"));
                }

                // No remote fallback configured: there is nothing to branch on,
                // so read every range straight from disk without ever stat-ing
                // for existence -- an open failure has nowhere else to fall
                // back to anyway and simply propagates. This is the common path
                // for local (non-remote) batch extraction.
                if !has_remote {
                    let mut out = Vec::with_capacity(ranges.len());
                    for &(start, end) in ranges {
                        if start == end {
                            out.push(Vec::new());
                        } else {
                            out.push(self.get_substring_from_disk(meta, start, end)?);
                        }
                    }
                    return Ok(out);
                }

                // A remote fallback exists, so we do need to know up front
                // whether the file is already cached locally in order to pick
                // between the on-disk and (possibly bulk-downloaded) remote
                // path below. This stat is paid once per batch call for this
                // sequence, not once per range.
                let mut local_seq_exists = self
                    .local_seq_path_parts(&meta.sha512t24u, meta.alphabet)
                    .map(|p| p.exists())
                    .unwrap_or(false);

                // Bulk remote extraction: one HTTP round-trip per range is slow
                // for large batches. Past REMOTE_BULK_FETCH_THRESHOLD, download the
                // whole stored record once (from the remote `.seq` or its span in
                // a remote pack; in parallel ranges, digest-checked) into the
                // shared cache as a `.seq`, then read every range from the cached
                // file. Single / small-batch reads fall through to pure byte-range.
                if !local_seq_exists
                    && self.seq_cache.is_some()
                    && ranges.len() >= REMOTE_BULK_FETCH_THRESHOLD
                {
                    self.read_stored_bytes(actual_key, true)?;
                    local_seq_exists = true;
                }

                let mut out = Vec::with_capacity(ranges.len());
                for &(start, end) in ranges {
                    // Zero-length range -> "" (see get_substring). Short-circuit
                    // before the partial-read paths, which reject start == end.
                    if start == end {
                        out.push(Vec::new());
                    } else if local_seq_exists {
                        out.push(self.get_substring_from_disk(meta, start, end)?);
                    } else {
                        out.push(self.get_substring_from_remote(meta, start, end)?);
                    }
                }
                Ok(out)
            }
        }
    }

    /// Read a substring directly from the on-disk `.seq` file, fetching only the
    /// bytes that cover `[start, end)`.
    ///
    /// `.seq` files are headerless raw byte arrays (bit-packed in Encoded mode,
    /// one ASCII byte per base in Raw mode), so byte offsets map directly to base
    /// positions. We compute the covering byte range with [`byte_range_for_bases`],
    /// read just that window with a positioned read, and decode it with
    /// [`decode_substring_from_bytes_at_offset`]. This is the partial-read path
    /// that makes random-access extraction cheap: it never materializes the whole
    /// sequence in memory.
    fn get_substring_from_disk(
        &self,
        metadata: &SeqCore,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>> {
        check_range(start, end, metadata.length)?;

        // Zstd: there is no base->byte mapping into a compressed frame, so read
        // the whole stored frame (a small CDC chunk), decompress once, and slice
        // the ASCII. Whole-frame decompress of a <=256 KiB chunk is microseconds.
        if self.mode == StorageMode::Zstd {
            let compressed = match self.pack_row(&metadata.sha512t24u, StoreSide::Local)? {
                Some(row) => self.read_pack_bytes_local(&row.pack, row.offset, row.offset + row.len)?,
                None => {
                    let file = self.get_cached_seq_file(metadata)?;
                    let file_len = file
                        .metadata()
                        .context("failed to stat zstd seq file")?
                        .len() as usize;
                    crate::posread::read_exact_window(&file, 0, file_len)?
                }
            };
            let ascii = zstd::decode_all(&compressed[..])
                .context("zstd decompression of seq file failed")?;
            return Ok(ascii[start..end].to_vec());
        }

        // Compute the covering byte window (identity under Raw). This math is
        // layout-independent; narrowing falls out of adding the pack offset to
        // an already-narrow window.
        let (byte_start, byte_end) = match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                byte_range_for_bases(start, end, alphabet.bits_per_symbol)
            }
            StorageMode::Raw => (start, end),
            StorageMode::Zstd => unreachable!("handled above"),
        };

        // Resolve where the bytes live and read only the covering window.
        let buf = match self.pack_row(&metadata.sha512t24u, StoreSide::Local)? {
            Some(row) => {
                // Add the pack offset to the already-narrow window; block-cache
                // read of the absolute range.
                let abs_start = row.offset + byte_start as u64;
                let abs_end = row.offset + byte_end as u64;
                self.read_pack_bytes_local(&row.pack, abs_start, abs_end)?
            }
            None => {
                // Fetch (or open + cache) the shared `.seq` file handle. The
                // critical section is tiny; the positioned read is outside the
                // lock, so concurrent reads are never serialized by the cache.
                let file = self.get_cached_seq_file(metadata)?;
                crate::posread::read_exact_window(&file, byte_start, byte_end - byte_start)?
            }
        };

        let decoded: Vec<u8> = match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                decode_substring_from_bytes_at_offset(&buf, byte_start, start, end, alphabet)
            }
            StorageMode::Raw => buf,
            StorageMode::Zstd => unreachable!("handled above"),
        };

        Ok(decoded)
    }

    /// Return a shared handle to the `.seq` file for `digest`, opening and
    /// caching it on a miss. Bounded LRU; the LRU handle is evicted (closed)
    /// when the cache is full.
    fn get_cached_seq_file(&self, meta: &SeqCore) -> Result<Arc<File>> {
        let full_path = self
            .local_seq_path_parts(&meta.sha512t24u, meta.alphabet)
            .ok_or_else(|| anyhow!("Partial read requires a local disk-backed store"))?;
        self.get_cached_file(&full_path)
    }

    /// Return a shared handle to the file at `full_path`, opening and caching it
    /// on a miss. The bounded LRU is keyed by path, so pack files and `.seq`
    /// files share one pool of open handles.
    pub(crate) fn get_cached_file(&self, full_path: &Path) -> Result<Arc<File>> {
        // Fast path: hit. Hold the lock only for the tiny map op.
        {
            let mut cache = self
                .seq_fd_cache
                .lock()
                .map_err(|_| anyhow!("seq fd cache mutex poisoned"))?;
            if let Some(f) = cache.get(full_path) {
                return Ok(f);
            }
        }

        // Miss: open OUTSIDE the lock (the open syscall must not serialize
        // other readers).
        let file = Arc::new(
            File::open(full_path)
                .with_context(|| format!("Failed to open file {}", full_path.display()))?,
        );

        // Insert under the lock. A concurrent miss may have inserted the same
        // path meanwhile; that is harmless (both handles read the same immutable
        // file) -- prefer the already-cached one to avoid duplicating.
        let mut cache = self
            .seq_fd_cache
            .lock()
            .map_err(|_| anyhow!("seq fd cache mutex poisoned"))?;
        if let Some(existing) = cache.get(full_path) {
            return Ok(existing);
        }
        cache.insert(full_path.to_path_buf(), Arc::clone(&file));
        Ok(file)
    }

    /// Read a substring from a remote `.seq` file via an HTTP byte-range
    /// request, fetching only the bytes that cover `[start, end)` and persisting
    /// nothing locally.
    ///
    /// This is the remote fallback for the partial-read path (flow 1): when a
    /// sequence is remote-only (no loaded or resident bytes and no cached local
    /// `.seq`), [`get_substring`](Self::get_substring) /
    /// [`get_substrings`](Self::get_substrings) fall through to here instead of
    /// erroring. It mirrors
    /// [`get_substring_from_disk`](Self::get_substring_from_disk) but fetches the
    /// covering byte window over HTTP `Range:` (the store's pooled,
    /// retrying [`Fetcher`](super::fetch::Fetcher), as [`stream_sequence`](Self::stream_sequence) uses) rather than
    /// reading from a local file. Decoded with
    /// [`decode_substring_from_bytes_at_offset`].
    ///
    /// Kept `&self` so the `Arc`-shared readonly store stays usable; this path
    /// never mutates the store and never writes to disk (use `load_sequence`,
    /// flow 3, for caching).
    fn get_substring_from_remote(
        &self,
        metadata: &SeqCore,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>> {
        check_range(start, end, metadata.length)?;

        // Zstd records are whole-frame compressed with no base->byte mapping, so
        // a covering byte window cannot be computed for a remote range request.
        // The novel/unplaced store this mode serves is always local, so this
        // path is unreachable in practice; fail loudly rather than fetch garbage.
        if self.mode == StorageMode::Zstd {
            return Err(anyhow!(
                "remote partial read is not supported for Zstd-mode stores; \
                 fetch and cache the whole compressed record first"
            ));
        }

        if self.remote_source.is_none() {
            return Err(anyhow!("Remote partial read requires a remote source"));
        }

        // Compute the covering byte window. Sequence bytes are headerless raw
        // byte arrays, so byte offsets map directly to base positions.
        let (byte_start, byte_end) = match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                byte_range_for_bases(start, end, alphabet.bits_per_symbol)
            }
            StorageMode::Raw => (start, end),
            StorageMode::Zstd => unreachable!("guarded above"),
        };

        // Resolve where the bytes live on the REMOTE, then range-fetch the
        // covering window. Packed reads go through packs/<pack>.pack; standalone
        // reads through the templated `.seq` path. The remote packed path issues
        // one direct HTTP range per read; its sidecar is fetched whole once and
        // kept in memory.
        let (relpath, abs_start, abs_end) =
            match self.pack_row(&metadata.sha512t24u, StoreSide::Remote)? {
                Some(row) => (
                    pack_relpath(&row.pack),
                    row.offset + byte_start as u64,
                    row.offset + byte_end as u64,
                ),
                None => (
                    self.sequence_relative_path(&metadata.sha512t24u)?,
                    byte_start as u64,
                    byte_end as u64,
                ),
            };
        let relpath_str = relpath.to_str().ok_or_else(|| anyhow!("Non-UTF8 sequence path"))?;
        let buf = self.remote_fetcher()?.get_range(relpath_str, abs_start, abs_end)?;

        let decoded: Vec<u8> = match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&metadata.alphabet);
                decode_substring_from_bytes_at_offset(&buf, byte_start, start, end, alphabet)
            }
            StorageMode::Raw => buf,
            StorageMode::Zstd => unreachable!("guarded above"),
        };

        // The Encoded decoder emits only ASCII alphabet bytes and Raw stores
        // ASCII sequence bytes, so this validation always passes; we keep it
        // checked (rather than skipping straight to bytes) because this is the
        // untrusted-remote path — bytes come from an HTTP server we don't
        // control — and the validation cost is nil against the network fetch.
        // Reuse the validated String buffer when returning bytes.
        String::from_utf8(decoded)
            .map(String::into_bytes)
            .map_err(|e| anyhow!("decoded remote sequence was not valid UTF-8: {}", e))
    }

    // =========================================================================
    // Streaming retrieval
    // =========================================================================

    /// The path of a sequence's stored `.seq` record relative to the store root
    /// (local directory or remote URL prefix), from the manifest's path template.
    ///
    /// The template is expanded with the digest and the result is checked for
    /// absolute paths and directory traversal, so it is safe to join onto a
    /// local root or a remote base URL.
    pub fn sequence_relative_path(&self, digest: &str) -> Result<PathBuf> {
        let template = self
            .seqdata_path_template
            .as_deref()
            .unwrap_or(DEFAULT_SEQDATA_PATH_TEMPLATE);
        let relpath = Self::expand_template(digest, template);
        let relpath_str = relpath
            .to_str()
            .ok_or_else(|| anyhow!("Non-UTF8 sequence path"))?;
        Self::sanitize_relative_path(relpath_str)?;
        Ok(relpath)
    }

    // =====================================================================
    // Packed layout: pack locations
    // =====================================================================

    /// The pack location of `digest` on `side`: the local side through the
    /// per-row pack columns, the remote side by binary search of its fetched
    /// sidecar. `Ok(None)` when that side is standalone; an error when it is
    /// packed and its sidecar lacks the digest.
    pub(crate) fn pack_row(
        &self,
        digest: &str,
        side: StoreSide,
    ) -> Result<Option<super::packidx::PackRow>> {
        let row = match (side, self.layout, self.remote_layout) {
            (StoreSide::Local, Layout::Packed, _) => match self.catalog.find_sha(digest) {
                Some(r) => self
                    .local_pack_cols()?
                    .get(r)
                    .map(|(pack, offset, len)| super::packidx::PackRow { pack, offset, len }),
                None => None,
            },
            (StoreSide::Remote, _, Layout::Packed) => {
                super::packidx::lookup(&self.remote_pack_index_bytes()?, digest)?
            }
            _ => return Ok(None),
        };
        row.map(Some).ok_or_else(|| anyhow!("sequence {} not in pack index", digest))
    }

    /// The local pack columns, building them from `sequences.pack.idx` on
    /// first use.
    fn local_pack_cols(&self) -> Result<Arc<super::catalog::PackColumns>> {
        if let Some(c) = self.pack_cols.lock().unwrap().as_ref() {
            return Ok(Arc::clone(c));
        }
        let local_path = self
            .local_path
            .as_ref()
            .ok_or_else(|| anyhow!("packed store requires a local path"))?;
        let bytes = super::packidx::load_pack_idx_or_empty(local_path)?;
        let cols = Arc::new(super::catalog::PackColumns::build(&self.catalog, &bytes)?);
        let mut g = self.pack_cols.lock().unwrap();
        if let Some(existing) = g.as_ref() {
            return Ok(Arc::clone(existing));
        }
        *g = Some(Arc::clone(&cols));
        Ok(cols)
    }

    /// The REMOTE side's sidecar bytes, fetched whole on first use.
    fn remote_pack_index_bytes(&self) -> Result<Arc<Vec<u8>>> {
        {
            let g = self.remote_pack_index.lock().unwrap();
            if let Some(h) = g.as_ref() {
                return Ok(Arc::clone(h));
            }
        }
        let handle = Arc::new(
            self.fetch_meta_remote(&["sequences.pack.idx".to_string()])?
                .pop()
                .flatten()
                .ok_or_else(|| anyhow!("Remote packed store has no sequences.pack.idx"))?,
        );
        let mut g = self.remote_pack_index.lock().unwrap();
        if let Some(existing) = g.as_ref() {
            return Ok(Arc::clone(existing));
        }
        *g = Some(Arc::clone(&handle));
        Ok(handle)
    }

    /// Invalidate the cached local pack index (call after a commit that changed
    /// the sidecar). The next resolution reloads from disk.
    pub(crate) fn invalidate_pack_index(&self) {
        *self.pack_cols.lock().unwrap() = None;
        self.pack_lengths.lock().unwrap().clear();
    }

    /// The byte length of a pack file, cached (the block cache needs it to
    /// clamp the final block of a read).
    fn pack_file_len(&self, pack: &Arc<str>) -> Result<u64> {
        if let Some(len) = self.pack_lengths.lock().unwrap().get(pack).copied() {
            return Ok(len);
        }
        let local_path = self
            .local_path
            .as_ref()
            .ok_or_else(|| anyhow!("packed read requires a local path"))?;
        let path = local_path.join("packs").join(format!("{}.pack", pack));
        let len = std::fs::metadata(&path)
            .with_context(|| format!("stat {}", path.display()))?
            .len();
        self.pack_lengths.lock().unwrap().insert(Arc::clone(pack), len);
        Ok(len)
    }

    /// Read `[abs_start, abs_end)` from a LOCAL pack file: one positional read
    /// on the cached file handle, or through the block cache when it is on.
    fn read_pack_bytes_local(&self, pack: &Arc<str>, abs_start: u64, abs_end: u64) -> Result<Vec<u8>> {
        let local_path = self
            .local_path
            .as_ref()
            .ok_or_else(|| anyhow!("packed read requires a local path"))?;
        let file = self.get_cached_file(&local_path.join(pack_relpath(pack)))?;
        let fetch = |off: u64, len: u64| -> Result<Vec<u8>> {
            crate::posread::read_exact_window(&file, off as usize, len as usize)
        };
        match &self.pack_cache {
            None => fetch(abs_start, abs_end - abs_start),
            Some(cache) => {
                let pack_len = self.pack_file_len(pack)?;
                cache.read_range(pack, pack_len, abs_start, abs_end, &fetch)
            }
        }
    }

    /// Set the per-pack size cap in bytes for new packs.
    pub fn set_pack_cap_bytes(&mut self, cap: u64) {
        self.pack_cap_bytes = cap.max(1);
    }

    /// Turn the pack block cache on with room for `mb` megabytes, or off with
    /// 0 (the default). Alpha. The cache helps clustered reads on a cold
    /// network file system and hurts scattered reads on stores larger than
    /// the cache. It fronts local pack reads only.
    pub fn set_pack_cache_mb(&mut self, mb: usize) {
        use super::packcache::{PackBlockCache, PACK_BLOCK_SIZE};
        self.pack_cache = (mb > 0)
            .then(|| PackBlockCache::new((mb << 20) / PACK_BLOCK_SIZE as usize, PACK_BLOCK_SIZE));
    }

    // =====================================================================
    // Packed layout: the session pack every packed write goes through
    // =====================================================================

    /// Append `bytes` (already in store-mode encoding) to the open session pack,
    /// sealing the current pack first if it would exceed the cap.
    pub(crate) fn session_pack_append(&self, digest: &str, bytes: &[u8]) -> Result<()> {
        use sha2::Digest as _;

        let local_path = self
            .local_path
            .as_ref()
            .ok_or_else(|| anyhow!("packed write requires a local path"))?;
        let packs_dir = local_path.join("packs");

        let mut guard = self.open_session_pack.lock().unwrap();
        if guard.is_none() {
            std::fs::create_dir_all(&packs_dir)?;
            let ordinal = 0usize;
            let token = super::atomic::unique_temp_token();
            let (tmp_path, file) = SessionPack::create_temp(&packs_dir, &token, ordinal)?;
            *guard = Some(SessionPack {
                packs_dir: packs_dir.clone(),
                token,
                ordinal,
                tmp_path,
                file,
                len: 0,
                hasher: sha2::Sha256::new(),
                rows: Vec::new(),
                sealed: Vec::new(),
                created: Vec::new(),
            });
        }

        let session = guard.as_mut().unwrap();
        let len = bytes.len() as u64;
        // Seal-before-place at cap (a sequence never straddles two packs).
        if session.len > 0 && session.len + len > self.pack_cap_bytes {
            session.seal()?;
        }
        // Write at the tracked offset, not at the file cursor: bytes a failed
        // write left behind are overwritten here or cut off at seal, so they
        // never shift this row or a later one.
        let offset = session.len;
        #[cfg(unix)]
        {
            use std::os::unix::fs::FileExt;
            session.file.write_all_at(bytes, offset)?;
        }
        #[cfg(not(unix))]
        {
            use std::io::{Seek, SeekFrom, Write};
            session.file.seek(SeekFrom::Start(offset))?;
            session.file.write_all(bytes)?;
        }
        session.hasher.update(bytes);
        session.len += len;
        session.rows.push((digest.to_string(), offset, len));
        Ok(())
    }

    /// Seal any open session pack and return its `(digest, PackRow)` rows for
    /// the sidecar merge. Called by the commit path under the lock.
    pub(crate) fn seal_open_session_pack(
        &self,
    ) -> Result<Vec<(String, super::packidx::PackRow)>> {
        let mut guard = self.open_session_pack.lock().unwrap();
        let Some(mut session) = guard.take() else {
            return Ok(Vec::new());
        };
        if let Err(e) = session.seal() {
            // Keep the session, with the rows of sub-packs it already
            // sealed, so a later commit can retry instead of losing them.
            *guard = Some(session);
            return Err(e);
        }
        // The rows go to the commit path now, so dropping the session must
        // not delete their packs.
        session.created.clear();
        Ok(std::mem::take(&mut session.sealed))
    }

    // =====================================================================
    // Packed deletion accounting + compaction
    // =====================================================================

    /// Dead (unreferenced) bytes inside pack files:
    /// `sum(pack sizes) - sum(live sidecar lengths)`. `Ok(0)` for a store with
    /// no packs. Computed, not tracked.
    pub fn packed_dead_bytes(&self) -> Result<u64> {
        let local_path = self
            .local_path
            .as_ref()
            .ok_or_else(|| anyhow!("packed_dead_bytes requires a local path"))?;
        let total_pack_bytes: u64 = Self::pack_file_sizes(&local_path.join("packs")).iter().sum();
        let mut live: u64 = 0;
        super::packidx::for_each_row(
            &super::packidx::load_pack_idx_or_empty(local_path)?,
            |_, _, _, len| {
                live += len;
                Ok(())
            },
        )?;
        Ok(total_pack_bytes.saturating_sub(live))
    }

    /// Rewrite the live spans of a packed store into fresh sealed packs,
    /// reclaiming dead bytes. Deterministic and locality-preserving (packs
    /// visited by name, rows by offset). Retrieval is byte-identical before and
    /// after. Same as [`Self::compact_with`] with default options.
    pub fn compact(&mut self) -> Result<CompactReport> {
        self.compact_with(CompactOptions::default())
    }

    /// [`Self::compact`] with options; see [`CompactOptions`].
    pub fn compact_with(&mut self, opts: CompactOptions) -> Result<CompactReport> {
        if self.layout != Layout::Packed {
            return Err(anyhow!("compact() only applies to a packed store"));
        }
        let local_path = self
            .local_path
            .clone()
            .ok_or_else(|| anyhow!("compact requires a local disk-backed store"))?;
        // A batch lock the caller already holds stays held.
        let held = self.holds_batch_lock();
        self.lock_for_batch("compact")?;
        let result = self.compact_inner(&local_path, opts);
        if !held {
            self.release_batch_lock();
        }
        result
    }

    fn compact_inner(&self, local_path: &Path, opts: CompactOptions) -> Result<CompactReport> {
        use super::packidx::{read_pack_idx_rows, write_pack_idx, PackRow};
        use std::collections::BTreeMap;

        let packs_dir = local_path.join("packs");
        let idx_path = local_path.join("sequences.pack.idx");

        let old_rows: BTreeMap<String, PackRow> =
            read_pack_idx_rows(&super::packidx::load_pack_idx_or_empty(local_path)?)?;

        let packs_before = Self::pack_file_sizes(&packs_dir).len();
        let dead_before = self.packed_dead_bytes()?;

        // Deterministic, locality-preserving order: packs by name, rows by offset.
        let mut by_pack: BTreeMap<Arc<str>, Vec<(String, u64, u64)>> = BTreeMap::new();
        for (digest, row) in &old_rows {
            by_pack
                .entry(Arc::clone(&row.pack))
                .or_default()
                .push((digest.clone(), row.offset, row.len));
        }
        for v in by_pack.values_mut() {
            v.sort_by_key(|(_, off, _)| *off);
        }

        // Stream every live span through the session-pack builder into fresh
        // sealed packs (reusing the incremental writer machinery).
        for (pack_name, rows) in &by_pack {
            let src_path = packs_dir.join(format!("{}.pack", pack_name));
            let src = File::open(&src_path)
                .with_context(|| format!("opening source pack {}", src_path.display()))?;
            for (digest, offset, len) in rows {
                let bytes =
                    crate::posread::read_exact_window(&src, *offset as usize, *len as usize)?;
                self.session_pack_append(digest, &bytes)?;
            }
        }
        let new_rows: BTreeMap<String, PackRow> =
            self.seal_open_session_pack()?.into_iter().collect();

        // Publish the replacement sidecar, then drop the old packs it no longer
        // uses. Only packs the old sidecar listed are dropped; see
        // `remove_fully_dead_packs`.
        write_pack_idx(&idx_path, &new_rows)?;
        self.invalidate_pack_index();
        let old_pack_names: HashSet<String> =
            old_rows.values().map(|r| r.pack.to_string()).collect();
        self.remove_fully_dead_packs(local_path, &old_pack_names)?;

        let new_pack_names: HashSet<String> =
            new_rows.values().map(|r| r.pack.to_string()).collect();
        let packs_pruned = if opts.prune_unreferenced {
            Self::remove_unreferenced_packs(&packs_dir, &new_pack_names)?
        } else {
            0
        };

        let dead_after = self.packed_dead_bytes()?;
        Ok(CompactReport {
            packs_before,
            packs_after: new_pack_names.len(),
            bytes_reclaimed: dead_before.saturating_sub(dead_after),
            packs_pruned,
        })
    }

    /// Unlink every sealed `*.pack` in `packs_dir` whose name is not in
    /// `referenced`, and return how many were unlinked. Temp files
    /// (`.rgstore.tmp.*`) are never touched. Called under the store lock
    /// right after the index that `referenced` describes was published.
    fn remove_unreferenced_packs(packs_dir: &Path, referenced: &HashSet<String>) -> Result<usize> {
        use std::io::ErrorKind::NotFound;
        let entries = match fs::read_dir(packs_dir) {
            Ok(entries) => entries,
            Err(e) if e.kind() == NotFound => return Ok(0),
            Err(e) => return Err(e).with_context(|| format!("listing {}", packs_dir.display())),
        };
        let mut removed = 0;
        for entry in entries {
            let entry = entry?;
            let file_name = entry.file_name();
            let Some(name) = file_name.to_str() else { continue };
            if name.starts_with(".rgstore.tmp.") {
                continue;
            }
            let Some(pack) = name.strip_suffix(".pack") else { continue };
            if referenced.contains(pack) || !entry.file_type()?.is_file() {
                continue;
            }
            match fs::remove_file(entry.path()) {
                Ok(()) => removed += 1,
                Err(e) if e.kind() == NotFound => {}
                Err(e) => {
                    return Err(e).with_context(|| {
                        format!("removing unreferenced pack {}", entry.path().display())
                    })
                }
            }
        }
        Ok(removed)
    }

    /// Names of the packs the on-disk sidecar references.
    fn sidecar_pack_names(local_path: &Path) -> Result<HashSet<String>> {
        let mut names = HashSet::new();
        super::packidx::for_each_row(
            &super::packidx::load_pack_idx_or_empty(local_path)?,
            |_, pack, _, _| {
                names.insert(pack.to_string());
                Ok(())
            },
        )?;
        Ok(names)
    }

    /// Unlink each pack in `candidates` that the current sidecar no longer
    /// references. Called under the store lock after the sidecar dropped rows.
    ///
    /// Only packs the caller saw listed in the sidecar, read earlier under the
    /// same lock, are candidates. A pack the sidecar has never listed may have
    /// been sealed by another writer that has not committed yet, so it is never
    /// deleted here, even though nothing references it.
    fn remove_fully_dead_packs(&self, local_path: &Path, candidates: &HashSet<String>) -> Result<()> {
        let live = Self::sidecar_pack_names(local_path)?;
        let packs_dir = local_path.join("packs");
        for name in candidates.difference(&live) {
            let _ = fs::remove_file(packs_dir.join(format!("{}.pack", name)));
        }
        Ok(())
    }

    /// Byte sizes of the sealed `*.pack` files (temp files excluded) in a
    /// `packs/` directory.
    fn pack_file_sizes(packs_dir: &Path) -> Vec<u64> {
        let Ok(entries) = fs::read_dir(packs_dir) else { return Vec::new() };
        entries
            .filter_map(|e| e.ok())
            .filter(|e| e.file_name().to_str().is_some_and(|n| n.ends_with(".pack")))
            .map(|e| e.metadata().map(|m| m.len()).unwrap_or(0))
            .collect()
    }

    /// Stream a (sub)sequence as decoded ASCII bases without loading the full
    /// sequence into memory.
    ///
    /// Opens a bounded byte window on the backing store (local file via
    /// `seek`+`take`, or remote via HTTP `Range: bytes=`), and wraps it in a
    /// [`StreamingDecoder`] when the store is in `Encoded` mode. In `Raw`
    /// mode the bounded byte source is returned directly, since bytes and
    /// bases are 1:1.
    ///
    /// Memory use is O(1) in the sequence length.
    pub fn stream_sequence<K: AsRef<[u8]>>(
        &self,
        sha512_digest: K,
        start: Option<u64>,
        end: Option<u64>,
    ) -> Result<Box<dyn Read + Send>> {
        // 1. Resolve digest (MD5 -> SHA512t24u fallback).
        let actual_key = self.seq_row(sha512_digest.as_ref()).ok_or_else(|| {
            anyhow!(
                "Sequence not found: {}",
                String::from_utf8_lossy(sha512_digest.as_ref())
            )
        })?;

        // 2. Fetch metadata only (do not load bytes).
        let metadata = &self.catalog.seq(actual_key).core();
        let owned = self.loaded_arc(actual_key);

        // 3. Resolve bounds.
        let length = metadata.length as u64;
        let start = start.unwrap_or(0);
        let end = end.unwrap_or(length);
        if start > end || end > length {
            return Err(anyhow!(
                "Invalid stream range: start={}, end={}, length={}",
                start,
                end,
                length
            ));
        }

        // 3b. Zstd: no base->byte mapping into a compressed frame, so the byte-
        // window streaming machinery below does not apply. Decompress the whole
        // record (small CDC chunk) once and stream the requested ASCII slice.
        if self.mode == StorageMode::Zstd {
            let (s, e) = (start as usize, end as usize);
            let stored = self.read_stored_bytes(actual_key, false)?;
            let ascii = zstd::decode_all(&stored[..])
                .context("zstd decompression of sequence record failed")?;
            return Ok(Box::new(std::io::Cursor::new(ascii[s..e].to_vec())));
        }

        // 4. Look up alphabet and compute bits-per-base.
        let alphabet = crate::digest::lookup_alphabet(&metadata.alphabet);
        let bps = match self.mode {
            StorageMode::Encoded => alphabet.bits_per_symbol as u64,
            StorageMode::Raw => 8,
            StorageMode::Zstd => unreachable!("handled above"),
        };

        // 5. Compute byte range + leading-skip bits.
        let start_bit = start * bps;
        let end_bit = end * bps;
        let byte_start = start_bit / 8;
        let byte_end = end_bit.div_ceil(8);
        let leading_skip_bits = (start_bit - byte_start * 8) as u8;
        let bases_to_emit = end - start;
        let byte_len = byte_end - byte_start;

        // Fast path for empty range.
        if bases_to_emit == 0 {
            return Ok(Box::new(std::io::empty()));
        }

        let sha = metadata.sha512t24u.to_string();

        // 6. Build the source reader. Loaded bytes are streamed
        // directly from their encoded byte buffer. Otherwise try the
        // on-disk path (local_path / remote_source).
        if let Some(sequence) = owned {
            // Share ownership of the already-resident buffer instead of
            // cloning the byte range. This keeps peak allocation O(1) in
            // the sequence size.
            let source: Box<dyn Read + Send> = Box::new(ArcSliceReader::new(
                sequence,
                byte_start as usize,
                byte_end as usize,
            ));
            return match self.mode {
                StorageMode::Encoded => Ok(Box::new(StreamingDecoder::new(
                    source,
                    alphabet,
                    leading_skip_bits,
                    bases_to_emit,
                ))),
                StorageMode::Raw => {
                    debug_assert_eq!(leading_skip_bits, 0);
                    Ok(source)
                }
                StorageMode::Zstd => unreachable!("handled above"),
            };
        }

        // Open a bounded window on a local file (seek + take), verifying it is
        // long enough. Shared by the standalone `.seq` and packed-pack paths.
        let open_local_window =
            |full: &Path, at: u64| -> Result<Box<dyn Read + Send>> {
                let mut file = File::open(full)
                    .with_context(|| format!("Failed to open local file: {}", full.display()))?;
                file.seek(SeekFrom::Start(at))
                    .context("Failed to seek in local file")?;
                let file_len = file
                    .metadata()
                    .with_context(|| format!("Failed to stat local file: {}", full.display()))?
                    .len();
                if file_len < at + byte_len {
                    return Err(anyhow!(
                        "Local file is too short to satisfy the requested byte range \
                         [{}..{}) (file has {} bytes): {}",
                        at,
                        at + byte_len,
                        file_len,
                        full.display()
                    ));
                }
                Ok(Box::new(BufReader::new(file).take(byte_len)))
            };

        // Prefer a locally-present file: the local store's own `.seq` or pack,
        // or (remote stores) the sequence's `.seq` in the shared cache, which
        // is standalone even when the remote is packed.
        let mut missing_local: Option<PathBuf> = None;
        let local_source: Option<Box<dyn Read + Send>> = match self.pack_row(&sha, StoreSide::Local)? {
            Some(row) => {
                let local = self
                    .local_path
                    .as_ref()
                    .ok_or_else(|| anyhow!("packed read requires a local path"))?;
                Some(open_local_window(&local.join(pack_relpath(&row.pack)), row.offset + byte_start)?)
            }
            None => match self.local_seq_path_parts(&sha, metadata.alphabet) {
                Some(full) if full.exists() => Some(open_local_window(&full, byte_start)?),
                other => {
                    missing_local = other;
                    None
                }
            },
        };

        let source: Box<dyn Read + Send> = if let Some(s) = local_source {
            s
        } else if self.remote_source.is_some() {
            let (relpath, s, e) = match self.pack_row(&sha, StoreSide::Remote)? {
                Some(row) => (pack_relpath(&row.pack), row.offset + byte_start, row.offset + byte_end),
                None => (self.sequence_relative_path(&sha)?, byte_start, byte_end),
            };
            let rel = relpath.to_str().ok_or_else(|| anyhow!("Non-UTF8 sequence path"))?;
            self.remote_fetcher()?.open_range(rel, s, e)?
        } else if let Some(full) = missing_local {
            return Err(anyhow!(
                "Sequence file missing locally and no remote source: {}",
                full.display()
            ));
        } else {
            return Err(anyhow!("No backing source configured for sequence {}", sha));
        };

        // 7. Dispatch on mode.
        match self.mode {
            StorageMode::Encoded => Ok(Box::new(StreamingDecoder::new(
                source,
                alphabet,
                leading_skip_bits,
                bases_to_emit,
            ))),
            StorageMode::Raw => {
                debug_assert_eq!(leading_skip_bits, 0);
                Ok(source)
            }
            StorageMode::Zstd => unreachable!("handled above"),
        }
    }

    // =========================================================================
    // Internal helpers
    // =========================================================================

    /// Expand a path template by substituting digest-based placeholders.
    pub(crate) fn expand_template(digest_str: &str, template: &str) -> PathBuf {
        debug_assert!(
            digest_str.len() >= 4,
            "Digest string must be at least 4 characters for template expansion, got {} chars",
            digest_str.len()
        );
        let path_str = template
            .replace("%s2", digest_str.get(0..2).unwrap_or(digest_str))
            .replace("%s4", digest_str.get(0..4).unwrap_or(digest_str))
            .replace("%s", digest_str);
        PathBuf::from(path_str)
    }

    /// Return the full filesystem path where a standalone `.seq` file for the
    /// given digest would be.
    ///
    /// For a remote store this is the file's place in the shared cache, which
    /// depends on the sequence's alphabet, so the digest must be known (`None`
    /// otherwise). For a local store it is the store's own template path. On a
    /// packed store no file exists there: the bytes sit in a pack, so readers
    /// should use [`stored_span`](Self::stored_span) or
    /// [`read_local_stored_bytes`](Self::read_local_stored_bytes). Returns
    /// `None` if the store has no local path or no seqdata path template.
    pub fn sequence_file_path(&self, seq_digest: &str) -> Option<PathBuf> {
        if self.seq_cache.is_some() {
            let view = self.sequence_view(seq_digest)?;
            return self.local_seq_path_parts(view.sha512t24u().as_str(), view.alphabet());
        }
        let local_path = self.local_path.as_ref()?;
        self.seqdata_path_template.as_ref()?;
        Some(local_path.join(self.sequence_relative_path(seq_digest).ok()?))
    }

    /// Where one sequence's stored bytes sit on this machine: a whole `.seq`
    /// file (standalone store, or a remote store's copy in the shared cache)
    /// or a span of a pack file (packed local store). The bytes are the
    /// store's own form (2/3/4-bit codes, ASCII, or one zstd frame); a pack
    /// adds no compression, so a span reads like a `.seq` file at an offset.
    ///
    /// The file is not checked for existence: a remote store's cache copy may
    /// not be pulled yet. A packed store needs its sequence index loaded.
    pub fn stored_span(&self, digest: &str) -> Result<StoredSpan> {
        let view = self
            .sequence_view(digest)
            .ok_or_else(|| anyhow!("Sequence not found: {}", digest))?;
        let sha = view.sha512t24u().to_string();
        let whole = |path: PathBuf| StoredSpan { path, offset: 0, len: None };
        if self.seq_cache.is_some() {
            return self
                .local_seq_path_parts(&sha, view.alphabet())
                .map(whole)
                .ok_or_else(|| anyhow!("No cache path for sequence {}", sha));
        }
        let local_path = self
            .local_path
            .as_ref()
            .ok_or_else(|| anyhow!("Sequence {} has no local file: the store is in memory", sha))?;
        match self.pack_row(&sha, StoreSide::Local)? {
            Some(row) => Ok(StoredSpan {
                path: local_path.join(pack_relpath(&row.pack)),
                offset: row.offset,
                len: Some(row.len),
            }),
            None => Ok(whole(local_path.join(self.sequence_relative_path(&sha)?))),
        }
    }

    /// The stored bytes (2/3/4-bit codes, ASCII, or one zstd frame) of one
    /// sequence, read from this machine in any layout: its pack span, its
    /// `.seq` file, or a remote store's copy in the shared cache. Never
    /// touches the network and bypasses the pack block cache; a pack's open
    /// handle comes from the store's fd cache.
    pub fn read_local_stored_bytes(&self, digest: &str) -> Result<Vec<u8>> {
        let span = self.stored_span(digest)?;
        match span.len {
            Some(len) => {
                let file = self.get_cached_file(&span.path)?;
                crate::posread::read_exact_window(&file, span.offset as usize, len as usize)
                    .with_context(|| {
                        format!("read {len} bytes at {} of {}", span.offset, span.path.display())
                    })
            }
            None => fs::read(&span.path)
                .with_context(|| format!("Failed to read seq file {}", span.path.display())),
        }
    }

    /// Return the full filesystem path to a collection RGSI file for the given digest.
    ///
    /// Returns `None` if the store has no local path.
    pub fn collection_file_path(&self, coll_digest: &str) -> Option<PathBuf> {
        let local_path = self.local_path.as_ref()?;
        Some(local_path.join(format!("collections/{}.rgsi", coll_digest)))
    }

    /// Validate a relative path to prevent directory traversal attacks.
    pub(crate) fn sanitize_relative_path(path: &str) -> Result<()> {
        if path.starts_with('/') || path.starts_with('\\') {
            return Err(anyhow!("Absolute paths not allowed: {}", path));
        }
        if path.contains("..") {
            return Err(anyhow!("Directory traversal not allowed: {}", path));
        }
        if path.contains('\0') {
            return Err(anyhow!("Null bytes not allowed in path"));
        }
        Ok(())
    }

    /// The remote store's fetcher, or an error for a local store.
    #[cfg(feature = "http")]
    pub(crate) fn remote_fetcher(&self) -> Result<&Arc<super::fetch::Fetcher>> {
        self.fetcher
            .as_ref()
            .ok_or_else(|| anyhow!("No remote source configured"))
    }

    #[cfg(not(feature = "http"))]
    fn remote_fetcher(&self) -> Result<&NoFetcher> {
        Err(anyhow!("Remote access requires the 'http' feature to be enabled"))
    }

    /// Several metadata files straight from the remote, in parallel, writing
    /// nothing. `None` for a file the remote doesn't have (or a local store).
    pub(crate) fn fetch_meta_remote(&self, rels: &[String]) -> Result<Vec<Option<Vec<u8>>>> {
        for rel in rels {
            Self::sanitize_relative_path(rel)?;
        }
        #[cfg(feature = "http")]
        if let Some(f) = &self.fetcher {
            return Ok(f
                .get_many(rels)?
                .into_iter()
                .map(|g| match g {
                    super::fetch::Got::Found(d) => Some(d),
                    super::fetch::Got::Missing { .. } => None,
                })
                .collect());
        }
        Ok(vec![None; rels.len()])
    }

    /// Several metadata files: cached ones from disk (unless `refresh`), the
    /// rest from the remote in one parallel batch, which are then cached.
    pub(crate) fn fetch_meta_many(&self, rels: &[String], refresh: bool) -> Result<Vec<Option<Vec<u8>>>> {
        let mut out: Vec<Option<Vec<u8>>> = Vec::with_capacity(rels.len());
        let mut missing: Vec<usize> = Vec::new();
        for (i, rel) in rels.iter().enumerate() {
            Self::sanitize_relative_path(rel)?;
            let cached = if refresh { None } else { self.read_cached_meta(rel)? };
            if cached.is_none() {
                missing.push(i);
            }
            out.push(cached);
        }
        if missing.is_empty() {
            return Ok(out);
        }
        let wanted: Vec<String> = missing.iter().map(|&i| rels[i].clone()).collect();
        for (i, data) in missing.into_iter().zip(self.fetch_meta_remote(&wanted)?) {
            if let Some(data) = data {
                self.write_cached_meta(&rels[i], &data)?;
                out[i] = Some(data);
            }
        }
        Ok(out)
    }

    /// A metadata file that must match a sha256 from the manifest. A cached
    /// copy was checked when it was downloaded. A download that doesn't match
    /// is fetched once more (the store may have been republished between the
    /// two GETs); a second mismatch is an error.
    pub(crate) fn fetch_meta_checked(&self, rel: &str, expected_sha256: Option<&str>) -> Result<Option<Vec<u8>>> {
        Self::sanitize_relative_path(rel)?;
        if let Some(data) = self.read_cached_meta(rel)? {
            return Ok(Some(data));
        }
        for attempt in 0..2 {
            let Some(data) = self.fetch_meta_remote(&[rel.to_string()])?.pop().flatten() else {
                return Ok(None);
            };
            let matches = expected_sha256.is_none_or(|want| {
                use sha2::{Digest, Sha256};
                format!("{:x}", Sha256::digest(&data)) == want
            });
            if matches {
                self.write_cached_meta(rel, &data)?;
                return Ok(Some(data));
            }
            if attempt == 1 {
                return Err(anyhow!(
                    "{} does not match rgstore.json; the store may be in the middle of an update",
                    rel
                ));
            }
        }
        unreachable!("the loop returns on its second pass")
    }

    /// Pull sidecar files (alias TSVs, FHR JSON) under `strategy`, in one
    /// parallel batch. Only a file the remote doesn't have counts as
    /// `not_found`; network errors propagate. Returns the bytes of each file
    /// that is now local (`KeepOurs`, `KeepTheirs`).
    pub(crate) fn sync_sidecars(
        &self,
        rels: &[String],
        strategy: &SyncStrategy,
    ) -> Result<(PullResult, Vec<Option<Vec<u8>>>)> {
        let mut result = PullResult::default();
        let local_exists: Vec<bool> = rels
            .iter()
            .map(|r| self.local_path.as_ref().is_some_and(|p| p.join(r).exists()))
            .collect();
        match strategy {
            SyncStrategy::KeepOurs | SyncStrategy::KeepTheirs => {
                let refresh = matches!(strategy, SyncStrategy::KeepTheirs);
                let datas = self.fetch_meta_many(rels, refresh)?;
                for (data, was_local) in datas.iter().zip(&local_exists) {
                    match data {
                        Some(_) if *was_local && !refresh => result.skipped += 1,
                        Some(_) => result.pulled += 1,
                        None => result.not_found += 1,
                    }
                }
                Ok((result, datas))
            }
            SyncStrategy::Notify => {
                let remote = self.fetch_meta_remote(rels)?;
                for ((rel, data), exists) in rels.iter().zip(remote).zip(local_exists) {
                    match data {
                        Some(remote_data) if exists => {
                            let local = self.local_path.as_ref().expect("exists implies a local path").join(rel);
                            if fs::read(&local)? != remote_data {
                                result.conflicts.push(rel.clone());
                            } else {
                                result.skipped += 1;
                            }
                        }
                        Some(_) => result.conflicts.push(rel.clone()),
                        None => result.not_found += 1,
                    }
                }
                Ok((result, vec![None; rels.len()]))
            }
        }
    }

    fn read_cached_meta(&self, rel: &str) -> Result<Option<Vec<u8>>> {
        let Some(lp) = &self.local_path else { return Ok(None) };
        let full = lp.join(rel);
        match fs::read(&full) {
            Ok(d) => Ok(Some(d)),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
            Err(e) => Err(anyhow!("Failed to read local file {}: {}", full.display(), e)),
        }
    }

    fn write_cached_meta(&self, rel: &str, data: &[u8]) -> Result<()> {
        if let Some(lp) = &self.local_path {
            let full = lp.join(rel);
            super::atomic::atomic_write_bytes(&full, data)
                .with_context(|| format!("Failed to cache file to: {}", full.display()))?;
        }
        Ok(())
    }

    /// A manifest field of this store's cached `rgstore.json`.
    pub(crate) fn manifest_digest<F: Fn(&StoreMetadata) -> Option<String>>(&self, f: F) -> Option<String> {
        let lp = self.local_path.as_ref()?;
        super::persistence::read_manifest(lp).ok().flatten().and_then(|m| f(&m))
    }

    /// The stored bytes (packed, raw or zstd, as on disk) of one sequence: its
    /// loaded or resident bytes if present, else the local pack or the local
    /// `.seq` file. A remote store reads its `.seq` in the shared cache, or
    /// downloads the remote `.seq` or pack span (in parallel ranges,
    /// digest-checked) through one prefetch. With `cache`, downloaded bytes
    /// are also published to the shared cache as a standalone `.seq`. Never
    /// adds the sequence to the loaded set.
    pub(crate) fn read_stored_bytes(&self, row: SeqRow, cache: bool) -> Result<Arc<Vec<u8>>> {
        if let Some(bytes) = self.loaded_arc(row).or_else(|| self.resident_lookup(row)) {
            return Ok(bytes);
        }
        let core = self.catalog.seq(row).core();
        let sha = core.sha512t24u.to_string();

        #[cfg(feature = "filesystem")]
        if self.remote_source.is_some() {
            let opts = PrefetchOptions { write_cache: cache, ..PrefetchOptions::default() };
            return Ok(self
                .prefetch_keys(&[sha.to_key()], opts, DEFAULT_MEMORY_BUDGET)?
                .next()
                .ok_or_else(|| anyhow!("Sequence data for {} was not delivered", sha))??
                .into_shared());
        }
        let _ = cache;

        if let Some(p) = self.pack_row(&sha, StoreSide::Local)? {
            return Ok(Arc::new(self.read_pack_bytes_local(&p.pack, p.offset, p.offset + p.len)?));
        }
        let path = self
            .local_seq_path_parts(&sha, core.alphabet)
            .ok_or_else(|| anyhow!("No sequence data path template configured"))?;
        fs::read(&path).map(Arc::new).map_err(|e| {
            if e.kind() == std::io::ErrorKind::NotFound && self.remote_source.is_none() {
                anyhow!(
                    "Sequence data for {} not found locally and no remote source configured",
                    sha
                )
            } else {
                anyhow!("Failed to read seq file {}: {}", path.display(), e)
            }
        })
    }

    /// Ensure a collection's member list is loaded into the catalog.
    pub(crate) fn ensure_collection_loaded(&mut self, collection_digest: &str) -> Result<()> {
        self.ensure_collections_loaded(std::slice::from_ref(&collection_digest))
    }

    /// Ensure several collections are loaded. Missing `collections/<d>.rgsi`
    /// files come down in one parallel batch (remote stores), then each is
    /// read into the catalog.
    pub(crate) fn ensure_collections_loaded<S: AsRef<str>>(&mut self, digests: &[S]) -> Result<()> {
        // (digest as given, digest as stored)
        let mut need: Vec<(String, String)> = Vec::new();
        let mut seen: HashSet<String> = HashSet::new();
        for d in digests {
            let d = d.as_ref();
            let existing = self.catalog.find_coll(d);
            if existing.is_some_and(|c| self.catalog.is_loaded(c)) {
                continue;
            }
            let digest_str = match existing {
                Some(c) => self.catalog.coll_meta(c).digest.clone(),
                None => d.to_string(),
            };
            if seen.insert(digest_str.clone()) {
                need.push((d.to_string(), digest_str));
            }
        }
        if need.is_empty() {
            return Ok(());
        }

        let rels: Vec<String> = need.iter().map(|(_, d)| format!("collections/{}.rgsi", d)).collect();
        if !self.quiet {
            let to_download = rels
                .iter()
                .filter(|r| !self.local_path.as_ref().is_some_and(|p| p.join(r).exists()))
                .count();
            if to_download > 0 && self.remote_source.is_some() {
                if to_download == 1 {
                    eprintln!("Downloading collection metadata {}...", need[0].1);
                } else {
                    eprintln!("Downloading {} collection metadata files...", to_download);
                }
            }
        }
        // Downloads are written to the cache, so every file is now on disk.
        let datas = self.fetch_meta_many(&rels, false)?;

        let local_path = self
            .local_path
            .clone()
            .ok_or_else(|| anyhow!("No local path configured"))?;
        for ((given, _), (rel, data)) in need.iter().zip(rels.iter().zip(datas)) {
            if data.is_none() {
                return Err(anyhow!(
                    "File not found locally and no remote source configured: {}",
                    rel
                ));
            }
            self.load_collection_rgsi(&local_path.join(rel), Some(given))?;
        }
        Ok(())
    }

    /// Ensure a sequence is loaded into memory. A remote store downloads it
    /// (in parallel ranges, digest-checked) into the shared cache.
    pub(crate) fn ensure_sequence_loaded(&mut self, row: SeqRow) -> Result<()> {
        if self.loaded.contains_key(&row) {
            return Ok(());
        }
        if !self.quiet && self.remote_source.is_some() {
            let core = self.catalog.seq(row).core();
            let cached = self
                .local_seq_path_parts(core.sha512t24u.as_str(), core.alphabet)
                .is_some_and(|p| p.exists());
            if !cached {
                eprintln!("Downloading sequence {}...", core.sha512t24u);
            }
        }
        let data = self.read_stored_bytes(row, self.persist_to_disk)?;
        self.loaded.insert(row, data);
        Ok(())
    }

    // =========================================================================
    // Write methods
    // =========================================================================

    /// Write the store using its configured paths.
    pub fn write(&self) -> Result<()> {
        if !self.persist_to_disk {
            return Err(anyhow!(
                "write() only works with disk-backed stores - use write_store_to_dir() instead"
            ));
        }
        self.write_index_files()
    }

    /// Write a RefgetStore object to a directory
    pub fn write_store_to_dir<P: AsRef<Path>>(
        &self,
        root_path: P,
        seqdata_path_template: Option<&str>,
    ) -> Result<()> {
        let root_path = root_path.as_ref();

        let template = seqdata_path_template
            .or(self.seqdata_path_template.as_deref())
            .unwrap_or(DEFAULT_SEQDATA_PATH_TEMPLATE);

        if !self.quiet {
            eprintln!(
                "Writing store to directory: {}; Using seqdata path template: {}",
                root_path.display(),
                template
            );
        }

        fs::create_dir_all(root_path)?;

        let sequences_dir = root_path.join("sequences");
        fs::create_dir_all(&sequences_dir)?;

        let collections_dir = root_path.join("collections");
        fs::create_dir_all(&collections_dir)?;

        for row in self.catalog.seq_rows() {
            let record = self.seq_record(row);
            match &record {
                SequenceRecord::Full { metadata, .. } => {
                    let rel_path = Self::expand_template(&metadata.sha512t24u, template);
                    let full_path = root_path.join(&rel_path);
                    record.to_file(full_path)?;
                }
                SequenceRecord::Stub(metadata) => {
                    // A PACKED source has no resident bytes; export each sequence
                    // as a STANDALONE `.seq` by reading its whole span from the
                    // pack. This is also the practical "convert a packed store"
                    // escape hatch.
                    if self.layout == Layout::Packed && self.local_path.is_some() {
                        let bytes = self.read_stored_bytes(row, false)?;
                        let rel_path = Self::expand_template(&metadata.sha512t24u, template);
                        let full_path = root_path.join(&rel_path);
                        if let Some(parent) = full_path.parent() {
                            fs::create_dir_all(parent)?;
                        }
                        super::atomic::publish_bytes_nosync(&full_path, &bytes)?;
                    }
                }
            }
        }

        for c in self.catalog.coll_rows() {
            let meta = self.catalog.coll_meta(c);
            let collection_file_path = root_path.join(format!("collections/{}.rgsi", meta.digest));
            crate::collection::write_collection_rgsi_rows(
                &collection_file_path,
                meta,
                self.catalog.members(c).map(|m| m.to_metadata()),
            )?;
        }

        let sequence_index_path = root_path.join("sequences.rgsi");
        self.write_sequences_rgsi(&sequence_index_path)?;

        let collection_index_path = root_path.join("collections.rgci");
        self.write_collections_rgci(&collection_index_path)?;

        let aliases_dir = root_path.join("aliases");
        self.aliases.write_to_dir(&aliases_dir)?;

        super::fhr_metadata::write_sidecars(&root_path.join("fhr"), &self.fhr_metadata)?;

        self.write_rgstore_json(root_path, template)?;

        Ok(())
    }

    /// Returns statistics about the store.
    ///
    /// The `*_in_memory` fields are a RAM-residency snapshot, not a record of
    /// what any import run did; see [`StoreStats`].
    pub fn stats(&self) -> StoreStats {
        let n_sequences = self.catalog.n_sequences();
        let n_sequences_in_memory = self.loaded.len();
        let n_collections = self.catalog.n_collections();
        let n_collections_in_memory = self
            .catalog
            .coll_rows()
            .filter(|&c| self.catalog.is_loaded(c))
            .count();
        let mode_str = match self.mode {
            StorageMode::Raw => "Raw",
            StorageMode::Encoded => "Encoded",
            StorageMode::Zstd => "Zstd",
        };
        let packed_dead_bytes = if self.layout == Layout::Packed {
            self.packed_dead_bytes().ok()
        } else {
            None
        };
        StoreStats {
            n_sequences,
            n_sequences_in_memory,
            n_collections,
            n_collections_in_memory,
            storage_mode: mode_str.to_string(),
            logical_sequence_bytes: self.logical_sequence_bytes() as u64,
            packed_dead_bytes,
        }
    }

    /// List alias namespaces available on this store (from manifest).
    pub fn available_alias_namespaces(&self) -> AvailableAliases<'_> {
        AvailableAliases {
            sequences: &self.available_sequence_alias_namespaces,
            collections: &self.available_collection_alias_namespaces,
        }
    }
}

impl Display for ReadonlyRefgetStore {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let total_size = self.logical_sequence_bytes();
        let size_str = format_bytes(total_size);
        writeln!(f, "ReadonlyRefgetStore object:")?;
        writeln!(f, "  Mode: {:?}", self.mode)?;
        writeln!(f, "  Disk size: {} ({} bytes)", size_str, total_size)?;
        writeln!(f, ">Sequences (n={}):", self.catalog.n_sequences())?;
        for (i, view) in self.sequence_views().take(10).enumerate() {
            let first_8_chars = match self.loaded_arc(view.row()) {
                None => "<stub>".to_string(),
                Some(seq) => match self.mode {
                    StorageMode::Encoded => {
                        let alphabet = lookup_alphabet(&view.alphabet());
                        let decoded =
                            decode_substring_from_bytes(&seq, 0, 8.min(view.length()), alphabet);
                        String::from_utf8(decoded).unwrap_or_else(|_| "???".to_string())
                    }
                    StorageMode::Raw => String::from_utf8(seq[0..8.min(seq.len())].to_vec())
                        .unwrap_or_else(|_| "???".to_string()),
                    StorageMode::Zstd => zstd::decode_all(&seq[..])
                        .ok()
                        .and_then(|ascii| {
                            String::from_utf8(ascii[0..8.min(ascii.len())].to_vec()).ok()
                        })
                        .unwrap_or_else(|| "???".to_string()),
                },
            };

            writeln!(
                f,
                "   - {}. {:02x?}, MD5: {:02x?}, Length: {}, Alphabet: {:?}, Start: {}",
                i + 1,
                view.sha512t24u().as_str(),
                view.md5().as_str(),
                view.length(),
                view.alphabet(),
                first_8_chars
            )?;
        }
        let loaded: Vec<CollRow> =
            self.catalog.coll_rows().filter(|&c| self.catalog.is_loaded(c)).collect();
        writeln!(f, ">Collections (n={:?}):", loaded.len())?;
        for (i, &c) in loaded.iter().enumerate() {
            let n = self.catalog.n_members(c);
            writeln!(
                f,
                "  {}. Collection Digest: {:02x?} ({} sequences)",
                i + 1,
                self.catalog.coll_meta(c).digest,
                n
            )?;
            for m in self.catalog.members(c).take(5) {
                writeln!(f, "   - Name: {}, SHA512: {:02x?}", m.name(), m.seq().sha512t24u().as_str())?;
            }
            if n > 5 {
                writeln!(f, "   - ... and {} more", n - 5)?;
            }
        }

        Ok(())
    }
}

/// A packed store served over HTTP reads back through every remote read path:
/// one range and a small batch (byte ranges into the remote pack), a bulk batch
/// (whole record fetched once and cached as a `.seq`), the bulk batch again
/// from that cache, `load_sequence`, and FASTA export.
#[cfg(all(test, feature = "http"))]
mod remote_packed_tests {
    use super::*;
    use crate::digest::digest_sequence;
    use crate::store::http_fixture::Fixture;
    use crate::store::{ExportOptions, ExportSelection, FetchConfig, LabelTemplate, RefgetStore};

    #[test]
    fn remote_packed_store_reads() {
        let short = "GGGGCCCCAAAATTTT";
        let long: String = (0..200usize).map(|i| b"ACGT"[(i * 7 + i / 5) % 4] as char).collect();
        let with_n = "NNNNACGTACGTTTGGCCAANNNNACGT";
        // chrLong is second, so its record starts at a non-zero pack offset.
        let fasta = format!(">chrA\n{short}\n>chrLong\n{long}\n>chrN\n{with_n}\n");
        let sha = |name: &str, bases: &str| {
            digest_sequence(name, bases.as_bytes()).metadata().sha512t24u.clone()
        };
        let (long_sha, n_sha) = (sha("chrLong", &long), sha("chrN", with_n));
        let wanted = |ranges: &[(usize, usize)]| -> Vec<String> {
            ranges.iter().map(|&(s, e)| long[s..e].to_string()).collect()
        };
        let few: Vec<(usize, usize)> = vec![(0, 9), (10, 11), (37, 80), (150, 199), (5, 5)];
        let bulk: Vec<(usize, usize)> = (0..20).map(|i| (i * 9 + 1, i * 9 + 8)).collect();
        assert!(bulk.len() >= REMOTE_BULK_FETCH_THRESHOLD);

        for mode in [StorageMode::Encoded, StorageMode::Raw, StorageMode::Zstd] {
            let work = tempfile::tempdir().unwrap();
            let fa = work.path().join("in.fa");
            std::fs::write(&fa, &fasta).unwrap();
            let served = work.path().join("served");
            let coll_digest = {
                let mut s = RefgetStore::on_disk_packed(&served).unwrap();
                s.set_quiet(true);
                s.set_encoding_mode(mode).unwrap();
                let (meta, _) = s
                    .add_sequence_collection_from_fasta(&fa, FastaImportOptions::new().jobs(1))
                    .unwrap();
                meta.digest
            };
            let pack_rows = super::super::packidx::read_pack_idx_rows(
                &std::fs::read(served.join("sequences.pack.idx")).unwrap(),
            )
            .unwrap();
            let long_row = &pack_rows[&long_sha];
            assert!(long_row.offset > 0, "{mode:?}: chrLong should not start its pack");
            let long_stored = std::fs::read(served.join(pack_relpath(&long_row.pack))).unwrap()
                [long_row.offset as usize..(long_row.offset + long_row.len) as usize]
                .to_vec();

            let fx = Fixture::serve(&served, &["store"]);
            let mut store =
                RefgetStore::open_remote(work.path().join("cache"), fx.url("store"), &FetchConfig::default())
                    .unwrap();
            store.set_quiet(true);
            store.load_sequence_index().unwrap();
            assert_eq!(store.remote_layout, Layout::Packed);

            // Zstd stores refuse remote partial reads by design.
            if mode != StorageMode::Zstd {
                let got = store.get_substring(&long_sha, 3, 50).unwrap();
                assert_eq!(got, long[3..50], "{mode:?}: one range");
                let got = store.get_substrings(&long_sha, &few).unwrap();
                assert_eq!(got, wanted(&few), "{mode:?}: small batch");
            }

            let cached = store.sequence_file_path(&long_sha).unwrap();
            assert!(!cached.exists(), "{mode:?}: nothing cached before the bulk batch");
            let got = store
                .get_substrings(&long_sha, &bulk)
                .unwrap_or_else(|e| panic!("{mode:?}: bulk batch failed: {e:#}"));
            assert_eq!(got, wanted(&bulk), "{mode:?}: bulk batch");
            assert_eq!(
                std::fs::read(&cached).unwrap(),
                long_stored,
                "{mode:?}: the cached .seq must be the record's bytes from the pack"
            );
            let got = store.get_substrings(&long_sha, &bulk).unwrap();
            assert_eq!(got, wanted(&bulk), "{mode:?}: bulk batch from the cache");

            store.load_sequence(&n_sha).unwrap();
            assert!(store.get_sequence(&n_sha).unwrap().is_loaded());
            let got = store.get_substring(&n_sha, 0, with_n.len()).unwrap();
            assert_eq!(got, with_n, "{mode:?}: loaded sequence");

            let out = work.path().join("out.fa");
            let opts = ExportOptions {
                label: LabelTemplate::parse("{name}").unwrap(),
                line_width: 0,
                ..Default::default()
            };
            let sel = ExportSelection::Collections { digests: vec![coll_digest], names: None };
            store.export_fasta(&sel, &out, &opts).unwrap();
            drop(fx);
            assert_eq!(std::fs::read_to_string(&out).unwrap(), fasta, "{mode:?}: export");
        }
    }
}
