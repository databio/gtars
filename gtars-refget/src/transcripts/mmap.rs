//! Native-only file-backed transcript store backends (gated behind `filesystem`).
//!
//! Provides two file-backed byte sources for a `.reftx` store:
//!
//! - **mmap** ([`TxBytes::Mmap`]): the whole file is memory-mapped. Best for
//!   the large, multi-process-shared (VRS fan-out) case where the OS page cache
//!   is shared. Sound only because the builder publishes files write-once via an
//!   atomic rename (see [`super::builder`]).
//! - **pread** ([`TxBytes::Pread`]): positioned reads on a shared `&File`, no
//!   `unsafe` mapping. The conservative default for the small / single-process
//!   case (a handful of scattered small reads per lookup).
//!
//! Both readonly backends reuse the WASM-safe lookup/parse logic in
//! [`super::store`] over the [`TxBytes`] abstraction — only the byte source
//! differs. The positioned read helper is COPIED/ADAPTED from
//! `gtars-refget/src/store/readonly.rs` (`read_exact_window`); it deliberately
//! does NOT call into the refget store.

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use anyhow::Result;
use memmap2::Mmap;

use crate::transcripts::models::Transcript;
use crate::transcripts::store::{
    fnv1a_64, linear_probe, parse_header, read_record, ReadonlyTxStore, TxBytes, INDEX_ENTRY_SIZE,
};

/// Which file-backed byte source to use when opening a `.reftx` file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TxBackend {
    /// Memory-map the whole file. For the large, multi-process-shared case.
    Mmap,
    /// Positioned (`pread`) reads. The conservative single-process default.
    Pread,
}

/// Positioned-read byte source over a single `.reftx` file.
///
/// Holds a shared `Arc<File>` plus the file length captured at open. Because
/// reftx opens ONE file (not many `.seq` files keyed by digest), a full
/// `FdCache` is unnecessary — the single `Arc<File>` IS the cache. `Arc<File>`
/// is auto `Send + Sync`, so no `unsafe` is required.
pub struct PreadSource {
    file: Arc<File>,
    len: usize,
}

impl PreadSource {
    /// Total file length (captured at open).
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    /// Read `len` bytes at `offset` via a positioned read; `None` if the
    /// requested window overflows the captured file length.
    pub(crate) fn read_at(&self, offset: usize, len: usize) -> Option<Vec<u8>> {
        let end = offset.checked_add(len)?;
        if end > self.len {
            return None;
        }
        read_exact_window(&self.file, offset, len).ok()
    }
}

/// Read exactly `len` bytes starting at `byte_start` from `file`.
///
/// COPIED/ADAPTED from `gtars-refget/src/store/readonly.rs::read_exact_window`
/// (do NOT call into that module). On unix this is a POSITIONED read (`pread`
/// via `FileExt::read_exact_at`), cursorless and safe to call concurrently on a
/// shared `&File`. On non-unix we `try_clone()` for an independent cursor, then
/// seek+read.
fn read_exact_window(file: &File, byte_start: usize, len: usize) -> Result<Vec<u8>> {
    let mut buf = vec![0u8; len];
    #[cfg(unix)]
    {
        use std::os::unix::fs::FileExt;
        file.read_exact_at(&mut buf, byte_start as u64)?;
    }
    #[cfg(not(unix))]
    {
        use std::io::{Read, Seek, SeekFrom};
        let mut handle = file.try_clone()?;
        handle.seek(SeekFrom::Start(byte_start as u64))?;
        handle.read_exact(&mut buf)?;
    }
    Ok(buf)
}

// ============================================================================
// TxStore: Mutable store for setup/build phase (mmap-backed)
// ============================================================================

/// Mutable transcript store for loading and setup, backed by a memory map.
///
/// Use this during initialization, then convert to `ReadonlyTxStore` for
/// concurrent access:
///
/// ```ignore
/// let store = TxStore::open("transcripts.reftx")?;   // or open_pread / open_mmap
/// let readonly = store.into_readonly();
/// let shared = Arc::new(readonly);
/// // Now safe for concurrent &self access
/// ```
///
/// Lookups run through a transient [`TxBytes::Mmap`] view rebuilt for each call.
/// Building a `TxBytes::Mmap` does not copy the file — it re-maps the same file
/// handle — but the build-phase store is primarily a staging area before
/// [`into_readonly`](TxStore::into_readonly); hot lookups go through the
/// readonly store.
pub struct TxStore {
    /// Wrapped in a `TxBytes::Mmap` so the backend-agnostic lookup core in
    /// `super::store` can be reused directly. The whole-file map is never copied.
    bytes: TxBytes,
    record_count: u64,
    index_offset: u64,
    mane_index_offset: u64,
    /// Pre-decoded transcript cache (populated during into_readonly).
    decoded_cache: HashMap<u64, Transcript>,
}

impl TxStore {
    /// Open a transcript store from disk via memory map.
    ///
    /// Equivalent to [`TxStore::open_mmap`]; kept as the historical name for
    /// the mutable build-phase store.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::open_mmap(path)
    }

    /// Open the mutable build-phase store via memory map.
    pub fn open_mmap<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())?;
        // SAFETY: the file is published write-once by the atomic-rename builder,
        // so no writer mutates the mapping under us.
        let mmap = unsafe { Mmap::map(&file)? };
        let bytes = TxBytes::Mmap(mmap);
        let header = parse_header(&bytes)?;
        Ok(Self {
            bytes,
            record_count: header.record_count,
            index_offset: header.index_offset,
            mane_index_offset: header.mane_index_offset,
            decoded_cache: HashMap::new(),
        })
    }

    /// Number of transcripts in store.
    #[inline]
    pub fn len(&self) -> u64 {
        self.record_count
    }

    /// Returns true if store is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.record_count == 0
    }

    /// Look up a transcript by accession (O(log n) binary search).
    pub fn lookup(&self, accession: &str) -> Option<Transcript> {
        let hash = fnv1a_64(accession.as_bytes());
        self.lookup_by_hash(accession, hash)
    }

    fn lookup_by_hash(&self, accession: &str, hash: u64) -> Option<Transcript> {
        let bytes = &self.bytes;
        let idx_start = self.index_offset as usize;
        let record_bound = self.index_offset as usize;
        let count = self.record_count as usize;

        let mut lo = 0usize;
        let mut hi = count;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let entry_offset = idx_start + mid * INDEX_ENTRY_SIZE;
            let entry_hash = bytes.read_u64_le(entry_offset)?;
            match entry_hash.cmp(&hash) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    let record_offset = bytes.read_u64_le(entry_offset + 8)? as usize;
                    let tx = read_record(bytes, record_offset, record_bound)?;
                    if tx.accession == accession {
                        return Some(tx);
                    }
                    return linear_probe(
                        bytes,
                        self.index_offset,
                        record_bound,
                        accession,
                        hash,
                        mid,
                        count,
                    );
                }
            }
        }
        None
    }

    /// Pre-decode specific transcripts into cache.
    pub fn ensure_decoded(&mut self, accession: &str) -> Result<()> {
        let hash = fnv1a_64(accession.as_bytes());
        if self.decoded_cache.contains_key(&hash) {
            return Ok(());
        }
        if let Some(tx) = self.lookup_by_hash(accession, hash) {
            self.decoded_cache.insert(hash, tx);
        }
        Ok(())
    }

    /// Pre-decode all transcripts matching a predicate.
    pub fn ensure_decoded_where<F>(&mut self, predicate: F) -> Result<usize>
    where
        F: Fn(&Transcript) -> bool,
    {
        let idx_start = self.index_offset as usize;
        let record_bound = self.index_offset as usize;
        let mut decoded: Vec<(u64, Transcript)> = Vec::new();

        for i in 0..self.record_count as usize {
            let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
            let hash = self
                .bytes
                .read_u64_le(entry_offset)
                .ok_or_else(|| anyhow::anyhow!("index entry out of bounds"))?;
            let record_offset = self
                .bytes
                .read_u64_le(entry_offset + 8)
                .ok_or_else(|| anyhow::anyhow!("index entry out of bounds"))?
                as usize;

            if let Some(tx) = read_record(&self.bytes, record_offset, record_bound)
                && predicate(&tx)
            {
                decoded.push((hash, tx));
            }
        }
        let count = decoded.len();
        for (hash, tx) in decoded {
            self.decoded_cache.insert(hash, tx);
        }
        Ok(count)
    }

    /// Convert to immutable store for concurrent access (mmap-backed).
    pub fn into_readonly(mut self) -> ReadonlyTxStore {
        if self.record_count < 500_000 && self.decoded_cache.is_empty() {
            let _ = self.ensure_decoded_where(|_| true);
        }
        self.build_readonly()
    }

    /// Convert to immutable store with lazy decoding (mmap-backed).
    pub fn into_readonly_lazy(self) -> ReadonlyTxStore {
        self.build_readonly()
    }

    /// Build the readonly store from the mmap WITHOUT copying the bytes
    /// (zero-copy): the `TxBytes::Mmap` moves directly into the readonly store.
    fn build_readonly(self) -> ReadonlyTxStore {
        ReadonlyTxStore::from_parts(
            self.bytes,
            self.record_count,
            self.index_offset,
            self.mane_index_offset,
            self.decoded_cache,
        )
    }
}

// ============================================================================
// Readonly file-backed constructors
// ============================================================================

impl ReadonlyTxStore {
    /// Open a readonly transcript store from disk using the given backend.
    ///
    /// Native-only. Callers pick the backend explicitly — there is no silent
    /// guesswork. Use [`TxBackend::Pread`] for the small / single-process case
    /// (the conservative default) and [`TxBackend::Mmap`] for the large,
    /// multi-process-shared (VRS fan-out) case.
    pub fn open_with_backend<P: AsRef<Path>>(path: P, backend: TxBackend) -> Result<Self> {
        match backend {
            TxBackend::Mmap => Self::open_mmap(path),
            TxBackend::Pread => Self::open_pread(path),
        }
    }

    /// Open a readonly transcript store backed by a whole-file memory map.
    pub fn open_mmap<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())?;
        // SAFETY: file is published write-once by the atomic-rename builder.
        let mmap = unsafe { Mmap::map(&file)? };
        Self::from_byte_source(TxBytes::Mmap(mmap))
    }

    /// Open a readonly transcript store backed by positioned (`pread`) reads.
    pub fn open_pread<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())?;
        let len = file.metadata()?.len() as usize;
        let source = PreadSource {
            file: Arc::new(file),
            len,
        };
        Self::from_byte_source(TxBytes::Pread(source))
    }
}
