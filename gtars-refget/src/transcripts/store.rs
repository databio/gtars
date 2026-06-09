//! Binary transcript store: WASM-safe lookup logic over an abstract byte source.
//!
//! This module is the CORE of the transcript store. It contains:
//! - the `.reftx` format constants and the header/record parse logic,
//! - the binary-search lookup over the sorted hash index,
//! - a tri-backend byte source ([`TxBytes`]) and `ReadonlyTxStore::from_bytes`,
//!   which build a working readonly transcript store from owned bytes.
//!
//! All parsing/search logic operates through the [`TxBytes`] abstraction
//! (`len` / `read_at` / `read_u64_le`), so it is completely agnostic to where
//! the bytes came from. The `InMemory` variant is always present and compiles
//! on ALL targets (including `wasm32-unknown-unknown`); the file-backed `Mmap`
//! and `Pread` variants are native-only, gated behind the `filesystem` feature
//! and provided by [`super::mmap`]. The on-disk `.reftx` layout is identical
//! across all three backends — only HOW bytes are read differs.

use std::collections::HashMap;
use std::sync::Arc;

use anyhow::{anyhow, Result};

use crate::transcripts::models::{Exon, ManeStatus, Strand, Transcript};

pub(crate) const MAGIC: &[u8; 4] = b"RFTX";
pub(crate) const VERSION: u32 = 2;
pub(crate) const HEADER_SIZE: usize = 40;
pub(crate) const INDEX_ENTRY_SIZE: usize = 16;
pub(crate) const NONE_SENTINEL: u32 = 0xFFFFFFFF;

// ============================================================================
// Byte source: the ONLY place the three backends differ.
// ============================================================================

/// The byte source backing one open `.reftx` store.
///
/// This is the single abstraction over which all lookup/parse logic operates.
/// Exactly one variant is selected at construction time:
///
/// - [`TxBytes::InMemory`] — WASM-safe. Owns the bytes in an `Arc<Vec<u8>>`.
///   Compiles on ALL targets; pulls in NO `memmap2` and NO `std::fs`. This is
///   the only backend when the `filesystem` feature is disabled (the wasm
///   config), and the analogue of `RefgetStore::in_memory()`.
/// - [`TxBytes::Mmap`] — native, behind `filesystem`. Memory-maps the whole
///   file; sound only because the builder publishes the file write-once via an
///   atomic rename (see [`super::builder`]). Best for the large, multi-process
///   shared (VRS fan-out) case where the OS page cache is shared.
///
///   SAFETY / FILE-INTEGRITY CONTRACT: the mapped `.reftx` file MUST be
///   immutable for the lifetime of the mapping. The atomic-rename builder
///   guarantees write-once publication, so no writer mutates the bytes under
///   us. If a `.reftx` is truncated or overwritten in place by an external
///   process while mapped, reads of the now-missing pages fault with SIGBUS
///   (a process abort Rust cannot catch). Callers that cannot guarantee an
///   immutable file should use [`TxBytes::Pread`] instead, which performs
///   bounds-checked positioned reads and surfaces a short read as a normal
///   error rather than a fault.
/// - [`TxBytes::Pread`] — native, behind `filesystem`. The conservative
///   default: positioned (`pread`) reads on a shared `&File`, no `unsafe`
///   mapping. Best for the small / single-process case (scattered small reads).
///
/// Every `match` on `TxBytes` keeps the `Mmap`/`Pread` arms gated so the enum
/// (and all lookup logic over it) compiles with ONLY the `InMemory` arm on wasm.
pub enum TxBytes {
    /// WASM-SAFE. Compiles on ALL targets. No memmap2, no std::fs.
    InMemory(Arc<Vec<u8>>),
    /// Native memory-mapped backend (whole-file mmap).
    #[cfg(feature = "filesystem")]
    Mmap(memmap2::Mmap),
    /// Native positioned-read backend (one `pread` per access).
    #[cfg(feature = "filesystem")]
    Pread(super::mmap::PreadSource),
}

impl TxBytes {
    /// Total length of the byte source, in bytes.
    #[inline]
    pub(crate) fn len(&self) -> usize {
        match self {
            TxBytes::InMemory(buf) => buf.len(),
            #[cfg(feature = "filesystem")]
            TxBytes::Mmap(mmap) => mmap.len(),
            #[cfg(feature = "filesystem")]
            TxBytes::Pread(p) => p.len(),
        }
    }

    /// Read `len` bytes starting at `offset`, or `None` if out of bounds.
    ///
    /// For `InMemory`/`Mmap` this is a (cheap) copy out of an in-memory slice;
    /// for `Pread` this is one positioned read.
    #[inline]
    pub(crate) fn read_at(&self, offset: usize, len: usize) -> Option<Vec<u8>> {
        match self {
            TxBytes::InMemory(buf) => buf.get(offset..offset.checked_add(len)?).map(<[u8]>::to_vec),
            #[cfg(feature = "filesystem")]
            TxBytes::Mmap(mmap) => mmap
                .get(offset..offset.checked_add(len)?)
                .map(<[u8]>::to_vec),
            #[cfg(feature = "filesystem")]
            TxBytes::Pread(p) => p.read_at(offset, len),
        }
    }

    /// Read a little-endian `u64` at `offset` (the hot index-probe path), or
    /// `None` if out of bounds. Avoids a heap allocation for `InMemory`/`Mmap`.
    #[inline]
    pub(crate) fn read_u64_le(&self, offset: usize) -> Option<u64> {
        match self {
            TxBytes::InMemory(buf) => {
                let s = buf.get(offset..offset.checked_add(8)?)?;
                Some(u64::from_le_bytes(s.try_into().ok()?))
            }
            #[cfg(feature = "filesystem")]
            TxBytes::Mmap(mmap) => {
                let s = mmap.get(offset..offset.checked_add(8)?)?;
                Some(u64::from_le_bytes(s.try_into().ok()?))
            }
            #[cfg(feature = "filesystem")]
            TxBytes::Pread(p) => {
                let buf = p.read_at(offset, 8)?;
                Some(u64::from_le_bytes(buf.as_slice().try_into().ok()?))
            }
        }
    }
}

/// Parsed `.reftx` header fields.
pub(crate) struct ReftxHeader {
    pub(crate) record_count: u64,
    pub(crate) index_offset: u64,
    pub(crate) mane_index_offset: u64,
}

/// Validate the `.reftx` MAGIC/VERSION and read the header fields.
///
/// Shared by the in-memory (`ReadonlyTxStore::from_bytes`) and native
/// (`super::mmap`) construction paths — the on-disk format is
/// target-independent and read entirely through the [`TxBytes`] abstraction.
pub(crate) fn parse_header(bytes: &TxBytes) -> Result<ReftxHeader> {
    if bytes.len() < HEADER_SIZE {
        return Err(anyhow!("File too small for header"));
    }
    let head = bytes
        .read_at(0, HEADER_SIZE)
        .ok_or_else(|| anyhow!("Failed to read header"))?;
    if &head[0..4] != MAGIC {
        return Err(anyhow!("Invalid magic number: expected RFTX"));
    }
    let version = u32::from_le_bytes(head[4..8].try_into()?);
    if version != VERSION {
        return Err(anyhow!("Unsupported format version: {}", version));
    }
    let record_count = u64::from_le_bytes(head[8..16].try_into()?);
    let index_offset = u64::from_le_bytes(head[16..24].try_into()?);
    let mane_index_offset = u64::from_le_bytes(head[24..32].try_into()?);
    Ok(ReftxHeader {
        record_count,
        index_offset,
        mane_index_offset,
    })
}

// ============================================================================
// ReadonlyTxStore: Immutable store for concurrent read access (WASM-safe)
// ============================================================================

/// Immutable transcript store for concurrent read access.
///
/// Backed by a [`TxBytes`] byte source: the WASM-safe in-memory buffer
/// ([`ReadonlyTxStore::from_bytes`]) on all targets, or the native mmap/pread
/// backends produced by `super::mmap::TxStore::into_readonly()` (gated behind
/// `filesystem`).
pub struct ReadonlyTxStore {
    bytes: TxBytes,
    record_count: u64,
    index_offset: u64,
    mane_index_offset: u64,
    transcript_cache: HashMap<u64, Transcript>,
}

// Soundness of `Send`/`Sync`:
//   - `TxBytes::InMemory(Arc<Vec<u8>>)`: auto `Send + Sync`.
//   - `TxBytes::Pread(PreadSource)` wraps `Arc<File>` + `usize`: auto `Send + Sync`.
//   - `TxBytes::Mmap(memmap2::Mmap)`: `Mmap` is read-only and is itself
//     `Send + Sync`; it is sound here only because the file it maps is published
//     write-once via an atomic rename (see `super::builder`), so no writer can
//     mutate the mapping under a reader.
// `HashMap<u64, Transcript>` is immutable after construction. Every field is
// therefore `Send + Sync` by auto-derivation, so NO manual `unsafe impl` is
// needed.

impl ReadonlyTxStore {
    /// Build a readonly transcript store from owned `.reftx` bytes (WASM-safe).
    ///
    /// Wraps the bytes in `TxBytes::InMemory`, runs the same MAGIC/header parse
    /// as the native path, then serves lookups from the owned buffer. Uses NO
    /// `std::fs` and NO `memmap2`, so it builds on `wasm32-unknown-unknown` and
    /// is the only backend when `filesystem` is disabled. Small stores are
    /// pre-decoded (per the native heuristic) so lookups never re-touch bytes.
    pub fn from_bytes(bytes: Vec<u8>) -> Result<Self> {
        let buf = TxBytes::InMemory(Arc::new(bytes));
        Self::from_byte_source(buf)
    }

    /// Construct from any [`TxBytes`] backend after header validation.
    ///
    /// Shared by `from_bytes` (in-memory) and the native mmap/pread
    /// constructors. Pre-decodes a small store into the transcript cache using
    /// the same `record_count < 500_000` heuristic the native `into_readonly`
    /// used, so the in-memory/wasm case never re-touches the buffer on lookup.
    pub(crate) fn from_byte_source(bytes: TxBytes) -> Result<Self> {
        let header = parse_header(&bytes)?;
        let mut store = Self {
            bytes,
            record_count: header.record_count,
            index_offset: header.index_offset,
            mane_index_offset: header.mane_index_offset,
            transcript_cache: HashMap::new(),
        };
        if store.record_count < 500_000 {
            store.predecode_all();
        }
        Ok(store)
    }

    /// Construct from a [`TxBytes`] backend with header fields and a
    /// pre-decoded cache already supplied (native `into_readonly` path).
    #[cfg(feature = "filesystem")]
    pub(crate) fn from_parts(
        bytes: TxBytes,
        record_count: u64,
        index_offset: u64,
        mane_index_offset: u64,
        transcript_cache: HashMap<u64, Transcript>,
    ) -> Self {
        Self {
            bytes,
            record_count,
            index_offset,
            mane_index_offset,
            transcript_cache,
        }
    }

    /// Pre-decode every record into the transcript cache (small-store path).
    fn predecode_all(&mut self) {
        let idx_start = self.index_offset as usize;
        for i in 0..self.record_count as usize {
            let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
            let hash = match self.bytes.read_u64_le(entry_offset) {
                Some(h) => h,
                None => break,
            };
            let record_offset = match self.bytes.read_u64_le(entry_offset + 8) {
                Some(o) => o as usize,
                None => break,
            };
            let next_bound = self.index_offset as usize;
            if let Some(tx) = read_record(&self.bytes, record_offset, next_bound) {
                self.transcript_cache.insert(hash, tx);
            }
        }
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

    /// Look up transcript by accession. O(log n) binary search.
    pub fn lookup(&self, accession: &str) -> Option<TranscriptRef<'_>> {
        let hash = fnv1a_64(accession.as_bytes());

        if let Some(tx) = self.transcript_cache.get(&hash)
            && tx.accession == accession
        {
            return Some(TranscriptRef::Cached(tx));
        }

        self.lookup_from_bytes(accession, hash)
            .map(TranscriptRef::Owned)
    }

    /// Look up the MANE Select transcript for a gene symbol (case-insensitive).
    ///
    /// O(log n) binary search over the MANE gene index. Returns `None` if the
    /// gene has no MANE Select transcript.
    pub fn lookup_mane(&self, gene: &str) -> Option<Transcript> {
        if self.mane_index_offset == 0 {
            return None;
        }
        let bytes = &self.bytes;
        let record_bound = self.index_offset as usize;
        let normalized = gene.to_uppercase();
        let hash = fnv1a_64(normalized.as_bytes());

        let idx_base = self.mane_index_offset as usize;
        let mane_count = bytes.read_u64_le(idx_base)? as usize;
        let idx_start = idx_base + 8;

        let mut lo = 0usize;
        let mut hi = mane_count;
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
                    if tx.gene.to_uppercase() == normalized {
                        return Some(tx);
                    }
                    return mane_linear_probe(
                        bytes,
                        idx_start,
                        mane_count,
                        record_bound,
                        &normalized,
                        hash,
                        mid,
                    );
                }
            }
        }
        None
    }

    /// Returns true if the store has any MANE index entries.
    #[inline]
    pub fn has_mane_index(&self) -> bool {
        self.mane_index_offset != 0
    }

    fn lookup_from_bytes(&self, accession: &str, hash: u64) -> Option<Transcript> {
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
}

/// Read a record from the byte source at the given byte offset.
///
/// The record length is not stored, so a bounded window is read up to
/// `record_bound` (the start of the index region) — a record read never spills
/// into the index. For `InMemory`/`Mmap` this is a slice; for `Pread` it is one
/// positioned read of the bounded window.
pub(crate) fn read_record(
    bytes: &TxBytes,
    offset: usize,
    record_bound: usize,
) -> Option<Transcript> {
    if offset > record_bound {
        return None;
    }
    let window = record_bound - offset;
    let data = bytes.read_at(offset, window)?;
    let data = data.as_slice();
    let mut pos = 0;

    let acc_len = *data.get(pos)? as usize;
    pos += 1;
    let accession = std::str::from_utf8(data.get(pos..pos + acc_len)?)
        .ok()?
        .to_string();
    pos += acc_len;

    let gene_len = *data.get(pos)? as usize;
    pos += 1;
    let gene = std::str::from_utf8(data.get(pos..pos + gene_len)?)
        .ok()?
        .to_string();
    pos += gene_len;

    let mut chrom_digest = [0u8; 24];
    chrom_digest.copy_from_slice(data.get(pos..pos + 24)?);
    pos += 24;

    let strand = Strand::from_i8(*data.get(pos)? as i8)?;
    pos += 1;

    // MANE flags byte (v2).
    let mane = ManeStatus::from_flags_byte(*data.get(pos)?);
    pos += 1;

    let cds_start_raw = u32::from_le_bytes(data.get(pos..pos + 4)?.try_into().ok()?);
    pos += 4;
    let cds_end_raw = u32::from_le_bytes(data.get(pos..pos + 4)?.try_into().ok()?);
    pos += 4;
    let cds_start = if cds_start_raw == NONE_SENTINEL {
        None
    } else {
        Some(cds_start_raw)
    };
    let cds_end = if cds_end_raw == NONE_SENTINEL {
        None
    } else {
        Some(cds_end_raw)
    };

    let exon_count = u16::from_le_bytes(data.get(pos..pos + 2)?.try_into().ok()?) as usize;
    pos += 2;
    let mut exons = Vec::with_capacity(exon_count);
    for _ in 0..exon_count {
        let start = u32::from_le_bytes(data.get(pos..pos + 4)?.try_into().ok()?);
        pos += 4;
        let end = u32::from_le_bytes(data.get(pos..pos + 4)?.try_into().ok()?);
        pos += 4;
        exons.push(Exon { start, end });
    }

    Some(Transcript {
        accession,
        gene,
        chrom_digest,
        strand,
        cds_start,
        cds_end,
        exons,
        mane,
    })
}

/// Linear probe over MANE index for hash collisions.
fn mane_linear_probe(
    bytes: &TxBytes,
    idx_start: usize,
    count: usize,
    record_bound: usize,
    gene_upper: &str,
    hash: u64,
    start_idx: usize,
) -> Option<Transcript> {
    for i in (start_idx + 1)..count {
        let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
        let entry_hash = bytes.read_u64_le(entry_offset)?;
        if entry_hash != hash {
            break;
        }
        let record_offset = bytes.read_u64_le(entry_offset + 8)? as usize;
        let tx = read_record(bytes, record_offset, record_bound)?;
        if tx.gene.to_uppercase() == gene_upper {
            return Some(tx);
        }
    }
    for i in (0..start_idx).rev() {
        let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
        let entry_hash = bytes.read_u64_le(entry_offset)?;
        if entry_hash != hash {
            break;
        }
        let record_offset = bytes.read_u64_le(entry_offset + 8)? as usize;
        let tx = read_record(bytes, record_offset, record_bound)?;
        if tx.gene.to_uppercase() == gene_upper {
            return Some(tx);
        }
    }
    None
}

/// Linear probe over the accession index for hash collisions.
pub(crate) fn linear_probe(
    bytes: &TxBytes,
    index_offset: u64,
    record_bound: usize,
    accession: &str,
    hash: u64,
    start_idx: usize,
    count: usize,
) -> Option<Transcript> {
    let idx_start = index_offset as usize;

    for i in (start_idx + 1)..count {
        let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
        let entry_hash = bytes.read_u64_le(entry_offset)?;
        if entry_hash != hash {
            break;
        }
        let record_offset = bytes.read_u64_le(entry_offset + 8)? as usize;
        let tx = read_record(bytes, record_offset, record_bound)?;
        if tx.accession == accession {
            return Some(tx);
        }
    }

    for i in (0..start_idx).rev() {
        let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
        let entry_hash = bytes.read_u64_le(entry_offset)?;
        if entry_hash != hash {
            break;
        }
        let record_offset = bytes.read_u64_le(entry_offset + 8)? as usize;
        let tx = read_record(bytes, record_offset, record_bound)?;
        if tx.accession == accession {
            return Some(tx);
        }
    }

    None
}

/// Reference to a transcript: either cached (zero-copy) or freshly decoded.
pub enum TranscriptRef<'a> {
    Cached(&'a Transcript),
    Owned(Transcript),
}

impl<'a> std::ops::Deref for TranscriptRef<'a> {
    type Target = Transcript;

    fn deref(&self) -> &Self::Target {
        match self {
            TranscriptRef::Cached(tx) => tx,
            TranscriptRef::Owned(tx) => tx,
        }
    }
}

/// FNV-1a 64-bit hash (stable, deterministic).
#[inline]
pub(crate) fn fnv1a_64(data: &[u8]) -> u64 {
    const FNV_OFFSET: u64 = 0xcbf29ce484222325;
    const FNV_PRIME: u64 = 0x100000001b3;
    let mut hash = FNV_OFFSET;
    for &byte in data {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

// ============================================================================
// The ONE `.reftx` encoder (per-record serializer + whole-file assembler).
//
// This is the single, authoritative writer for the `.reftx` format; it lives
// next to the decoder (`read_record`/`parse_header`) and the format constants
// so a layout change touches encoder + decoder + constants together. It is
// WASM-safe (depends only on `models`, the constants, and `fnv1a_64` — no
// `std::fs`, no `memmap2`), so it is available to the core (no-`filesystem`)
// test config, to the native `builder`, and — via the public
// `build_reftx_bytes_in_memory` wrapper — to the `tests/` integration files.
// ============================================================================

/// Serialize one transcript record into `buf` (clears it first).
///
/// Mirrors [`read_record`] exactly. Fails if the accession or gene exceeds the
/// u8 length field (255 bytes) rather than silently truncating via `as u8`.
/// This is the ONLY length-validation point in the encoder.
pub(crate) fn serialize_record_into(buf: &mut Vec<u8>, tx: &Transcript) -> Result<()> {
    buf.clear();
    if tx.accession.len() > u8::MAX as usize {
        return Err(anyhow!(
            "accession {:?} is {} bytes; exceeds 255-byte .reftx limit",
            tx.accession,
            tx.accession.len()
        ));
    }
    if tx.gene.len() > u8::MAX as usize {
        return Err(anyhow!(
            "gene {:?} is {} bytes; exceeds 255-byte .reftx limit",
            tx.gene,
            tx.gene.len()
        ));
    }
    buf.push(tx.accession.len() as u8);
    buf.extend_from_slice(tx.accession.as_bytes());
    buf.push(tx.gene.len() as u8);
    buf.extend_from_slice(tx.gene.as_bytes());
    buf.extend_from_slice(&tx.chrom_digest);
    buf.push(tx.strand.to_byte());
    buf.push(tx.mane.to_flags_byte());
    buf.extend_from_slice(&tx.cds_start.unwrap_or(NONE_SENTINEL).to_le_bytes());
    buf.extend_from_slice(&tx.cds_end.unwrap_or(NONE_SENTINEL).to_le_bytes());
    let exon_count: u16 = tx.exons.len().try_into().map_err(|_| {
        anyhow!(
            "transcript {:?} has {} exons; exceeds {}-exon .reftx limit",
            tx.accession,
            tx.exons.len(),
            u16::MAX
        )
    })?;
    buf.extend_from_slice(&exon_count.to_le_bytes());
    for exon in &tx.exons {
        buf.extend_from_slice(&exon.start.to_le_bytes());
        buf.extend_from_slice(&exon.end.to_le_bytes());
    }
    Ok(())
}

/// Build a complete `.reftx` byte image entirely in memory.
///
/// Reproduces, byte-for-byte, the layout `read_record`/`parse_header` parse:
/// a 40-byte header, then the records sorted by `fnv1a_64(accession)`, then the
/// sorted accession index, then (only if any `mane_select` transcript exists)
/// the sorted MANE gene index. Sorts a borrowed view rather than mutating the
/// caller's slice. Returns `Err` if any record fails the 255-byte length
/// validation in [`serialize_record_into`].
pub(crate) fn build_reftx_bytes(transcripts: &[Transcript]) -> Result<Vec<u8>> {
    let mut sorted: Vec<&Transcript> = transcripts.iter().collect();
    sorted.sort_by_key(|t| fnv1a_64(t.accession.as_bytes()));

    let mut out = vec![0u8; HEADER_SIZE];
    let mut index_entries: Vec<(u64, u64)> = Vec::with_capacity(sorted.len());
    let mut mane_entries: Vec<(u64, u64)> = Vec::new();
    let mut current_offset = HEADER_SIZE as u64;
    let mut record_buf = Vec::with_capacity(512);

    for tx in &sorted {
        let hash = fnv1a_64(tx.accession.as_bytes());
        index_entries.push((hash, current_offset));
        if tx.mane.mane_select {
            let gene_hash = fnv1a_64(tx.gene.to_uppercase().as_bytes());
            mane_entries.push((gene_hash, current_offset));
        }
        serialize_record_into(&mut record_buf, tx)?;
        current_offset += record_buf.len() as u64;
        out.extend_from_slice(&record_buf);
    }

    let index_offset = current_offset;
    for (hash, offset) in &index_entries {
        out.extend_from_slice(&hash.to_le_bytes());
        out.extend_from_slice(&offset.to_le_bytes());
    }

    // MANE gene index: 0 if empty (sentinel for "no MANE info").
    let mane_index_offset = if mane_entries.is_empty() {
        0u64
    } else {
        mane_entries.sort_by_key(|(h, _)| *h);
        let off = out.len() as u64;
        out.extend_from_slice(&(mane_entries.len() as u64).to_le_bytes());
        for (hash, offset) in &mane_entries {
            out.extend_from_slice(&hash.to_le_bytes());
            out.extend_from_slice(&offset.to_le_bytes());
        }
        off
    };

    out[0..4].copy_from_slice(MAGIC);
    out[4..8].copy_from_slice(&VERSION.to_le_bytes());
    out[8..16].copy_from_slice(&(sorted.len() as u64).to_le_bytes());
    out[16..24].copy_from_slice(&index_offset.to_le_bytes());
    out[24..32].copy_from_slice(&mane_index_offset.to_le_bytes());
    // bytes [32..40] reserved (already zero)
    Ok(out)
}

/// Public, WASM-safe wrapper over [`build_reftx_bytes`] for the in-memory build
/// path and the `tests/` integration files (which can only see the crate's
/// public API). Mirrors `RefgetStore::in_memory` / `ReadonlyTxStore::from_bytes`.
pub fn build_reftx_bytes_in_memory(transcripts: &[Transcript]) -> Result<Vec<u8>> {
    build_reftx_bytes(transcripts)
}

// ============================================================================
// Core tests: in-memory round-trip (run under `transcripts` without `filesystem`)
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fnv1a_64_deterministic() {
        let hash1 = fnv1a_64(b"NM_004333.6");
        let hash2 = fnv1a_64(b"NM_004333.6");
        assert_eq!(hash1, hash2);

        let hash3 = fnv1a_64(b"NM_000546.6");
        assert_ne!(hash1, hash3);
    }

    /// Convenience: build via the ONE shared encoder, unwrapping the `Result`.
    fn build_reftx_bytes_ok(transcripts: &[Transcript]) -> Vec<u8> {
        super::build_reftx_bytes(transcripts).unwrap()
    }

    fn sample_transcript() -> Transcript {
        Transcript {
            accession: "NM_004333.6".to_string(),
            gene: "BRAF".to_string(),
            chrom_digest: [1u8; 24],
            strand: Strand::Forward,
            cds_start: Some(100),
            cds_end: Some(400),
            exons: vec![Exon { start: 50, end: 500 }],
            mane: ManeStatus {
                mane_select: true,
                mane_clinical: false,
            },
        }
    }

    #[test]
    fn test_from_bytes_round_trip_lookup() {
        let bytes = build_reftx_bytes_ok(&[sample_transcript()]);
        let store = ReadonlyTxStore::from_bytes(bytes).unwrap();
        assert_eq!(store.len(), 1);

        let tx = store.lookup("NM_004333.6").expect("must find transcript");
        assert_eq!(tx.gene, "BRAF");
        assert_eq!(tx.exons.len(), 1);
        assert!(matches!(tx.strand, Strand::Forward));

        assert!(store.lookup("NM_NONEXISTENT.1").is_none());
    }

    #[test]
    fn test_from_bytes_mane_lookup() {
        let mut other = sample_transcript();
        other.accession = "NM_OTHER.1".to_string();
        other.gene = "OTHER".to_string();
        other.mane = Default::default();

        let bytes = build_reftx_bytes_ok(&[sample_transcript(), other]);
        let store = ReadonlyTxStore::from_bytes(bytes).unwrap();
        assert!(store.has_mane_index());

        let tx = store.lookup_mane("BRAF").expect("MANE BRAF must resolve");
        assert_eq!(tx.accession, "NM_004333.6");
        assert!(tx.mane.mane_select);
        // Case-insensitive.
        assert_eq!(store.lookup_mane("braf").unwrap().accession, "NM_004333.6");
        // Non-MANE / missing genes.
        assert!(store.lookup_mane("OTHER").is_none());
        assert!(store.lookup_mane("MISSING").is_none());
    }

    #[test]
    fn test_from_bytes_invalid_magic() {
        let mut bytes = build_reftx_bytes_ok(&[sample_transcript()]);
        bytes[0] = b'X';
        assert!(ReadonlyTxStore::from_bytes(bytes).is_err());
    }

    #[test]
    fn test_n_to_g_zero_returns_error_not_panic() {
        use crate::transcripts::mapper::{CoordinateMapper, MappingError};

        let bytes = build_reftx_bytes_ok(&[sample_transcript()]);
        let store = ReadonlyTxStore::from_bytes(bytes).unwrap();

        let mapper = CoordinateMapper::new(&store);
        // n. positions are 1-based; n.0 is invalid and must not underflow/panic.
        let result = mapper.n_to_g("NM_004333.6", 0);
        assert!(matches!(result, Err(MappingError::OutsideTranscript(0))));
    }

    // ------------------------------------------------------------------------
    // Length validation: accession/gene must fit the u8 length field.
    // ------------------------------------------------------------------------

    #[test]
    fn test_serialize_record_rejects_oversized_accession() {
        let mut tx = sample_transcript();
        tx.accession = "A".repeat(256); // 256 bytes > 255-byte field
        let mut buf = Vec::new();
        assert!(serialize_record_into(&mut buf, &tx).is_err());
    }

    #[test]
    fn test_serialize_record_rejects_oversized_gene() {
        let mut tx = sample_transcript();
        tx.gene = "G".repeat(256);
        let mut buf = Vec::new();
        assert!(serialize_record_into(&mut buf, &tx).is_err());
    }

    #[test]
    fn test_serialize_record_accepts_max_length() {
        // Exactly 255 bytes must succeed for both fields.
        let mut tx = sample_transcript();
        tx.accession = "A".repeat(255);
        tx.gene = "G".repeat(255);
        let mut buf = Vec::new();
        assert!(serialize_record_into(&mut buf, &tx).is_ok());
    }

    #[test]
    fn test_build_reftx_bytes_propagates_length_error() {
        let mut tx = sample_transcript();
        tx.accession = "A".repeat(300);
        assert!(super::build_reftx_bytes(&[tx]).is_err());
    }

    // ------------------------------------------------------------------------
    // Golden fixture: a hand-built `.reftx` byte image NOT produced by the
    // crate's own encoder. This is the ONE test that fails if the on-disk
    // layout silently drifts (every other test round-trips through the shared
    // encoder, so a matched encoder+decoder bug would pass them). The bytes are
    // annotated with their field offsets; see the `golden.py` derivation in the
    // consolidation plan.
    // ------------------------------------------------------------------------

    #[test]
    fn test_decode_external_golden_fixture() {
        // Single transcript: accession "NM_1.1", gene "G", chrom_digest [0x07;24],
        // strand Forward (1), MANE select (flags 0x01), cds 100..200, one exon
        // 10..20. Header (40B) + record (53B) + accession index (16B) + MANE
        // index (8B count + 16B entry) = 133 bytes total.
        #[rustfmt::skip]
        let fixture: &[u8] = &[
            // header [0..40]
            0x52, 0x46, 0x54, 0x58, // MAGIC "RFTX"
            0x02, 0x00, 0x00, 0x00, // VERSION = 2
            0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // record_count = 1
            0x5d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // index_offset = 93
            0x6d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // mane_index_offset = 109
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // reserved [32..40]
            // record [40..93]
            0x06, 0x4e, 0x4d, 0x5f, 0x31, 0x2e, 0x31, // acc_len=6 "NM_1.1"
            0x01, 0x47, // gene_len=1 "G"
            0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, // chrom_digest[0..8]
            0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, // chrom_digest[8..16]
            0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, // chrom_digest[16..24]
            0x01, // strand Forward
            0x01, // MANE flags: mane_select
            0x64, 0x00, 0x00, 0x00, // cds_start = 100
            0xc8, 0x00, 0x00, 0x00, // cds_end = 200
            0x01, 0x00, // exon_count = 1
            0x0a, 0x00, 0x00, 0x00, // exon start = 10
            0x14, 0x00, 0x00, 0x00, // exon end = 20
            // accession index [93..109]: hash(NM_1.1) + offset 40
            0xa1, 0x1a, 0xf1, 0xdd, 0xe4, 0x0a, 0xd6, 0x3f,
            0x28, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            // MANE index [109..133]: count 1, then hash(upper "G") + offset 40
            0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x86, 0x1f, 0x02, 0x86, 0x4c, 0xfa, 0x63, 0xaf,
            0x28, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let store = ReadonlyTxStore::from_bytes(fixture.to_vec())
            .expect("hand-built fixture must decode");
        assert_eq!(store.len(), 1);

        let tx = store.lookup("NM_1.1").expect("must find NM_1.1");
        assert_eq!(tx.gene, "G");
        assert_eq!(tx.chrom_digest, [0x07u8; 24]);
        assert!(matches!(tx.strand, Strand::Forward));
        assert!(tx.mane.mane_select);
        assert_eq!(tx.cds_start, Some(100));
        assert_eq!(tx.cds_end, Some(200));
        assert_eq!(tx.exons.len(), 1);
        assert_eq!(tx.exons[0].start, 10);
        assert_eq!(tx.exons[0].end, 20);

        // MANE lookup by gene (case-insensitive).
        assert_eq!(store.lookup_mane("g").unwrap().accession, "NM_1.1");
    }

    #[test]
    fn test_shared_encoder_matches_golden_layout() {
        // Byte-identity: the shared encoder, fed the SAME transcript the golden
        // fixture encodes, must reproduce the fixture byte-for-byte. This ties
        // the (encoder-produced) and (hand-built) layouts together, proving the
        // refactor did not change the on-disk bytes.
        let tx = Transcript {
            accession: "NM_1.1".to_string(),
            gene: "G".to_string(),
            chrom_digest: [0x07u8; 24],
            strand: Strand::Forward,
            cds_start: Some(100),
            cds_end: Some(200),
            exons: vec![Exon { start: 10, end: 20 }],
            mane: ManeStatus {
                mane_select: true,
                mane_clinical: false,
            },
        };
        let produced = super::build_reftx_bytes(&[tx]).unwrap();

        #[rustfmt::skip]
        let fixture: Vec<u8> = vec![
            0x52, 0x46, 0x54, 0x58, 0x02, 0x00, 0x00, 0x00,
            0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x5d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x6d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x06, 0x4e, 0x4d, 0x5f, 0x31, 0x2e, 0x31, 0x01,
            0x47, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
            0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
            0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07, 0x07,
            0x07, 0x01, 0x01, 0x64, 0x00, 0x00, 0x00, 0xc8,
            0x00, 0x00, 0x00, 0x01, 0x00, 0x0a, 0x00, 0x00,
            0x00, 0x14, 0x00, 0x00, 0x00, 0xa1, 0x1a, 0xf1,
            0xdd, 0xe4, 0x0a, 0xd6, 0x3f, 0x28, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x86, 0x1f, 0x02,
            0x86, 0x4c, 0xfa, 0x63, 0xaf, 0x28, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00,
        ];
        assert_eq!(produced, fixture, "shared encoder drifted from golden layout");
    }
}
