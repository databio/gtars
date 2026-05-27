//! Binary transcript store with mmap and concurrent access.

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use anyhow::{anyhow, Result};
use memmap2::Mmap;

use crate::models::{Exon, ManeStatus, Strand, Transcript};

pub(crate) const MAGIC: &[u8; 4] = b"RFTX";
pub(crate) const VERSION: u32 = 2;
pub(crate) const HEADER_SIZE: usize = 40;
pub(crate) const INDEX_ENTRY_SIZE: usize = 16;
pub(crate) const NONE_SENTINEL: u32 = 0xFFFFFFFF;

// ============================================================================
// TxStore: Mutable store for setup/build phase
// ============================================================================

/// Mutable transcript store for loading and setup.
///
/// Use this during initialization, then convert to `ReadonlyTxStore` for
/// concurrent access:
///
/// ```ignore
/// let store = TxStore::open("transcripts.reftx")?;
/// let readonly = store.into_readonly();
/// let shared = Arc::new(readonly);
/// // Now safe for concurrent &self access
/// ```
pub struct TxStore {
    mmap: Mmap,
    record_count: u64,
    index_offset: u64,
    mane_index_offset: u64,
    /// Pre-decoded transcript cache (populated during into_readonly).
    decoded_cache: HashMap<u64, Transcript>,
}

impl TxStore {
    /// Open a transcript store from disk (mmap).
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())?;
        let mmap = unsafe { Mmap::map(&file)? };

        // Validate header
        if mmap.len() < HEADER_SIZE {
            return Err(anyhow!("File too small for header"));
        }
        if &mmap[0..4] != MAGIC {
            return Err(anyhow!("Invalid magic number: expected RFTX"));
        }
        let version = u32::from_le_bytes(mmap[4..8].try_into()?);
        if version != VERSION {
            return Err(anyhow!("Unsupported format version: {}", version));
        }
        let record_count = u64::from_le_bytes(mmap[8..16].try_into()?);
        let index_offset = u64::from_le_bytes(mmap[16..24].try_into()?);
        let mane_index_offset = u64::from_le_bytes(mmap[24..32].try_into()?);

        Ok(Self {
            mmap,
            record_count,
            index_offset,
            mane_index_offset,
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
        let idx_start = self.index_offset as usize;
        let count = self.record_count as usize;

        let mut lo = 0usize;
        let mut hi = count;

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let entry_offset = idx_start + mid * INDEX_ENTRY_SIZE;
            let entry_hash = u64::from_le_bytes(
                self.mmap[entry_offset..entry_offset + 8].try_into().ok()?,
            );

            match entry_hash.cmp(&hash) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    let record_offset = u64::from_le_bytes(
                        self.mmap[entry_offset + 8..entry_offset + 16].try_into().ok()?,
                    ) as usize;
                    let tx = read_record(&self.mmap, record_offset)?;
                    if tx.accession == accession {
                        return Some(tx);
                    }
                    return linear_probe(&self.mmap, self.index_offset, accession, hash, mid, count);
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
        let mut count = 0;

        for i in 0..self.record_count as usize {
            let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
            let hash = u64::from_le_bytes(
                self.mmap[entry_offset..entry_offset + 8].try_into()?,
            );
            let record_offset = u64::from_le_bytes(
                self.mmap[entry_offset + 8..entry_offset + 16].try_into()?,
            ) as usize;

            if let Some(tx) = read_record(&self.mmap, record_offset) {
                if predicate(&tx) {
                    self.decoded_cache.insert(hash, tx);
                    count += 1;
                }
            }
        }
        Ok(count)
    }

    /// Convert to immutable store for concurrent access.
    pub fn into_readonly(mut self) -> ReadonlyTxStore {
        if self.record_count < 500_000 && self.decoded_cache.is_empty() {
            let _ = self.ensure_decoded_where(|_| true);
        }

        ReadonlyTxStore {
            mmap: self.mmap,
            record_count: self.record_count,
            index_offset: self.index_offset,
            mane_index_offset: self.mane_index_offset,
            transcript_cache: self.decoded_cache,
        }
    }

    /// Convert to immutable store with lazy decoding.
    pub fn into_readonly_lazy(self) -> ReadonlyTxStore {
        ReadonlyTxStore {
            mmap: self.mmap,
            record_count: self.record_count,
            index_offset: self.index_offset,
            mane_index_offset: self.mane_index_offset,
            transcript_cache: self.decoded_cache,
        }
    }
}

// ============================================================================
// ReadonlyTxStore: Immutable store for concurrent access
// ============================================================================

/// Immutable transcript store for concurrent read access.
pub struct ReadonlyTxStore {
    mmap: Mmap,
    record_count: u64,
    index_offset: u64,
    mane_index_offset: u64,
    transcript_cache: HashMap<u64, Transcript>,
}

// SAFETY: Mmap is read-only after construction, HashMap is immutable.
unsafe impl Send for ReadonlyTxStore {}
unsafe impl Sync for ReadonlyTxStore {}

impl ReadonlyTxStore {
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

        if let Some(tx) = self.transcript_cache.get(&hash) {
            if tx.accession == accession {
                return Some(TranscriptRef::Cached(tx));
            }
        }

        self.lookup_from_mmap(accession, hash)
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
        let normalized = gene.to_uppercase();
        let hash = fnv1a_64(normalized.as_bytes());

        let idx_base = self.mane_index_offset as usize;
        if idx_base + 8 > self.mmap.len() {
            return None;
        }
        let mane_count = u64::from_le_bytes(
            self.mmap[idx_base..idx_base + 8].try_into().ok()?,
        ) as usize;
        let idx_start = idx_base + 8;

        let mut lo = 0usize;
        let mut hi = mane_count;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let entry_offset = idx_start + mid * INDEX_ENTRY_SIZE;
            let entry_hash = u64::from_le_bytes(
                self.mmap[entry_offset..entry_offset + 8].try_into().ok()?,
            );
            match entry_hash.cmp(&hash) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    let record_offset = u64::from_le_bytes(
                        self.mmap[entry_offset + 8..entry_offset + 16].try_into().ok()?,
                    ) as usize;
                    let tx = read_record(&self.mmap, record_offset)?;
                    if tx.gene.to_uppercase() == normalized {
                        return Some(tx);
                    }
                    return mane_linear_probe(
                        &self.mmap,
                        idx_start,
                        mane_count,
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

    fn lookup_from_mmap(&self, accession: &str, hash: u64) -> Option<Transcript> {
        let idx_start = self.index_offset as usize;
        let count = self.record_count as usize;

        let mut lo = 0usize;
        let mut hi = count;

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let entry_offset = idx_start + mid * INDEX_ENTRY_SIZE;
            let entry_hash = u64::from_le_bytes(
                self.mmap[entry_offset..entry_offset + 8].try_into().ok()?,
            );

            match entry_hash.cmp(&hash) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    let record_offset = u64::from_le_bytes(
                        self.mmap[entry_offset + 8..entry_offset + 16].try_into().ok()?,
                    ) as usize;
                    let tx = read_record(&self.mmap, record_offset)?;
                    if tx.accession == accession {
                        return Some(tx);
                    }
                    return linear_probe(&self.mmap, self.index_offset, accession, hash, mid, count);
                }
            }
        }
        None
    }
}

/// Read a record from the mmap at the given byte offset.
pub(crate) fn read_record(mmap: &[u8], offset: usize) -> Option<Transcript> {
    let data = &mmap[offset..];
    let mut pos = 0;

    let acc_len = *data.get(pos)? as usize;
    pos += 1;
    let accession = std::str::from_utf8(data.get(pos..pos + acc_len)?).ok()?.to_string();
    pos += acc_len;

    let gene_len = *data.get(pos)? as usize;
    pos += 1;
    let gene = std::str::from_utf8(data.get(pos..pos + gene_len)?).ok()?.to_string();
    pos += gene_len;

    let mut chrom_digest = [0u8; 24];
    chrom_digest.copy_from_slice(data.get(pos..pos + 24)?);
    pos += 24;

    let strand = Strand::from_i8(data[pos] as i8)?;
    pos += 1;

    // MANE flags byte (v2). For v1 stores this would be missing, but we always
    // read v2 now per the bumped VERSION constant.
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
    mmap: &[u8],
    idx_start: usize,
    count: usize,
    gene_upper: &str,
    hash: u64,
    start_idx: usize,
) -> Option<Transcript> {
    for i in (start_idx + 1)..count {
        let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
        let entry_hash =
            u64::from_le_bytes(mmap[entry_offset..entry_offset + 8].try_into().ok()?);
        if entry_hash != hash {
            break;
        }
        let record_offset = u64::from_le_bytes(
            mmap[entry_offset + 8..entry_offset + 16].try_into().ok()?,
        ) as usize;
        let tx = read_record(mmap, record_offset)?;
        if tx.gene.to_uppercase() == gene_upper {
            return Some(tx);
        }
    }
    for i in (0..start_idx).rev() {
        let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
        let entry_hash =
            u64::from_le_bytes(mmap[entry_offset..entry_offset + 8].try_into().ok()?);
        if entry_hash != hash {
            break;
        }
        let record_offset = u64::from_le_bytes(
            mmap[entry_offset + 8..entry_offset + 16].try_into().ok()?,
        ) as usize;
        let tx = read_record(mmap, record_offset)?;
        if tx.gene.to_uppercase() == gene_upper {
            return Some(tx);
        }
    }
    None
}

fn linear_probe(
    mmap: &[u8],
    index_offset: u64,
    accession: &str,
    hash: u64,
    start_idx: usize,
    count: usize,
) -> Option<Transcript> {
    let idx_start = index_offset as usize;

    for i in (start_idx + 1)..count {
        let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
        let entry_hash =
            u64::from_le_bytes(mmap[entry_offset..entry_offset + 8].try_into().ok()?);
        if entry_hash != hash {
            break;
        }
        let record_offset = u64::from_le_bytes(
            mmap[entry_offset + 8..entry_offset + 16].try_into().ok()?,
        ) as usize;
        let tx = read_record(mmap, record_offset)?;
        if tx.accession == accession {
            return Some(tx);
        }
    }

    for i in (0..start_idx).rev() {
        let entry_offset = idx_start + i * INDEX_ENTRY_SIZE;
        let entry_hash =
            u64::from_le_bytes(mmap[entry_offset..entry_offset + 8].try_into().ok()?);
        if entry_hash != hash {
            break;
        }
        let record_offset = u64::from_le_bytes(
            mmap[entry_offset + 8..entry_offset + 16].try_into().ok()?,
        ) as usize;
        let tx = read_record(mmap, record_offset)?;
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
