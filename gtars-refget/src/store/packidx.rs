//! The `sequences.pack.idx` sidecar: a sorted, fixed-width offset index that
//! maps a sequence digest to a byte span inside a sealed pack file.
//!
//! # Format
//!
//! ```text
//! byte 0..64    header, exactly 64 bytes, space-padded, '\n' at byte 63:
//!               "#gtars-pack-idx\tv1\trowbytes=100\trows=%016x"
//! byte 64..     rows, 100 bytes each, sorted bytewise by digest:
//!               digest(32) '\t' pack(32-hex) '\t' offset(16-hex) '\t' len(16-hex) '\n'
//! ```
//!
//! NOTE: the plan specified a 43-char digest / 111-byte row. In this codebase a
//! sha512t24u digest is base64url of 24 bytes = exactly **32** chars, so the row
//! is 100 bytes. The format is otherwise as specified (fixed width, digest-
//! sorted, binary-searchable over range requests).
//!
//! `digest` = sha512t24u (always exactly 32 base64url chars); `pack` = the
//! content-hash pack name (32 hex chars; the file is `packs/<pack>.pack`);
//! `offset`/`len` = zero-padded lowercase hex. Row *i* lives at byte
//! `64 + i*100`, which is the binary-search invariant a remote client relies on
//! to bisect the file over HTTP range requests without downloading it whole.
//!
//! The sidecar is a SHARED artifact. It is merged and republished only under the
//! store lock inside `write_index_files`, exactly like `sequences.rgsi` — never
//! rewritten from one handle's in-memory view.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Arc;

use anyhow::{anyhow, Context, Result};

use super::atomic::atomic_write;

/// Fixed width of one sidecar row, in bytes.
pub(crate) const ROW_BYTES: usize = 100;

/// Fixed width of the sidecar header, in bytes.
pub(crate) const HEADER_BYTES: usize = 64;

/// A digest length (sha512t24u = base64url of 24 bytes) in chars.
const DIGEST_LEN: usize = 32;

/// A resolved sidecar row: where a sequence's bytes live.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PackRow {
    /// The content-hash pack name (32 hex chars); the file is `packs/<pack>.pack`.
    pub(crate) pack: Arc<str>,
    /// Byte offset of the sequence within that pack.
    pub(crate) offset: u64,
    /// Byte length of the sequence within that pack.
    pub(crate) len: u64,
}

/// Build the 64-byte header for `rows` rows.
fn header_bytes(rows: u64) -> [u8; HEADER_BYTES] {
    let s = format!("#gtars-pack-idx\tv1\trowbytes={}\trows={:016x}", ROW_BYTES, rows);
    let mut buf = [b' '; HEADER_BYTES];
    let bytes = s.as_bytes();
    debug_assert!(bytes.len() < HEADER_BYTES);
    buf[..bytes.len()].copy_from_slice(bytes);
    buf[HEADER_BYTES - 1] = b'\n';
    buf
}

/// Encode one row into its fixed-width form.
fn encode_row(digest: &str, row: &PackRow) -> [u8; ROW_BYTES] {
    debug_assert_eq!(digest.len(), DIGEST_LEN, "digest must be 32 chars");
    debug_assert_eq!(row.pack.len(), 32, "pack name must be 32 hex chars");
    let line = format!(
        "{}\t{}\t{:016x}\t{:016x}\n",
        digest, row.pack, row.offset, row.len
    );
    let bytes = line.as_bytes();
    debug_assert_eq!(bytes.len(), ROW_BYTES, "encoded row must be {} bytes", ROW_BYTES);
    let mut buf = [0u8; ROW_BYTES];
    buf.copy_from_slice(bytes);
    buf
}

/// Parse one fixed-width row into `(digest, PackRow)`.
fn parse_row(row: &[u8]) -> Result<(String, PackRow)> {
    if row.len() != ROW_BYTES {
        return Err(anyhow!("pack index row is {} bytes, expected {}", row.len(), ROW_BYTES));
    }
    let digest = std::str::from_utf8(&row[0..DIGEST_LEN])
        .context("pack index row digest not UTF-8")?
        .to_string();
    // Layout: digest(32) \t pack(32) \t offset(16) \t len(16) \n
    let pack = std::str::from_utf8(&row[33..65]).context("pack name not UTF-8")?;
    let offset = parse_hex(&row[66..82])?;
    let len = parse_hex(&row[83..99])?;
    Ok((
        digest,
        PackRow {
            pack: Arc::from(pack),
            offset,
            len,
        },
    ))
}

fn parse_hex(bytes: &[u8]) -> Result<u64> {
    let s = std::str::from_utf8(bytes).context("hex field not UTF-8")?;
    u64::from_str_radix(s, 16).with_context(|| format!("invalid hex field '{}'", s))
}

/// Atomically publish `sequences.pack.idx` from a digest-sorted row map.
///
/// The `BTreeMap` key ordering (bytewise over the digest string) is exactly the
/// on-disk row order, so two writes of the same rows are byte-identical — the
/// determinism guarantee that mirrors `write_sequences_rgsi_rows`.
pub(crate) fn write_pack_idx(path: &Path, rows: &BTreeMap<String, PackRow>) -> Result<()> {
    let header = header_bytes(rows.len() as u64);
    atomic_write(path, |w| {
        w.write_all(&header)?;
        for (digest, row) in rows {
            w.write_all(&encode_row(digest, row))?;
        }
        Ok(())
    })
}

/// Parse a whole sidecar file's bytes into a digest-keyed row map.
pub(crate) fn read_pack_idx_rows(bytes: &[u8]) -> Result<BTreeMap<String, PackRow>> {
    let mut out = BTreeMap::new();
    if bytes.is_empty() {
        return Ok(out);
    }
    if bytes.len() < HEADER_BYTES {
        return Err(anyhow!("pack index too small: {} bytes", bytes.len()));
    }
    let body = &bytes[HEADER_BYTES..];
    if !body.len().is_multiple_of(ROW_BYTES) {
        return Err(anyhow!(
            "pack index body {} bytes is not a multiple of row width {}",
            body.len(),
            ROW_BYTES
        ));
    }
    for chunk in body.chunks_exact(ROW_BYTES) {
        let (digest, row) = parse_row(chunk)?;
        out.insert(digest, row);
    }
    Ok(out)
}

/// A backend-agnostic sorted, fixed-width row store. One binary search serves
/// every backend (in-memory today; a ranged-HTTP backend is trivial to add).
pub(crate) trait RowSource: Send + Sync {
    fn rows(&self) -> u64;
    fn read_row(&self, i: u64) -> Result<[u8; ROW_BYTES]>;
}

/// In-memory sidecar backend: the whole file's bytes held in RAM.
#[derive(Debug)]
pub(crate) struct InMemory {
    bytes: Vec<u8>,
    rows: u64,
}

impl InMemory {
    /// Build from the raw sidecar file bytes.
    pub(crate) fn from_bytes(bytes: Vec<u8>) -> Result<Self> {
        let rows = if bytes.is_empty() {
            0
        } else {
            if bytes.len() < HEADER_BYTES {
                return Err(anyhow!("pack index too small: {} bytes", bytes.len()));
            }
            let body = bytes.len() - HEADER_BYTES;
            if !body.is_multiple_of(ROW_BYTES) {
                return Err(anyhow!(
                    "pack index body {} bytes is not a multiple of row width {}",
                    body,
                    ROW_BYTES
                ));
            }
            (body / ROW_BYTES) as u64
        };
        Ok(Self { bytes, rows })
    }
}

impl RowSource for InMemory {
    fn rows(&self) -> u64 {
        self.rows
    }
    fn read_row(&self, i: u64) -> Result<[u8; ROW_BYTES]> {
        let start = HEADER_BYTES + (i as usize) * ROW_BYTES;
        let end = start + ROW_BYTES;
        if end > self.bytes.len() {
            return Err(anyhow!("pack index row {} out of range", i));
        }
        let mut buf = [0u8; ROW_BYTES];
        buf.copy_from_slice(&self.bytes[start..end]);
        Ok(buf)
    }
}

/// A sorted fixed-width row store that resolves a digest to a `PackRow` via one
/// binary search. THE single point of pack-offset resolution for a backend.
#[derive(Debug)]
pub(crate) struct PackIndex<S: RowSource> {
    source: S,
}

impl<S: RowSource> PackIndex<S> {
    pub(crate) fn new(source: S) -> Self {
        Self { source }
    }

    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn rows(&self) -> u64 {
        self.source.rows()
    }

    /// Binary-search the sidecar for `digest`. `Ok(None)` on a genuine miss.
    pub(crate) fn lookup(&self, digest: &str) -> Result<Option<PackRow>> {
        let needle = digest.as_bytes();
        let (mut lo, mut hi) = (0u64, self.source.rows());
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let row = self.source.read_row(mid)?;
            match row[..DIGEST_LEN].cmp(needle) {
                Ordering::Less => lo = mid + 1,
                Ordering::Greater => hi = mid,
                Ordering::Equal => {
                    let (_d, parsed) = parse_row(&row)?;
                    return Ok(Some(parsed));
                }
            }
        }
        Ok(None)
    }
}

/// Convenience: an in-memory `PackIndex` built from raw sidecar bytes.
pub(crate) fn in_memory_index(bytes: Vec<u8>) -> Result<PackIndex<InMemory>> {
    Ok(PackIndex::new(InMemory::from_bytes(bytes)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

    fn row(pack: &str, offset: u64, len: u64) -> PackRow {
        PackRow { pack: Arc::from(pack), offset, len }
    }

    fn sample_rows() -> BTreeMap<String, PackRow> {
        let mut m = BTreeMap::new();
        // 32-char digests (base64url-ish). Use distinct sortable prefixes.
        m.insert("A".repeat(32), row(&"a".repeat(32), 0, 100));
        m.insert("M".repeat(32), row(&"b".repeat(32), 100, 250));
        m.insert("Z".repeat(32), row(&"c".repeat(32), 0, 42));
        m
    }

    #[test]
    fn header_is_64_bytes_newline_terminated() {
        let h = header_bytes(3);
        assert_eq!(h.len(), HEADER_BYTES);
        assert_eq!(h[HEADER_BYTES - 1], b'\n');
        let s = std::str::from_utf8(&h[..HEADER_BYTES - 1]).unwrap();
        assert!(s.starts_with("#gtars-pack-idx\tv1\trowbytes=100\trows=0000000000000003"));
    }

    #[test]
    fn round_trip_write_read() {
        let rows = sample_rows();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sequences.pack.idx");
        write_pack_idx(&path, &rows).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        // File length is header + N rows.
        assert_eq!(bytes.len(), HEADER_BYTES + rows.len() * ROW_BYTES);
        let parsed = read_pack_idx_rows(&bytes).unwrap();
        assert_eq!(parsed, rows);
    }

    #[test]
    fn two_writes_are_byte_identical() {
        let rows = sample_rows();
        let dir = tempfile::tempdir().unwrap();
        let p1 = dir.path().join("a.idx");
        let p2 = dir.path().join("b.idx");
        write_pack_idx(&p1, &rows).unwrap();
        write_pack_idx(&p2, &rows).unwrap();
        assert_eq!(std::fs::read(&p1).unwrap(), std::fs::read(&p2).unwrap());
    }

    #[test]
    fn lookup_hit_first_last_and_miss() {
        let rows = sample_rows();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sequences.pack.idx");
        write_pack_idx(&path, &rows).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        let idx = in_memory_index(bytes).unwrap();

        let first = "A".repeat(32);
        let last = "Z".repeat(32);
        let mid = "M".repeat(32);
        assert_eq!(idx.lookup(&first).unwrap(), Some(row(&"a".repeat(32), 0, 100)));
        assert_eq!(idx.lookup(&mid).unwrap(), Some(row(&"b".repeat(32), 100, 250)));
        assert_eq!(idx.lookup(&last).unwrap(), Some(row(&"c".repeat(32), 0, 42)));
        assert_eq!(idx.lookup(&"B".repeat(32)).unwrap(), None);
        assert_eq!(idx.lookup(&"0".repeat(32)).unwrap(), None);
    }

    #[test]
    fn empty_index_lookup_is_none() {
        let idx = in_memory_index(Vec::new()).unwrap();
        assert_eq!(idx.rows(), 0);
        assert_eq!(idx.lookup(&"A".repeat(32)).unwrap(), None);
    }

    /// A `RowSource` that counts `read_row` calls, to bound binary-search cost.
    struct CountingSource {
        inner: InMemory,
        reads: AtomicU64,
    }
    impl RowSource for CountingSource {
        fn rows(&self) -> u64 {
            self.inner.rows()
        }
        fn read_row(&self, i: u64) -> Result<[u8; ROW_BYTES]> {
            self.reads.fetch_add(1, AtomicOrdering::SeqCst);
            self.inner.read_row(i)
        }
    }

    #[test]
    fn lookup_reads_are_logarithmic() {
        // Build 1000 rows with sortable 43-char digests.
        let mut m = BTreeMap::new();
        for i in 0..1000u32 {
            let d = format!("{:032}", i);
            m.insert(d, row(&"a".repeat(32), i as u64, 1));
        }
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("s.idx");
        write_pack_idx(&path, &m).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        let source = CountingSource {
            inner: InMemory::from_bytes(bytes).unwrap(),
            reads: AtomicU64::new(0),
        };
        let idx = PackIndex::new(source);
        let target = format!("{:032}", 777u32);
        assert!(idx.lookup(&target).unwrap().is_some());
        // ceil(log2(1000)) == 10; allow a small slack.
        let reads = idx.source.reads.load(AtomicOrdering::SeqCst);
        assert!(reads <= 12, "binary search used {} reads (> ~log2 n)", reads);
    }
}
