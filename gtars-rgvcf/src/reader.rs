//! rgvcf SPEC-v1 reader. Mmap-only; no copies for variant bytes.

use std::fs::File;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use memmap2::Mmap;
use smallvec::SmallVec;

use crate::{DIGEST_LEN, MAGIC_FOOTER, MAGIC_HEADER, VERSION};

/// Per-chromosome header-table entry.
pub struct ChromEntry {
    pub name: String,
    /// Raw base64url sha512t24u digest (32 ASCII bytes; no `SQ.` prefix).
    /// A 32-byte sequence of `0x00` indicates an unresolved chromosome in
    /// the source collection — the encoder emitted this as a sentinel.
    pub seq_digest: String,
    pub variant_count: u32,
    pub stream_byte_offset: u64,
    pub stream_byte_len: u64,
}

/// A decoded rgvcf record. Inline SNVs carry only their ALT (REF is
/// recovered via the refget store); sidecar entries carry both the
/// explicit REF and the list of ALTs.
pub enum Record<'a> {
    Snv { chrom_idx: u16, pos: u64, alt: u8 },
    Complex {
        chrom_idx: u16,
        pos: u64,
        ref_: &'a [u8],
        alts: SmallVec<[&'a [u8]; 2]>,
    },
}

/// Read-only rgvcf reader. Owns the `Mmap` handle and exposes
/// zero-copy views into it.
pub struct RgvcfReader {
    mmap: Mmap,
    collection_digest: String,
    chroms: Vec<ChromEntry>,
    sidecar_entry_offsets: Vec<u64>,
}

impl RgvcfReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file =
            File::open(path.as_ref()).with_context(|| format!("open {:?}", path.as_ref()))?;
        let mmap = unsafe { Mmap::map(&file).context("mmap rgvcf")? };
        Self::from_mmap(mmap)
    }

    pub fn from_mmap(mmap: Mmap) -> Result<Self> {
        let buf: &[u8] = &mmap[..];
        // Header.
        if buf.len() < 12 || &buf[0..4] != MAGIC_HEADER {
            bail!("bad rgvcf header magic");
        }
        if buf[4] != VERSION {
            bail!("unsupported rgvcf version: {}", buf[4]);
        }
        // flags @ buf[5], reserved buf[6..8]; probe-level checks skipped here.

        let mut off = 8usize;
        let digest_len = read_u8(buf, &mut off)? as usize;
        if digest_len != DIGEST_LEN {
            bail!("collection_digest_len = {}, expected {}", digest_len, DIGEST_LEN);
        }
        let collection_digest = std::str::from_utf8(&buf[off..off + DIGEST_LEN])
            .context("collection digest utf-8")?
            .to_owned();
        off += DIGEST_LEN;

        // Chrom table.
        let num_chroms = read_u16_le(buf, &mut off)? as usize;
        let mut chroms = Vec::with_capacity(num_chroms);
        for _ in 0..num_chroms {
            let name_len = read_u8(buf, &mut off)? as usize;
            let name = std::str::from_utf8(&buf[off..off + name_len])
                .context("chrom name utf-8")?
                .to_owned();
            off += name_len;
            let seq_digest = std::str::from_utf8(&buf[off..off + DIGEST_LEN])
                .context("seq digest utf-8")?
                .to_owned();
            off += DIGEST_LEN;
            let variant_count = read_u32_le(buf, &mut off)?;
            let stream_byte_offset = read_u64_le(buf, &mut off)?;
            let stream_byte_len = read_u64_le(buf, &mut off)?;
            chroms.push(ChromEntry {
                name,
                seq_digest,
                variant_count,
                stream_byte_offset,
                stream_byte_len,
            });
        }

        // Footer.
        let fsz = buf.len();
        if fsz < 12 || &buf[fsz - 4..fsz] != MAGIC_FOOTER {
            bail!("bad rgvcf footer magic");
        }
        let sidecar_byte_offset =
            u64::from_le_bytes(buf[fsz - 12..fsz - 4].try_into().unwrap()) as usize;

        // Sidecar index.
        let mut scur = sidecar_byte_offset;
        let num_sidecar = read_u32_le(buf, &mut scur)? as usize;
        let mut sidecar_entry_offsets = Vec::with_capacity(num_sidecar);
        for _ in 0..num_sidecar {
            sidecar_entry_offsets.push(read_u64_le(buf, &mut scur)?);
        }

        Ok(Self {
            mmap,
            collection_digest,
            chroms,
            sidecar_entry_offsets,
        })
    }

    pub fn collection_digest(&self) -> &str {
        &self.collection_digest
    }

    pub fn chromosomes(&self) -> &[ChromEntry] {
        &self.chroms
    }

    pub fn sidecar_entry_offsets(&self) -> &[u64] {
        &self.sidecar_entry_offsets
    }

    pub fn mmap_slice(&self) -> &[u8] {
        &self.mmap[..]
    }

    /// Iterate all records for the given chromosome name.
    pub fn iter_chrom(&self, chrom: &str) -> Result<RgvcfRecordIter<'_>> {
        let (idx, entry) = self
            .chroms
            .iter()
            .enumerate()
            .find(|(_, c)| c.name == chrom)
            .ok_or_else(|| anyhow!("chromosome {} not in rgvcf", chrom))?;
        Ok(self.iter_chrom_slice(idx, entry, 0, entry.stream_byte_len as usize, 0))
    }

    /// Iterate a sub-range of a chromosome's stream. `start_off` / `end_off`
    /// are relative to the stream's start; `start_prev_pos` is the cumulative
    /// absolute pos at the first record in the sub-range (0 if starting at
    /// record 0).
    ///
    /// This is the entry point for parallel chunked readers; pair with a
    /// call to [`scan_breakpoints`].
    pub fn iter_chrom_slice<'a>(
        &'a self,
        chrom_idx: usize,
        entry: &'a ChromEntry,
        start_off: usize,
        end_off: usize,
        start_prev_pos: u64,
    ) -> RgvcfRecordIter<'a> {
        let base = entry.stream_byte_offset as usize;
        let s_end = base + end_off;
        let stream = &self.mmap[base + start_off..s_end];
        RgvcfRecordIter {
            chrom_idx: chrom_idx as u16,
            mmap: &self.mmap[..],
            sidecar_entry_offsets: &self.sidecar_entry_offsets,
            stream,
            cursor: 0,
            prev_pos: start_prev_pos,
        }
    }

    /// Scan a chromosome's stream once, returning `(byte_offset_in_stream,
    /// absolute_pos_before_record)` at each of `num_chunks` equally-spaced
    /// record indices.
    pub fn scan_breakpoints(&self, chrom_idx: usize, num_chunks: usize) -> Result<Vec<(usize, u64)>> {
        let entry = &self.chroms[chrom_idx];
        let base = entry.stream_byte_offset as usize;
        let stream = &self.mmap[base..base + entry.stream_byte_len as usize];

        // Count records.
        let mut n = 0usize;
        let mut cursor = 0usize;
        while cursor < stream.len() {
            while cursor < stream.len() && stream[cursor] & 0x80 != 0 {
                cursor += 1;
            }
            cursor += 1;
            if cursor > stream.len() {
                bail!("malformed stream");
            }
            let code = stream[cursor - 1 + 1] & 0x07; // redundant; recompute below
            let _ = code;
            let code_byte = stream[cursor];
            cursor += 1;
            if code_byte & 0x07 == 4 {
                cursor += 4;
            }
            n += 1;
        }
        if n == 0 {
            return Ok(vec![(0, 0)]);
        }
        let mut targets: Vec<usize> = (0..num_chunks).map(|i| i * n / num_chunks).collect();
        targets.dedup();

        let mut bps = Vec::with_capacity(targets.len());
        let mut ti = 0usize;
        let mut cursor = 0usize;
        let mut prev_pos = 0u64;
        let mut rec_idx = 0usize;
        while cursor < stream.len() && ti < targets.len() {
            if rec_idx == targets[ti] {
                bps.push((cursor, prev_pos));
                ti += 1;
            }
            let mut val: u64 = 0;
            let mut shift = 0u32;
            loop {
                let b = stream[cursor];
                cursor += 1;
                val |= ((b & 0x7f) as u64) << shift;
                if b & 0x80 == 0 {
                    break;
                }
                shift += 7;
            }
            let code = stream[cursor] & 0x07;
            cursor += 1;
            if code == 4 {
                cursor += 4;
            }
            prev_pos += val;
            rec_idx += 1;
        }
        Ok(bps)
    }
}

/// Iterator over rgvcf records in a (range of a) chromosome stream.
pub struct RgvcfRecordIter<'a> {
    chrom_idx: u16,
    mmap: &'a [u8],
    sidecar_entry_offsets: &'a [u64],
    stream: &'a [u8],
    cursor: usize,
    prev_pos: u64,
}

impl<'a> Iterator for RgvcfRecordIter<'a> {
    type Item = Result<Record<'a>>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor >= self.stream.len() {
            return None;
        }
        // LEB128 pos_delta.
        let mut val: u64 = 0;
        let mut shift = 0u32;
        loop {
            if self.cursor >= self.stream.len() {
                return Some(Err(anyhow!("EOF reading LEB128")));
            }
            let b = self.stream[self.cursor];
            self.cursor += 1;
            val |= ((b & 0x7f) as u64) << shift;
            if b & 0x80 == 0 {
                break;
            }
            shift += 7;
            if shift >= 64 {
                return Some(Err(anyhow!("LEB128 overflow")));
            }
        }
        let pos = self.prev_pos + val;
        self.prev_pos = pos;

        if self.cursor >= self.stream.len() {
            return Some(Err(anyhow!("EOF reading code byte")));
        }
        let code = self.stream[self.cursor];
        self.cursor += 1;
        if code & 0xF8 != 0 {
            return Some(Err(anyhow!("reserved upper bits set in code byte: {}", code)));
        }
        let code_low = code & 0x07;
        if code_low < 4 {
            let alt = b"ACGT"[code_low as usize];
            Some(Ok(Record::Snv {
                chrom_idx: self.chrom_idx,
                pos,
                alt,
            }))
        } else if code_low == 4 {
            if self.cursor + 4 > self.stream.len() {
                return Some(Err(anyhow!("EOF reading sidecar_index")));
            }
            let sidecar_index = u32::from_le_bytes(
                self.stream[self.cursor..self.cursor + 4].try_into().unwrap(),
            ) as usize;
            self.cursor += 4;
            if sidecar_index >= self.sidecar_entry_offsets.len() {
                return Some(Err(anyhow!("sidecar_index out of bounds")));
            }
            let entry_off = self.sidecar_entry_offsets[sidecar_index] as usize;
            match parse_sidecar_entry(self.mmap, entry_off) {
                Ok((ci, ep, r, a)) => Some(Ok(Record::Complex {
                    chrom_idx: ci,
                    pos: ep,
                    ref_: r,
                    alts: a,
                })),
                Err(e) => Some(Err(e)),
            }
        } else {
            Some(Err(anyhow!("reserved code {}", code_low)))
        }
    }
}

fn parse_sidecar_entry(
    mmap: &[u8],
    off: usize,
) -> Result<(u16, u64, &[u8], SmallVec<[&[u8]; 2]>)> {
    let mut cur = off;
    let chrom_index = read_u16_le(mmap, &mut cur)?;
    let pos = read_u64_le(mmap, &mut cur)?;
    let ref_len = read_u16_le(mmap, &mut cur)? as usize;
    let ref_ = &mmap[cur..cur + ref_len];
    cur += ref_len;
    let num_alts = read_u8(mmap, &mut cur)? as usize;
    let mut alts: SmallVec<[&[u8]; 2]> = SmallVec::with_capacity(num_alts);
    for _ in 0..num_alts {
        let alt_len = read_u16_le(mmap, &mut cur)? as usize;
        alts.push(&mmap[cur..cur + alt_len]);
        cur += alt_len;
    }
    Ok((chrom_index, pos, ref_, alts))
}

// ---- low-level byte readers --------------------------------------------

fn read_u8(buf: &[u8], off: &mut usize) -> Result<u8> {
    if *off >= buf.len() {
        bail!("EOF u8");
    }
    let v = buf[*off];
    *off += 1;
    Ok(v)
}
fn read_u16_le(buf: &[u8], off: &mut usize) -> Result<u16> {
    if *off + 2 > buf.len() {
        bail!("EOF u16");
    }
    let v = u16::from_le_bytes(buf[*off..*off + 2].try_into().unwrap());
    *off += 2;
    Ok(v)
}
fn read_u32_le(buf: &[u8], off: &mut usize) -> Result<u32> {
    if *off + 4 > buf.len() {
        bail!("EOF u32");
    }
    let v = u32::from_le_bytes(buf[*off..*off + 4].try_into().unwrap());
    *off += 4;
    Ok(v)
}
fn read_u64_le(buf: &[u8], off: &mut usize) -> Result<u64> {
    if *off + 8 > buf.len() {
        bail!("EOF u64");
    }
    let v = u64::from_le_bytes(buf[*off..*off + 8].try_into().unwrap());
    *off += 8;
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    // Minimal SPEC-v1 file:
    //   header 'RGVC' v=1 flags=1 reserved
    //   anchor digest_len=32 digest="A"*32
    //   chrom table: 1 chrom "chr1", seq_digest="B"*32, variant_count=2
    //   stream:  pos=5 code=2 (SNV->G);  pos=10 code=0 (SNV->A)
    //   sidecar: 0 entries
    //   footer: sidecar_byte_offset + CFGR
    fn build_minimal() -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(b"RGVC");
        out.extend_from_slice(&[1u8, 1, 0, 0]);
        out.push(32); // digest_len
        out.extend_from_slice(&[b'A'; 32]);
        out.extend_from_slice(&(1u16).to_le_bytes()); // num_chroms
        // chrom 0: name "chr1"
        out.push(4);
        out.extend_from_slice(b"chr1");
        out.extend_from_slice(&[b'B'; 32]);
        out.extend_from_slice(&(2u32).to_le_bytes()); // variant_count
        // stream_byte_offset placeholder; fill after.
        let stream_offset_pos = out.len();
        out.extend_from_slice(&(0u64).to_le_bytes());
        let stream_len_pos = out.len();
        out.extend_from_slice(&(0u64).to_le_bytes());

        let stream_start = out.len();
        // record 1: pos=5 -> pos_delta=5 (LEB128 0x05), code=2
        out.push(5);
        out.push(2);
        // record 2: pos=10 -> pos_delta=5, code=0
        out.push(5);
        out.push(0);
        let stream_end = out.len();

        let sidecar_offset = out.len();
        out.extend_from_slice(&(0u32).to_le_bytes());

        // footer
        out.extend_from_slice(&(sidecar_offset as u64).to_le_bytes());
        out.extend_from_slice(b"CFGR");

        // patch stream offsets/lens.
        let off = stream_start as u64;
        let len = (stream_end - stream_start) as u64;
        out[stream_offset_pos..stream_offset_pos + 8]
            .copy_from_slice(&off.to_le_bytes());
        out[stream_len_pos..stream_len_pos + 8].copy_from_slice(&len.to_le_bytes());
        out
    }

    #[test]
    fn header_and_chrom_parse() {
        let bytes = build_minimal();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.as_file().write_all(&bytes).unwrap();
        let reader = RgvcfReader::open(tmp.path()).unwrap();
        assert_eq!(reader.collection_digest(), &"A".repeat(32));
        let c = &reader.chromosomes()[0];
        assert_eq!(c.name, "chr1");
        assert_eq!(c.seq_digest, "B".repeat(32));
        assert_eq!(c.variant_count, 2);
    }

    #[test]
    fn iter_chrom_snv() {
        let bytes = build_minimal();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.as_file().write_all(&bytes).unwrap();
        let reader = RgvcfReader::open(tmp.path()).unwrap();
        let recs: Vec<_> = reader
            .iter_chrom("chr1")
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();
        assert_eq!(recs.len(), 2);
        match &recs[0] {
            Record::Snv { pos, alt, .. } => {
                assert_eq!(*pos, 5);
                assert_eq!(*alt, b'G');
            }
            _ => panic!(),
        }
        match &recs[1] {
            Record::Snv { pos, alt, .. } => {
                assert_eq!(*pos, 10);
                assert_eq!(*alt, b'A');
            }
            _ => panic!(),
        }
    }

    #[test]
    fn scan_breakpoints_2() {
        let bytes = build_minimal();
        let tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.as_file().write_all(&bytes).unwrap();
        let reader = RgvcfReader::open(tmp.path()).unwrap();
        let bps = reader.scan_breakpoints(0, 2).unwrap();
        assert_eq!(bps.len(), 2);
        assert_eq!(bps[0], (0, 0));
        // second breakpoint at second record (index 1)
        // pos_delta=5 LEB128 is 1 byte (0x05), code=1 byte -> cursor=2
        assert_eq!(bps[1], (2, 5));
    }
}
