//! Query FASTA reading and the concatenated query block.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/blastinput/blast_fasta_input.cpp` (defaults: id = first
//! whitespace-delimited token, `Query_N` when the defline is empty),
//! `algo/blast/api/blast_setup_cxx.cpp` (`SetupQueryInfo_OMF`, `s_QueryInfo_SetContext`,
//! `SetupQueries_OMF`, `s_AddMask`: one blastna buffer, a sentinel before and after every
//! context, context 2i = plus strand, 2i+1 = reverse complement), `algo/blast/core/blast_filter.c`
//! (`BlastSeqLocCombine`, `BLAST_ComplementMaskLocations`, `BlastSetUp_MaskQuery`,
//! `Blast_MaskTheResidues`) and `blast_util.c::BlastCompressBlastnaSequence`.

use super::dust;
use super::encoding::{ASCII_TO_BLASTNA, NUCL_MASK, NUCL_SENTINEL, reverse_complement};
use super::options::{BlastOptions, Strand};
use crate::error::Error;
use serde::{Deserialize, Serialize};
use std::path::Path;

/// One query as read from FASTA.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Query {
    pub id: String,
    /// Full defline text after '>'.
    pub defline: String,
    /// blastna letters.
    pub seq: Vec<u8>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct QuerySet {
    pub queries: Vec<Query>,
}

impl QuerySet {
    pub fn from_fasta_path(p: &Path) -> Result<Self, Error> {
        let bytes = std::fs::read(p).map_err(|e| Error::io(p, e))?;
        Self::from_fasta_bytes(&bytes)
    }

    /// Parse FASTA. Lowercase is accepted (no masking), `U` reads as `T`, any other byte that
    /// is not an IUPAC nucleotide letter is `Error::Query` naming the position.
    pub fn from_fasta_bytes(b: &[u8]) -> Result<Self, Error> {
        let mut queries = Vec::new();
        let mut cur: Option<Query> = None;
        let mut n = 0usize;
        let mut pos_in_seq = 0usize;
        let mut line_no = 0usize;
        // CFastaReader numbers every query (`Query_N`); the id printed is the title's first
        // token, or `Query_N` when the defline is empty
        let new_query = |text: &str, n: usize| -> Query {
            let id = text
                .split_whitespace()
                .next()
                .map(|s| s.to_string())
                .unwrap_or_else(|| format!("Query_{n}"));
            Query {
                id,
                defline: text.to_string(),
                seq: Vec::new(),
            }
        };
        for line in b.split(|&c| c == b'\n') {
            line_no += 1;
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            if let Some(rest) = line.strip_prefix(b">") {
                if let Some(q) = cur.take() {
                    queries.push(q);
                }
                n += 1;
                let text = String::from_utf8_lossy(rest).trim().to_string();
                cur = Some(new_query(&text, n));
                pos_in_seq = 0;
                continue;
            }
            if line.iter().all(|c| c.is_ascii_whitespace()) {
                continue;
            }
            if cur.is_none() {
                // Sequence without a defline: blast names it Query_1.
                n += 1;
                cur = Some(new_query("", n));
                pos_in_seq = 0;
            }
            let q = cur.as_mut().unwrap();
            for &c in line {
                if c.is_ascii_whitespace() {
                    continue;
                }
                let v = ASCII_TO_BLASTNA[c as usize];
                if v == 255 || v == 15 {
                    return Err(Error::Query {
                        id: q.id.clone(),
                        pos: pos_in_seq,
                        msg: format!(
                            "byte {:?} (0x{c:02x}) on line {line_no} is not a nucleotide letter",
                            c as char
                        ),
                    });
                }
                q.seq.push(v);
                pos_in_seq += 1;
            }
        }
        if let Some(q) = cur.take() {
            queries.push(q);
        }
        if queries.is_empty() {
            return Err(Error::Query {
                id: String::new(),
                pos: 0,
                msg: "no sequences found in the query file".into(),
            });
        }
        Ok(QuerySet { queries })
    }
}

/// `BlastContextInfo` (the layout part; statistics live in `options::ContextStats`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct ContextInfo {
    pub query_offset: i32,
    pub query_length: i32,
    pub query_index: u32,
    pub frame: i32,
    pub is_valid: bool,
}

/// `BLAST_SequenceBlk` + `BlastQueryInfo` for the concatenated queries.
#[derive(Debug, Clone)]
pub struct QueryBlock {
    /// `sequence_start`: sentinel, ctx0, sentinel, ctx1, sentinel, ...
    pub buf: Vec<u8>,
    /// `query->length`: last offset + last length.
    pub length: i32,
    pub contexts: Vec<ContextInfo>,
    pub first_context: usize,
    pub max_length: i32,
    pub min_length: i32,
    /// Mask intervals per query, plus-strand coordinates, inclusive, sorted and combined.
    pub masks: Vec<Vec<(i32, i32)>>,
    /// Query letters before hard masking (`sequence_nomask`), when it was applied.
    pub nomask: Option<Vec<u8>>,
    /// `compressed_nuc_seq` (SmallNa tables), index `-3..len` stored at `compressed[i+3]`.
    pub compressed: Option<Vec<u8>>,
    pub strand: Strand,
}

impl QueryBlock {
    /// `sequence` (sentinel-free view: `buf[1..]`).
    #[inline]
    pub fn sequence(&self) -> &[u8] {
        &self.buf[1..]
    }

    #[inline]
    pub fn context_seq(&self, c: usize) -> &[u8] {
        let ci = &self.contexts[c];
        let o = ci.query_offset as usize + 1;
        &self.buf[o..o + ci.query_length as usize]
    }

    /// `context_seq` plus the sentinel after it: `BlastGetStartForGappedAlignmentNucl` can
    /// put the gapped start at `query_length` and `ALIGN_EX` then reads that sentinel.
    #[inline]
    pub fn context_seq_with_sentinel(&self, c: usize) -> &[u8] {
        let ci = &self.contexts[c];
        let o = ci.query_offset as usize + 1;
        &self.buf[o..o + ci.query_length as usize + 1]
    }

    pub fn num_queries(&self) -> usize {
        self.contexts.len() / 2
    }

    /// Build the block: dust (per `o.dust`), contexts per `o.strand`, hard masking unless
    /// `o.mask_at_hash`.
    pub fn build(qs: &QuerySet, o: &BlastOptions) -> QueryBlock {
        let nq = qs.queries.len();
        // masks (plus-strand coordinates)
        let mut masks: Vec<Vec<(i32, i32)>> = Vec::with_capacity(nq);
        for q in &qs.queries {
            let m = match o.dust {
                Some(d) => {
                    let iupac: Vec<u8> = q
                        .seq
                        .iter()
                        .map(|&b| super::encoding::BLASTNA_TO_IUPAC[b as usize])
                        .collect();
                    dust::dust_query(&iupac, d.level as u32, d.window as usize, d.linker as usize)
                }
                None => Vec::new(),
            };
            masks.push(combine(m));
        }
        // contexts (SetupQueryInfo_OMF)
        let mut contexts = Vec::with_capacity(nq * 2);
        let mut max_length = 0i32;
        let mut min_length = i32::MAX;
        let set_ctx = |contexts: &mut Vec<ContextInfo>, qi: u32, frame: i32, len: i32| {
            let off = match contexts.last() {
                Some(prev) => {
                    let shift = if prev.query_length != 0 {
                        prev.query_length + 1
                    } else {
                        0
                    };
                    prev.query_offset + shift
                }
                None => 0,
            };
            contexts.push(ContextInfo {
                query_offset: off,
                query_length: len,
                query_index: qi,
                frame,
                is_valid: len != 0,
            });
        };
        for (qi, q) in qs.queries.iter().enumerate() {
            let len = q.seq.len() as i32;
            max_length = max_length.max(len);
            min_length = min_length.min(len);
            match o.strand {
                Strand::Plus => {
                    set_ctx(&mut contexts, qi as u32, 1, len);
                    set_ctx(&mut contexts, qi as u32, -1, 0);
                    min_length = 0;
                }
                Strand::Minus => {
                    set_ctx(&mut contexts, qi as u32, 1, 0);
                    set_ctx(&mut contexts, qi as u32, -1, len);
                    min_length = 0;
                }
                Strand::Both => {
                    set_ctx(&mut contexts, qi as u32, 1, len);
                    set_ctx(&mut contexts, qi as u32, -1, len);
                }
            }
        }
        let last = contexts.last().unwrap();
        let length = last.query_offset + last.query_length;
        // buffer (SetupQueries_OMF): calloc(buflen+1) with buflen = length + 2
        let mut buf = vec![NUCL_SENTINEL; length as usize + 2];
        for c in &contexts {
            if c.query_length == 0 {
                continue;
            }
            let q = &qs.queries[c.query_index as usize];
            let o0 = c.query_offset as usize + 1;
            if c.frame == 1 {
                buf[o0..o0 + q.seq.len()].copy_from_slice(&q.seq);
            } else {
                let rc = reverse_complement(&q.seq);
                buf[o0..o0 + rc.len()].copy_from_slice(&rc);
            }
        }
        let first_context = if o.strand == Strand::Minus { 1 } else { 0 };
        let mut qb = QueryBlock {
            buf,
            length,
            contexts,
            first_context,
            max_length,
            min_length,
            masks,
            nomask: None,
            compressed: None,
            strand: o.strand,
        };
        if !o.mask_at_hash {
            qb.mask_residues();
        }
        qb
    }

    /// `BlastSetUp_MaskQuery`: replace masked letters with N (14) in every valid context.
    fn mask_residues(&mut self) {
        if self.masks.iter().all(|m| m.is_empty()) {
            return;
        }
        self.nomask = Some(self.buf.clone());
        for c in 0..self.contexts.len() {
            let ci = self.contexts[c];
            if !ci.is_valid {
                continue;
            }
            let len = ci.query_length;
            let base = ci.query_offset as usize + 1;
            let reverse = ci.frame < 0;
            for &(l, r) in &self.masks[ci.query_index as usize] {
                let (start, stop) = if reverse {
                    (len - 1 - r, len - 1 - l)
                } else {
                    (l, r)
                };
                for i in start..=stop {
                    self.buf[base + i as usize] = NUCL_MASK;
                }
            }
        }
    }

    /// `BLAST_ComplementMaskLocations`: the unmasked intervals of every valid context, in
    /// `sequence` coordinates, inclusive, in the order the C builds them.
    pub fn lookup_segments(&self) -> Vec<(i32, i32)> {
        let mut out = Vec::new();
        for c in self.first_context..self.contexts.len() {
            let ci = &self.contexts[c];
            if !ci.is_valid {
                continue;
            }
            let start_offset = ci.query_offset;
            let end_offset = ci.query_length + start_offset - 1;
            let m = &self.masks[ci.query_index as usize];
            if m.is_empty() {
                out.push((start_offset, end_offset));
                continue;
            }
            let reverse = ci.frame < 0;
            let ranges: Vec<(i32, i32)> = if reverse {
                m.iter()
                    .rev()
                    .map(|&(l, r)| (end_offset - r, end_offset - l))
                    .collect()
            } else {
                m.iter()
                    .map(|&(l, r)| (start_offset + l, start_offset + r))
                    .collect()
            };
            let mut first = true;
            let mut last_interval_open = true;
            let mut left = 0;
            for (filter_start, filter_end) in ranges {
                if first {
                    last_interval_open = true;
                    first = false;
                    if filter_start > start_offset {
                        left = start_offset;
                    } else {
                        left = filter_end + 1;
                        continue;
                    }
                }
                let right = filter_start - 1;
                out.push((left, right));
                if filter_end >= end_offset {
                    last_interval_open = false;
                    break;
                } else {
                    left = filter_end + 1;
                }
            }
            if last_interval_open {
                out.push((left, end_offset));
            }
        }
        out
    }

    /// `EstimateNumTableEntries`: (num_entries, max_q_off).
    pub fn estimate_table_entries(segments: &[(i32, i32)]) -> (i32, i32) {
        let mut n = 0i32;
        let mut mx = 0i32;
        for &(l, r) in segments {
            n += r - l;
            mx = mx.max(r);
        }
        (n, mx)
    }

    /// `BlastCompressBlastnaSequence`: `compressed_nuc_seq_start`, with `compressed[i+3]` =
    /// the byte for position `i` of `sequence` (i from -3).
    pub fn compress(&mut self) {
        let len = self.length as usize;
        let old = self.sequence();
        let mut new = vec![0u8; len + 3];
        let max_start = len.min(3);
        let mut curr: u32 = 0;
        let mut i = 0;
        while i < max_start {
            curr = (curr << 2) | (old[i] & 3) as u32;
            new[3 + i - max_start] = curr as u8;
            i += 1;
        }
        while i < len {
            curr = (curr << 2) | (old[i] & 3) as u32;
            new[3 + i - max_start] = curr as u8;
            i += 1;
        }
        for i in 0..max_start {
            curr <<= 2;
            new[3 + len - (max_start - i)] = curr as u8;
        }
        self.compressed = Some(new);
    }
}

/// `BlastSeqLocCombine(.., 0)`: sort by start and merge overlapping intervals.
fn combine(mut m: Vec<(i32, i32)>) -> Vec<(i32, i32)> {
    if m.is_empty() {
        return m;
    }
    m.sort_by_key(|x| x.0);
    let mut out: Vec<(i32, i32)> = Vec::with_capacity(m.len());
    out.push(m[0]);
    for &(l, r) in &m[1..] {
        let last = out.last_mut().unwrap();
        if last.1 > l {
            last.1 = last.1.max(r);
        } else {
            out.push((l, r));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fasta_ids_and_letters() {
        let qs = QuerySet::from_fasta_bytes(b">q1 some title\nACGTn\nacgu\n>\nTTTT\n").unwrap();
        assert_eq!(qs.queries.len(), 2);
        assert_eq!(qs.queries[0].id, "q1");
        assert_eq!(qs.queries[0].defline, "q1 some title");
        assert_eq!(qs.queries[0].seq, vec![0, 1, 2, 3, 14, 0, 1, 2, 3]);
        assert_eq!(qs.queries[1].id, "Query_2");
        let e = QuerySet::from_fasta_bytes(b">q\nACGX\n").unwrap_err();
        assert!(matches!(e, Error::Query { pos: 3, .. }), "{e}");
    }

    #[test]
    fn block_layout_both_strands() {
        let qs = QuerySet::from_fasta_bytes(b">a\nAACG\n>b\nTT\n").unwrap();
        let mut o = BlastOptions::megablast();
        o.dust = None;
        let qb = QueryBlock::build(&qs, &o);
        assert_eq!(qb.contexts.len(), 4);
        assert_eq!(qb.contexts[0].query_offset, 0);
        assert_eq!(qb.contexts[1].query_offset, 5);
        assert_eq!(qb.contexts[2].query_offset, 10);
        assert_eq!(qb.contexts[3].query_offset, 13);
        assert_eq!(qb.length, 15);
        assert_eq!(qb.buf.len(), 17);
        // ctx1 = revcomp(AACG) = CGTT
        assert_eq!(qb.context_seq(1), &[1, 2, 3, 3]);
        assert_eq!(qb.buf[0], 15);
        assert_eq!(qb.buf[5], 15);
        assert_eq!(
            qb.lookup_segments(),
            vec![(0, 3), (5, 8), (10, 11), (13, 14)]
        );
    }

    #[test]
    fn strand_plus_and_masks() {
        let qs = QuerySet::from_fasta_bytes(b">a\nAACGTTGCA\n").unwrap();
        let mut o = BlastOptions::megablast();
        o.dust = None;
        o.strand = Strand::Plus;
        let mut qb = QueryBlock::build(&qs, &o);
        assert_eq!(qb.contexts[1].query_length, 0);
        assert!(!qb.contexts[1].is_valid);
        assert_eq!(qb.contexts[1].query_offset, 10);
        qb.masks[0] = vec![(2, 4)];
        assert_eq!(qb.lookup_segments(), vec![(0, 1), (5, 8)]);
        // minus strand of the same mask
        o.strand = Strand::Both;
        let mut qb = QueryBlock::build(&qs, &o);
        qb.masks[0] = vec![(0, 2), (6, 8)];
        // ctx0 [0,8]: masks [0,2],[6,8] -> segments (3,5); ctx1 offset 10: reversed masks
        // (10+8-8=10, 10+8-6=12), (16, 18) -> (13,15)
        assert_eq!(qb.lookup_segments(), vec![(3, 5), (13, 15)]);
        o.mask_at_hash = false;
        o.dust = Some(super::super::options::DustSpec::default());
        let qs2 = QuerySet::from_fasta_bytes(
            b">a\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGT\n",
        )
        .unwrap();
        let qb = QueryBlock::build(&qs2, &o);
        let masked = qb
            .context_seq(0)
            .iter()
            .take(64)
            .filter(|&&b| b == 14)
            .count();
        assert!(masked >= 50, "{masked} of the A run masked");
        assert!(qb.nomask.is_some());
    }

    #[test]
    fn compress_matches_c() {
        let qs = QuerySet::from_fasta_bytes(b">a\nACGTAC\n").unwrap();
        let mut o = BlastOptions::megablast();
        o.dust = None;
        o.strand = Strand::Plus;
        let mut qb = QueryBlock::build(&qs, &o);
        qb.compress();
        let c = qb.compressed.as_ref().unwrap();
        // compressed[i] holds letters i..i+3 (the C's new_seq[i], stored at c[i+3]); the
        // query length counts the trailing sentinel (15 & 3 = 3).
        assert_eq!(qb.length, 7);
        assert_eq!(c[3], 0b0001_1011); // ACGT
        assert_eq!(c[0], 0b0000_0000); // A (letters -3..0)
        assert_eq!(c[1], 0b0000_0001); // AC
        assert_eq!(c[2], 0b0000_0110); // ACG
        assert_eq!(c[3 + 3], 0b1100_0111); // T A C sentinel(3)
    }
}
