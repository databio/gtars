//! The subjects a search reads (`BlastSeqSrc` in BLAST+, see
//! `algo/blast/api/seqsrc_seqdb.cpp`): the foundation's `SubjectPlan` over the store's mapped
//! files (`store::access::MappedFiles`).
//!
//! `for_each_chunk` decodes a subject to ncbi2na (makeblastdb's layout with the random fill
//! for ambiguous bases, `convert::Pieces`) one chunk at a time, padded so the scanners may
//! over-read four bytes past a chunk as the C does. `blastna_window` decodes letters of one
//! subject, with a trailing sentinel, for the traceback stage and the aligned-sequence output;
//! `blastna_sparse` only the letters a long subject's traceback reads.

use super::prelim::subject_chunks;
use crate::convert::Pieces;
use crate::error::Error;
use crate::plan::SubjectPlan;
use crate::store::access::MappedFiles;
use crate::store::decode;
use gtars_refget::store::StorageMode;
use memmap2::MmapMut;
use std::ops::Range;
use std::sync::Arc;

/// Padding after each subject in the arena (the MB scanner reads up to 4 bytes past the word).
pub const ARENA_PAD: usize = 8;

/// A search thread's buffers for `StoreSubjects::for_each_chunk`: the decoded chunk and the
/// stored bytes read for it.
#[derive(Default)]
pub struct ChunkBuf {
    seq: Vec<u8>,
    stored: Vec<u8>,
}

/// `STATS_NSEQ` / `STATS_TOTLEN` of a wrapper alias.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct DbStatsOverride {
    pub num_seqs: u32,
    pub total_length: u64,
}

/// The plan + mapped files the engine searches.
pub struct StoreSubjects {
    pub plan: Arc<SubjectPlan>,
    access: MappedFiles,
    mode: StorageMode,
    stats: Option<DbStatsOverride>,
}

impl StoreSubjects {
    pub fn new(
        plan: Arc<SubjectPlan>,
        access: MappedFiles,
        stats: Option<DbStatsOverride>,
    ) -> StoreSubjects {
        let mode = plan.store_mode;
        StoreSubjects {
            plan,
            access,
            mode,
            stats,
        }
    }

    pub fn num_seqs(&self) -> i32 {
        self.plan.totals.num_seqs as i32
    }
    /// The real total length (what the query batching reads, whatever the STATS override).
    pub fn total_length(&self) -> i64 {
        self.plan.totals.total_length as i64
    }
    /// `BlastSeqSrcGetNumSeqsStats` / `BlastSeqSrcGetTotLenStats`: 0 when there is no override.
    pub fn num_seqs_stats(&self) -> i32 {
        self.stats.map(|s| s.num_seqs as i32).unwrap_or(0)
    }
    pub fn total_length_stats(&self) -> i64 {
        self.stats.map(|s| s.total_length as i64).unwrap_or(0)
    }
    pub fn max_len(&self) -> i32 {
        self.plan.totals.max_len as i32
    }
    /// `s_SeqDbGetAvgLength`: total length over the sequence count.
    pub fn avg_len(&self) -> i32 {
        (self.total_length() / (self.num_seqs().max(1) as i64)) as i32
    }
    pub fn seq_len(&self, oid: i32) -> i32 {
        self.plan.subjects[oid as usize].length as i32
    }
    /// Every subject's length, in OID order (`engine::search::build_work_items`).
    pub fn lengths(&self) -> Vec<u32> {
        self.plan.subjects.iter().map(|s| s.length).collect()
    }

    /// A subject chunk by chunk (`prelim::subject_chunks`; one chunk up to `MAX_DBSEQ_LEN`):
    /// `f` gets each chunk's ncbi2na bytes as they sit in the whole subject's (the chunk, then
    /// the `ARENA_PAD` bytes that follow it there; zeros after the last), its offset and its
    /// length. Each chunk is decoded just before `f` reads it, the random fill carrying over,
    /// so a thread holds one chunk, not a long subject. A large packed or raw record is read
    /// piece by piece too (`MappedFiles::read_at`); a small or zstd one whole.
    pub fn for_each_chunk(
        &self,
        oid: u32,
        buf: &mut ChunkBuf,
        mut f: impl FnMut(&[u8], i32, i32),
    ) -> Result<(), Error> {
        let s = &self.plan.subjects[oid as usize];
        let ChunkBuf { seq, stored } = buf;
        let mut run = || match (self.mode, s.stored_len) {
            (StorageMode::Encoded | StorageMode::Raw, Some(n)) if n >= MappedFiles::MMAP_MIN => {
                let alphabet = (self.mode == StorageMode::Encoded).then_some(s.alphabet);
                let mut pieces = Pieces::new(alphabet, n as usize, s.length)?;
                let bits = pieces.bits;
                let mut file = self.access.open(oid)?;
                chunks(s.length, seq, &mut f, |r, seq| {
                    let b = r.start as u64 * bits / 8;
                    stored.resize(((r.end as u64 * bits).div_ceil(8) - b) as usize, 0);
                    self.access.read_at(oid, &mut file, b, stored)?;
                    pieces.next(stored, b as usize, r.end, seq)
                })
            }
            _ => {
                let bytes = self.access.bytes(oid)?;
                let (mut pieces, src) =
                    decode::to_blast_pieces(self.mode, s.alphabet, bytes.as_slice(), s.length)?;
                chunks(s.length, seq, &mut f, |r, seq| {
                    pieces.next(&src, 0, r.end, seq)
                })
            }
        };
        run().map_err(|e| e.with_digest(&s.digest))
    }

    /// Letters `range` of the subject in blastna followed by one sentinel byte (15), into
    /// `out` (cleared first; callers reuse one buffer across windows).
    pub fn blastna_window(
        &self,
        oid: i32,
        range: Range<i32>,
        out: &mut Vec<u8>,
    ) -> Result<(), Error> {
        let s = &self.plan.subjects[oid as usize];
        let bytes = self.access.bytes(oid as u32)?;
        decode::to_blastna_window(
            self.mode,
            s.alphabet,
            bytes.as_slice(),
            s.length,
            range.start.max(0) as usize,
            range.end.clamp(0, s.length as i32) as usize,
            out,
        )
        .map_err(|e| e.with_digest(&s.digest))
    }

    /// The whole subject in blastna, as SeqDB's partial fetch hands blastn's traceback a long
    /// subject: only the letters of `ranges` are decoded, the rest of the buffer is anonymous
    /// memory never written (no memory until touched; small pages, or each window would take
    /// a 2 MB huge page), and a sentinel (15) follows the last letter. `tmp` is a reusable
    /// decode buffer.
    pub fn blastna_sparse(
        &self,
        oid: i32,
        ranges: &[Range<i32>],
        tmp: &mut Vec<u8>,
    ) -> Result<MmapMut, Error> {
        let len = self.seq_len(oid) as usize;
        let mut buf = MmapMut::map_anon(len + 1).map_err(|e| Error::io("anonymous memory", e))?;
        #[cfg(target_os = "linux")]
        let _ = buf.advise(memmap2::Advice::NoHugePage);
        for r in ranges {
            self.blastna_window(oid, r.clone(), tmp)?;
            let (a, b) = (r.start as usize, r.end as usize);
            buf[a..b].copy_from_slice(&tmp[..b - a]);
        }
        buf[len] = 15;
        Ok(buf)
    }
}

/// The chunk loop of `StoreSubjects::for_each_chunk`: `decode(range, buf)` appends the bytes
/// of letters `range` (consecutive ranges from 0) to `buf`.
fn chunks(
    len: u32,
    buf: &mut Vec<u8>,
    f: &mut impl FnMut(&[u8], i32, i32),
    mut decode: impl FnMut(Range<u32>, &mut Vec<u8>) -> Result<(), Error>,
) -> Result<(), Error> {
    buf.clear();
    // buf holds letters start..done (chunk offsets and ends are multiples of 4)
    let (mut start, mut done) = (0u32, 0u32);
    for (offset, chunk_len) in subject_chunks(len as i32) {
        buf.drain(..(offset as u32 - start) as usize / 4);
        start = offset as u32;
        let want = ((offset + chunk_len) as u32 + 4 * ARENA_PAD as u32).min(len);
        if want > done {
            decode(done..want, buf)?;
            if want == len {
                buf.resize(buf.len() + ARENA_PAD, 0);
            }
            done = want;
        }
        f(buf, offset, chunk_len);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_refget::digest::{AlphabetType, encode_sequence, lookup_alphabet};

    /// Each chunk's bytes are the whole subject's at the chunk's offset, with at least the
    /// `ARENA_PAD` bytes that follow the chunk there: 3-bit subjects with N runs (the random
    /// fill carries over; odd chunks start mid-group) across chunk boundaries, and 2-bit ones.
    #[test]
    fn chunks_are_slices_of_the_whole_subject() {
        let step = 5_000_000 - 100;
        let mut x = 7u64;
        for (len, ns) in [
            (10_000_003usize, true),
            (step + 132, true),
            (step + 133, false),
            (4_999_999, true),
        ] {
            let mut s: Vec<u8> = (0..len)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"ACGT"[(x % 4) as usize]
                })
                .collect();
            if ns {
                for (at, n) in [
                    (10, 5),
                    (step - 37, 60),
                    (step + 2_001, 5_000),
                    (len - 3, 3),
                ] {
                    s[at.min(len)..(at + n).min(len)].fill(b'N');
                }
            }
            let alphabet = if ns {
                AlphabetType::Dna3bit
            } else {
                AlphabetType::Dna2bit
            };
            let packed = encode_sequence(&s, lookup_alphabet(&alphabet)).unwrap();
            let (mut whole, mut amb) = (Vec::new(), Vec::new());
            crate::convert::packed_to_blast(alphabet, &packed, len as u32, &mut whole, &mut amb)
                .unwrap();
            whole.resize(whole.len() + ARENA_PAD, 0);
            let mut seen = Vec::new();
            let mut pieces = Pieces::new(Some(alphabet), packed.len(), len as u32).unwrap();
            let bits = pieces.bits as usize;
            let mut check = |seq: &[u8], offset: i32, clen: i32| {
                let at = offset as usize / 4;
                let end = (offset + clen) as usize / 4 + ARENA_PAD;
                assert!(seq.len() >= end - at, "len {len} chunk at {offset}");
                assert_eq!(
                    seq,
                    &whole[at..at + seq.len()],
                    "len {len} chunk at {offset}"
                );
                seen.push((offset, clen));
            };
            // each range from just the stored bytes it needs, as `for_each_chunk` reads them
            chunks(len as u32, &mut Vec::new(), &mut check, |r, buf| {
                let (b0, b1) = (
                    r.start as usize * bits / 8,
                    (r.end as usize * bits).div_ceil(8),
                );
                pieces.next(&packed[b0..b1], b0, r.end, buf)
            })
            .unwrap();
            assert_eq!(seen, subject_chunks(len as i32), "len {len}");
        }
    }
}
