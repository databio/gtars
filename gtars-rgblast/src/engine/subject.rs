//! The subjects a search reads (`BlastSeqSrc` in BLAST+, see
//! `algo/blast/api/seqsrc_seqdb.cpp`): the foundation's `SubjectPlan` over the store's stored
//! records, each a `.seq` file or a pack span (`store::access::MappedFiles`).
//!
//! `for_each_chunk` decodes a subject to ncbi2na (makeblastdb's layout with the random fill
//! for ambiguous bases, `convert::Pieces`) one chunk at a time, padded so the scanners may
//! over-read four bytes past a chunk as the C does; `DecodedSubjects` keeps each subject
//! decoded whole for the later query batches of a search. `blastna_window` decodes letters of one
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
use std::sync::{Arc, OnceLock};

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

/// The plan + mapped records the engine searches.
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
        let ChunkBuf { seq, stored } = buf;
        let len = self.plan.subjects[oid as usize].length;
        self.with_decoder(oid, stored, |decode| chunks(len, seq, &mut f, decode))
    }

    /// The whole subject in ncbi2na, `ARENA_PAD` zeros after it: the bytes `for_each_chunk`
    /// hands out, in one buffer (`DecodedSubjects`). Decoded in pieces of `WHOLE_STEP`
    /// letters, so a large packed record is read piece by piece.
    pub fn decode_whole(&self, oid: u32, stored: &mut Vec<u8>) -> Result<Box<[u8]>, Error> {
        const WHOLE_STEP: u32 = 1 << 24;
        let len = self.plan.subjects[oid as usize].length;
        self.with_decoder(oid, stored, |decode| {
            let mut out = Vec::with_capacity(len as usize / 4 + 1 + ARENA_PAD);
            let mut at = 0;
            while at < len {
                let hi = len.min(at + WHOLE_STEP);
                decode(at..hi, &mut out)?;
                at = hi;
            }
            out.resize(out.len() + ARENA_PAD, 0);
            Ok(out.into_boxed_slice())
        })
    }

    /// Run `body` with the subject's decoder: `decode(range, buf)` appends the ncbi2na bytes
    /// of letters `range` (consecutive ranges from 0) to `buf`. A large packed or raw record
    /// is read piece by piece (`MappedFiles::read_at`, into `stored`); a small or zstd one
    /// whole.
    fn with_decoder<R>(
        &self,
        oid: u32,
        stored: &mut Vec<u8>,
        body: impl FnOnce(
            &mut dyn FnMut(Range<u32>, &mut Vec<u8>) -> Result<(), Error>,
        ) -> Result<R, Error>,
    ) -> Result<R, Error> {
        let s = &self.plan.subjects[oid as usize];
        let run = || match (self.mode, s.stored_len) {
            (StorageMode::Encoded | StorageMode::Raw, Some(n)) if n >= MappedFiles::MMAP_MIN => {
                let alphabet = (self.mode == StorageMode::Encoded).then_some(s.alphabet);
                let mut pieces = Pieces::new(alphabet, n as usize, s.length)?;
                let bits = pieces.bits;
                let reader = self.access.reader(oid)?;
                body(&mut |r, seq| {
                    let b = r.start as u64 * bits / 8;
                    stored.resize(((r.end as u64 * bits).div_ceil(8) - b) as usize, 0);
                    self.access.read_at(oid, &reader, b, stored)?;
                    pieces.next(stored, b as usize, r.end, seq)
                })
            }
            _ => {
                let bytes = self.access.bytes(oid)?;
                let (mut pieces, src) =
                    decode::to_blast_pieces(self.mode, s.alphabet, bytes.as_slice(), s.length)?;
                body(&mut |r, seq| pieces.next(&src, 0, r.end, seq))
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

/// The subjects of a search with several query batches, each decoded once (on first use, by
/// the thread that searches it) and read by every later batch: the bytes
/// `StoreSubjects::decode_whole` gives, so the chunks are those `for_each_chunk` hands out.
pub struct DecodedSubjects {
    seqs: Vec<OnceLock<Box<[u8]>>>,
}

impl DecodedSubjects {
    /// Room for every subject of `subjects`; nothing decoded yet.
    pub fn new(subjects: &StoreSubjects) -> DecodedSubjects {
        DecodedSubjects {
            seqs: (0..subjects.plan.subjects.len())
                .map(|_| OnceLock::new())
                .collect(),
        }
    }

    /// The bytes all subjects take decoded: `decode_whole`'s length summed.
    pub fn bytes_needed(subjects: &StoreSubjects) -> u64 {
        let per = |l: u32| l as u64 / 4 + 1 + ARENA_PAD as u64;
        subjects.plan.subjects.iter().map(|s| per(s.length)).sum()
    }

    /// `StoreSubjects::for_each_chunk` on the decoded subject (decoding it first if no batch
    /// has).
    pub fn for_each_chunk(
        &self,
        subjects: &StoreSubjects,
        oid: u32,
        buf: &mut ChunkBuf,
        mut f: impl FnMut(&[u8], i32, i32),
    ) -> Result<(), Error> {
        let cell = &self.seqs[oid as usize];
        let whole = match cell.get() {
            Some(w) => w,
            None => {
                let w = subjects.decode_whole(oid, &mut buf.stored)?;
                cell.get_or_init(|| w)
            }
        };
        whole_chunks(subjects.plan.subjects[oid as usize].length, whole, &mut f);
        Ok(())
    }
}

/// `chunks` on a subject decoded whole (`StoreSubjects::decode_whole`): `f` gets the same
/// slices, the chunk through the letters `chunks` would have decoded for it.
fn whole_chunks(len: u32, whole: &[u8], f: &mut impl FnMut(&[u8], i32, i32)) {
    for (offset, chunk_len) in subject_chunks(len as i32) {
        let want = ((offset + chunk_len) as u32 + 4 * ARENA_PAD as u32).min(len);
        let end = if want == len {
            whole.len()
        } else {
            want as usize / 4
        };
        f(&whole[offset as usize / 4..end], offset, chunk_len);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_refget::digest::{AlphabetType, encode_sequence, lookup_alphabet};

    /// Each chunk's bytes are the whole subject's at the chunk's offset, with at least the
    /// `ARENA_PAD` bytes that follow the chunk there, and `whole_chunks` hands out the same
    /// slices: 3-bit subjects with N runs (the random fill carries over; odd chunks start
    /// mid-group) across chunk boundaries, and 2-bit ones.
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
            let (mut seen, mut slices) = (Vec::new(), Vec::new());
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
                slices.push(seq.to_vec());
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
            // a subject decoded whole (`DecodedSubjects`) hands out the very same slices
            let mut k = 0;
            whole_chunks(len as u32, &whole, &mut |seq, offset, clen| {
                assert_eq!((offset, clen), seen[k], "len {len}");
                assert!(seq == slices[k].as_slice(), "len {len} chunk at {offset}");
                k += 1;
            });
            assert_eq!(k, seen.len(), "len {len}");
        }
    }
}
