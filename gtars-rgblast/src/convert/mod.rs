//! Store sequence bytes -> BLAST nucleotide volume bytes.
//!
//! Output matches NCBI WriteDB (writedb_convert.cpp) exactly:
//! * ncbi2na packing (A=0 C=1 G=2 T=3, first base in the high bits),
//!   `len/4 + 1` bytes, the low 2 bits of the last byte = `len % 4`;
//! * ambiguous bases are replaced by bases drawn from NCBI's CRandom seeded
//!   with the sequence length, one draw per ambiguous base, in order;
//! * an ambiguity block (runs of one ncbi4na value, at most 4095 long) follows
//!   unless the sequence is pure ACGT.
//!
//! Inputs are gtars-refget packed records (`Dna2bit`, `Dna3bit`, `DnaIupac`) or
//! ASCII. The code tables are derived from gtars' alphabet definitions at
//! compile time, so a change there breaks the build instead of the output.

pub mod blastna;
pub mod random;

use crate::error::Error;
use gtars_refget::digest::{
    Alphabet, AlphabetType, DNA_2BIT_ALPHABET, DNA_3BIT_ALPHABET, DNA_IUPAC_ALPHABET,
};
pub use random::NcbiRandom;

/// Longest run one ambiguity entry may describe (CAmbigDataBuilder).
pub const MAX_RUN: u32 = 0xFFF;
const NO: u8 = 0xFF;
/// ncbi4na value -> ncbi2na code, or NO for ambiguity codes.
const LOG2: [u8; 16] = [NO, 0, 1, NO, 2, NO, NO, NO, 3, NO, NO, NO, NO, NO, NO, NO];
const HIGH: [u8; 4] = [0, 0xC0, 0xF0, 0xFC];

/// IUPAC letter -> ncbi4na value (A1 C2 M3 G4 R5 S6 V7 T8 U8 W9 Y10 H11 K12 D13 B14 N15,
/// lowercase the same); 0 for any other byte. BLAST has no U, so U is stored as T.
pub const ASCII_TO_NA4: [u8; 256] = {
    let mut t = [0u8; 256];
    let pairs: [(u8, u8); 16] = [
        (b'A', 1),
        (b'C', 2),
        (b'M', 3),
        (b'G', 4),
        (b'R', 5),
        (b'S', 6),
        (b'V', 7),
        (b'T', 8),
        (b'U', 8),
        (b'W', 9),
        (b'Y', 10),
        (b'H', 11),
        (b'K', 12),
        (b'D', 13),
        (b'B', 14),
        (b'N', 15),
    ];
    let mut i = 0;
    while i < pairs.len() {
        t[pairs[i].0 as usize] = pairs[i].1;
        t[pairs[i].0.to_ascii_lowercase() as usize] = pairs[i].1;
        i += 1;
    }
    t
};

/// ncbi4na value -> uppercase IUPAC letter (the reader's inverse; 0 is a gap, written '-').
pub const NA4_TO_IUPAC: [u8; 16] = *b"-ACMGRSVTWYHKDBN";

/// gtars packed code -> ncbi4na, for the first `n` codes of an alphabet.
const fn na4_table(a: &Alphabet, n: usize) -> [u8; 16] {
    let mut t = [0u8; 16];
    let mut code = 0;
    while code < n {
        t[code] = ASCII_TO_NA4[a.decoding_array[code] as usize];
        assert!(t[code] != 0, "gtars alphabet symbol has no ncbi4na value");
        code += 1;
    }
    t
}

/// Dna3bit code -> ncbi4na (gtars order A C G T N R Y U).
pub const NA4_FROM_3: [u8; 16] = na4_table(&DNA_3BIT_ALPHABET, 8);
/// DnaIupac code -> ncbi4na (gtars order N A C M G R S K T W Y U B D H V).
pub const NA4_FROM_IUPAC: [u8; 16] = na4_table(&DNA_IUPAC_ALPHABET, 16);

// `map2` turns gtars' 2-bit layout (T=00 C=01 A=10 G=11) into ncbi2na with bit
// arithmetic, so it depends on that exact order.
const _: () = assert!(
    DNA_2BIT_ALPHABET.decoding_array[0] == b'T'
        && DNA_2BIT_ALPHABET.decoding_array[1] == b'C'
        && DNA_2BIT_ALPHABET.decoding_array[2] == b'A'
        && DNA_2BIT_ALPHABET.decoding_array[3] == b'G'
        && DNA_2BIT_ALPHABET.bits_per_symbol == 2
        && DNA_3BIT_ALPHABET.bits_per_symbol == 3
        && DNA_IUPAC_ALPHABET.bits_per_symbol == 4
);
// The 3-bit all-N fast path matches the word 0x924924 (eight code-4 symbols).
const _: () = assert!(DNA_3BIT_ALPHABET.decoding_array[4] == b'N');
// The IUPAC all-N fast path matches zero bytes (code 0 = N).
const _: () = assert!(DNA_IUPAC_ALPHABET.decoding_array[0] == b'N');

/// Bits per base of a nucleotide alphabet; `None` for alphabets BLAST can't hold.
pub fn nucleotide_bits(a: AlphabetType) -> Option<u64> {
    match a {
        AlphabetType::Dna2bit => Some(2),
        AlphabetType::Dna3bit => Some(3),
        AlphabetType::DnaIupac => Some(4),
        _ => None,
    }
}

/// 12-bit Dna3bit chunk (4 symbols) -> ncbi2na byte, bit 8 set if any is ambiguous.
static LUT3: [u16; 4096] = {
    let mut t = [0u16; 4096];
    let mut x = 0;
    while x < 4096 {
        let mut out = 0u16;
        let mut k = 0;
        while k < 4 {
            let c = (x >> (9 - 3 * k)) & 7;
            let l = LOG2[NA4_FROM_3[c] as usize];
            if l == NO {
                out |= 0x100;
            } else {
                out |= (l as u16) << (6 - 2 * k);
            }
            k += 1;
        }
        t[x] = out;
        x += 1;
    }
    t
};

/// DnaIupac byte (2 symbols) -> 4 bits of ncbi2na, bit 7 set if either is ambiguous.
static LUT4: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut b = 0;
    while b < 256 {
        let l1 = LOG2[NA4_FROM_IUPAC[b >> 4] as usize];
        let l2 = LOG2[NA4_FROM_IUPAC[b & 15] as usize];
        t[b] = if l1 == NO || l2 == NO {
            0x80
        } else {
            (l1 << 2) | l2
        };
        b += 1;
    }
    t
};

/// Base counts a conversion produces as a by-product: ambiguous bases (anything but
/// A, C, G, T, U) and, among them, N.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct BaseCounts {
    pub ambiguous: u64,
    pub n: u64,
}

/// The converted subjects of one batch of consecutive OIDs in one arena: `.nsq` sequence bytes
/// and ambiguity bytes back to back, with offset tables, and the base counts of each subject.
/// Workers reuse these arenas from batch to batch.
#[derive(Debug, Default)]
pub struct ConvertedBatch {
    /// The batch's first subject (an OID, or an index into the subjects being written).
    pub first_oid: u32,
    pub seq: Vec<u8>,
    /// n + 1 offsets into `seq`.
    pub seq_offsets: Vec<u32>,
    pub amb: Vec<u8>,
    /// n + 1 offsets into `amb`.
    pub amb_offsets: Vec<u32>,
    pub counts: Vec<BaseCounts>,
}

impl ConvertedBatch {
    /// Empty the arena (keeping its capacity) for a batch starting at `first_oid`.
    pub fn reset(&mut self, first_oid: u32) {
        self.first_oid = first_oid;
        self.seq.clear();
        self.amb.clear();
        self.counts.clear();
        self.seq_offsets.clear();
        self.seq_offsets.push(0);
        self.amb_offsets.clear();
        self.amb_offsets.push(0);
    }

    /// Append the next subject.
    pub fn push(&mut self, seq: &[u8], amb: &[u8], counts: BaseCounts) -> Result<(), Error> {
        self.seq.extend_from_slice(seq);
        self.amb.extend_from_slice(amb);
        let off = |n: usize| {
            u32::try_from(n).map_err(|_| {
                Error::InvalidOptions(
                    "a batch's converted bytes exceed 4 GiB; lower --batch-mb".into(),
                )
            })
        };
        self.seq_offsets.push(off(self.seq.len())?);
        self.amb_offsets.push(off(self.amb.len())?);
        self.counts.push(counts);
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// Sequence bytes, ambiguity bytes and base counts of the batch's `i`-th subject.
    pub fn subject(&self, i: usize) -> (&[u8], &[u8], BaseCounts) {
        let (s0, s1) = (
            self.seq_offsets[i] as usize,
            self.seq_offsets[i + 1] as usize,
        );
        let (a0, a1) = (
            self.amb_offsets[i] as usize,
            self.amb_offsets[i + 1] as usize,
        );
        (&self.seq[s0..s1], &self.amb[a0..a1], self.counts[i])
    }
}

/// Ambiguity accumulator: CAmbigDataBuilder.
pub(crate) struct Amb {
    regions: Vec<(u8, u32, u32)>, // (ncbi4na value, start, length)
    rng: Option<NcbiRandom>,
    seed: u32,
    pub(crate) count: u64,
    pub(crate) n: u64,
}

impl Amb {
    fn new(len: u32) -> Self {
        Amb {
            regions: Vec::new(),
            rng: None,
            seed: len,
            count: 0,
            n: 0,
        }
    }
    #[inline(always)]
    fn rand(&mut self) -> u32 {
        if self.rng.is_none() {
            self.rng = Some(NcbiRandom::new(self.seed));
        }
        self.rng.as_mut().unwrap().next()
    }
    /// Record ambiguous value `v` at `off`; return the random base to store.
    #[inline(always)]
    fn add(&mut self, v: u8, off: u32) -> u8 {
        self.count += 1;
        self.n += (v == 15) as u64;
        match self.regions.last_mut() {
            Some(r) if r.0 == v && r.1 + r.2 == off && r.2 < MAX_RUN => r.2 += 1,
            _ => self.regions.push((v, off, 1)),
        }
        if v == 15 || v == 0 {
            return (self.rand() & 3) as u8;
        }
        let bc = (v as u32).count_ones();
        let mut pick = self.rand() % bc;
        for i in 0..4u8 {
            if v & (1 << i) == 0 {
                continue;
            }
            if pick == 0 {
                return i;
            }
            pick -= 1;
        }
        unreachable!()
    }
    /// `n` consecutive N (value 15) starting at `off`; returns packed codes, first in high bits.
    #[inline(always)]
    fn add_n(&mut self, off: u32, n: u32) -> u32 {
        let fast = matches!(self.regions.last(), Some(r) if r.0 == 15 && r.1 + r.2 == off && r.2 + n <= MAX_RUN);
        let mut out = 0u32;
        if fast {
            self.regions.last_mut().unwrap().2 += n;
            self.count += n as u64;
            self.n += n as u64;
            for _ in 0..n {
                out = (out << 2) | (self.rand() & 3);
            }
        } else {
            for k in 0..n {
                out = (out << 2) | self.add(15, off + k) as u32;
            }
        }
        out
    }
    fn encode(&self, len: u32, out: &mut Vec<u8>) {
        if self.regions.is_empty() {
            return;
        }
        let newf = len > 0xFF_FFFF || self.regions.iter().any(|r| r.2 > 0xF);
        let n = self.regions.len() as u32;
        let words = if newf { 0x8000_0000 | (n * 2) } else { n };
        out.reserve(4 + self.regions.len() * if newf { 8 } else { 4 });
        out.extend_from_slice(&words.to_be_bytes());
        for &(v, start, l) in &self.regions {
            let lm1 = l - 1;
            if newf {
                out.extend_from_slice(&[(v << 4) | (lm1 >> 8) as u8, (lm1 & 0xFF) as u8, 0, 0]);
                out.extend_from_slice(&start.to_be_bytes());
            } else {
                out.extend_from_slice(&[
                    (v << 4) | lm1 as u8,
                    (start >> 16) as u8,
                    (start >> 8) as u8,
                    start as u8,
                ]);
            }
        }
    }
}

/// Letter `p` of a piece that starts at letter `lo` (a multiple of 4) at `seq[0]`.
#[inline(always)]
fn slow(seq: &mut [u8], a: &mut Amb, p: usize, lo: usize, na4: u8) {
    let l = LOG2[na4 as usize];
    let c = if l != NO { l } else { a.add(na4, p as u32) };
    seq[(p - lo) >> 2] |= c << (6 - 2 * (p & 3));
}

/// Dna2bit (T C A G) -> ncbi2na (A C G T) on whole bytes: field (h,l) -> (!(h^l), !h).
#[inline(always)]
fn map2(x: u64) -> u64 {
    let h = x & 0xAAAA_AAAA_AAAA_AAAA;
    let l = x & 0x5555_5555_5555_5555;
    (!(h ^ (l << 1)) & 0xAAAA_AAAA_AAAA_AAAA) | ((!h >> 1) & 0x5555_5555_5555_5555)
}

// The conv functions write letters `lo..hi` (`lo` a multiple of 4) to `seq` from `seq[0]`,
// reading a record whose byte `base` is `packed[0]`.
fn conv2(packed: &[u8], base: usize, lo: usize, hi: usize, seq: &mut [u8]) {
    let full = (hi - lo) / 4;
    let packed = &packed[lo / 4 - base..];
    let (src, dst) = (&packed[..full], &mut seq[..full]);
    let mut sc = src.chunks_exact(8);
    let mut dc = dst.chunks_exact_mut(8);
    for (s, d) in (&mut sc).zip(&mut dc) {
        let x = u64::from_ne_bytes(s.try_into().unwrap());
        d.copy_from_slice(&map2(x).to_ne_bytes());
    }
    for (s, d) in sc.remainder().iter().zip(dc.into_remainder()) {
        *d = map2(*s as u64) as u8;
    }
    let r = hi & 3;
    if r > 0 {
        seq[full] = map2(packed[full] as u64) as u8 & HIGH[r];
    }
}

fn conv3(packed: &[u8], base: usize, lo: usize, hi: usize, seq: &mut [u8], a: &mut Amb) {
    let code = |p: usize| {
        let bit = 3 * p;
        let i = bit / 8 - base;
        let w = (packed[i] as usize) << 8 | *packed.get(i + 1).unwrap_or(&0) as usize;
        NA4_FROM_3[(w >> (13 - bit % 8)) & 7]
    };
    // whole 8-letter groups (3 bytes); letters before the first and after the last one by
    // one, in order (the random fill draws in letter order)
    let (g0, g1) = (lo.div_ceil(8), (hi / 8).max(lo.div_ceil(8)));
    for p in lo..(8 * g0).min(hi) {
        slow(seq, a, p, lo, code(p));
    }
    // (a piece inside one group has no whole group, and its bytes may end before 3 * g0)
    let (src, dst): (&[u8], &mut [u8]) = if g0 < g1 {
        let d = &mut seq[(8 * g0 - lo) / 4..(8 * g1 - lo) / 4];
        (&packed[3 * g0 - base..3 * g1 - base], d)
    } else {
        (&[], &mut [])
    };
    for (g, (b, d)) in (g0..).zip(src.chunks_exact(3).zip(dst.chunks_exact_mut(2))) {
        let w = (b[0] as u32) << 16 | (b[1] as u32) << 8 | b[2] as u32;
        let e1 = LUT3[(w >> 12) as usize];
        let e2 = LUT3[(w & 0xFFF) as usize];
        if (e1 | e2) & 0x100 == 0 {
            d[0] = e1 as u8;
            d[1] = e2 as u8;
        } else if w == 0x92_4924 {
            let n = a.add_n((8 * g) as u32, 8);
            d[0] = (n >> 8) as u8;
            d[1] = n as u8;
        } else {
            for k in 0..8 {
                let c = ((w >> (21 - 3 * k)) & 7) as usize;
                slow(d, a, 8 * g + k, 8 * g, NA4_FROM_3[c]);
            }
        }
    }
    for p in (8 * g1).max(lo)..hi {
        slow(seq, a, p, lo, code(p));
    }
}

fn conv4(packed: &[u8], base: usize, lo: usize, hi: usize, seq: &mut [u8], a: &mut Amb) {
    for i in lo / 4..hi / 4 {
        let (b0, b1) = (packed[2 * i - base], packed[2 * i + 1 - base]);
        let (e0, e1) = (LUT4[b0 as usize], LUT4[b1 as usize]);
        let o = i - lo / 4;
        if (e0 | e1) & 0x80 == 0 {
            seq[o] = (e0 << 4) | e1;
        } else if b0 == 0 && b1 == 0 {
            seq[o] = a.add_n((4 * i) as u32, 4) as u8;
        } else {
            let w = (b0 as usize) << 8 | b1 as usize;
            for k in 0..4 {
                slow(
                    seq,
                    a,
                    4 * i + k,
                    lo,
                    NA4_FROM_IUPAC[(w >> (12 - 4 * k)) & 15],
                );
            }
        }
    }
    for p in (hi & !3).max(lo)..hi {
        let b = packed[p / 2 - base];
        let c = if p % 2 == 0 { b >> 4 } else { b & 15 };
        slow(seq, a, p, lo, NA4_FROM_IUPAC[c as usize]);
    }
}

/// One record converted to `.nsq` sequence bytes piece by piece, in order: the random fill
/// carries over from piece to piece, so the pieces join into the bytes of a whole conversion
/// (the search decodes a long subject one chunk at a time this way).
pub struct Pieces {
    /// The packed alphabet, or `None` for ASCII letters.
    alphabet: Option<AlphabetType>,
    /// Bits per letter: 2, 3 or 4 packed, 8 for ASCII.
    pub bits: u64,
    len: u32,
    pos: u32,
    amb: Amb,
}

impl Pieces {
    /// Check that a record of `stored` bytes (packed codes, or ASCII letters when `alphabet` is
    /// `None`) holds `len` letters.
    pub fn new(alphabet: Option<AlphabetType>, stored: usize, len: u32) -> Result<Pieces, Error> {
        let bits = match alphabet {
            Some(a) => nucleotide_bits(a).ok_or_else(|| {
                Error::sequence("", format!("alphabet {a} is not a nucleotide alphabet"))
            })?,
            None => 8,
        };
        let need = (len as u64 * bits).div_ceil(8);
        if stored as u64 != need {
            let of = alphabet.map_or(String::new(), |a| format!(" for {len} bp of {a}"));
            return Err(Error::sequence(
                "",
                format!("stored record is {stored} bytes, expected {need}{of}"),
            ));
        }
        Ok(Pieces {
            alphabet,
            bits,
            len,
            pos: 0,
            amb: Amb::new(len),
        })
    }

    /// Append the bytes of the next letters, up to `hi` (a multiple of 4, or the length, where
    /// the last byte also gets the length code `len % 4`), to `seq`. `src` holds the record
    /// from byte `base` on, at least through the bytes of letter `hi - 1`.
    pub fn next(
        &mut self,
        src: &[u8],
        base: usize,
        hi: u32,
        seq: &mut Vec<u8>,
    ) -> Result<(), Error> {
        let (lo, h) = (self.pos as usize, hi as usize);
        let at = seq.len();
        seq.resize(at + (h - lo).div_ceil(4), 0);
        let out = &mut seq[at..];
        let a = &mut self.amb;
        match self.alphabet {
            Some(AlphabetType::Dna2bit) => conv2(src, base, lo, h, out),
            Some(AlphabetType::Dna3bit) => conv3(src, base, lo, h, out, a),
            Some(_) => conv4(src, base, lo, h, out, a),
            None => {
                for p in lo..h {
                    let b = src[p - base];
                    let na4 = ASCII_TO_NA4[b as usize];
                    if na4 == 0 {
                        return Err(Error::sequence(
                            "",
                            format!(
                                "byte {:?} (0x{b:02x}) at position {p} is not a nucleotide letter",
                                b as char
                            ),
                        ));
                    }
                    slow(out, a, p, lo, na4);
                }
            }
        }
        if hi == self.len {
            if hi.is_multiple_of(4) {
                seq.push(0);
            }
            *seq.last_mut().unwrap() |= (hi & 3) as u8;
        }
        self.pos = hi;
        Ok(())
    }

    /// The whole record into `seq` and its ambiguity table into `amb`.
    pub fn whole(
        mut self,
        src: &[u8],
        seq: &mut Vec<u8>,
        amb: &mut Vec<u8>,
    ) -> Result<BaseCounts, Error> {
        seq.clear();
        self.next(src, 0, self.len, seq)?;
        amb.clear();
        self.amb.encode(self.len, amb);
        Ok(BaseCounts {
            ambiguous: self.amb.count,
            n: self.amb.n,
        })
    }
}

/// Packed gtars record -> `.nsq` sequence bytes and ambiguity bytes; returns the ambiguous
/// and N base counts. The record size must be exactly `ceil(len * bits / 8)`.
pub fn packed_to_blast(
    alphabet: AlphabetType,
    packed: &[u8],
    len: u32,
    seq: &mut Vec<u8>,
    amb: &mut Vec<u8>,
) -> Result<BaseCounts, Error> {
    Pieces::new(Some(alphabet), packed.len(), len)?.whole(packed, seq, amb)
}

/// ASCII (IUPAC letters, either case) -> the same output as [`packed_to_blast`].
/// Any byte that is not an IUPAC nucleotide letter is an error naming its position.
pub fn ascii_to_blast(
    ascii: &[u8],
    seq: &mut Vec<u8>,
    amb: &mut Vec<u8>,
) -> Result<BaseCounts, Error> {
    let len: u32 = ascii
        .len()
        .try_into()
        .map_err(|_| Error::sequence("", format!("{} bp is too long for BLAST", ascii.len())))?;
    Pieces::new(None, ascii.len(), len)?.whole(ascii, seq, amb)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtars_refget::digest::{encode_sequence, lookup_alphabet};

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    fn blast(ascii: &[u8]) -> (Vec<u8>, Vec<u8>, u64) {
        let (mut s, mut a) = (Vec::new(), Vec::new());
        let n = ascii_to_blast(ascii, &mut s, &mut a).unwrap();
        (s, a, n.ambiguous)
    }

    #[test]
    fn map2_matches_table() {
        let m = [3u8, 1, 0, 2];
        for b in 0..256u32 {
            let mut e = 0u8;
            for k in 0..4 {
                e |= m[((b >> (6 - 2 * k)) & 3) as usize] << (6 - 2 * k);
            }
            assert_eq!(map2(b as u64) as u8, e, "byte {b:08b}");
        }
    }

    #[test]
    fn ncbi_random_first_outputs() {
        let first = |seed: u32| {
            let mut r = NcbiRandom::new(seed);
            (0..16).map(|_| r.next()).collect::<Vec<u32>>()
        };
        assert_eq!(first(0), PINNED_SEED_0);
        assert_eq!(first(1), PINNED_SEED_1);
        assert_eq!(first(248_956_422), PINNED_SEED_CHR1);
    }

    // Pinned from this implementation, which was verified byte-identical to
    // makeblastdb 2.17 on the whole jungle store (11.6 G ambiguous bases).
    // 248,956,422 is the length of hg38 chr1.
    const PINNED_SEED_0: [u32; 16] = [
        268829939, 1366391869, 263636910, 1101374817, 229910017, 957137640, 1502329388, 615052060,
        1217726492, 1546006862, 838525884, 1503548984, 523556860, 1988938449, 532144258,
        1444840934,
    ];
    const PINNED_SEED_1: [u32; 16] = [
        73902710, 1005518751, 421776705, 756398104, 1668674573, 167634074, 942066300, 2040517292,
        244848787, 416061467, 968830145, 2016625440, 759469809, 1706259361, 1285316677, 1581700453,
    ];
    const PINNED_SEED_CHR1: [u32; 16] = [
        1564515080, 209779078, 1424096803, 340671656, 1661969740, 93272340, 44825612, 1438355580,
        158470627, 552413212, 145218269, 1927016613, 672321466, 2102461108, 1467633169, 448175328,
    ];

    #[test]
    fn ambiguity_old_format() {
        // short sequence, runs <= 15: 4-byte entries, plain word count
        let (_, a, n) = blast(b"ACGTNNNACGTRACGT");
        assert_eq!(n, 4);
        assert_eq!(&a[..4], &2u32.to_be_bytes());
        assert_eq!(a.len(), 4 + 2 * 4);
        assert_eq!(a[4], (15 << 4) | 2); // N, length 3
        assert_eq!(&a[5..8], &[0, 0, 4]);
        assert_eq!(a[8], 5 << 4); // R, length 1
        assert_eq!(&a[9..12], &[0, 0, 11]);
    }

    #[test]
    fn ambiguity_new_format_long_run() {
        let mut s = b"ACGT".to_vec();
        s.extend(std::iter::repeat_n(b'N', 16));
        s.extend_from_slice(b"ACGT");
        let (_, a, n) = blast(&s);
        assert_eq!(n, 16);
        assert_eq!(
            u32::from_be_bytes(a[..4].try_into().unwrap()),
            0x8000_0000 | 2
        );
        assert_eq!(a.len(), 4 + 8);
        assert_eq!(a[4], 15 << 4);
        assert_eq!(a[5], 15); // length - 1
        assert_eq!(&a[8..12], &4u32.to_be_bytes());
    }

    #[test]
    fn ambiguity_new_format_long_sequence() {
        // longer than 0xFFFFFF bases: new format even with short runs
        let len = 16_777_300usize;
        let mut s = vec![b'A'; len];
        s[len - 10] = b'N';
        let (seq, a, n) = blast(&s);
        assert_eq!(n, 1);
        assert_eq!(seq.len(), len / 4 + 1);
        assert_eq!(
            u32::from_be_bytes(a[..4].try_into().unwrap()),
            0x8000_0000 | 2
        );
        assert_eq!(a[4], 15 << 4);
        assert_eq!(a[5], 0);
        assert_eq!(&a[8..12], &((len - 10) as u32).to_be_bytes());
    }

    #[test]
    fn max_run_splits_at_4095() {
        for (run, entries) in [
            (4095usize, 1usize),
            (4096, 2),
            (4097, 2),
            (8190, 2),
            (8191, 3),
        ] {
            let mut s = b"AC".to_vec();
            s.extend(std::iter::repeat_n(b'N', run));
            let (_, a, n) = blast(&s);
            assert_eq!(n as usize, run);
            let words = u32::from_be_bytes(a[..4].try_into().unwrap());
            assert_eq!(words, 0x8000_0000 | (2 * entries as u32), "run {run}");
            // first entry is a full 4095 run (or the whole run)
            let lm1 = ((a[4] & 15) as u32) << 8 | a[5] as u32;
            assert_eq!(lm1 + 1, (run as u32).min(MAX_RUN));
        }
    }

    #[test]
    fn ascii_path_equals_packed_path() {
        let mut rng = Rng(0x1234_5678_9abc_def1);
        let alphabets: [(AlphabetType, &[u8]); 3] = [
            (AlphabetType::Dna2bit, b"ACGT"),
            (AlphabetType::Dna3bit, b"ACGTNRYU"),
            (AlphabetType::DnaIupac, b"ACGTNMRSKWYUBDHV"),
        ];
        for (at, symbols) in alphabets {
            for trial in 0..60 {
                let len = rng.below(3000) as usize + trial;
                let mut s = Vec::with_capacity(len);
                while s.len() < len {
                    if symbols.contains(&b'N') && rng.below(20) == 0 {
                        let run = rng.below(40) as usize + 1;
                        s.extend(std::iter::repeat_n(b'N', run));
                    } else {
                        s.push(symbols[rng.below(symbols.len() as u64) as usize]);
                    }
                }
                s.truncate(len);
                let packed = encode_sequence(&s, lookup_alphabet(&at)).unwrap();
                let (mut s1, mut a1) = (Vec::new(), Vec::new());
                let n1 = packed_to_blast(at, &packed, len as u32, &mut s1, &mut a1).unwrap();
                let (mut s2, mut a2) = (Vec::new(), Vec::new());
                let n2 = ascii_to_blast(&s, &mut s2, &mut a2).unwrap();
                assert_eq!((s1, a1, n1), (s2, a2, n2), "{at} trial {trial} len {len}");
                let n_count = s.iter().filter(|&&b| b == b'N').count() as u64;
                assert_eq!(n1.n, n_count, "{at} trial {trial}: N count");
                // lowercase gives the same result
                let lower: Vec<u8> = s.iter().map(|b| b.to_ascii_lowercase()).collect();
                assert_eq!(blast(&lower), blast(&s));
            }
        }
    }

    /// Pieces cut at random multiples of 4 (so 3-bit pieces start mid-group) join into the
    /// whole conversion: the random fill carries over from piece to piece.
    #[test]
    fn pieces_join_into_the_whole_record() {
        let mut rng = Rng(0x0dd5_eed5_0fc0_ffee);
        let alphabets: [(AlphabetType, &[u8]); 3] = [
            (AlphabetType::Dna2bit, b"ACGT"),
            (AlphabetType::Dna3bit, b"ACGTNRYU"),
            (AlphabetType::DnaIupac, b"ACGTNMRSKWYUBDHV"),
        ];
        for (at, symbols) in alphabets {
            for trial in 0..40 {
                let len = rng.below(2000) as usize + trial;
                let s: Vec<u8> = (0..len)
                    .map(|_| match rng.below(4) {
                        0 if symbols.contains(&b'N') => b'N',
                        _ => symbols[rng.below(symbols.len() as u64) as usize],
                    })
                    .collect();
                let packed = encode_sequence(&s, lookup_alphabet(&at)).unwrap();
                let (mut whole, mut a) = (Vec::new(), Vec::new());
                packed_to_blast(at, &packed, len as u32, &mut whole, &mut a).unwrap();
                // each piece from just the bytes it needs, as a long record is read
                for (alphabet, src) in [(Some(at), &packed), (None, &s)] {
                    let mut p = Pieces::new(alphabet, src.len(), len as u32).unwrap();
                    let bits = p.bits as usize;
                    let mut joined = Vec::new();
                    let mut at_letter = 0;
                    loop {
                        let next = (at_letter + 4 * rng.below(40) as usize).min(len);
                        let (b0, b1) = (at_letter * bits / 8, (next * bits).div_ceil(8));
                        p.next(&src[b0..b1], b0, next as u32, &mut joined).unwrap();
                        at_letter = next;
                        if at_letter == len {
                            break;
                        }
                    }
                    assert_eq!(joined, whole, "{at} trial {trial} len {len}");
                }
            }
        }
    }

    #[test]
    fn tables_derived_from_gtars() {
        assert_eq!(&NA4_FROM_3[..8], &[1, 2, 4, 8, 15, 5, 10, 8]);
        assert_eq!(
            NA4_FROM_IUPAC,
            [15, 1, 2, 3, 4, 5, 6, 12, 8, 9, 10, 8, 14, 13, 11, 7]
        );
        for (code, &v) in NA4_FROM_IUPAC.iter().enumerate() {
            let letter = DNA_IUPAC_ALPHABET.decoding_array[code];
            let back = NA4_TO_IUPAC[v as usize];
            assert!(
                back == letter || (letter == b'U' && back == b'T'),
                "code {code}"
            );
        }
    }

    #[test]
    fn packed_size_is_checked() {
        let (mut s, mut a) = (Vec::new(), Vec::new());
        let e = packed_to_blast(AlphabetType::Dna2bit, &[0u8; 3], 13, &mut s, &mut a).unwrap_err();
        assert!(e.to_string().contains("expected 4"), "{e}");
        let e = ascii_to_blast(b"ACGX", &mut s, &mut a).unwrap_err();
        assert!(e.to_string().contains("position 3"), "{e}");
    }
}
