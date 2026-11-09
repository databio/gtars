//! Nucleotide lookup tables.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/core/blast_nalookup.c` (`BlastChooseNaLookupTable`,
//! `BlastMBLookupTableNew`, `s_FillContigMBTable`, `s_SeqLocListInvert`,
//! `BlastSmallNaLookupTableNew`, `s_BlastSmallNaLookupFinalize`, `BlastNaLookupTableNew`,
//! `s_BlastNaLookupFinalize`), `blast_lookup.c` (`BlastLookupIndexQueryExactMatches`,
//! `BlastLookupAddWordHit`), `lookup_wrap.c` (`LookupTableWrapInit` with its SmallNa to NA
//! fallback, `GetOffsetArraySize`), `lookup_util.c` (`ilog2`, `EstimateNumTableEntries`) and
//! `na_ungapped.c` (`s_MBLookup`, `s_SmallNaLookup`, `s_NaLookup`, `BlastChooseNaExtend`).
//!
//! `s_RemovePolyAWords` runs only for the mapping program (magicblast), not blastn, so it is
//! not ported. The NA (thick-backbone) table is used when `BlastChooseNaLookupTable` says so
//! (query offsets past 32767 or too many entries for 15-bit offsets) or when the SmallNa
//! table's overflow array would reach 32768 cells.

use super::query::QueryBlock;
use crate::error::Error;

pub const BITS_PER_NUC: u32 = 2;
pub const BLAST2NA_MASK: u8 = 0xfc;
pub const OFFSET_ARRAY_SIZE: usize = 4096;
const PV_ARRAY_BTS: u32 = 5;
const PV_ARRAY_MASK: u64 = 31;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LutType {
    Mb,
    SmallNa,
    Na,
}

/// Which `extend_callback` `BlastChooseNaExtend` picks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Extender {
    /// `s_BlastNaExtendDirect`
    Direct,
    /// `s_BlastNaExtend`
    Na,
    /// `s_BlastNaExtendAligned`
    NaAligned,
    /// `s_BlastSmallNaExtend`
    SmallNa,
    /// `s_BlastSmallNaExtendAlignedOneByte`
    SmallNaAlignedOneByte,
}

/// `BlastMBLookupTable` (contiguous words).
#[derive(Debug, Clone)]
pub struct MbTable {
    pub word_length: i32,
    pub lut_word_length: i32,
    pub hashsize: u64,
    pub scan_step: i32,
    /// `hashtable[word]` = 1 + offset of the word's start in `sequence`, or 0.
    pub hashtable: Vec<i32>,
    /// `next_pos[q_off + 1]`: next occurrence (same +1 convention), 0 ends the chain.
    pub next_pos: Vec<i32>,
    pub pv_array: Vec<u32>,
    pub pv_array_bts: u32,
    pub longest_chain: u32,
}

impl MbTable {
    #[inline(always)]
    pub fn has_hits(&self, index: u64) -> bool {
        (self.pv_array[(index >> self.pv_array_bts) as usize] & (1u32 << (index & PV_ARRAY_MASK)))
            != 0
    }

    /// `s_MBLookup`: is query position `q_pos` (word start) in the chain for `index`?
    #[inline]
    pub fn lookup(&self, index: u64, q_pos: i32) -> bool {
        let index = index & (self.hashsize - 1);
        let q_pos = q_pos + 1;
        if !self.has_hits(index) {
            return false;
        }
        let mut q_off = self.hashtable[index as usize];
        while q_off != 0 {
            if q_off == q_pos {
                return true;
            }
            q_off = self.next_pos[q_off as usize];
        }
        false
    }
}

/// `BlastSmallNaLookupTable`.
#[derive(Debug, Clone)]
pub struct SmallNaTable {
    pub mask: i32,
    pub word_length: i32,
    pub lut_word_length: i32,
    pub scan_step: i32,
    pub backbone_size: i32,
    pub longest_chain: u32,
    /// -1 empty; >= 0 one offset; < -1: -(index into overflow).
    pub final_backbone: Vec<i16>,
    pub overflow: Vec<i16>,
}

impl SmallNaTable {
    /// `s_SmallNaLookup`.
    #[inline]
    pub fn lookup(&self, index: i32, q_pos: i32) -> bool {
        let mut index = self.final_backbone[(index & self.mask) as usize] as i32;
        if index == q_pos {
            return true;
        }
        if index == -1 || index >= 0 {
            return false;
        }
        let mut src_off = (-index) as usize;
        index = self.overflow[src_off] as i32;
        src_off += 1;
        loop {
            if index == q_pos {
                return true;
            }
            index = self.overflow[src_off] as i32;
            src_off += 1;
            if index < 0 {
                return false;
            }
        }
    }
}

/// `NA_HITS_PER_CELL`: offsets that fit in a thick-backbone cell before spilling over.
pub const NA_HITS_PER_CELL: usize = 3;

/// `BlastNaLookupTable`: thick backbone (`num_used` + up to three offsets per cell, or an
/// overflow cursor), overflow array and presence vector.
#[derive(Debug, Clone)]
pub struct NaTable {
    pub mask: i32,
    pub word_length: i32,
    pub lut_word_length: i32,
    pub scan_step: i32,
    pub backbone_size: i32,
    pub longest_chain: u32,
    /// `thick_backbone[i].num_used`
    pub num_used: Vec<i32>,
    /// `payload.entries` when `num_used <= NA_HITS_PER_CELL`, else `payload.overflow_cursor`
    /// in `[0]`.
    pub payload: Vec<[i32; NA_HITS_PER_CELL]>,
    pub overflow: Vec<i32>,
    pub pv: Vec<u32>,
}

impl NaTable {
    /// `PV_TEST`.
    #[inline(always)]
    pub fn has_hits(&self, index: i32) -> bool {
        (self.pv[(index >> PV_ARRAY_BTS) as usize] & (1u32 << (index as u64 & PV_ARRAY_MASK))) != 0
    }

    /// `s_BlastLookupGetNumHits`.
    #[inline(always)]
    pub fn num_hits(&self, index: i32) -> i32 {
        if self.has_hits(index) {
            self.num_used[index as usize]
        } else {
            0
        }
    }

    /// The query offsets stored for `index` (`s_BlastLookupRetrieve`'s `lookup_pos`).
    #[inline]
    pub fn positions(&self, index: i32) -> &[i32] {
        let n = self.num_used[index as usize] as usize;
        if n <= NA_HITS_PER_CELL {
            &self.payload[index as usize][..n]
        } else {
            let c = self.payload[index as usize][0] as usize;
            &self.overflow[c..c + n]
        }
    }

    /// `s_NaLookup`.
    #[inline]
    pub fn lookup(&self, index: i32, q_pos: i32) -> bool {
        let index = index & self.mask;
        if !self.has_hits(index) {
            return false;
        }
        self.positions(index).contains(&q_pos)
    }
}

#[derive(Debug, Clone)]
pub enum Table {
    Mb(MbTable),
    SmallNa(SmallNaTable),
    Na(NaTable),
}

/// `LookupTableWrap` for nucleotide searches.
#[derive(Debug, Clone)]
pub struct Lookup {
    pub table: Table,
    pub lut_type: LutType,
    pub word_length: i32,
    pub lut_word_length: i32,
    pub scan_step: i32,
    pub extender: Extender,
    /// `masked_locations`: the inverse of the lookup segments (regions the full word may not
    /// cross with soft masking), `None` when `word_length == lut_word_length` or no masking.
    pub masked_locations: Option<Vec<(i32, i32)>>,
    pub longest_chain: u32,
}

/// `BlastChooseNaLookupTable`: (type, width).
pub fn choose_na_lookup_table(
    word_size: i32,
    approx_table_entries: i32,
    max_q_off: i32,
) -> (LutType, i32) {
    let e = approx_table_entries;
    let (mut t, w) = match word_size {
        4..=6 => (LutType::SmallNa, word_size),
        7 => (LutType::SmallNa, if e < 250 { 6 } else { 7 }),
        8 => (LutType::SmallNa, if e < 8500 { 7 } else { 8 }),
        9 => {
            if e < 1250 {
                (LutType::SmallNa, 7)
            } else if e < 21000 {
                (LutType::SmallNa, 8)
            } else {
                (LutType::Mb, 9)
            }
        }
        10 => {
            if e < 1250 {
                (LutType::SmallNa, 7)
            } else if e < 8500 {
                (LutType::SmallNa, 8)
            } else if e < 18000 {
                (LutType::Mb, 9)
            } else {
                (LutType::Mb, 10)
            }
        }
        11 => {
            if e < 12000 {
                (LutType::SmallNa, 8)
            } else if e < 180000 {
                (LutType::Mb, 10)
            } else {
                (LutType::Mb, 11)
            }
        }
        12 => {
            if e < 8500 {
                (LutType::SmallNa, 8)
            } else if e < 18000 {
                (LutType::Mb, 9)
            } else if e < 60000 {
                (LutType::Mb, 10)
            } else if e < 900000 {
                (LutType::Mb, 11)
            } else {
                (LutType::Mb, 12)
            }
        }
        _ => {
            if e < 8500 {
                (LutType::SmallNa, 8)
            } else if e < 300000 {
                (LutType::Mb, 11)
            } else {
                (LutType::Mb, 12)
            }
        }
    };
    if t == LutType::SmallNa && (e >= 32767 || max_q_off >= 32768) {
        t = LutType::Na;
    }
    (t, w)
}

pub fn ilog2(mut x: u64) -> u32 {
    let mut lg = 0;
    if x == 0 {
        return 0;
    }
    loop {
        x >>= 1;
        if x == 0 {
            break;
        }
        lg += 1;
    }
    lg
}

/// `s_SeqLocListInvert(locations, length)`.
fn invert_locations(locations: &[(i32, i32)], length: i32) -> Vec<(i32, i32)> {
    let mut out = Vec::new();
    if locations.is_empty() {
        return out;
    }
    let start = 0;
    let stop = (locations[0].0 - 1).max(0);
    if stop - start > 2 {
        out.push((start, stop));
    }
    for (i, loc) in locations.iter().enumerate() {
        let start = loc.1 + 1;
        let stop = if i + 1 < locations.len() {
            locations[i + 1].0 - 1
        } else {
            length - 1
        };
        out.push((start, stop));
    }
    out
}

pub(crate) fn build_mb(
    qb: &QueryBlock,
    segments: &[(i32, i32)],
    word_size: i32,
    lut_width: i32,
    approx_table_entries: i32,
) -> MbTable {
    let k_target_pv_size: i64 = 131072;
    let k_small_query_cutoff = 15000;
    let k_large_query_cutoff = 800000;
    let hashsize: u64 = 1u64 << (BITS_PER_NUC * lut_width as u32);
    let mut pv_size: i64 = if lut_width <= 12 {
        if (hashsize as i64) <= 8 * k_target_pv_size {
            (hashsize >> PV_ARRAY_BTS) as i64
        } else {
            k_target_pv_size / 4
        }
    } else {
        k_target_pv_size * 64 / 4
    };
    if approx_table_entries <= k_small_query_cutoff || approx_table_entries >= k_large_query_cutoff
    {
        pv_size /= 2;
    }
    let pv_array_bts = ilog2(hashsize / pv_size as u64);
    let mut t = MbTable {
        word_length: word_size,
        lut_word_length: lut_width,
        hashsize,
        scan_step: word_size - lut_width + 1,
        hashtable: vec![0; hashsize as usize],
        next_pos: vec![0; qb.length as usize + 1],
        pv_array: vec![0u32; pv_size as usize],
        pv_array_bts,
        longest_chain: 2,
    };
    // s_FillContigMBTable
    let k_compression_factor: usize = 2048;
    let mut helper = vec![0u32; (hashsize as usize) / k_compression_factor];
    let lut_mask = hashsize - 1;
    let seq_start = &qb.buf; // sequence_start
    for &(left, right) in segments {
        let from = left;
        let to = right - lut_width;
        if word_size > right - left + 1 {
            continue;
        }
        let mut ecode: u64 = 0;
        // seq points at sequence_start + from; pos = seq + lut_width
        let mut seq = from as isize; // index into seq_start
        let mut pos = seq + lut_width as isize;
        let from = from - (lut_width - 2);
        let last_offset = to + 2;
        let mut index = from;
        while index <= last_offset {
            seq += 1;
            let val = seq_start[seq as usize];
            if val & BLAST2NA_MASK != 0 {
                ecode = 0;
                pos = seq + lut_width as isize;
                index += 1;
                continue;
            }
            ecode = ((ecode << BITS_PER_NUC) & lut_mask) + val as u64;
            if seq < pos {
                index += 1;
                continue;
            }
            let e = ecode as usize;
            if t.hashtable[e] == 0 {
                t.pv_array[e >> pv_array_bts] |= 1u32 << (ecode & PV_ARRAY_MASK);
            } else {
                helper[e / k_compression_factor] += 1;
            }
            t.next_pos[index as usize] = t.hashtable[e];
            t.hashtable[e] = index;
            index += 1;
        }
    }
    for &h in &helper {
        t.longest_chain = t.longest_chain.max(h);
    }
    t
}

/// `BlastLookupIndexQueryExactMatches` into a thin backbone of offset lists.
fn index_exact_matches(
    qb: &QueryBlock,
    segments: &[(i32, i32)],
    word_length: i32,
    lut_word_length: i32,
) -> Vec<Vec<i32>> {
    let backbone_size = 1usize << (BITS_PER_NUC * lut_word_length as u32);
    let mut thin: Vec<Vec<i32>> = vec![Vec::new(); backbone_size];
    let sequence = qb.sequence();
    let add = |thin: &mut Vec<Vec<i32>>, word: &[u8], q_off: i32| {
        let mut idx = 0usize;
        for &w in word {
            idx = (idx << BITS_PER_NUC) | w as usize;
        }
        thin[idx].push(q_off);
    };
    for &(from, to) in segments {
        if word_length > to - from + 1 {
            continue;
        }
        let mut seq = from; // index into sequence
        let mut word_target = seq + lut_word_length;
        let mut offset = from;
        while offset <= to {
            if seq >= word_target {
                let s = (seq - lut_word_length) as usize;
                add(
                    &mut thin,
                    &sequence[s..s + lut_word_length as usize],
                    offset - lut_word_length,
                );
            }
            if sequence[seq as usize] & BLAST2NA_MASK != 0 {
                word_target = seq + lut_word_length + 1;
            }
            offset += 1;
            seq += 1;
        }
        if seq >= word_target {
            let s = (seq - lut_word_length) as usize;
            add(
                &mut thin,
                &sequence[s..s + lut_word_length as usize],
                offset - lut_word_length,
            );
        }
    }
    thin
}

/// `BlastSmallNaLookupTableNew` + `s_BlastSmallNaLookupFinalize`; `None` when the overflow
/// array would exceed 32768 cells (the C then falls back to the NA table).
pub(crate) fn build_small_na(
    qb: &mut QueryBlock,
    segments: &[(i32, i32)],
    word_size: i32,
    lut_width: i32,
) -> Option<SmallNaTable> {
    let thin = index_exact_matches(qb, segments, word_size, lut_width);
    let backbone_size = thin.len();
    let mut overflow_cells_needed = 2usize;
    let mut longest_chain = 0usize;
    for chain in &thin {
        let n = chain.len();
        if n > 1 {
            overflow_cells_needed += n + 1;
        }
        longest_chain = longest_chain.max(n);
    }
    if overflow_cells_needed >= 32768 {
        return None;
    }
    qb.compress();
    let mut final_backbone = vec![-1i16; backbone_size];
    let mut overflow = vec![0i16; overflow_cells_needed];
    let mut cursor = 2usize;
    for (i, chain) in thin.iter().enumerate() {
        match chain.len() {
            0 => {}
            1 => final_backbone[i] = chain[0] as i16,
            _ => {
                final_backbone[i] = -(cursor as i16);
                for &q in chain {
                    overflow[cursor] = q as i16;
                    cursor += 1;
                }
                overflow[cursor] = -1;
                cursor += 1;
            }
        }
    }
    Some(SmallNaTable {
        mask: backbone_size as i32 - 1,
        word_length: word_size,
        lut_word_length: lut_width,
        scan_step: word_size - lut_width + 1,
        backbone_size: backbone_size as i32,
        longest_chain: longest_chain as u32,
        final_backbone,
        overflow,
    })
}

/// `BlastNaLookupTableNew` + `s_BlastNaLookupFinalize`.
pub(crate) fn build_na(
    qb: &QueryBlock,
    segments: &[(i32, i32)],
    word_size: i32,
    lut_width: i32,
) -> NaTable {
    let thin = index_exact_matches(qb, segments, word_size, lut_width);
    let backbone_size = thin.len();
    let mut overflow_cells_needed = 0usize;
    let mut longest_chain = 0usize;
    for chain in &thin {
        let n = chain.len();
        if n > NA_HITS_PER_CELL {
            overflow_cells_needed += n;
        }
        longest_chain = longest_chain.max(n);
    }
    let mut t = NaTable {
        mask: backbone_size as i32 - 1,
        word_length: word_size,
        lut_word_length: lut_width,
        scan_step: word_size - lut_width + 1,
        backbone_size: backbone_size as i32,
        longest_chain: longest_chain as u32,
        num_used: vec![0; backbone_size],
        payload: vec![[0; NA_HITS_PER_CELL]; backbone_size],
        overflow: Vec::with_capacity(overflow_cells_needed),
        pv: vec![0u32; (backbone_size >> PV_ARRAY_BTS) + 1],
    };
    for (i, chain) in thin.iter().enumerate() {
        if chain.is_empty() {
            continue;
        }
        t.num_used[i] = chain.len() as i32;
        t.pv[i >> PV_ARRAY_BTS] |= 1u32 << (i as u64 & PV_ARRAY_MASK);
        if chain.len() <= NA_HITS_PER_CELL {
            for (j, &q) in chain.iter().enumerate() {
                t.payload[i][j] = q;
            }
        } else {
            t.payload[i][0] = t.overflow.len() as i32;
            t.overflow.extend_from_slice(chain);
        }
    }
    t
}

impl Lookup {
    fn from_na(t: NaTable, lut_type: LutType, masked_locations: Option<Vec<(i32, i32)>>) -> Lookup {
        let extender = if t.lut_word_length == t.word_length {
            Extender::Direct
        } else if t.lut_word_length % 4 == 0 && t.scan_step % 4 == 0 {
            Extender::NaAligned
        } else {
            Extender::Na
        };
        Lookup {
            lut_type,
            word_length: t.word_length,
            lut_word_length: t.lut_word_length,
            scan_step: t.scan_step,
            extender,
            masked_locations,
            longest_chain: t.longest_chain,
            table: Table::Na(t),
        }
    }

    /// `LookupTableWrapInit` for `eMBLookupTable`/`eNaLookupTable` requests, then
    /// `BlastChooseNaExtend`. `mask_at_hash` decides whether `masked_locations` is kept.
    pub fn build(qb: &mut QueryBlock, word_size: i32, mask_at_hash: bool) -> Result<Lookup, Error> {
        let segments = qb.lookup_segments();
        let (entries, max_q_off) = QueryBlock::estimate_table_entries(&segments);
        let (lut_type, width) = choose_na_lookup_table(word_size, entries, max_q_off);
        let masked_locations = if word_size > width && mask_at_hash && !segments.is_empty() {
            Some(invert_locations(&segments, qb.length))
        } else {
            None
        };
        match lut_type {
            LutType::Mb => {
                let t = build_mb(qb, &segments, word_size, width, entries);
                let extender = if t.lut_word_length == t.word_length {
                    Extender::Direct
                } else if t.lut_word_length % 4 == 0 && t.scan_step % 4 == 0 {
                    Extender::NaAligned
                } else {
                    Extender::Na
                };
                Ok(Lookup {
                    lut_type,
                    word_length: t.word_length,
                    lut_word_length: t.lut_word_length,
                    scan_step: t.scan_step,
                    extender,
                    masked_locations,
                    longest_chain: t.longest_chain,
                    table: Table::Mb(t),
                })
            }
            LutType::SmallNa => {
                let Some(t) = build_small_na(qb, &segments, word_size, width) else {
                    // the overflow array would need 32768 cells: the C falls back to the
                    // NA table of the same width
                    let t = build_na(qb, &segments, word_size, width);
                    return Ok(Lookup::from_na(t, LutType::Na, masked_locations));
                };
                let extender = if t.lut_word_length == t.word_length {
                    Extender::Direct
                } else if t.lut_word_length % 4 == 0
                    && t.scan_step % 4 == 0
                    && t.word_length - t.lut_word_length <= 4
                {
                    Extender::SmallNaAlignedOneByte
                } else {
                    Extender::SmallNa
                };
                // masked_locations for SmallNa does not depend on mask_at_hash in the C.
                let ml = if word_size > width && !segments.is_empty() {
                    Some(invert_locations(&segments, qb.length))
                } else {
                    None
                };
                Ok(Lookup {
                    lut_type,
                    word_length: t.word_length,
                    lut_word_length: t.lut_word_length,
                    scan_step: t.scan_step,
                    extender,
                    masked_locations: ml,
                    longest_chain: t.longest_chain,
                    table: Table::SmallNa(t),
                })
            }
            LutType::Na => {
                let t = build_na(qb, &segments, word_size, width);
                Ok(Lookup::from_na(t, LutType::Na, masked_locations))
            }
        }
    }

    /// `GetOffsetArraySize`.
    pub fn offset_array_size(&self) -> usize {
        OFFSET_ARRAY_SIZE + self.longest_chain as usize
    }

    /// The `lookup_callback`: is `q_pos` one of the offsets stored for word `index`?
    #[inline]
    pub fn contains(&self, index: u64, q_pos: i32) -> bool {
        match &self.table {
            Table::Mb(t) => t.lookup(index, q_pos),
            Table::SmallNa(t) => t.lookup(index as i32, q_pos),
            Table::Na(t) => t.lookup(index as i32, q_pos),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::options::{BlastOptions, Strand};
    use crate::engine::query::QuerySet;

    #[test]
    fn chooser() {
        assert_eq!(
            choose_na_lookup_table(28, 2000, 1000),
            (LutType::SmallNa, 8)
        );
        assert_eq!(choose_na_lookup_table(28, 20000, 20000), (LutType::Mb, 11));
        assert_eq!(
            choose_na_lookup_table(28, 400000, 400000),
            (LutType::Mb, 12)
        );
        assert_eq!(
            choose_na_lookup_table(11, 2000, 1000),
            (LutType::SmallNa, 8)
        );
        assert_eq!(choose_na_lookup_table(11, 20000, 20000), (LutType::Mb, 10));
        assert_eq!(choose_na_lookup_table(28, 2000, 40000), (LutType::Na, 8));
        assert_eq!(ilog2(4096), 12);
    }

    #[test]
    fn small_na_contents() {
        // ACGTACGTAC with word 4: lut width 4 == word -> Direct; ACGT at offsets 0 and 4.
        let qs = QuerySet::from_fasta_bytes(b">a\nACGTACGTAC\n").unwrap();
        let mut o = BlastOptions::megablast();
        o.dust = None;
        o.strand = Strand::Plus;
        o.word_size = 4;
        let mut qb = QueryBlock::build(&qs, &o);
        let lk = Lookup::build(&mut qb, 4, true).unwrap();
        assert_eq!(lk.lut_type, LutType::SmallNa);
        assert_eq!(lk.extender, Extender::Direct);
        let Table::SmallNa(t) = &lk.table else {
            panic!()
        };
        let acgt = 0b00_01_10_11usize;
        assert!(t.final_backbone[acgt] < -1, "chain goes to overflow");
        assert!(lk.contains(acgt as u64, 0));
        assert!(lk.contains(acgt as u64, 4));
        assert!(!lk.contains(acgt as u64, 1));
        let cgta = 0b01_10_11_00usize;
        assert!(t.final_backbone[cgta] < -1);
        assert!(lk.contains(cgta as u64, 1) && lk.contains(cgta as u64, 5));
        assert!(qb.compressed.is_some());
        // word 28 on this short query: width 8, SmallNa extender
        let mut qb = QueryBlock::build(&qs, &o);
        let lk = Lookup::build(&mut qb, 28, true).unwrap();
        assert_eq!((lk.lut_word_length, lk.extender), (8, Extender::SmallNa));
    }

    #[test]
    fn na_table_when_offsets_exceed_15_bits() {
        // 40 kb of low complexity (dust masks it) followed by a 300-base random slice: few
        // words but query offsets past 32768, so BlastChooseNaLookupTable picks eNaLookupTable
        let mut fa = b">a\n".to_vec();
        for _ in 0..20000 {
            fa.extend_from_slice(b"AT");
        }
        let mut x = 12345u64;
        for _ in 0..300 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            fa.push(b"ACGT"[((x >> 33) % 4) as usize]);
        }
        let qs = QuerySet::from_fasta_bytes(&fa).unwrap();
        let o = BlastOptions::blastn();
        let mut qb = QueryBlock::build(&qs, &o);
        let lk = Lookup::build(&mut qb, 11, true).unwrap();
        assert_eq!(
            (lk.lut_type, lk.lut_word_length, lk.scan_step),
            (LutType::Na, 8, 4)
        );
        assert_eq!(lk.extender, Extender::NaAligned);
        let Table::Na(t) = &lk.table else { panic!() };
        // every indexed word is found again through the callback
        let seq = qb.sequence();
        let segs = qb.lookup_segments();
        let mut n = 0;
        for &(from, to) in &segs {
            let mut off = from;
            while off + 8 <= to + 1 {
                let w = &seq[off as usize..off as usize + 8];
                if w.iter().any(|&c| c > 3) {
                    off += 1;
                    continue;
                }
                let idx = w.iter().fold(0i32, |a, &c| (a << 2) | c as i32);
                assert!(lk.contains(idx as u64, off), "word at {off} missing");
                n += 1;
                off += 1;
            }
        }
        assert!(n > 500, "{n} words indexed");
        assert!(t.longest_chain >= 1 && !t.overflow.is_empty() || t.longest_chain <= 3);
        assert!(!lk.contains(0, 1_000_000));
        // the SmallNa table is chosen for the short query in small_na_contents; the fallback
        // for overflow needs >= 32768 overflow cells and is exercised by blastn-short with
        // many queries in the golden tests
    }

    #[test]
    fn mb_chains() {
        // enough letters for the MB table: 20 kb of repeating ACGT -> entries ~ 40000 (both strands)
        let s: Vec<u8> = (0..20000).map(|i| b"ACGT"[i % 4]).collect();
        let mut fa = b">a\n".to_vec();
        fa.extend_from_slice(&s);
        let qs = QuerySet::from_fasta_bytes(&fa).unwrap();
        let mut o = BlastOptions::megablast();
        o.dust = None;
        let mut qb = QueryBlock::build(&qs, &o);
        let lk = Lookup::build(&mut qb, 28, true).unwrap();
        assert_eq!(
            (lk.lut_type, lk.lut_word_length, lk.scan_step),
            (LutType::Mb, 11, 18)
        );
        assert_eq!(lk.extender, Extender::Na);
        let Table::Mb(t) = &lk.table else { panic!() };
        // word ACGTACGTACG (11 letters) occurs at offsets 0,4,8,... on the plus strand
        let mut w = 0u64;
        for &c in b"ACGTACGTACG" {
            w = (w << 2) | super::super::encoding::ASCII_TO_BLASTNA[c as usize] as u64;
        }
        assert!(t.has_hits(w));
        assert!(lk.contains(w, 0));
        assert!(lk.contains(w, 4));
        assert!(!lk.contains(w, 1));
        assert!(t.longest_chain > 2);
        assert!(lk.masked_locations.is_some());
    }
}
