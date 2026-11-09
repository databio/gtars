//! Word finder: scan pieces of a subject, verify seeds to full words, track diagonals and
//! run the ungapped extension.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/core/na_ungapped.c` (`BlastNaWordFinder`,
//! `s_DetermineScanningOffsets` (masksubj.inl), `s_BlastNaExtendDirect`, `s_BlastNaExtend`,
//! `s_BlastNaExtendAligned`, `s_BlastSmallNaExtend`, `s_BlastSmallNaExtendAlignedOneByte`,
//! `s_TypeOfWord`, `s_IsSeedMasked`, `s_BlastnDiagTableExtendInitialHit`,
//! `s_BlastnDiagHashExtendInitialHit`, `s_BlastDiagHashRetrieve/Insert`,
//! `s_NuclUngappedExtend`, `s_NuclUngappedExtendExact`), `blast_extend.c`
//! (`s_BlastDiagTableNew`, `BlastExtendWordNew`, `Blast_ExtendWordExit`, `s_BlastDiagClear`,
//! `BLAST_SaveInitialHit`, `score_compare_match`, `Blast_InitHitListSortByScore`) and
//! `blast_query_info.c::BSearchContextInfo`.

use super::lookup::{Extender, Lookup, Table};
use super::options::{DiagContainer, Parameters};
use super::query::QueryBlock;
use super::scan::{OffsetPair, mb_scan_any, na_scan_8_4, na_scan_any, small_na_scan_any};
use std::cmp::Ordering;

const DIAGHASH_NUM_BUCKETS: u32 = 512;
const DIAGHASH_CHAIN_LENGTH: usize = 256;

/// One subject (or chunk of one) in ncbi2na: `seq[0]` holds the first four bases of the
/// chunk; `seq_ranges` are `[left, right)` in chunk coordinates (`right` = length when the
/// subject is unmasked, as `s_BackupSubject`/`s_GetNextSubjectChunk` set them).
#[derive(Debug, Clone, Copy)]
pub struct SubjectBlk<'a> {
    pub seq: &'a [u8],
    pub length: i32,
    pub seq_ranges: &'a [(i32, i32)],
}

#[inline(always)]
pub fn unpack_base(byte: u8, base: i32) -> usize {
    ((byte >> (2 * base)) & 3) as usize
}

/// `BlastUngappedData`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct UngappedData {
    pub q_start: i32,
    pub s_start: i32,
    pub length: i32,
    pub score: i32,
}

/// `BlastInitHSP` (ungapped data is always present for blastn).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InitHsp {
    pub q_off: i32,
    pub s_off: i32,
    pub ungapped: UngappedData,
}

/// `score_compare_match`.
pub fn score_compare_match(a: &InitHsp, b: &InitHsp) -> Ordering {
    let (u1, u2) = (&a.ungapped, &b.ungapped);
    u2.score
        .cmp(&u1.score)
        .then(u1.s_start.cmp(&u2.s_start))
        .then(u2.length.cmp(&u1.length))
        .then(u1.q_start.cmp(&u2.q_start))
        .then(u2.length.cmp(&u1.length))
}

/// `BLAST_DiagTable` for one-hit searches (`window_size` 0, every remaining task): only
/// `DiagStruct::last_hit` is ever read, so the flags and the `hit_len_array` are not kept.
#[derive(Debug, Clone)]
pub struct DiagTable {
    diag_array_length: i32,
    diag_mask: i32,
    offset: i32,
    last_hit: Vec<i32>,
}

/// A `BLAST_DiagHash` cell, without the two-hit fields (`hit_len`, `hit_saved`).
#[derive(Debug, Clone, Copy, Default)]
struct DiagHashCell {
    diag: i32,
    level: i32,
    next: u32,
}

/// `BLAST_DiagHash`.
#[derive(Debug, Clone)]
pub struct DiagHash {
    occupancy: usize,
    backbone: Vec<u32>,
    chain: Vec<DiagHashCell>,
    offset: i32,
}

/// `Blast_ExtendWord`.
#[derive(Debug, Clone)]
pub enum ExtendWord {
    Table(DiagTable),
    Hash(DiagHash),
}

impl ExtendWord {
    /// `BlastExtendWordNew` with `window_size` 0.
    pub fn new(query_length: i32, params: &Parameters) -> ExtendWord {
        match params.container_type {
            DiagContainer::DiagHash => ExtendWord::Hash(DiagHash {
                occupancy: 1,
                backbone: vec![0; DIAGHASH_NUM_BUCKETS as usize],
                chain: vec![DiagHashCell::default(); DIAGHASH_CHAIN_LENGTH],
                offset: 0,
            }),
            DiagContainer::DiagArray => {
                let mut diag_array_length = 1;
                while diag_array_length < query_length {
                    diag_array_length <<= 1;
                }
                ExtendWord::Table(DiagTable {
                    diag_array_length,
                    diag_mask: diag_array_length - 1,
                    offset: 0,
                    last_hit: vec![0; diag_array_length as usize],
                })
            }
        }
    }

    /// `Blast_ExtendWordExit` (window 0).
    pub fn exit(&mut self, subject_length: i32) {
        match self {
            ExtendWord::Table(d) => {
                if d.offset >= i32::MAX / 4 {
                    d.offset = 0;
                    d.last_hit.fill(0);
                } else {
                    d.offset += subject_length;
                }
            }
            ExtendWord::Hash(h) => {
                if h.offset >= i32::MAX / 4 {
                    h.occupancy = 1;
                    h.offset = 0;
                    h.backbone.fill(0);
                } else {
                    h.offset += subject_length;
                }
            }
        }
    }
}

impl DiagHash {
    #[inline]
    fn bucket(diag: i32) -> usize {
        ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS) as usize
    }

    fn retrieve(&self, diag: i32) -> Option<i32> {
        let mut index = self.backbone[Self::bucket(diag)];
        while index != 0 {
            let c = &self.chain[index as usize];
            if c.diag == diag {
                return Some(c.level);
            }
            index = c.next;
        }
        None
    }

    fn insert(&mut self, diag: i32, level: i32, s_off: i32, window_size: i32) {
        let bucket = Self::bucket(diag);
        let mut index = self.backbone[bucket];
        while index != 0 {
            let c = &mut self.chain[index as usize];
            if c.diag == diag {
                c.level = level;
                return;
            } else if s_off - c.level > window_size {
                c.diag = diag;
                c.level = level;
                return;
            }
            index = c.next;
        }
        if self.occupancy == self.chain.len() {
            self.chain
                .resize(self.chain.len() * 2, DiagHashCell::default());
        }
        let cell = DiagHashCell {
            diag,
            level,
            next: self.backbone[bucket],
        };
        self.chain[self.occupancy] = cell;
        self.backbone[bucket] = self.occupancy as u32;
        self.occupancy += 1;
    }
}

/// `BSearchContextInfo`.
pub fn bsearch_context(qb: &QueryBlock, n: i32) -> usize {
    let size = qb.contexts.len() as i32;
    let (mut b, mut e);
    if qb.min_length > 0 && qb.max_length > 0 && qb.first_context == 0 {
        b = (n / (qb.max_length + 1)).min(size - 1);
        e = (n / (qb.min_length + 1) + 1).min(size);
    } else {
        b = 0;
        e = size;
    }
    while b < e - 1 {
        let m = (b + e) / 2;
        if qb.contexts[m as usize].query_offset > n {
            e = m;
        } else {
            b = m;
        }
    }
    b as usize
}

/// Per-thread state for the word finder.
pub struct WordFinderScratch {
    pub offset_pairs: Vec<OffsetPair>,
    pub ewp: ExtendWord,
    pub init_hitlist: Vec<InitHsp>,
}

impl WordFinderScratch {
    pub fn new(qb: &QueryBlock, lookup: &Lookup, params: &Parameters) -> Self {
        WordFinderScratch {
            offset_pairs: vec![OffsetPair::default(); lookup.offset_array_size()],
            ewp: ExtendWord::new(qb.length, params),
            init_hitlist: Vec::with_capacity(1024),
        }
    }
}

/// Everything the extension needs that is constant across subjects.
pub struct WordFinderContext<'a> {
    pub qb: &'a QueryBlock,
    pub lookup: &'a Lookup,
    pub matrix: &'a [[i32; 16]; 16],
    pub params: &'a Parameters,
}

/// `s_NuclUngappedExtendExact`.
fn ungapped_extend_exact(
    ctx: &WordFinderContext,
    subject: &SubjectBlk,
    q_off: i32,
    s_off: i32,
    x: i32,
) -> UngappedData {
    let query = ctx.qb.sequence();
    let subj = subject.seq;
    let matrix = ctx.matrix;
    let mut base: i32 = 3 - (s_off % 4);
    let q_avail = ctx.qb.length - q_off;
    let s_avail = subject.length - s_off;
    let mut q = q_off as isize; // index into query (sequence); q-1 reaches the sentinel at buf[0]
    let mut q_beg = q;
    let mut q_end = q;
    let mut s = (s_off / 4) as isize;
    let (start, remainder): (isize, i32) = if q_off < s_off {
        (((s_off - q_off) / 4) as isize, 3 - ((s_off - q_off) % 4))
    } else {
        (0, 3)
    };
    let mut score = 0i32;
    let mut sum = 0i32;
    let buf = &ctx.qb.buf; // sequence_start; query[i] == buf[i+1]
    while s > start || (s == start && base < remainder) {
        if base == 3 {
            s -= 1;
            base = 0;
        } else {
            base += 1;
        }
        let ch = subj[s as usize];
        q -= 1;
        let qc = buf[(q + 1) as usize] as usize;
        sum += matrix[qc][unpack_base(ch, base)];
        if sum > 0 {
            q_beg = q;
            score += sum;
            sum = 0;
        } else if sum < x {
            break;
        }
    }
    let q_start = q_beg as i32;
    let s_start = s_off - (q_off - q_start);
    let (sf, remainder): (isize, i32) = if q_avail < s_avail {
        (
            ((s_off + q_avail) / 4) as isize,
            3 - ((s_off + q_avail) % 4),
        )
    } else {
        ((subject.length / 4) as isize, 3 - (subject.length % 4))
    };
    let mut q = q_off as isize;
    let mut s = (s_off / 4) as isize;
    sum = 0;
    base = 3 - (s_off % 4);
    let mut x_current = x;
    while s < sf || (s == sf && base > remainder) {
        let ch = subj[s as usize];
        let qc = query[q as usize] as usize;
        q += 1;
        sum += matrix[qc][unpack_base(ch, base)];
        if sum > 0 {
            q_end = q;
            score += sum;
            x_current = if -score > x { -score } else { x };
            sum = 0;
        } else if sum < x_current {
            break;
        }
        if base == 0 {
            base = 3;
            s += 1;
        } else {
            base -= 1;
        }
    }
    UngappedData {
        q_start,
        s_start,
        length: (q_end - q_beg) as i32,
        score,
    }
}

/// `s_NuclUngappedExtend`.
#[allow(clippy::too_many_arguments)]
fn ungapped_extend(
    ctx: &WordFinderContext,
    subject: &SubjectBlk,
    q_off: i32,
    s_match_end: i32,
    s_off: i32,
    x: i32,
    reduced_cutoff: i32,
) -> UngappedData {
    let query = ctx.qb.sequence();
    let subj = subject.seq;
    let table = &ctx.params.nucl_score_table;
    let mut len = (4 - (s_off % 4)) % 4;
    let q_ext = q_off + len;
    let s_ext = s_off + len;
    let mut q = q_ext as isize;
    let mut s = (s_ext / 4) as isize;
    len = q_ext.min(s_ext) / 4;
    let mut score = 0i32;
    let mut sum = 0i32;
    let mut new_q = q;
    for _ in 0..len {
        let s_byte = subj[(s - 1) as usize];
        // C: `Uint1 q_byte = (q[-4] << 6) | (q[-3] << 4) | (q[-2] << 2) | q[-1]`, OR'd as
        // ints then truncated; the terms overlap when a query letter is an ambiguity code.
        let q_byte = (((query[(q - 4) as usize] as i32) << 6)
            | ((query[(q - 3) as usize] as i32) << 4)
            | ((query[(q - 2) as usize] as i32) << 2)
            | (query[(q - 1) as usize] as i32)) as u8;
        sum += table[(q_byte ^ s_byte) as usize];
        if sum > 0 {
            new_q = q - 4;
            score += sum;
            sum = 0;
        }
        if sum < x {
            break;
        }
        s -= 1;
        q -= 4;
    }
    let q_start = new_q as i32;
    let s_start = s_ext - (q_ext - q_start);
    let mut q = q_ext as isize;
    let mut s = (s_ext / 4) as isize;
    len = (ctx.qb.length - q_ext).min(subject.length - s_ext) / 4;
    sum = 0;
    new_q = q;
    for _ in 0..len {
        let s_byte = subj[s as usize];
        let q_byte = (((query[q as usize] as i32) << 6)
            | ((query[(q + 1) as usize] as i32) << 4)
            | ((query[(q + 2) as usize] as i32) << 2)
            | (query[(q + 3) as usize] as i32)) as u8;
        sum += table[(q_byte ^ s_byte) as usize];
        if sum > 0 {
            new_q = q + 3;
            score += sum;
            sum = 0;
        }
        if sum < x {
            break;
        }
        s += 1;
        q += 4;
    }
    if score >= reduced_cutoff {
        ungapped_extend_exact(ctx, subject, q_off, s_off, x)
    } else {
        UngappedData {
            q_start,
            s_start,
            score,
            length: (s_match_end - s_start).max(new_q as i32 - q_start + 1),
        }
    }
}

/// `s_IsSeedMasked`: true when the lut word at `s_off` is NOT in the table for `q_pos`.
#[inline]
fn is_seed_masked(
    lookup: &Lookup,
    subject: &SubjectBlk,
    s_off: i32,
    lut_word_length: i32,
    q_pos: i32,
) -> bool {
    let s = &subject.seq[(s_off / 4) as usize..];
    let shift = 2 * (16 - s_off % 4 - lut_word_length);
    let index: u64 = match shift {
        8 | 10 | 12 | 14 => {
            ((s[0] as u32) << 24 | (s[1] as u32) << 16 | (s[2] as u32) << 8) as u64 >> shift
        }
        16 | 18 | 20 | 22 => ((s[0] as u32) << 24 | (s[1] as u32) << 16) as u64 >> shift,
        24 => s[0] as u64,
        _ => {
            ((s[0] as u32) << 24 | (s[1] as u32) << 16 | (s[2] as u32) << 8 | s[3] as u32) as u64
                >> shift
        }
    };
    !lookup.contains(index, q_pos)
}

/// `s_TypeOfWord` with `check_double` FALSE (one-hit): false rejects the word; adjusts
/// `q_off`/`s_off` and `extended`.
#[allow(clippy::too_many_arguments)]
fn type_of_word(
    ctx: &WordFinderContext,
    subject: &SubjectBlk,
    q_off: &mut i32,
    s_off: &mut i32,
    has_locations: bool,
    s_range: i32,
    word_length: i32,
    lut_word_length: i32,
    extended: &mut i32,
) -> bool {
    *extended = 0;
    if word_length == lut_word_length {
        return true;
    }
    let mut q_end = *q_off + word_length;
    let mut s_end = *s_off + word_length;
    let context = bsearch_context(ctx.qb, q_end);
    let c = &ctx.qb.contexts[context];
    let q_range = c.query_offset + c.query_length;
    let lookup = ctx.lookup;
    if has_locations {
        if is_seed_masked(
            lookup,
            subject,
            s_end - lut_word_length,
            lut_word_length,
            q_end - lut_word_length,
        ) {
            return false;
        }
        loop {
            if !is_seed_masked(lookup, subject, *s_off, lut_word_length, *q_off) {
                break;
            }
            *s_off += 1;
            *q_off += 1;
        }
    }
    let ext_to = word_length - (q_end - *q_off);
    let ext_max = (q_range - q_end).min(s_range - s_end);
    if ext_to != 0 || has_locations {
        if ext_to > ext_max {
            return false;
        }
        q_end += ext_to;
        s_end += ext_to;
        let mut s_pos = s_end - lut_word_length;
        let mut q_pos = q_end - lut_word_length;
        while s_pos > *s_off {
            if is_seed_masked(lookup, subject, s_pos, lut_word_length, q_pos) {
                return false;
            }
            s_pos -= lut_word_length;
            q_pos -= lut_word_length;
        }
        *extended = ext_to;
    }
    true
}

/// `s_BlastnDiagTableExtendInitialHit` / `s_BlastnDiagHashExtendInitialHit` for one-hit
/// searches (`window_size` 0; one body, the containers differ only in how `last_hit` is
/// stored). Returns hit_ready.
#[allow(clippy::too_many_arguments)]
fn extend_initial_hit(
    ctx: &WordFinderContext,
    subject: &SubjectBlk,
    mut q_off: i32,
    mut s_off: i32,
    has_locations: bool,
    s_range: i32,
    word_length: i32,
    lut_word_length: i32,
    ewp: &mut ExtendWord,
    init_hitlist: &mut Vec<InitHsp>,
    check_masks: bool,
) -> i32 {
    let mut s_end = s_off + word_length;
    let (table_offset, last_hit, real_diag, diag) = match ewp {
        ExtendWord::Table(d) => {
            let diag = s_off + d.diag_array_length - q_off;
            let real_diag = diag & d.diag_mask;
            (d.offset, d.last_hit[real_diag as usize], real_diag, diag)
        }
        ExtendWord::Hash(h) => {
            let diag = s_off - q_off;
            (h.offset, h.retrieve(diag).unwrap_or(0), 0, diag)
        }
    };
    let s_off_pos = s_off + table_offset;
    let mut s_end_pos = s_end + table_offset;
    if s_off_pos < last_hit {
        return 0;
    }
    if check_masks {
        let mut extended = 0;
        if !type_of_word(
            ctx,
            subject,
            &mut q_off,
            &mut s_off,
            has_locations,
            s_range,
            word_length,
            lut_word_length,
            &mut extended,
        ) {
            return 0;
        }
        s_end += extended;
        s_end_pos += extended;
    }
    let context = bsearch_context(ctx.qb, q_off);
    let cutoffs = &ctx.params.word_cutoffs[context];
    let ud = if word_length < 11 {
        ungapped_extend_exact(ctx, subject, q_off, s_off, -cutoffs.x_dropoff)
    } else {
        ungapped_extend(
            ctx,
            subject,
            q_off,
            s_end,
            s_off,
            -cutoffs.x_dropoff,
            cutoffs.reduced_nucl_cutoff_score,
        )
    };
    let hit_ready = if ud.score >= cutoffs.cutoff_score {
        init_hitlist.push(InitHsp {
            q_off,
            s_off,
            ungapped: ud,
        });
        s_end_pos = ud.length + ud.s_start + table_offset;
        1
    } else {
        0
    };
    match ewp {
        ExtendWord::Table(d) => d.last_hit[real_diag as usize] = s_end_pos,
        // the C passes `window_size + delta + 1` with delta = MIN(scan_range, window_size -
        // word_length) = -word_length: a cell is reused once `s_off - level > 1 - word_length`
        ExtendWord::Hash(h) => h.insert(diag, s_end_pos, s_off_pos, 1 - word_length),
    }
    hit_ready
}

const fn exact_match_extend_left() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        let b = i as u8;
        t[i] = if b == 0 {
            4
        } else if b & 0x3f == 0 {
            3
        } else if b & 0x0f == 0 {
            2
        } else if b & 0x03 == 0 {
            1
        } else {
            0
        };
        i += 1;
    }
    t
}

const fn exact_match_extend_right() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        let b = i as u8;
        t[i] = if b == 0 {
            4
        } else if b & 0xfc == 0 {
            3
        } else if b & 0xf0 == 0 {
            2
        } else if b & 0xc0 == 0 {
            1
        } else {
            0
        };
        i += 1;
    }
    t
}

const EXT_LEFT: [u8; 256] = exact_match_extend_left();
const EXT_RIGHT: [u8; 256] = exact_match_extend_right();

/// `BlastNaWordFinder`: scan the subject in pieces, extend seeds, then sort the init hits.
pub fn na_word_finder(
    ctx: &WordFinderContext,
    subject: &SubjectBlk,
    scratch: &mut WordFinderScratch,
) {
    let lookup = ctx.lookup;
    let word_length = lookup.word_length;
    let lut_word_length = lookup.lut_word_length;
    let max_hits = scratch.offset_pairs.len() as i32;
    let mut range_idx = 0usize;
    let mut scan_range = [0i32, subject.length - lut_word_length];
    scratch.init_hitlist.clear();
    loop {
        // s_DetermineScanningOffsets
        let mut ok = true;
        while scan_range[0] > scan_range[1] {
            range_idx += 1;
            if range_idx >= subject.seq_ranges.len() {
                ok = false;
                break;
            }
            scan_range[0] = subject.seq_ranges[range_idx].0 + word_length - lut_word_length;
            scan_range[1] = subject.seq_ranges[range_idx].1 - lut_word_length;
        }
        if !ok {
            break;
        }
        let hits = match &lookup.table {
            Table::Mb(t) => mb_scan_any(
                t,
                subject.seq,
                &mut scratch.offset_pairs,
                max_hits,
                &mut scan_range,
            ),
            Table::SmallNa(t) => small_na_scan_any(
                t,
                subject.seq,
                &mut scratch.offset_pairs,
                max_hits,
                &mut scan_range,
            ),
            Table::Na(t) => {
                if t.lut_word_length == 8 && t.scan_step == 4 {
                    na_scan_8_4(
                        t,
                        subject.seq,
                        &mut scratch.offset_pairs,
                        max_hits,
                        &mut scan_range,
                    )
                } else {
                    na_scan_any(
                        t,
                        subject.seq,
                        &mut scratch.offset_pairs,
                        max_hits,
                        &mut scan_range,
                    )
                }
            }
        };
        if hits == 0 {
            continue;
        }
        let s_range = scan_range[1] + lut_word_length;
        extend_hits(ctx, subject, hits as usize, s_range, scratch);
    }
    scratch.ewp.exit(subject.length);
    scratch.init_hitlist.sort_by(score_compare_match);
    debug_assert!(
        scratch
            .init_hitlist
            .windows(2)
            .all(|w| score_compare_match(&w[0], &w[1]) != Ordering::Equal),
        "init hit list has ties"
    );
}

/// The `extend_callback`: `s_BlastNaExtendDirect`, `s_BlastNaExtend`, `s_BlastNaExtendAligned`,
/// `s_BlastSmallNaExtend` or `s_BlastSmallNaExtendAlignedOneByte` per `lookup.extender`.
fn extend_hits(
    ctx: &WordFinderContext,
    subject: &SubjectBlk,
    num_hits: usize,
    s_range: i32,
    scratch: &mut WordFinderScratch,
) {
    let lookup = ctx.lookup;
    let word_length = lookup.word_length;
    let lut_word_length = lookup.lut_word_length;
    let has_locations = lookup.masked_locations.is_some();
    let query = ctx.qb.sequence();
    // `query[-1]` is the leading sentinel in the C (sequence_start[0]); index through buf.
    let buf = &ctx.qb.buf;
    let qat = |i: isize| -> u8 { buf[(i + 1) as usize] };
    let subj = subject.seq;
    let qlen = ctx.qb.length;
    for i in 0..num_hits {
        let OffsetPair { q_off, s_off } = scratch.offset_pairs[i];
        let (mut q_offset, mut s_offset) = (q_off, s_off);
        match lookup.extender {
            Extender::Direct => {
                extend_initial_hit(
                    ctx,
                    subject,
                    q_offset,
                    s_offset,
                    false,
                    s_range,
                    word_length,
                    word_length,
                    &mut scratch.ewp,
                    &mut scratch.init_hitlist,
                    true,
                );
                continue;
            }
            Extender::Na => {
                let ext_to = word_length - lut_word_length;
                let mut ext_left = 0;
                let mut so = s_offset;
                let mut q = q_offset as isize;
                let mut s = (so / 4) as isize;
                let lim = ext_to.min(s_offset);
                while ext_left < lim {
                    so -= 1;
                    q -= 1;
                    if so % 4 == 3 {
                        s -= 1;
                    }
                    if ((subj[s as usize] << (2 * (so % 4))) >> 6) != qat(q) {
                        break;
                    }
                    ext_left += 1;
                }
                if ext_left < ext_to {
                    let mut ext_right = 0;
                    let mut so = s_offset + lut_word_length;
                    if so + ext_to - ext_left > s_range {
                        continue;
                    }
                    let mut q = (q_offset + lut_word_length) as isize;
                    let mut s = (so / 4) as isize;
                    while ext_right < ext_to - ext_left {
                        if ((subj[s as usize] << (2 * (so % 4))) >> 6) != query[q as usize] {
                            break;
                        }
                        so += 1;
                        q += 1;
                        if so % 4 == 0 {
                            s += 1;
                        }
                        ext_right += 1;
                    }
                    if ext_left + ext_right < ext_to {
                        continue;
                    }
                }
                q_offset -= ext_left;
                s_offset -= ext_left;
            }
            Extender::NaAligned => {
                let ext_to = word_length - lut_word_length;
                let mut ext_left = 0;
                let mut ext_max = ext_to.min(s_offset);
                let mut q = q_offset as isize;
                let mut s = (s_offset / 4) as isize;
                'left: while ext_left < ext_max {
                    let byte = subj[(s - 1) as usize];
                    let qq = |k: isize| qat(q - k);
                    if (byte & 3) != qq(1) || {
                        ext_left += 1;
                        ext_left == ext_max
                    } {
                        break 'left;
                    }
                    if ((byte >> 2) & 3) != qq(2) || {
                        ext_left += 1;
                        ext_left == ext_max
                    } {
                        break 'left;
                    }
                    if ((byte >> 4) & 3) != qq(3) || {
                        ext_left += 1;
                        ext_left == ext_max
                    } {
                        break 'left;
                    }
                    if (byte >> 6) != qq(4) {
                        break 'left;
                    }
                    s -= 1;
                    q -= 4;
                    ext_left += 1;
                }
                if ext_left < ext_to {
                    let mut ext_right = 0;
                    ext_max = ext_to - ext_left;
                    if s_offset + lut_word_length + ext_max > s_range {
                        continue;
                    }
                    let mut q = (q_offset + lut_word_length) as isize;
                    let mut s = ((s_offset + lut_word_length) / 4) as isize;
                    'right: while ext_right < ext_max {
                        let byte = subj[s as usize];
                        let qq = |k: isize| query[(q + k) as usize];
                        if (byte >> 6) != qq(0) || {
                            ext_right += 1;
                            ext_right == ext_max
                        } {
                            break 'right;
                        }
                        if ((byte >> 4) & 3) != qq(1) || {
                            ext_right += 1;
                            ext_right == ext_max
                        } {
                            break 'right;
                        }
                        if ((byte >> 2) & 3) != qq(2) || {
                            ext_right += 1;
                            ext_right == ext_max
                        } {
                            break 'right;
                        }
                        if (byte & 3) != qq(3) {
                            break 'right;
                        }
                        s += 1;
                        q += 4;
                        ext_right += 1;
                    }
                    if ext_left + ext_right < ext_to {
                        continue;
                    }
                }
                q_offset -= ext_left;
                s_offset -= ext_left;
            }
            Extender::SmallNaAlignedOneByte => {
                let comp = ctx.qb.compressed.as_ref().expect("compressed query");
                let ext_to = word_length - lut_word_length;
                let mut ext_left = 0;
                let context = bsearch_context(ctx.qb, q_offset);
                let c = &ctx.qb.contexts[context];
                let q_start = c.query_offset;
                let q_range = q_start + c.query_length;
                if s_offset > 0 && q_offset > 0 {
                    let q_byte = comp[(q_offset - 4 + 3) as usize];
                    let s_byte = subj[(s_offset / 4 - 1) as usize];
                    ext_left = EXT_LEFT[(q_byte ^ s_byte) as usize] as i32;
                    ext_left = ext_left.min(ext_to).min(q_offset - q_start);
                }
                if ext_left < ext_to && (q_offset + lut_word_length) < qlen {
                    let q_byte = comp[(q_offset + lut_word_length + 3) as usize];
                    let s_byte = subj[((s_offset + lut_word_length) / 4) as usize];
                    let mut ext_right = EXT_RIGHT[(q_byte ^ s_byte) as usize] as i32;
                    ext_right = ext_right
                        .min(s_range - (s_offset + lut_word_length))
                        .min(q_range - (q_offset + lut_word_length));
                    if ext_left + ext_right < ext_to {
                        continue;
                    }
                }
                q_offset -= ext_left;
                s_offset -= ext_left;
            }
            Extender::SmallNa => {
                let comp = ctx.qb.compressed.as_ref().expect("compressed query");
                let context = bsearch_context(ctx.qb, q_offset);
                let c = &ctx.qb.contexts[context];
                let q_start = c.query_offset;
                let q_range = q_start + c.query_length;
                let mut ext_max = (word_length - lut_word_length)
                    .min(s_offset)
                    .min(q_offset - q_start);
                let rsdl = 4 - (s_offset % 4);
                s_offset += rsdl;
                q_offset += rsdl;
                ext_max += rsdl;
                let mut so = s_offset;
                let mut qo = q_offset;
                let mut ext_left = 0;
                while ext_left < ext_max {
                    let q_byte = comp[(qo - 4 + 3) as usize];
                    let s_byte = subj[(so / 4 - 1) as usize];
                    let bases = EXT_LEFT[(q_byte ^ s_byte) as usize] as i32;
                    ext_left += bases;
                    if bases < 4 {
                        break;
                    }
                    qo -= 4;
                    so -= 4;
                }
                ext_left = ext_left.min(ext_max);
                so = s_offset;
                qo = q_offset;
                ext_max = (word_length - ext_left).min(s_range - so).min(q_range - qo);
                let mut ext_right = 0;
                while ext_right < ext_max {
                    let q_byte = comp[(qo + 3) as usize];
                    let s_byte = subj[(so / 4) as usize];
                    let bases = EXT_RIGHT[(q_byte ^ s_byte) as usize] as i32;
                    ext_right += bases;
                    if bases < 4 {
                        break;
                    }
                    qo += 4;
                    so += 4;
                }
                ext_right = ext_right.min(ext_max);
                if ext_left + ext_right < word_length {
                    continue;
                }
                q_offset -= ext_left;
                s_offset -= ext_left;
            }
        }
        extend_initial_hit(
            ctx,
            subject,
            q_offset,
            s_offset,
            has_locations,
            s_range,
            word_length,
            lut_word_length,
            &mut scratch.ewp,
            &mut scratch.init_hitlist,
            true,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extend_tables_match_c() {
        // from na_ungapped.c: s_ExactMatchExtendLeft starts 4,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,2,...
        assert_eq!(
            &EXT_LEFT[..17],
            &[4, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2]
        );
        assert_eq!(EXT_LEFT[64], 3);
        assert_eq!(EXT_LEFT[128], 3);
        assert_eq!(
            &EXT_RIGHT[..16],
            &[4, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        );
        assert_eq!(EXT_RIGHT[16], 1);
        assert_eq!(EXT_RIGHT[63], 1);
        assert_eq!(EXT_RIGHT[64], 0);
    }
}
