//! Dynamic-programming gapped extension: the `-task blastn` (`eDynProgScoreOnly` /
//! `eDynProgTbck`) path.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/core/blast_gapalign.c`: `BLAST_GapAlignStructNew` (the
//! `dp_mem` part), `s_GapGetState`/`s_GapPurgeState` (the traceback state arena, kept here as
//! one growable byte array with a start offset per row), `ALIGN_EX`,
//! `s_BlastAlignPackedNucl`, `s_BlastDynProgNtGappedAlignment`,
//! `BLAST_GappedAlignmentWithTraceback` and `AdjustSubjectRange`.
//!
//! `Blast_SemiGappedAlign` (score-only on two blastna sequences) is not reached by blastn:
//! the preliminary stage scores on the packed subject and the traceback always calls
//! `ALIGN_EX`, so it is not ported.

use super::greedy::{PrelimEditBlock, prelim_to_edit_script};
use super::hsp::{EditOp, EditScript};

/// `MININT` of blast_gapalign.c.
const MININT: i32 = i32::MIN / 2;
pub(crate) const MAX_SUBJECT_OFFSET: i32 = 90000;
pub(crate) const MAX_TOTAL_GAPS: i32 = 3000;

const SCRIPT_SUB: u8 = 3;
const SCRIPT_GAP_IN_A: u8 = 0;
const SCRIPT_GAP_IN_B: u8 = 6;
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;

/// `BlastGapDP`.
#[derive(Debug, Clone, Copy, Default)]
struct GapDp {
    best: i32,
    best_gap: i32,
}

/// What `s_BlastDynProgNtGappedAlignment` / `BLAST_GappedAlignmentWithTraceback` leave in
/// `gap_align`.
#[derive(Debug, Clone, Default)]
pub struct DpResult {
    pub query_start: i32,
    pub query_stop: i32,
    pub subject_start: i32,
    pub subject_stop: i32,
    pub score: i32,
    pub edit_script: Option<EditScript>,
}

/// `BlastGapAlignStruct` for the dynamic-programming path.
pub struct DpAligner {
    pub gap_x_dropoff: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    /// `dp_mem` (`dp_mem_alloc` = its length).
    dp_mem: Vec<GapDp>,
    /// `state_struct`: every traceback row's cells, back to back.
    state: Vec<u8>,
    /// `edit_script[a]`: where row `a` starts in `state`.
    row_base: Vec<usize>,
    /// `edit_start_offset[a]`.
    edit_start_offset: Vec<i32>,
    pub fwd: PrelimEditBlock,
    pub rev: PrelimEditBlock,
}

#[inline(always)]
fn unpack(byte: u8, n: i32) -> usize {
    ((byte >> (2 * n)) & 3) as usize
}

/// A letter of `seq` at `base + i`, 15 (the sentinel the C would read) before the start.
#[inline(always)]
fn letter(seq: &[u8], base: isize, i: i32) -> usize {
    let idx = base + i as isize;
    if idx < 0 {
        15
    } else {
        seq[idx as usize] as usize
    }
}

/// `AdjustSubjectRange`: (new subject offset, new subject length, start_shift).
pub fn adjust_subject_range(
    subject_offset: i32,
    subject_length: i32,
    query_offset: i32,
    query_length: i32,
) -> (i32, i32, i32) {
    if subject_length < MAX_SUBJECT_OFFSET {
        return (subject_offset, subject_length, 0);
    }
    let s_offset = subject_offset;
    let max_extension_left = query_offset + MAX_TOTAL_GAPS;
    let max_extension_right = query_length - query_offset + MAX_TOTAL_GAPS;
    let (start_shift, new_offset) = if s_offset <= max_extension_left {
        (0, s_offset)
    } else {
        (s_offset - max_extension_left, max_extension_left)
    };
    let new_length = subject_length.min(s_offset + max_extension_right) - start_shift;
    (new_offset, new_length, start_shift)
}

impl DpAligner {
    /// `BLAST_GapAlignStructNew` for `eDynProgScoreOnly`: `dp_mem_alloc = 1000`.
    pub fn new(gap_open: i32, gap_extend: i32, gap_x_dropoff: i32) -> DpAligner {
        DpAligner {
            gap_x_dropoff,
            gap_open,
            gap_extend,
            dp_mem: vec![GapDp::default(); 1000],
            state: Vec::new(),
            row_base: Vec::new(),
            edit_start_offset: Vec::new(),
            fwd: PrelimEditBlock::default(),
            rev: PrelimEditBlock::default(),
        }
    }

    /// `s_BlastAlignPackedNucl`: score-only extension of `b` (blastna letters, `b[base_b + i]`)
    /// against the packed subject `a` (`a[base_a + i]` bytes), `n` letters of B and `m` of A.
    /// Returns (best score, a_offset, b_offset).
    #[rustfmt::skip]
    fn align_packed_nucl(
        &mut self,
        b: &[u8], base_b: isize,
        a: &[u8], base_a: usize,
        n: i32, m: i32,
        matrix: &[[i32; 16]; 16],
        reverse_sequence: bool,
        mut x_dropoff: i32,
    ) -> (i32, i32, i32) {
        let gap_open_extend = self.gap_open + self.gap_extend;
        let gap_extend = self.gap_extend;
        if x_dropoff < gap_open_extend {
            x_dropoff = gap_open_extend;
        }
        if n <= 0 || m <= 0 {
            return (0, 0, 0);
        }
        let num_extra_cells = if gap_extend > 0 { x_dropoff / gap_extend + 3 } else { n + 3 };
        if num_extra_cells as usize > self.dp_mem.len() {
            let want = (num_extra_cells as usize + 100).max(2 * self.dp_mem.len());
            self.dp_mem.resize(want, GapDp::default());
        }
        let score_array = &mut self.dp_mem;
        let mut score = -gap_open_extend;
        score_array[0].best = 0;
        score_array[0].best_gap = -gap_open_extend;
        let mut i = 1;
        while i <= n {
            if score < -x_dropoff {
                break;
            }
            score_array[i as usize].best = score;
            score_array[i as usize].best_gap = score - gap_open_extend;
            score -= gap_extend;
            i += 1;
        }
        let mut b_size = i;
        let mut best_score = 0;
        let mut a_offset = 0;
        let mut b_offset = 0;
        let mut first_b_index = 0;
        let b_increment: i32 = if reverse_sequence { -1 } else { 1 };
        let mut a_index = 1;
        while a_index <= m {
            let a_base_pair = if reverse_sequence {
                unpack(a[base_a + ((m - a_index) / 4) as usize], (a_index - 1) % 4)
            } else {
                unpack(a[base_a + 1 + ((a_index - 1) / 4) as usize], 3 - ((a_index - 1) % 4))
            };
            let matrix_row = &matrix[a_base_pair];
            let mut b_ptr: i32 = if reverse_sequence { n - first_b_index } else { first_b_index };
            score = MININT;
            let mut score_gap_row = MININT;
            let mut last_b_index = first_b_index;
            let mut b_index = first_b_index;
            while b_index < b_size {
                b_ptr += b_increment;
                let score_gap_col = score_array[b_index as usize].best_gap;
                let next_score = score_array[b_index as usize].best + matrix_row[letter(b, base_b, b_ptr)];
                if score < score_gap_col {
                    score = score_gap_col;
                }
                if score < score_gap_row {
                    score = score_gap_row;
                }
                if best_score - score > x_dropoff {
                    if b_index == first_b_index {
                        first_b_index += 1;
                    } else {
                        score_array[b_index as usize].best = MININT;
                    }
                } else {
                    last_b_index = b_index;
                    if score > best_score {
                        best_score = score;
                        a_offset = a_index;
                        b_offset = b_index;
                    }
                    score_gap_row -= gap_extend;
                    let score_gap_col = score_gap_col - gap_extend;
                    score_array[b_index as usize].best_gap = (score - gap_open_extend).max(score_gap_col);
                    score_gap_row = (score - gap_open_extend).max(score_gap_row);
                    score_array[b_index as usize].best = score;
                }
                score = next_score;
                b_index += 1;
            }
            if first_b_index == b_size {
                break;
            }
            if (last_b_index + num_extra_cells + 3) as usize >= score_array.len() {
                let want = ((last_b_index + num_extra_cells + 100) as usize).max(2 * score_array.len());
                score_array.resize(want, GapDp::default());
            }
            if last_b_index < b_size - 1 {
                b_size = last_b_index + 1;
            } else {
                while score_gap_row >= (best_score - x_dropoff) && b_size <= n {
                    score_array[b_size as usize].best = score_gap_row;
                    score_array[b_size as usize].best_gap = score_gap_row - gap_open_extend;
                    score_gap_row -= gap_extend;
                    b_size += 1;
                }
            }
            if b_size <= n {
                score_array[b_size as usize].best = MININT;
                score_array[b_size as usize].best_gap = MININT;
                b_size += 1;
            }
            a_index += 1;
        }
        (best_score, a_offset, b_offset)
    }

    /// `ALIGN_EX`: gapped extension with traceback of `a` (`a[base_a + i]`) against `b`
    /// (`b[base_b + i]`), both blastna, `m` letters of A and `n` of B. The edit operations
    /// go to `self.rev` when `reverse_sequence`, else `self.fwd`. Returns (best score,
    /// a_offset, b_offset).
    #[rustfmt::skip]
    fn align_ex(
        &mut self,
        a: &[u8], base_a: isize,
        b: &[u8], base_b: isize,
        m: i32, n: i32,
        matrix: &[[i32; 16]; 16],
        reverse_sequence: bool,
    ) -> (i32, i32, i32) {
        let gap_open_extend = self.gap_open + self.gap_extend;
        let gap_extend = self.gap_extend;
        let mut x_dropoff = self.gap_x_dropoff;
        if x_dropoff < gap_open_extend {
            x_dropoff = gap_open_extend;
        }
        if n <= 0 || m <= 0 {
            return (0, 0, 0);
        }
        let DpAligner { dp_mem, state, row_base, edit_start_offset, fwd, rev, .. } = self;
        let edit_block = if reverse_sequence { rev } else { fwd };
        // s_GapPurgeState
        state.clear();
        row_base.clear();
        edit_start_offset.clear();
        let num_extra_cells = if gap_extend > 0 { x_dropoff / gap_extend + 3 } else { n + 3 };
        if num_extra_cells as usize > dp_mem.len() {
            let want = (num_extra_cells as usize + 100).max(2 * dp_mem.len());
            dp_mem.resize(want, GapDp::default());
        }
        let score_array = dp_mem;
        // row 0
        state.resize(num_extra_cells as usize + 1, 0);
        row_base.push(0);
        edit_start_offset.push(0);
        let mut score = -gap_open_extend;
        score_array[0].best = 0;
        score_array[0].best_gap = -gap_open_extend;
        let mut i = 1;
        while i <= n {
            if score < -x_dropoff {
                break;
            }
            score_array[i as usize].best = score;
            score_array[i as usize].best_gap = score - gap_open_extend;
            score -= gap_extend;
            state[i as usize] = SCRIPT_GAP_IN_A;
            i += 1;
        }
        state.truncate(i as usize + 1);
        let mut b_size = i;
        let mut best_score = 0;
        let mut a_offset = 0;
        let mut b_offset = 0;
        let mut first_b_index = 0;
        let b_increment: i32 = if reverse_sequence { -1 } else { 1 };
        let mut a_index = 1;
        while a_index <= m {
            // s_GapGetState for this row
            let row_len = if gap_extend > 0 { b_size - first_b_index + num_extra_cells } else { n + 3 - first_b_index };
            let base = state.len();
            state.resize(base + row_len as usize + 1, 0);
            row_base.push(base);
            edit_start_offset.push(first_b_index);
            let orig_b_index = first_b_index;
            let matrix_row = if reverse_sequence {
                &matrix[letter(a, base_a, m - a_index)]
            } else {
                &matrix[letter(a, base_a, a_index)]
            };
            let mut b_ptr: i32 = if reverse_sequence { n - first_b_index } else { first_b_index };
            score = MININT;
            let mut score_gap_row = MININT;
            let mut last_b_index = first_b_index;
            let mut b_index = first_b_index;
            while b_index < b_size {
                b_ptr += b_increment;
                let score_gap_col = score_array[b_index as usize].best_gap;
                let next_score = score_array[b_index as usize].best + matrix_row[letter(b, base_b, b_ptr)];
                let mut script = SCRIPT_SUB;
                let script_col = SCRIPT_EXTEND_GAP_B;
                let script_row = SCRIPT_EXTEND_GAP_A;
                if score < score_gap_col {
                    script = SCRIPT_GAP_IN_B;
                    score = score_gap_col;
                }
                if score < score_gap_row {
                    script = SCRIPT_GAP_IN_A;
                    score = score_gap_row;
                }
                if best_score - score > x_dropoff {
                    if first_b_index == b_index {
                        first_b_index += 1;
                    } else {
                        score_array[b_index as usize].best = MININT;
                    }
                } else {
                    last_b_index = b_index;
                    if score > best_score {
                        best_score = score;
                        a_offset = a_index;
                        b_offset = b_index;
                    }
                    score_gap_row -= gap_extend;
                    let score_gap_col = score_gap_col - gap_extend;
                    if score_gap_col < (score - gap_open_extend) {
                        score_array[b_index as usize].best_gap = score - gap_open_extend;
                    } else {
                        score_array[b_index as usize].best_gap = score_gap_col;
                        script += script_col;
                    }
                    if score_gap_row < (score - gap_open_extend) {
                        score_gap_row = score - gap_open_extend;
                    } else {
                        script += script_row;
                    }
                    score_array[b_index as usize].best = score;
                }
                score = next_score;
                state[base + (b_index - orig_b_index) as usize] = script;
                b_index += 1;
            }
            if first_b_index == b_size {
                break;
            }
            if (last_b_index + num_extra_cells + 3) as usize >= score_array.len() {
                let want = ((last_b_index + num_extra_cells + 100) as usize).max(2 * score_array.len());
                score_array.resize(want, GapDp::default());
            }
            if last_b_index < b_size - 1 {
                b_size = last_b_index + 1;
            } else {
                while score_gap_row >= (best_score - x_dropoff) && b_size <= n {
                    score_array[b_size as usize].best = score_gap_row;
                    score_array[b_size as usize].best_gap = score_gap_row - gap_open_extend;
                    score_gap_row -= gap_extend;
                    state[base + (b_size - orig_b_index) as usize] = SCRIPT_GAP_IN_A;
                    b_size += 1;
                }
            }
            // state_struct->used += MAX(b_index, b_size) - orig_b_index + 1
            state.truncate(base + (b_index.max(b_size) - orig_b_index) as usize + 1);
            if b_size <= n {
                score_array[b_size as usize].best = MININT;
                score_array[b_size as usize].best_gap = MININT;
                b_size += 1;
            }
            a_index += 1;
        }
        // the optimal path through edit_script[][]
        let mut a_index = a_offset;
        let mut b_index = b_offset;
        let mut script = SCRIPT_SUB;
        while a_index > 0 || b_index > 0 {
            let next_script = state[row_base[a_index as usize] + (b_index - edit_start_offset[a_index as usize]) as usize];
            match script {
                SCRIPT_GAP_IN_A => {
                    script = next_script & SCRIPT_OP_MASK;
                    if next_script & SCRIPT_EXTEND_GAP_A != 0 {
                        script = SCRIPT_GAP_IN_A;
                    }
                }
                SCRIPT_GAP_IN_B => {
                    script = next_script & SCRIPT_OP_MASK;
                    if next_script & SCRIPT_EXTEND_GAP_B != 0 {
                        script = SCRIPT_GAP_IN_B;
                    }
                }
                _ => {
                    script = next_script & SCRIPT_OP_MASK;
                }
            }
            let op = if script == SCRIPT_GAP_IN_A {
                b_index -= 1;
                EditOp::Del
            } else if script == SCRIPT_GAP_IN_B {
                a_index -= 1;
                EditOp::Ins
            } else {
                a_index -= 1;
                b_index -= 1;
                EditOp::Sub
            };
            edit_block.add(op, 1);
        }
        (best_score, a_offset, b_offset)
    }

    /// `s_BlastDynProgNtGappedAlignment`. `qbuf[q_base + i]` is letter `i` of the context
    /// (`q_base >= 1`, so the sentinel before it is readable), `subject` is the packed
    /// (ncbi2na) chunk of `subject_length` bases, `(q_off, s_off)` the initial HSP's
    /// offsets and `ungapped_score` its ungapped score.
    pub fn dyn_prog_nt_gapped_alignment(
        &mut self,
        qbuf: &[u8],
        q_base: usize,
        query_length: i32,
        subject: &[u8],
        subject_length: i32,
        q_off: i32,
        s_off: i32,
        ungapped_score: i32,
        matrix: &[[i32; 16]; 16],
    ) -> DpResult {
        let mut x_dropoff = self.gap_x_dropoff;
        if ungapped_score < x_dropoff {
            x_dropoff = ungapped_score;
        }
        let offset_adjustment = 4 - (s_off % 4);
        let mut q_length = q_off + offset_adjustment;
        let mut s_length = s_off + offset_adjustment;
        if q_length > query_length || s_length > subject_length {
            q_length -= 4;
            s_length -= 4;
        }
        // s_BlastAlignPackedNucl(B = query, A = subject, ..., b_offset = &private_q_start,
        // a_offset = &private_s_start): the a_offset is the subject's
        let (score_left, private_s_start, private_q_start) = self.align_packed_nucl(
            qbuf,
            q_base as isize,
            subject,
            0,
            q_length,
            s_length,
            matrix,
            true,
            x_dropoff,
        );
        let mut out = DpResult {
            query_start: q_length - private_q_start,
            subject_start: s_length - private_s_start,
            ..Default::default()
        };
        let mut score_right = 0;
        if q_length < query_length && s_length < subject_length {
            let (sr, ss, qs) = self.align_packed_nucl(
                qbuf,
                q_base as isize + q_length as isize - 1,
                subject,
                ((s_length + 3) / 4 - 1) as usize,
                query_length - q_length,
                subject_length - s_length,
                matrix,
                false,
                x_dropoff,
            );
            score_right = sr;
            out.query_stop = qs + q_length;
            out.subject_stop = ss + s_length;
        } else {
            out.query_stop = q_length;
            out.subject_stop = s_length;
        }
        out.score = score_right + score_left;
        out
    }

    /// `BLAST_GappedAlignmentWithTraceback` for blastn (no out-of-frame). `query` is the
    /// context's letters (`query.len() >= query_length + 1`), `subject` the blastna subject
    /// (`subject[s] for s in 0..=subject_length` readable).
    pub fn gapped_alignment_with_traceback(
        &mut self,
        query: &[u8],
        subject: &[u8],
        q_start: i32,
        s_start: i32,
        query_length: i32,
        subject_length: i32,
        matrix: &[[i32; 16]; 16],
    ) -> DpResult {
        self.fwd.reset();
        self.rev.reset();
        let q_length = query_length;
        let s_length = subject_length;
        // left extension includes [q_start, s_start]
        let (mut score_left, private_q_length, private_s_length) =
            self.align_ex(query, 0, subject, 0, q_start + 1, s_start + 1, matrix, true);
        let mut out = DpResult {
            query_start: q_start - private_q_length + 1,
            subject_start: s_start - private_s_length + 1,
            ..Default::default()
        };
        let mut score_right = 0;
        let mut found_end = false;
        if q_start < q_length && s_start < s_length {
            found_end = true;
            let (sr, pq, ps) = self.align_ex(
                query,
                q_start as isize,
                subject,
                s_start as isize,
                q_length - q_start - 1,
                s_length - s_start - 1,
                matrix,
                false,
            );
            score_right = sr;
            out.query_stop = q_start + pq + 1;
            out.subject_stop = s_start + ps + 1;
        }
        if !found_end {
            out.query_stop = q_start - 1;
            out.subject_stop = s_start - 1;
        }
        let mut esp = prelim_to_edit_script(&self.rev, &self.fwd);
        // prune leading and trailing gaps
        while !esp.ops.is_empty() && esp.ops[0].0 != EditOp::Sub {
            score_left += self.gap_open + esp.ops[0].1 * self.gap_extend;
            if esp.ops[0].0 == EditOp::Del {
                out.subject_start += esp.ops[0].1;
            } else {
                out.query_start += esp.ops[0].1;
            }
            esp.ops.remove(0);
        }
        while !esp.ops.is_empty() && esp.ops[esp.ops.len() - 1].0 != EditOp::Sub {
            let (op, num) = esp.ops[esp.ops.len() - 1];
            score_right += self.gap_open + num * self.gap_extend;
            if op == EditOp::Del {
                out.subject_stop -= num;
            } else {
                out.query_stop -= num;
            }
            esp.ops.pop();
        }
        out.edit_script = Some(esp);
        out.score = score_right + score_left;
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::stat::nucl_matrix;

    fn letters(s: &[u8]) -> Vec<u8> {
        s.iter()
            .map(|&c| crate::engine::encoding::ASCII_TO_BLASTNA[c as usize])
            .collect()
    }

    fn pack(l: &[u8]) -> Vec<u8> {
        let mut out = vec![0u8; l.len() / 4 + 1];
        for (i, &b) in l.iter().enumerate() {
            out[i / 4] |= (b & 3) << (2 * (3 - (i % 4)));
        }
        out
    }

    #[test]
    fn identical_sequences_score_fully() {
        let m = nucl_matrix(2, -3);
        let q = letters(b"ACGTTGCAAGGCTTAACCGGTTAAGGCC");
        let mut qbuf = vec![15u8];
        qbuf.extend_from_slice(&q);
        qbuf.push(15);
        let s = pack(&q);
        let mut dp = DpAligner::new(5, 2, 50);
        let r = dp.dyn_prog_nt_gapped_alignment(
            &qbuf,
            1,
            q.len() as i32,
            &s,
            q.len() as i32,
            10,
            10,
            20,
            &m,
        );
        assert_eq!(r.score, 2 * q.len() as i32);
        assert_eq!((r.query_start, r.query_stop), (0, q.len() as i32));
        assert_eq!((r.subject_start, r.subject_stop), (0, q.len() as i32));
        // traceback: one substitution block
        let mut sb = q.clone();
        sb.push(15);
        let t =
            dp.gapped_alignment_with_traceback(&q, &sb, 10, 10, q.len() as i32, q.len() as i32, &m);
        assert_eq!(t.score, 2 * q.len() as i32);
        assert_eq!(
            t.edit_script.unwrap().ops,
            vec![(EditOp::Sub, q.len() as i32)]
        );
    }

    #[test]
    fn mismatch_and_gap_are_scored() {
        let m = nucl_matrix(2, -3);
        // subject = query with one deleted base at 14 and a mismatch at 5
        let q = letters(b"ACGTTGCAAGGCTTAACCGGTTAAGGCCATGC");
        let mut s_letters = q.clone();
        s_letters[5] = (s_letters[5] + 1) % 4;
        s_letters.remove(14);
        let mut qbuf = vec![15u8];
        qbuf.extend_from_slice(&q);
        qbuf.push(15);
        let s = pack(&s_letters);
        let mut dp = DpAligner::new(5, 2, 50);
        let r = dp.dyn_prog_nt_gapped_alignment(
            &qbuf,
            1,
            q.len() as i32,
            &s,
            s_letters.len() as i32,
            20,
            19,
            20,
            &m,
        );
        // 31 matches, one mismatch, one gap of length 1 in the subject
        assert_eq!(r.score, 2 * 30 - 3 - 7);
        // the box: the query is one letter longer than the subject (the a/b offsets of
        // s_BlastAlignPackedNucl are the subject's and the query's, in that order)
        assert_eq!((r.query_start, r.query_stop), (0, q.len() as i32));
        assert_eq!(
            (r.subject_start, r.subject_stop),
            (0, s_letters.len() as i32)
        );
        let mut sb = s_letters.clone();
        sb.push(15);
        let t = dp.gapped_alignment_with_traceback(
            &q,
            &sb,
            20,
            19,
            q.len() as i32,
            s_letters.len() as i32,
            &m,
        );
        assert_eq!(t.score, r.score);
        let esp = t.edit_script.unwrap();
        let gaps: i32 = esp
            .ops
            .iter()
            .filter(|(o, _)| *o == EditOp::Ins)
            .map(|(_, n)| *n)
            .sum();
        assert_eq!(gaps, 1);
        assert_eq!((t.query_start, t.query_stop), (0, q.len() as i32));
    }

    #[test]
    fn subject_range_window() {
        assert_eq!(adjust_subject_range(1000, 50000, 10, 100), (1000, 50000, 0));
        let (o, l, shift) = adjust_subject_range(500_000, 1_000_000, 10, 100);
        assert_eq!(shift, 500_000 - 3010);
        assert_eq!(o, 3010);
        assert_eq!(l, 500_000 + 90 + 3000 - shift);
    }
}
