//! Greedy gapped extension (linear gap costs), score-only and with traceback.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/core/greedy_align.c` (`BLAST_GreedyAlign`,
//! `s_FindFirstMismatch`, `s_GetNextNonAffineTback`), `blast_gapalign.c`
//! (`s_BlastGreedyAlignMemAlloc`,
//! `BLAST_GapAlignStructNew`'s `max_dist` rule, `BLAST_GreedyGappedAlignment`,
//! `Blast_PrelimEditBlockToGapEditScript`, `s_UpdateEditScript`, `s_RebuildEditScript`,
//! `s_ReduceGaps`) and `gapinfo.c` (`GapPrelimEditBlockAdd`).
//!
//! Megablast's gap costs are 0/0, so the aligner tracks the edit distance (`greedy_align`);
//! `BLAST_AffineGreedyAlign` (non-zero greedy gap costs, `-gapopen`/`-gapextend` with
//! megablast) is not ported.

use super::hsp::{EditOp, EditScript};

const K_INVALID_OFFSET: i32 = -2;
pub const GREEDY_MAX_COST: i32 = 1000;
pub const GREEDY_MAX_COST_FRACTION: i32 = 2;

/// `SGreedySeed`.
#[derive(Debug, Clone, Copy, Default)]
pub struct GreedySeed {
    pub start_q: i32,
    pub start_s: i32,
    pub match_length: i32,
}

/// `GapPrelimEditBlock`.
#[derive(Debug, Clone, Default)]
pub struct PrelimEditBlock {
    pub ops: Vec<(EditOp, i32)>,
}

impl PrelimEditBlock {
    pub fn reset(&mut self) {
        self.ops.clear();
    }
    /// `GapPrelimEditBlockAdd`.
    pub fn add(&mut self, op: EditOp, num: i32) {
        if num == 0 {
            return;
        }
        if let Some(last) = self.ops.last_mut()
            && last.0 == op
        {
            last.1 += num;
            return;
        }
        self.ops.push((op, num));
    }
}

/// A sequence the greedy aligner walks: blastna letters (`rem == 4`) or packed ncbi2na
/// with `rem` leading bases to skip in the first byte.
#[derive(Clone, Copy)]
pub struct Seq2<'a> {
    pub bytes: &'a [u8],
    pub rem: u8,
}

#[inline(always)]
fn unpack(bytes: &[u8], pos: usize) -> u8 {
    (bytes[pos / 4] >> (2 * (3 - (pos % 4)))) & 3
}

/// `s_FindFirstMismatch`: letters in seq2 at logical position p are `bytes[(p+rem)/4]`'s base
/// `(p+rem)%4` when packed; for the reverse direction the C indexes from the end.
#[inline]
fn find_first_mismatch(
    seq1: &[u8],
    seq2: Seq2,
    len1: i32,
    len2: i32,
    mut i1: i32,
    mut i2: i32,
    reverse: bool,
) -> i32 {
    let start = i1;
    if seq2.rem == 4 {
        // Both sequences are one letter per byte: compare eight at a time. A byte pair
        // matches when the letters are equal and the query letter is a base (< 4, so its
        // top six bits are clear); the first offending byte of the word is the mismatch.
        const AMBIG: u64 = 0xFCFC_FCFC_FCFC_FCFC;
        let n = (len1 as isize - i1 as isize)
            .min(len2 as isize - i2 as isize)
            .min(seq1.len() as isize - i1 as isize)
            .min(seq2.bytes.len() as isize - i2 as isize)
            .max(0) as usize;
        if n == 0 {
            return 0;
        }
        let matched = if reverse {
            let (e1, e2) = ((len1 - i1) as usize, (len2 - i2) as usize);
            let a = &seq1[e1 - n..e1];
            let b = &seq2.bytes[e2 - n..e2];
            let mut k = 0usize;
            loop {
                if k + 8 <= n {
                    let x = u64::from_be_bytes(a[n - k - 8..n - k].try_into().unwrap());
                    let y = u64::from_be_bytes(b[n - k - 8..n - k].try_into().unwrap());
                    let m = (x ^ y) | (x & AMBIG);
                    if m != 0 {
                        break k + (m.trailing_zeros() / 8) as usize;
                    }
                    k += 8;
                } else if k < n && a[n - 1 - k] < 4 && a[n - 1 - k] == b[n - 1 - k] {
                    k += 1;
                } else {
                    break k;
                }
            }
        } else {
            let a = &seq1[i1 as usize..i1 as usize + n];
            let b = &seq2.bytes[i2 as usize..i2 as usize + n];
            let mut k = 0usize;
            loop {
                if k + 8 <= n {
                    let x = u64::from_le_bytes(a[k..k + 8].try_into().unwrap());
                    let y = u64::from_le_bytes(b[k..k + 8].try_into().unwrap());
                    let m = (x ^ y) | (x & AMBIG);
                    if m != 0 {
                        break k + (m.trailing_zeros() / 8) as usize;
                    }
                    k += 8;
                } else if k < n && a[k] < 4 && a[k] == b[k] {
                    k += 1;
                } else {
                    break k;
                }
            }
        };
        i1 += matched as i32;
    } else if reverse {
        while i1 < len1 && i2 < len2 {
            let p = (len2 - 1 - i2) as usize;
            if seq1[(len1 - 1 - i1) as usize] == unpack(seq2.bytes, p) {
                i1 += 1;
                i2 += 1;
            } else {
                break;
            }
        }
    } else {
        let rem = seq2.rem as i32;
        while i1 < len1 && i2 < len2 {
            let p = (i2 + rem) as usize;
            if seq1[i1 as usize] == unpack(seq2.bytes, p) {
                i1 += 1;
                i2 += 1;
            } else {
                break;
            }
        }
    }
    i1 - start
}

/// One row of `last_seq2_off`: values for diagonal `k` at `vals[k - base]`.
#[derive(Debug, Clone, Default)]
struct Row {
    base: i32,
    vals: Vec<i32>,
}

impl Row {
    #[inline(always)]
    fn get(&self, k: i32) -> i32 {
        let i = k - self.base;
        if i < 0 || i as usize >= self.vals.len() {
            K_INVALID_OFFSET
        } else {
            self.vals[i as usize]
        }
    }
    #[inline(always)]
    fn set(&mut self, k: i32, v: i32) {
        let i = k - self.base;
        debug_assert!(
            i >= 0 && (i as usize) < self.vals.len(),
            "greedy row write out of range"
        );
        if i >= 0 && (i as usize) < self.vals.len() {
            self.vals[i as usize] = v;
        }
    }
}

/// `SGreedyAlignMem`: the greedy aligner's scratch.
pub struct GreedyAlignMem {
    pub max_dist: i32,
    pub xdrop: i32,
    max_score: Vec<i32>,
    /// rows[d]; rows 0 and 1 are full width, later rows compact (traceback) or reuse 0/1.
    rows: Vec<Row>,
}

impl GreedyAlignMem {
    /// `s_BlastGreedyAlignMemAlloc` for zero gap costs.
    pub fn new(reward: i32, penalty: i32, max_d: i32, xdrop: i32) -> GreedyAlignMem {
        let (reward2, penalty2) = if reward % 2 == 1 {
            (2 * reward, -2 * penalty)
        } else {
            (reward, -penalty)
        };
        let width = (max_d + max_d + 6) as usize;
        let d_diff = (xdrop + reward2 / 2) / (penalty2 + reward2) + 1;
        GreedyAlignMem {
            max_dist: max_d,
            xdrop,
            max_score: vec![0; (max_d + 1 + d_diff) as usize],
            rows: vec![
                Row {
                    base: 0,
                    vals: vec![0; width],
                },
                Row {
                    base: 0,
                    vals: vec![0; width],
                },
            ],
        }
    }

    /// The `max_subject_length` rule of `BLAST_GapAlignStructNew` -> initial `max_dist`.
    pub fn initial_max_dist(max_subject_length: i32) -> i32 {
        let m = max_subject_length.min(super::options::MAX_DBSEQ_LEN);
        GREEDY_MAX_COST.min(m / GREEDY_MAX_COST_FRACTION + 1)
    }
}

/// `BLAST_GreedyAlign`. Returns `best_dist` (or -1 when the aligner did not converge within
/// `max_dist`), and the aligned lengths of seq1/seq2. `match_cost`/`mismatch_cost` are the
/// (possibly doubled) scores.
#[allow(clippy::too_many_arguments)]
pub fn greedy_align(
    seq1: &[u8],
    len1: i32,
    seq2: Seq2,
    len2: i32,
    reverse: bool,
    xdrop_threshold: i32,
    match_cost: i32,
    mismatch_cost: i32,
    aux: &mut GreedyAlignMem,
    mut edit_block: Option<&mut PrelimEditBlock>,
    seed: &mut GreedySeed,
) -> (i32, i32, i32) {
    let max_dist = aux.max_dist;
    let diag_origin = max_dist + 2;
    let xdrop_offset = (xdrop_threshold + match_cost / 2) / (match_cost + mismatch_cost) + 1;
    let index = find_first_mismatch(seq1, seq2, len1, len2, 0, 0, reverse);
    let mut seq1_align_len = index;
    let mut seq2_align_len = index;
    let mut seq1_index = index;
    seed.start_q = 0;
    seed.start_s = 0;
    seed.match_length = index;
    let mut longest_match_run = index;
    if index == len1 || index == len2 {
        if let Some(eb) = edit_block.as_deref_mut() {
            eb.add(EditOp::Sub, index);
        }
        return (0, seq1_align_len, seq2_align_len);
    }
    let traceback = edit_block.is_some();
    // max_score[d] lives at index d + xdrop_offset
    let ms = &mut aux.max_score;
    let need = (max_dist + 1 + xdrop_offset) as usize;
    if ms.len() < need {
        ms.resize(need, 0);
    }
    for v in ms.iter_mut().take(xdrop_offset as usize) {
        *v = 0;
    }
    let xo = xdrop_offset as usize;
    // rows
    if aux.rows.len() < 2 {
        aux.rows.resize_with(2, Row::default);
    }
    let width = (max_dist + max_dist + 6) as usize;
    for r in aux.rows.iter_mut().take(2) {
        r.base = 0;
        if r.vals.len() < width {
            r.vals.resize(width, 0);
        }
    }
    aux.rows[0].set(diag_origin, seq1_index);
    ms[xo] = seq1_index * match_cost;
    let mut diag_lower = diag_origin - 1;
    let mut diag_upper = diag_origin + 1;
    let mut end1_reached = false;
    let mut end2_reached = false;
    let mut best_dist = 0;
    let mut best_diag = 0;
    let mut converged = false;
    // row index mapping: in score-only mode row d uses buffer d % 2
    let row_idx = |d: i32| -> usize {
        if traceback {
            d as usize
        } else {
            (d % 2) as usize
        }
    };
    let mut d = 1;
    while d <= max_dist {
        let mut curr_extent = 0;
        let mut curr_seq2_index = 0;
        let mut curr_diag = 0;
        let tmp_diag_lower = diag_lower;
        let tmp_diag_upper = diag_upper;
        let prev = row_idx(d - 1);
        let cur = row_idx(d);
        if traceback && aux.rows.len() <= cur {
            aux.rows.resize_with(cur + 1, Row::default);
        }
        if traceback && d >= 2 {
            // row d was allocated at the end of iteration d-1 with the bounds then current
            // (the C's s_GetMBSpace); make sure it covers [diag_lower-3, diag_upper+3].
            let r = &mut aux.rows[cur];
            let base = tmp_diag_lower - 3;
            let len = (tmp_diag_upper - tmp_diag_lower + 9) as usize;
            if r.vals.len() < len {
                r.vals.resize(len, K_INVALID_OFFSET);
            }
            r.base = base;
        }
        {
            let r = &mut aux.rows[prev];
            r.set(diag_lower - 1, K_INVALID_OFFSET);
            r.set(diag_lower, K_INVALID_OFFSET);
            r.set(diag_upper, K_INVALID_OFFSET);
            r.set(diag_upper + 1, K_INVALID_OFFSET);
        }
        let mut xdrop_score =
            ms[(d as usize + xo) - xo] + (match_cost + mismatch_cost) * d - xdrop_threshold;
        xdrop_score = (xdrop_score as f64 / (match_cost / 2) as f64).ceil() as i32;
        let mut k = tmp_diag_lower;
        while k <= tmp_diag_upper {
            let (p_k1, p_k, p_km1) = {
                let r = &aux.rows[prev];
                (r.get(k + 1), r.get(k), r.get(k - 1))
            };
            let mut seq2_index = p_k1.max(p_k) + 1;
            seq2_index = seq2_index.max(p_km1);
            seq1_index = seq2_index + k - diag_origin;
            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == diag_lower {
                    diag_lower += 1;
                } else {
                    aux.rows[cur].set(k, K_INVALID_OFFSET);
                }
                k += 1;
                continue;
            }
            diag_upper = k;
            let index =
                find_first_mismatch(seq1, seq2, len1, len2, seq1_index, seq2_index, reverse);
            if index > longest_match_run {
                seed.start_q = seq1_index;
                seed.start_s = seq2_index;
                longest_match_run = index;
                seed.match_length = index;
            }
            seq1_index += index;
            seq2_index += index;
            aux.rows[cur].set(k, seq2_index);
            if seq1_index + seq2_index > curr_extent {
                curr_extent = seq1_index + seq2_index;
                curr_seq2_index = seq2_index;
                curr_diag = k;
            }
            if seq2_index == len2 {
                diag_lower = k + 1;
                end2_reached = true;
            }
            if seq1_index == len1 {
                diag_upper = k - 1;
                end1_reached = true;
            }
            k += 1;
        }
        let curr_score = curr_extent * (match_cost / 2) - d * (match_cost + mismatch_cost);
        if curr_score >= ms[(d - 1) as usize + xo] {
            ms[d as usize + xo] = curr_score;
            best_dist = d;
            best_diag = curr_diag;
            seq2_align_len = curr_seq2_index;
            seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        } else {
            ms[d as usize + xo] = ms[(d - 1) as usize + xo];
        }
        if diag_lower > diag_upper {
            converged = true;
            break;
        }
        if !end2_reached {
            diag_lower -= 1;
        }
        if !end1_reached {
            diag_upper += 1;
        }
        if traceback {
            // allocate row d+1 (the C's s_GetMBSpace), sized for the new bounds
            let nxt = (d + 1) as usize;
            if aux.rows.len() <= nxt {
                aux.rows.resize_with(nxt + 1, Row::default);
            }
            let r = &mut aux.rows[nxt];
            let base = diag_lower - 3;
            let len = (diag_upper - diag_lower + 9) as usize;
            r.vals.clear();
            r.vals.resize(len, K_INVALID_OFFSET);
            r.base = base;
        }
        d += 1;
    }
    if !converged {
        return (-1, seq1_align_len, seq2_align_len);
    }
    let Some(eb) = edit_block else {
        return (best_dist, seq1_align_len, seq2_align_len);
    };
    let mut d = best_dist;
    let mut seq2_index = seq2_align_len;
    while d > 0 {
        // s_GetNextNonAffineTback
        let r = &aux.rows[(d - 1) as usize];
        let (a, b, c) = (r.get(best_diag - 1), r.get(best_diag), r.get(best_diag + 1));
        let (new_diag, new_seq2_index) = if a > b.max(c) {
            (best_diag - 1, a)
        } else if b > c {
            (best_diag, b)
        } else {
            (best_diag + 1, c)
        };
        if new_diag == best_diag {
            if seq2_index - new_seq2_index > 0 {
                eb.add(EditOp::Sub, seq2_index - new_seq2_index);
            }
        } else if new_diag < best_diag {
            if seq2_index - new_seq2_index > 0 {
                eb.add(EditOp::Sub, seq2_index - new_seq2_index);
            }
            eb.add(EditOp::Ins, 1);
        } else {
            if seq2_index - new_seq2_index - 1 > 0 {
                eb.add(EditOp::Sub, seq2_index - new_seq2_index - 1);
            }
            eb.add(EditOp::Del, 1);
        }
        d -= 1;
        best_diag = new_diag;
        seq2_index = new_seq2_index;
    }
    eb.add(EditOp::Sub, aux.rows[0].get(diag_origin));
    (best_dist, seq1_align_len, seq2_align_len)
}

/// `Blast_PrelimEditBlockToGapEditScript`.
pub fn prelim_to_edit_script(rev: &PrelimEditBlock, fwd: &PrelimEditBlock) -> EditScript {
    let merge_ops = !fwd.ops.is_empty()
        && !rev.ops.is_empty()
        && fwd.ops[fwd.ops.len() - 1].0 == rev.ops[rev.ops.len() - 1].0;
    let mut ops: Vec<(EditOp, i32)> = Vec::with_capacity(rev.ops.len() + fwd.ops.len());
    ops.extend_from_slice(&rev.ops);
    if fwd.ops.is_empty() {
        return EditScript { ops };
    }
    let mut i = fwd.ops.len() as i32 - 1;
    if merge_ops {
        let last = ops.len() - 1;
        ops[last].1 += fwd.ops[fwd.ops.len() - 1].1;
        i = fwd.ops.len() as i32 - 2;
    }
    while i >= 0 {
        ops.push(fwd.ops[i as usize]);
        i -= 1;
    }
    EditScript { ops }
}

/// `s_UpdateEditScript`.
fn update_edit_script(esp: &mut EditScript, pos: usize, bf: i32, af: i32) {
    let n = esp.ops.len() as i32;
    if bf > 0 {
        let mut op = pos as i32;
        let (mut qd, mut sd) = (bf, bf);
        loop {
            op -= 1;
            if op < 0 {
                return;
            }
            let (t, num) = esp.ops[op as usize];
            match t {
                EditOp::Sub => {
                    qd -= num;
                    sd -= num;
                }
                EditOp::Ins => qd -= num,
                EditOp::Del => sd -= num,
            }
            if !(qd > 0 || sd > 0) {
                break;
            }
        }
        esp.ops[op as usize].1 = -qd.max(sd);
        esp.ops[op as usize].0 = EditOp::Sub;
        op += 1;
        while op < pos as i32 - 1 {
            esp.ops[op as usize].1 = 0;
            op += 1;
        }
        esp.ops[pos].1 += bf;
        qd -= sd;
        esp.ops[pos - 1].0 = if qd > 0 { EditOp::Del } else { EditOp::Ins };
        esp.ops[pos - 1].1 = if qd > 0 { qd } else { -qd };
    }
    if af > 0 {
        let mut op = pos as i32;
        let (mut qd, mut sd) = (af, af);
        loop {
            op += 1;
            if op >= n {
                return;
            }
            let (t, num) = esp.ops[op as usize];
            match t {
                EditOp::Sub => {
                    qd -= num;
                    sd -= num;
                }
                EditOp::Ins => qd -= num,
                EditOp::Del => sd -= num,
            }
            if !(qd > 0 || sd > 0) {
                break;
            }
        }
        esp.ops[op as usize].1 = -qd.max(sd);
        esp.ops[op as usize].0 = EditOp::Sub;
        op -= 1;
        while op > pos as i32 + 1 {
            esp.ops[op as usize].1 = 0;
            op -= 1;
        }
        esp.ops[pos].1 += af;
        qd -= sd;
        esp.ops[pos + 1].0 = if qd > 0 { EditOp::Del } else { EditOp::Ins };
        esp.ops[pos + 1].1 = if qd > 0 { qd } else { -qd };
    }
}

/// `s_RebuildEditScript`.
fn rebuild_edit_script(esp: &mut EditScript) {
    let n = esp.ops.len();
    let mut j: i32 = -1;
    for i in 0..n {
        if esp.ops[i].1 == 0 {
            continue;
        }
        if j >= 0 && esp.ops[i].0 == esp.ops[j as usize].0 {
            esp.ops[j as usize].1 += esp.ops[i].1;
        } else if j == -1 || esp.ops[i].0 == EditOp::Sub || esp.ops[j as usize].0 == EditOp::Sub {
            j += 1;
            esp.ops[j as usize] = esp.ops[i];
        } else {
            let d = esp.ops[j as usize].1 - esp.ops[i].1;
            if d > 0 {
                esp.ops[(j - 1) as usize].1 += esp.ops[i].1;
                esp.ops[j as usize].1 = d;
            } else if d < 0 {
                if j == 0 && i as i32 - j > 0 {
                    esp.ops[j as usize].0 = EditOp::Sub;
                    j += 1;
                } else {
                    esp.ops[(j - 1) as usize].1 += esp.ops[j as usize].1;
                }
                esp.ops[j as usize].1 = -d;
                esp.ops[j as usize].0 = esp.ops[i].0;
            } else {
                esp.ops[(j - 1) as usize].1 += esp.ops[j as usize].1;
                j -= 1;
            }
        }
    }
    esp.ops.truncate((j + 1) as usize);
}

/// `s_ReduceGaps`: `q`/`s` are the aligned query and subject letters from the alignment
/// start; `qf`/`sf` their ends (exclusive), as slices of the same buffers.
pub fn reduce_gaps(
    esp: &mut EditScript,
    q: &[u8],
    s: &[u8],
    q_start: usize,
    s_start: usize,
    qf: usize,
    sf: usize,
) {
    // pointer arithmetic on indices into q / s
    let mut q1 = q_start as isize;
    let mut s1 = s_start as isize;
    let n = esp.ops.len();
    for i in 0..n {
        let (op, num) = esp.ops[i];
        if num == 0 {
            continue;
        }
        match op {
            EditOp::Sub => {
                if num >= 12 {
                    let mut nm1 = 1isize;
                    if i > 0 {
                        while q1 - nm1 >= q_start as isize
                            && q[(q1 - nm1) as usize] == s[(s1 - nm1) as usize]
                        {
                            nm1 += 1;
                        }
                    }
                    q1 += num as isize;
                    s1 += num as isize;
                    let mut nm2 = 0i32;
                    if i < n - 1 {
                        while (q1 + 1) < qf as isize
                            && (s1 + 1) < sf as isize
                            && q[q1 as usize] == s[s1 as usize]
                        {
                            q1 += 1;
                            s1 += 1;
                            nm2 += 1;
                        }
                        // the C's `(*(q1++) == *(s1++))` advances once more on the failing compare
                        if (q1 + 1) < qf as isize && (s1 + 1) < sf as isize {
                            q1 += 1;
                            s1 += 1;
                        }
                    }
                    if nm1 > 1 || nm2 > 0 {
                        update_edit_script(esp, i, nm1 as i32 - 1, nm2);
                    }
                    q1 -= 1;
                    s1 -= 1;
                } else {
                    q1 += num as isize;
                    s1 += num as isize;
                }
            }
            EditOp::Ins => q1 += num as isize,
            EditOp::Del => s1 += num as isize,
        }
    }
    rebuild_edit_script(esp);
    let mut qp = q_start as isize;
    let mut sp = s_start as isize;
    let n = esp.ops.len();
    let mut i = 0;
    while i < n {
        let (op, num) = esp.ops[i];
        if op == EditOp::Sub {
            qp += num as isize;
            sp += num as isize;
            i += 1;
            continue;
        }
        if i > 1 && esp.ops[i].0 != esp.ops[i - 2].0 && esp.ops[i - 2].1 > 0 {
            let d = esp.ops[i].1 + esp.ops[i - 1].1 + esp.ops[i - 2].1;
            if d == 3 {
                esp.ops[i - 2].1 = 0;
                esp.ops[i - 1].1 = 2;
                esp.ops[i].1 = 0;
                if esp.ops[i].0 == EditOp::Ins {
                    qp += 1;
                } else {
                    sp += 1;
                }
            } else if d < 12 {
                let mut nm1 = 0;
                let mut nm2 = 0;
                let d = esp.ops[i].1.min(esp.ops[i - 2].1);
                qp -= esp.ops[i - 1].1 as isize;
                sp -= esp.ops[i - 1].1 as isize;
                let q1 = qp;
                let s1 = sp;
                if esp.ops[i].0 == EditOp::Ins {
                    sp -= d as isize;
                } else {
                    qp -= d as isize;
                }
                let (mut qq1, mut ss1) = (q1, s1);
                for _ in 0..esp.ops[i - 1].1 {
                    if q[qq1 as usize] == s[ss1 as usize] {
                        nm1 += 1;
                    }
                    if q[qp as usize] == s[sp as usize] {
                        nm2 += 1;
                    }
                    qq1 += 1;
                    ss1 += 1;
                    qp += 1;
                    sp += 1;
                }
                for _ in 0..d {
                    if q[qp as usize] == s[sp as usize] {
                        nm2 += 1;
                    }
                    qp += 1;
                    sp += 1;
                }
                if nm2 >= nm1 - d {
                    esp.ops[i - 2].1 -= d;
                    esp.ops[i - 1].1 += d;
                    esp.ops[i].1 -= d;
                } else {
                    // the C's q1/s1 were advanced by the j loop: back to the entry position
                    qp = qq1;
                    sp = ss1;
                }
            }
        }
        let (op, num) = esp.ops[i];
        if op == EditOp::Ins {
            qp += num as isize;
        } else {
            sp += num as isize;
        }
        i += 1;
    }
    rebuild_edit_script(esp);
}

/// `BlastGapAlignStruct` for the greedy path plus the scoring it needs.
pub struct GreedyAligner {
    pub mem: GreedyAlignMem,
    pub gap_x_dropoff: i32,
    pub reward: i32,
    pub penalty: i32,
    pub fwd: PrelimEditBlock,
    pub rev: PrelimEditBlock,
}

/// What `BLAST_GreedyGappedAlignment` fills in `gap_align`.
#[derive(Debug, Clone, Default)]
pub struct GreedyResult {
    pub query_start: i32,
    pub query_stop: i32,
    pub subject_start: i32,
    pub subject_stop: i32,
    pub seed_q: i32,
    pub seed_s: i32,
    pub score: i32,
    pub edit_script: Option<EditScript>,
}

impl GreedyAligner {
    /// `BLAST_GapAlignStructNew` + `s_BlastGreedyAlignMemAlloc(score_params, ext_params, ..)`.
    pub fn new(
        reward: i32,
        penalty: i32,
        gap_x_dropoff: i32,
        gap_x_dropoff_final: i32,
        max_subject_length: i32,
    ) -> GreedyAligner {
        let xdrop = if reward % 2 == 1 {
            2 * gap_x_dropoff.max(gap_x_dropoff_final)
        } else {
            gap_x_dropoff.max(gap_x_dropoff_final)
        };
        let max_d = GreedyAlignMem::initial_max_dist(max_subject_length);
        GreedyAligner {
            mem: GreedyAlignMem::new(reward, penalty, max_d, xdrop),
            gap_x_dropoff,
            reward,
            penalty,
            fwd: PrelimEditBlock::default(),
            rev: PrelimEditBlock::default(),
        }
    }

    fn realloc_mem(&mut self) {
        // BLAST_GreedyGappedAlignment: double the max distance and try again
        let new_dist = self.mem.max_dist * 2;
        self.mem = GreedyAlignMem::new(self.reward, self.penalty, new_dist, self.mem.xdrop);
    }

    /// One direction of `BLAST_GreedyAlign`, with the doubled scores for an odd reward.
    #[allow(clippy::too_many_arguments)]
    fn extend(
        &mut self,
        seq1: &[u8],
        len1: i32,
        seq2: Seq2,
        len2: i32,
        reverse: bool,
        do_traceback: bool,
        seed: &mut GreedySeed,
    ) -> (i32, i32, i32) {
        let x = self.gap_x_dropoff;
        let (mc, mmc, xd) = if self.reward % 2 == 1 {
            (2 * self.reward, -2 * self.penalty, 2 * x)
        } else {
            (self.reward, -self.penalty, x)
        };
        let eb = if do_traceback {
            Some(if reverse {
                &mut self.rev
            } else {
                &mut self.fwd
            })
        } else {
            None
        };
        greedy_align(
            seq1,
            len1,
            seq2,
            len2,
            reverse,
            xd,
            mc,
            mmc,
            &mut self.mem,
            eb,
            seed,
        )
    }

    /// `BLAST_GreedyGappedAlignment`. `subject` is packed ncbi2na when `compressed`, else
    /// blastna; `query` is one context's letters.
    #[allow(clippy::too_many_arguments)]
    pub fn align(
        &mut self,
        query: &[u8],
        subject: &[u8],
        query_length: i32,
        subject_length: i32,
        q_off: i32,
        s_off: i32,
        compressed: bool,
        do_traceback: bool,
    ) -> GreedyResult {
        let q_avail = query_length - q_off;
        let s_avail = subject_length - s_off;
        let (s2_fwd, rem) = if compressed {
            (
                Seq2 {
                    bytes: &subject[(s_off / 4) as usize..],
                    rem: (s_off % 4) as u8,
                },
                (s_off % 4) as u8,
            )
        } else {
            (
                Seq2 {
                    bytes: &subject[s_off as usize..],
                    rem: 4,
                },
                4u8,
            )
        };
        if do_traceback {
            self.fwd.reset();
            self.rev.reset();
        }
        let mut fwd_seed = GreedySeed::default();
        let mut rev_seed = GreedySeed::default();
        let q_fwd = &query[q_off as usize..];
        let (mut score, mut q_ext_r, mut s_ext_r);
        loop {
            let r = self.extend(
                q_fwd,
                q_avail,
                s2_fwd,
                s_avail,
                false,
                do_traceback,
                &mut fwd_seed,
            );
            score = r.0;
            q_ext_r = r.1;
            s_ext_r = r.2;
            if score >= 0 {
                break;
            }
            self.realloc_mem();
        }
        let rem_rev = if compressed { 0 } else { rem };
        let s2_rev = Seq2 {
            bytes: subject,
            rem: rem_rev,
        };
        let (mut q_ext_l, mut s_ext_l);
        loop {
            let r = self.extend(
                query,
                q_off,
                s2_rev,
                s_off,
                true,
                do_traceback,
                &mut rev_seed,
            );
            q_ext_l = r.1;
            s_ext_l = r.2;
            if r.0 >= 0 {
                score += r.0;
                break;
            }
            self.realloc_mem();
        }
        // distance -> score (linear gaps)
        score = (q_ext_r + s_ext_r + q_ext_l + s_ext_l) * self.reward / 2
            - score * (self.reward - self.penalty);
        let mut out = GreedyResult {
            query_start: q_off - q_ext_l,
            query_stop: q_off + q_ext_r,
            subject_start: s_off - s_ext_l,
            subject_stop: s_off + s_ext_r,
            seed_q: q_off,
            seed_s: s_off,
            score,
            edit_script: None,
        };
        if do_traceback {
            let mut esp = prelim_to_edit_script(&self.rev, &self.fwd);
            if !esp.ops.is_empty() {
                // s_ReduceGaps on the uncompressed subject
                debug_assert!(!compressed, "traceback needs blastna subject");
                reduce_gaps(
                    &mut esp,
                    query,
                    subject,
                    (q_off - q_ext_l) as usize,
                    (s_off - s_ext_l) as usize,
                    (q_off + q_ext_r) as usize,
                    (s_off + s_ext_r) as usize,
                );
            }
            out.edit_script = Some(esp);
        } else {
            let q_box_l = q_off - q_ext_l;
            let s_box_l = s_off - s_ext_l;
            let q_box_r = q_off + q_ext_r;
            let s_box_r = s_off + s_ext_r;
            let mut q_seed_start_l = q_off - rev_seed.start_q;
            let mut s_seed_start_l = s_off - rev_seed.start_s;
            let mut q_seed_start_r = q_off + fwd_seed.start_q;
            let mut s_seed_start_r = s_off + fwd_seed.start_s;
            let mut valid_seed_len_l = 0;
            let mut valid_seed_len_r = 0;
            if q_seed_start_r < q_box_r && s_seed_start_r < s_box_r {
                valid_seed_len_r = (q_box_r - q_seed_start_r).min(s_box_r - s_seed_start_r);
                valid_seed_len_r = valid_seed_len_r.min(fwd_seed.match_length) / 2;
            } else {
                q_seed_start_r = q_off;
                s_seed_start_r = s_off;
            }
            if q_seed_start_l > q_box_l && s_seed_start_l > s_box_l {
                valid_seed_len_l = (q_seed_start_l - q_box_l).min(s_seed_start_l - s_box_l);
                valid_seed_len_l = valid_seed_len_l.min(rev_seed.match_length) / 2;
            } else {
                q_seed_start_l = q_off;
                s_seed_start_l = s_off;
            }
            if valid_seed_len_r > valid_seed_len_l {
                out.seed_q = q_seed_start_r + valid_seed_len_r;
                out.seed_s = s_seed_start_r + valid_seed_len_r;
            } else {
                out.seed_q = q_seed_start_l - valid_seed_len_l;
                out.seed_s = s_seed_start_l - valid_seed_len_l;
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pack(letters: &[u8]) -> Vec<u8> {
        let mut v = vec![0u8; letters.len() / 4 + 1 + 4];
        for (i, &b) in letters.iter().enumerate() {
            v[i / 4] |= b << (2 * (3 - (i % 4)));
        }
        v
    }

    #[test]
    fn identical_sequences_extend_fully() {
        let q: Vec<u8> = (0..200u32).map(|i| ((i * 7 + i / 3) % 4) as u8).collect();
        let s = pack(&q);
        let mut g = GreedyAligner::new(1, -2, 13, 54, 1000);
        let r = g.align(&q, &s, 200, 200, 100, 100, true, false);
        assert_eq!(
            (r.query_start, r.query_stop, r.subject_start, r.subject_stop),
            (0, 200, 0, 200)
        );
        assert_eq!(r.score, 200);
        // traceback on the blastna subject
        let r = g.align(&q, &q, 200, 200, 100, 100, false, true);
        assert_eq!(r.score, 200);
        assert_eq!(r.edit_script.unwrap().ops, vec![(EditOp::Sub, 200)]);
    }

    #[test]
    fn mismatch_and_gap_are_scored() {
        let mut q: Vec<u8> = (0..120u32).map(|i| ((i * 5 + i / 7) % 4) as u8).collect();
        let s_letters = q.clone();
        q[60] = (q[60] + 1) % 4; // one mismatch
        let mut g = GreedyAligner::new(1, -2, 13, 54, 1000);
        let r = g.align(&q, &s_letters, 120, 120, 20, 20, false, true);
        assert_eq!(r.score, 119 - 2);
        let esp = r.edit_script.unwrap();
        assert_eq!(esp.ops.iter().map(|o| o.1).sum::<i32>(), 120);
        // a deletion in the query
        let mut q2 = s_letters.clone();
        q2.remove(60);
        let r = g.align(&q2, &s_letters, 119, 120, 20, 20, false, true);
        let esp = r.edit_script.unwrap();
        assert!(
            esp.ops.iter().any(|o| o.0 == EditOp::Del && o.1 == 1),
            "{esp:?}"
        );
        // score: 119 matches, one gap of length 1: linear cost = reward/2 + penalty... = (119+120)/2 - 1*(1+2) = 119 - 3 = 116
        assert_eq!(r.score, 116);
    }

    #[test]
    fn xdrop_stops_extension() {
        let q: Vec<u8> = (0..100u32).map(|i| (i % 4) as u8).collect();
        let mut s = q.clone();
        let mut x: u32 = 12345;
        for b in s.iter_mut().skip(70) {
            x = x.wrapping_mul(1103515245).wrapping_add(12345);
            *b = (*b + 1 + ((x >> 16) % 3) as u8) % 4; // 30 mismatches at the end, aperiodic
        }
        let mut g = GreedyAligner::new(1, -2, 13, 54, 1000);
        let r = g.align(&q, &s, 100, 100, 10, 10, false, false);
        assert_eq!(r.query_stop, 70);
        assert_eq!(r.query_start, 0);
        assert_eq!(r.score, 70);
    }

    /// `find_first_mismatch` on unpacked sequences (the traceback's case) equals the C's
    /// byte loop: equal letters, query letter a base, in both directions, at offsets and
    /// lengths that cross the eight-byte words.
    #[test]
    fn find_first_mismatch_words_equal_the_byte_loop() {
        fn naive(s1: &[u8], s2: &[u8], l1: i32, l2: i32, i1: i32, i2: i32, rev: bool) -> i32 {
            let (mut a, mut b) = (i1, i2);
            while a < l1 && b < l2 {
                let (x, y) = if rev {
                    (s1[(l1 - 1 - a) as usize], s2[(l2 - 1 - b) as usize])
                } else {
                    (s1[a as usize], s2[b as usize])
                };
                if x < 4 && x == y {
                    a += 1;
                    b += 1;
                } else {
                    break;
                }
            }
            a - i1
        }
        let mut x = 12345u64;
        let mut next = || {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (x >> 33) as usize
        };
        for _ in 0..400 {
            let l1 = 1 + next() % 70;
            let l2 = 1 + next() % 70;
            // mostly matching letters, a few mismatches and ambiguity codes (N = 14)
            let base: Vec<u8> = (0..80).map(|_| (next() % 4) as u8).collect();
            let s1: Vec<u8> = (0..l1)
                .map(|i| match next() % 23 {
                    0 => 14,
                    1 => (base[i] + 1) % 4,
                    _ => base[i],
                })
                .collect();
            let s2: Vec<u8> = (0..l2)
                .map(|i| {
                    if next() % 29 == 0 {
                        (base[i] + 2) % 4
                    } else {
                        base[i]
                    }
                })
                .collect();
            for rev in [false, true] {
                for i1 in 0..l1.min(10) {
                    for i2 in 0..l2.min(10) {
                        let want = naive(&s1, &s2, l1 as i32, l2 as i32, i1 as i32, i2 as i32, rev);
                        let got = find_first_mismatch(
                            &s1,
                            Seq2 { bytes: &s2, rem: 4 },
                            l1 as i32,
                            l2 as i32,
                            i1 as i32,
                            i2 as i32,
                            rev,
                        );
                        assert_eq!(got, want, "l1 {l1} l2 {l2} i1 {i1} i2 {i2} rev {rev}");
                    }
                }
            }
        }
    }
}
