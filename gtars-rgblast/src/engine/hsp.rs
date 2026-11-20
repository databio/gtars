//! HSPs, HSP lists, hit lists and results, with every comparator ported verbatim.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/core/blast_hits.c` (`Blast_HSPInit`, `ScoreCompareHSPs`,
//! `Blast_HSPListSortByScore`, `s_EvalueComp`, `s_EvalueCompareHSPs`,
//! `Blast_HSPListSortByEvalue`, `s_BlastMergeTwoHSPs`, `Blast_HSPListNew`, `s_Heapify`,
//! `s_CreateHeap`, `s_BlastHSPListInsertHSPInHeap`, `Blast_HSPListSaveHSP`,
//! `Blast_HSPListGetEvalues`, `Blast_HSPListGetBitScores`, `Blast_HSPListReapByEvalue`,
//! `Blast_HSPListReapByQueryCoverage`, `Blast_TrimHSPListByMaxHsps`,
//! `Blast_HSPListPurgeNullHSPs`, `s_QueryOffsetCompareHSPs`, `s_QueryEndCompareHSPs`,
//! `s_CutOffGapEditScript`, `Blast_HSPListPurgeHSPsWithCommonEndpoints`,
//! `s_BlastHSPListsCombineByScore`, `Blast_HSPListAppend`, `Blast_HSPListsMerge`,
//! `Blast_HSPListAdjustOffsets`, `s_EvalueCompareHSPLists`, `Blast_HitListNew`,
//! `s_BlastHitListPurge`, `s_BlastHitListInsertHSPListInHeap`, `Blast_HitListUpdate`,
//! `Blast_HSPResultsSortByEvalue`, `Blast_HSPResultsInsertHSPList`,
//! `Blast_HSPReevaluateWithAmbiguitiesGapped`, `s_UpdateReevaluatedHSP`,
//! `s_Blast_HSPGetNumIdentitiesAndPositives`, `Blast_HSPTestIdentityAndLength`,
//! `Blast_HSPGetQueryCoverage`, `Blast_HSPCalcLengthAndGaps`, `GetPrelimHitlistSize`) and
//! `blast_engine.c::s_Blast_HSPListReapByPrelimEvalue`.
//!
//! The C sorts with `qsort`; Rust uses stable sorts and asserts (in debug builds) that no two
//! distinct elements compare equal, so the order cannot depend on the sort algorithm.

use super::options::{NCBIMATH_LN2, Parameters};
use super::stat::{KarlinBlk, ScoreBlk, karlin_stoe_simple};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EditOp {
    /// `eGapAlignDel`: gap in the query (letters only in the subject).
    Del = 0,
    /// `eGapAlignSub`
    Sub = 3,
    /// `eGapAlignIns`: gap in the subject (letters only in the query).
    Ins = 6,
}

/// `GapEditScript`: parallel op/num arrays.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct EditScript {
    pub ops: Vec<(EditOp, i32)>,
}

impl EditScript {
    pub fn size(&self) -> usize {
        self.ops.len()
    }
}

/// `BlastSeg`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Seg {
    pub frame: i32,
    pub offset: i32,
    pub end: i32,
    pub gapped_start: i32,
}

/// `BlastHSP`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Hsp {
    pub score: i32,
    pub num_ident: i32,
    pub bit_score: f64,
    pub evalue: f64,
    pub query: Seg,
    pub subject: Seg,
    pub context: i32,
    pub gap_info: Option<EditScript>,
}

impl Hsp {
    /// `Blast_HSPInit`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        query_start: i32,
        query_end: i32,
        subject_start: i32,
        subject_end: i32,
        query_gapped_start: i32,
        subject_gapped_start: i32,
        context: i32,
        query_frame: i32,
        subject_frame: i32,
        score: i32,
        gap_info: Option<EditScript>,
    ) -> Hsp {
        Hsp {
            score,
            num_ident: 0,
            bit_score: 0.0,
            evalue: 0.0,
            query: Seg {
                frame: query_frame,
                offset: query_start,
                end: query_end,
                gapped_start: query_gapped_start,
            },
            subject: Seg {
                frame: subject_frame,
                offset: subject_start,
                end: subject_end,
                gapped_start: subject_gapped_start,
            },
            context,
            gap_info,
        }
    }
}

/// `ScoreCompareHSPs`.
pub fn score_compare_hsps(a: &Hsp, b: &Hsp) -> Ordering {
    b.score
        .cmp(&a.score)
        .then(a.subject.offset.cmp(&b.subject.offset))
        .then(b.subject.end.cmp(&a.subject.end))
        .then(a.query.offset.cmp(&b.query.offset))
        .then(b.query.end.cmp(&a.query.end))
}

/// `s_EvalueComp`.
pub fn evalue_comp(e1: f64, e2: f64) -> Ordering {
    const EPSILON: f64 = 1.0e-180;
    if e1 < EPSILON && e2 < EPSILON {
        return Ordering::Equal;
    }
    if e1 < e2 {
        Ordering::Less
    } else if e1 > e2 {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

/// `s_EvalueCompareHSPs`.
pub fn evalue_compare_hsps(a: &Hsp, b: &Hsp) -> Ordering {
    evalue_comp(a.evalue, b.evalue).then_with(|| score_compare_hsps(a, b))
}

/// `s_QueryOffsetCompareHSPs`.
fn query_offset_compare(a: &Hsp, b: &Hsp) -> Ordering {
    a.context
        .cmp(&b.context)
        .then(a.query.offset.cmp(&b.query.offset))
        .then(a.subject.offset.cmp(&b.subject.offset))
        .then(b.score.cmp(&a.score))
        .then(b.query.end.cmp(&a.query.end))
        .then(b.subject.end.cmp(&a.subject.end))
}

/// `s_QueryEndCompareHSPs`.
fn query_end_compare(a: &Hsp, b: &Hsp) -> Ordering {
    a.context
        .cmp(&b.context)
        .then(a.query.end.cmp(&b.query.end))
        .then(a.subject.end.cmp(&b.subject.end))
        .then(b.score.cmp(&a.score))
        .then(b.query.offset.cmp(&a.query.offset))
        .then(b.subject.offset.cmp(&a.subject.offset))
}

/// `BlastHSPList`. `hsps` holds only live HSPs (the C's NULL slots are purged eagerly where
/// the C would purge them at the next step that matters).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct HspList {
    pub oid: i32,
    pub query_index: i32,
    pub hsps: Vec<Hsp>,
    pub best_evalue: f64,
}

impl HspList {
    pub fn hspcnt(&self) -> i32 {
        self.hsps.len() as i32
    }

    /// `Blast_HSPListSortByScore`.
    pub fn sort_by_score(&mut self) {
        if self.hsps.len() <= 1 {
            return;
        }
        if !self
            .hsps
            .windows(2)
            .all(|w| score_compare_hsps(&w[0], &w[1]) != Ordering::Greater)
        {
            self.hsps.sort_by(score_compare_hsps);
        }
    }

    /// `Blast_HSPListSortByEvalue`.
    pub fn sort_by_evalue(&mut self) {
        if self.hsps.len() <= 1 {
            return;
        }
        if !self
            .hsps
            .windows(2)
            .all(|w| evalue_compare_hsps(&w[0], &w[1]) != Ordering::Greater)
        {
            self.hsps.sort_by(evalue_compare_hsps);
        }
    }

    /// `s_BlastGetBestEvalue`.
    pub fn best_evalue(&self) -> f64 {
        let mut best = i32::MAX as f64;
        for h in &self.hsps {
            best = best.min(h.evalue);
        }
        best
    }

    /// `Blast_HSPListGetEvalues` (no RPS, no gap decay) + `best_evalue`.
    pub fn get_evalues(
        &mut self,
        kbp: &[Option<KarlinBlk>],
        eff_searchsp: &[i64],
        round_down: bool,
    ) {
        for h in self.hsps.iter_mut() {
            let k = kbp[h.context as usize]
                .as_ref()
                .expect("Karlin block for context");
            let mut score = h.score;
            if round_down {
                score &= !1;
            }
            h.evalue = karlin_stoe_simple(score, k, eff_searchsp[h.context as usize]);
        }
        self.best_evalue = self.best_evalue();
    }

    /// `Blast_HSPListGetBitScores`.
    pub fn get_bit_scores(&mut self, kbp: &[Option<KarlinBlk>]) {
        for h in self.hsps.iter_mut() {
            let k = kbp[h.context as usize]
                .as_ref()
                .expect("Karlin block for context");
            h.bit_score = (h.score as f64 * k.lambda - k.log_k) / NCBIMATH_LN2;
        }
    }

    /// `Blast_HSPListReapByEvalue` / `s_Blast_HSPListReapByPrelimEvalue`.
    pub fn reap_by_evalue(&mut self, cutoff: f64) {
        self.hsps.retain(|h| h.evalue <= cutoff);
    }

    /// `Blast_TrimHSPListByMaxHsps`.
    pub fn trim_by_max_hsps(&mut self, max_hsps_per_subject: i32) {
        if max_hsps_per_subject == 0 || self.hspcnt() <= max_hsps_per_subject {
            return;
        }
        self.hsps.truncate(max_hsps_per_subject as usize);
    }

    /// `Blast_HSPListAdjustOffsets`.
    pub fn adjust_offsets(&mut self, offset: i32) {
        if offset == 0 {
            return;
        }
        for h in self.hsps.iter_mut() {
            h.subject.offset += offset;
            h.subject.end += offset;
            h.subject.gapped_start += offset;
        }
    }

    /// `Blast_HSPListPurgeHSPsWithCommonEndpoints`. With `purge`, duplicates are removed and
    /// the live count is returned; without it (blastn traceback), the later HSP of a pair with
    /// a common start or end is cut back instead and moved after the first `n` entries, where
    /// `n` is returned.
    pub fn purge_common_endpoints(&mut self, purge: bool) -> usize {
        if self.hsps.is_empty() {
            return 0;
        }
        self.hsps.sort_by(query_offset_compare);
        let mut dropped: Vec<Hsp> = Vec::new(); // the cut HSPs kept after the live range
        let mut hsp_count = self.hsps.len();
        let mut i = 0;
        while i < hsp_count {
            let j = 1;
            while i + j < hsp_count
                && self.hsps[i].context == self.hsps[i + j].context
                && self.hsps[i].query.offset == self.hsps[i + j].query.offset
                && self.hsps[i].subject.offset == self.hsps[i + j].subject.offset
                && self.hsps[i].subject.frame == self.hsps[i + j].subject.frame
            {
                hsp_count -= 1;
                let mut hsp = self.hsps.remove(i + j);
                if !purge && hsp.query.end > self.hsps[i].query.end {
                    cut_off_edit_script(
                        &mut hsp,
                        self.hsps[i].query.end,
                        self.hsps[i].subject.end,
                        true,
                    );
                    dropped.push(hsp);
                }
            }
            i += j;
        }
        self.hsps.sort_by(query_end_compare);
        i = 0;
        while i < hsp_count {
            let j = 1;
            while i + j < hsp_count
                && self.hsps[i].context == self.hsps[i + j].context
                && self.hsps[i].query.end == self.hsps[i + j].query.end
                && self.hsps[i].subject.end == self.hsps[i + j].subject.end
                && self.hsps[i].subject.frame == self.hsps[i + j].subject.frame
            {
                hsp_count -= 1;
                let mut hsp = self.hsps.remove(i + j);
                if !purge && hsp.query.offset < self.hsps[i].query.offset {
                    cut_off_edit_script(
                        &mut hsp,
                        self.hsps[i].query.offset,
                        self.hsps[i].subject.offset,
                        false,
                    );
                    dropped.push(hsp);
                }
            }
            i += j;
        }
        if !purge {
            // The C leaves the cut HSPs after the live range: each is written to
            // `hsp_array[hsp_count]` after the count was decremented, so it lands just after
            // the live range and *below* the ones moved before it. The moved HSPs therefore
            // sit in reverse order of moving (the last moved first). The order matters when
            // two cut HSPs re-evaluate to the same box and score: the later stable sorts keep
            // the first and the TRUE purge drops the other.
            self.hsps.extend(dropped.into_iter().rev());
        }
        hsp_count
    }
}

/// `s_CutOffGapEditScript`.
fn cut_off_edit_script(hsp: &mut Hsp, q_cut: i32, s_cut: i32, cut_begin: bool) {
    let Some(esp) = hsp.gap_info.as_mut() else {
        return;
    };
    let q_cut = q_cut - hsp.query.offset;
    let s_cut = s_cut - hsp.subject.offset;
    let (mut qid, mut sid, mut opid) = (0, 0, 0);
    let mut found = false;
    let mut index = 0usize;
    while index < esp.ops.len() {
        let (op, num) = esp.ops[index];
        opid = 0;
        while opid < num {
            match op {
                EditOp::Sub => {
                    qid += 1;
                    sid += 1;
                    opid += 1;
                }
                EditOp::Del => {
                    sid += num;
                    opid += num;
                }
                EditOp::Ins => {
                    qid += num;
                    opid += num;
                }
            }
            if qid >= q_cut && sid >= s_cut {
                found = true;
            }
            if found {
                break;
            }
        }
        if found {
            break;
        }
        index += 1;
    }
    if !found {
        return;
    }
    if cut_begin {
        let mut new_ops = Vec::with_capacity(esp.ops.len());
        if opid < esp.ops[index].1 {
            new_ops.push((esp.ops[index].0, esp.ops[index].1 - opid));
        }
        new_ops.extend_from_slice(&esp.ops[index + 1..]);
        esp.ops = new_ops;
        hsp.query.offset += qid;
        hsp.subject.offset += sid;
    } else {
        if opid < esp.ops[index].1 {
            esp.ops[index].1 = opid;
        }
        esp.ops.truncate(index + 1);
        hsp.query.end = hsp.query.offset + qid;
        hsp.subject.end = hsp.subject.offset + sid;
    }
}

/// `s_Heapify` on `v[base]`, with `lim` and `last` as indices.
fn heapify<T>(v: &mut [T], mut base: usize, lim: usize, last: usize, cmp: fn(&T, &T) -> Ordering) {
    let mut left_son = 2 * base + 1;
    while base <= lim {
        let large_son = if left_son == last {
            left_son
        } else if cmp(&v[left_son], &v[left_son + 1]) != Ordering::Less {
            left_son
        } else {
            left_son + 1
        };
        if cmp(&v[base], &v[large_son]) == Ordering::Less {
            v.swap(base, large_son);
            base = large_son;
            left_son = 2 * base + 1;
        } else {
            break;
        }
    }
}

/// `s_CreateHeap`.
fn create_heap<T>(v: &mut [T], cmp: fn(&T, &T) -> Ordering) {
    let nel = v.len();
    if nel < 2 {
        return;
    }
    let lim = (nel - 2) / 2;
    let last = nel - 1;
    let mut i = nel / 2;
    while i > 0 {
        heapify(v, i - 1, lim, last, cmp);
        i -= 1;
    }
}

const OVERLAP_DIAG_CLOSE: i32 = 10;

/// `s_BlastMergeTwoHSPs`.
fn merge_two_hsps(hsp1: &mut Hsp, hsp2: &Hsp, allow_gap: bool) -> bool {
    if !allow_gap
        && (hsp1.subject.offset - hsp2.subject.offset - hsp1.query.offset + hsp2.query.offset) != 0
    {
        return false;
    }
    if hsp1.subject.frame != hsp2.subject.frame {
        return false;
    }
    let contained =
        |a: i32, b: i32, c: i32, d: i32, e: i32, f: i32| (a <= c && b >= c) && (d <= f && e >= f);
    if contained(
        hsp1.query.offset,
        hsp1.query.end,
        hsp2.query.offset,
        hsp1.subject.offset,
        hsp1.subject.end,
        hsp2.subject.offset,
    ) || contained(
        hsp1.query.offset,
        hsp1.query.end,
        hsp2.query.end,
        hsp1.subject.offset,
        hsp1.subject.end,
        hsp2.subject.end,
    ) {
        let score_density = (hsp1.score + hsp2.score) as f64 * 1.0
            / ((hsp1.query.end - hsp1.query.offset) + (hsp2.query.end - hsp2.query.offset)) as f64;
        hsp1.query.offset = hsp1.query.offset.min(hsp2.query.offset);
        hsp1.subject.offset = hsp1.subject.offset.min(hsp2.subject.offset);
        hsp1.query.end = hsp1.query.end.max(hsp2.query.end);
        hsp1.subject.end = hsp1.subject.end.max(hsp2.subject.end);
        if hsp2.score > hsp1.score {
            hsp1.query.gapped_start = hsp2.query.gapped_start;
            hsp1.subject.gapped_start = hsp2.subject.gapped_start;
            hsp1.score = hsp2.score;
        }
        hsp1.score =
            ((score_density * (hsp1.query.end - hsp1.query.offset) as f64) as i32).max(hsp1.score);
        return true;
    }
    false
}

/// `Blast_HSPListsMerge` for a split subject (`contexts_per_query < 0`): merges HSPs across
/// the chunk overlap and combines the lists.
pub fn merge_chunk_lists(
    hsp_list: HspList,
    combined: &mut Option<HspList>,
    split_offset: i32,
    chunk_overlap_size: i32,
    allow_gap: bool,
) {
    if hsp_list.hsps.is_empty() {
        return;
    }
    let Some(c) = combined else {
        *combined = Some(hsp_list);
        return;
    };
    let mut hsp_list = hsp_list;
    // Put HSPs that intersect the overlap region first, preserving the C's swap order.
    let mut hspcnt1 = 0usize;
    for index1 in 0..c.hsps.len() {
        if c.hsps[index1].subject.end > split_offset {
            c.hsps.swap(hspcnt1, index1);
            hspcnt1 += 1;
        }
    }
    let mut hspcnt2 = 0usize;
    for index2 in 0..hsp_list.hsps.len() {
        if hsp_list.hsps[index2].subject.offset < split_offset + chunk_overlap_size {
            hsp_list.hsps.swap(hspcnt2, index2);
            hspcnt2 += 1;
        }
    }
    if hspcnt1 > 0 && hspcnt2 > 0 {
        let mut dead = vec![false; hsp_list.hsps.len()];
        for index1 in 0..hspcnt1 {
            for index2 in 0..hspcnt2 {
                if dead[index2] || c.hsps[index1].context != hsp_list.hsps[index2].context {
                    continue;
                }
                let hsp1 = &c.hsps[index1];
                let hsp2 = &hsp_list.hsps[index2];
                let end_diag = hsp1.query.end - hsp1.subject.end;
                let start_diag = hsp2.query.offset - hsp2.subject.offset;
                if (end_diag - start_diag).abs() < OVERLAP_DIAG_CLOSE {
                    let h2 = hsp_list.hsps[index2].clone();
                    if merge_two_hsps(&mut c.hsps[index1], &h2, allow_gap) {
                        dead[index2] = true;
                    }
                }
            }
        }
        let mut k = 0;
        hsp_list.hsps.retain(|_| {
            let keep = !dead[k];
            k += 1;
            keep
        });
    }
    // s_BlastHSPListsCombineByScore with room for every HSP (blastn never sets
    // hsp_num_max, so the combined list is never capped)
    c.hsps.append(&mut hsp_list.hsps);
    c.sort_by_score();
}

/// `s_EvalueCompareHSPLists`.
pub fn evalue_compare_hsp_lists(a: &HspList, b: &HspList) -> Ordering {
    if a.hsps.is_empty() && b.hsps.is_empty() {
        return Ordering::Equal;
    }
    if a.hsps.is_empty() {
        return Ordering::Greater;
    }
    if b.hsps.is_empty() {
        return Ordering::Less;
    }
    let r = evalue_comp(a.best_evalue, b.best_evalue);
    if r != Ordering::Equal {
        return r;
    }
    if a.hsps[0].score > b.hsps[0].score {
        return Ordering::Less;
    }
    if a.hsps[0].score < b.hsps[0].score {
        return Ordering::Greater;
    }
    b.oid.cmp(&a.oid)
}

/// `BlastHitList`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct HitList {
    pub hsplist_max: i32,
    pub worst_evalue: f64,
    pub low_score: i32,
    pub heapified: bool,
    pub hsplists: Vec<HspList>,
}

impl HitList {
    /// `Blast_HitListNew`.
    pub fn new(hitlist_size: i32) -> HitList {
        HitList {
            hsplist_max: hitlist_size,
            worst_evalue: 0.0,
            low_score: i32::MAX,
            heapified: false,
            hsplists: Vec::new(),
        }
    }

    /// `Blast_HitListUpdate`.
    pub fn update(&mut self, mut hsp_list: HspList) {
        hsp_list.best_evalue = hsp_list.best_evalue();
        if (self.hsplists.len() as i32) < self.hsplist_max {
            self.worst_evalue = hsp_list.best_evalue.max(self.worst_evalue);
            self.low_score = hsp_list.hsps[0].score.min(self.low_score);
            self.hsplists.push(hsp_list);
        } else {
            if !self.heapified {
                for l in self.hsplists.iter_mut() {
                    l.sort_by_evalue();
                    l.best_evalue = l.best_evalue();
                }
                create_heap(&mut self.hsplists, evalue_compare_hsp_lists);
                self.heapified = true;
            }
            hsp_list.sort_by_evalue();
            hsp_list.best_evalue = hsp_list.best_evalue();
            if evalue_compare_hsp_lists(&self.hsplists[0], &hsp_list) == Ordering::Less {
                // the new list is worse than the worst kept: discard
            } else {
                // s_BlastHitListInsertHSPListInHeap
                self.hsplists[0] = hsp_list;
                let n = self.hsplists.len();
                if n >= 2 {
                    heapify(
                        &mut self.hsplists,
                        0,
                        n / 2 - 1,
                        n - 1,
                        evalue_compare_hsp_lists,
                    );
                }
                self.worst_evalue = self.hsplists[0].best_evalue;
                self.low_score = self.hsplists[0].hsps[0].score;
            }
        }
    }

    /// `Blast_HSPResultsSortByEvalue` for one hit list: sort, then `s_BlastHitListPurge`.
    pub fn sort_by_evalue_and_purge(&mut self) {
        if self.hsplists.len() > 1 {
            self.hsplists.sort_by(evalue_compare_hsp_lists);
            debug_assert!(
                self.hsplists
                    .windows(2)
                    .all(
                        |w| evalue_compare_hsp_lists(&w[0], &w[1]) != Ordering::Equal
                            || (w[0].hsps.is_empty() && w[1].hsps.is_empty())
                    ),
                "two distinct HSP lists compare equal"
            );
        }
        let n = self
            .hsplists
            .iter()
            .take_while(|l| !l.hsps.is_empty())
            .count();
        self.hsplists.truncate(n);
    }
}

/// `BlastHSPResults`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct HspResults {
    pub hitlists: Vec<Option<HitList>>,
}

impl HspResults {
    pub fn new(num_queries: usize) -> HspResults {
        HspResults {
            hitlists: vec![None; num_queries],
        }
    }

    /// `Blast_HSPResultsInsertHSPList`.
    pub fn insert_hsp_list(&mut self, hsp_list: HspList, hitlist_size: i32) {
        if hsp_list.hsps.is_empty() {
            return;
        }
        let qi = hsp_list.query_index as usize;
        let hl = self.hitlists[qi].get_or_insert_with(|| HitList::new(hitlist_size));
        hl.update(hsp_list);
    }

    /// `Blast_HSPResultsSortByEvalue`.
    pub fn sort_by_evalue(&mut self) {
        for hl in self.hitlists.iter_mut().flatten() {
            hl.sort_by_evalue_and_purge();
        }
    }
}

/// `GetPrelimHitlistSize` without composition-based statistics.
pub fn prelim_hitlist_size(hitlist_size: i32, gapped: bool) -> i32 {
    if gapped {
        (2 * hitlist_size).max(10).min(hitlist_size + 50)
    } else {
        hitlist_size
    }
}

/// `Blast_HSPCalcLengthAndGaps`: (length, gaps, gap_opens).
pub fn calc_length_and_gaps(h: &Hsp) -> (i32, i32, i32) {
    let mut length = h.query.end - h.query.offset;
    let s_length = h.subject.end - h.subject.offset;
    let mut gap_opens = 0;
    let mut gaps = 0;
    if let Some(esp) = &h.gap_info {
        for &(op, num) in &esp.ops {
            match op {
                EditOp::Del => {
                    length += num;
                    gaps += num;
                    gap_opens += 1;
                }
                EditOp::Ins => {
                    gap_opens += 1;
                    gaps += num;
                }
                EditOp::Sub => {}
            }
        }
    } else if s_length > length {
        length = s_length;
    }
    (length, gaps, gap_opens)
}

/// `s_Blast_HSPGetNumIdentitiesAndPositives` for nucleotides: (num_ident, align_length).
/// `query`/`subject` are the context's letters and the subject's blastna.
pub fn num_identities(query: &[u8], subject: &[u8], h: &Hsp) -> (i32, i32) {
    let mut q = h.query.offset as usize;
    let mut s = h.subject.offset as usize;
    let mut num_ident = 0;
    let mut align_length = 0;
    match &h.gap_info {
        None => {
            let ql = h.query.end - h.query.offset;
            align_length = ql;
            for _ in 0..ql {
                if query[q] == subject[s] {
                    num_ident += 1;
                }
                q += 1;
                s += 1;
            }
        }
        Some(esp) => {
            for &(op, num) in &esp.ops {
                align_length += num;
                match op {
                    EditOp::Sub => {
                        for _ in 0..num {
                            if query[q] == subject[s] {
                                num_ident += 1;
                            }
                            q += 1;
                            s += 1;
                        }
                    }
                    EditOp::Del => s += num as usize,
                    EditOp::Ins => q += num as usize,
                }
            }
        }
    }
    (num_ident, align_length)
}

/// `Blast_HSPTestIdentityAndLength`: sets `num_ident`, returns true when the HSP fails the
/// `-perc_identity` test (`min_hit_length` is 0 for blastn).
pub fn test_identity_and_length(
    h: &mut Hsp,
    query: &[u8],
    subject: &[u8],
    percent_identity: f64,
) -> bool {
    let (ni, al) = num_identities(query, subject, h);
    h.num_ident = ni;
    (ni as f64 * 100.0) < al as f64 * percent_identity
}

/// `Blast_HSPReevaluateWithAmbiguitiesGapped`: rescore along the edit script with the real
/// subject letters; returns true when the HSP should be deleted. `q`/`s` are the context's
/// query letters and the subject's blastna (full length `qlen`/`slen`).
pub fn reevaluate_with_ambiguities_gapped(
    h: &mut Hsp,
    q: &[u8],
    qlen: i32,
    s: &[u8],
    slen: i32,
    params: &Parameters,
    sbp: &ScoreBlk,
) -> bool {
    let cutoff_score = params.hit_cutoffs[h.context as usize].cutoff_score;
    let matrix = &sbp.matrix;
    let (factor, gap_open, gap_extend) = if params.gap_open == 0 && params.gap_extend == 0 {
        let factor = if params.reward % 2 == 1 { 2 } else { 1 };
        (factor, 0, (params.reward - 2 * params.penalty) * factor / 2)
    } else {
        (1, params.gap_open, params.gap_extend)
    };
    let Some(esp) = h.gap_info.as_mut() else {
        return true;
    };
    let mut query = h.query.offset as isize;
    let mut subject = h.subject.offset as isize;
    let mut score = 0i32;
    let mut sum = 0i32;
    let (mut best_q_start, mut best_q_end, mut current_q_start) = (query, query, query);
    let (mut best_s_start, mut best_s_end, mut current_s_start) = (subject, subject, subject);
    let mut best_start_esp_index = 0usize;
    let mut best_end_esp_index = 0usize;
    let mut current_start_esp_index = 0usize;
    let mut best_end_esp_num: i32 = -1;
    let size = esp.ops.len();
    for index in 0..size {
        let mut op_index = 0i32;
        while op_index < esp.ops[index].1 {
            let (op, num) = esp.ops[index];
            match op {
                EditOp::Sub => {
                    sum += factor
                        * matrix[(q[query as usize] & 0x0f) as usize][s[subject as usize] as usize];
                    query += 1;
                    subject += 1;
                    op_index += 1;
                }
                EditOp::Del => {
                    sum -= gap_open + gap_extend * num;
                    subject += num as isize;
                    op_index += num;
                }
                EditOp::Ins => {
                    sum -= gap_open + gap_extend * num;
                    query += num as isize;
                    op_index += num;
                }
            }
            if sum < 0 {
                if op_index < esp.ops[index].1 {
                    esp.ops[index].1 -= op_index;
                    current_start_esp_index = index;
                    op_index = 0;
                } else {
                    current_start_esp_index = index + 1;
                }
                sum = 0;
                current_q_start = query;
                current_s_start = subject;
                if score < cutoff_score {
                    best_q_start = query;
                    best_s_start = subject;
                    score = 0;
                    best_start_esp_index = current_start_esp_index;
                    best_end_esp_index = current_start_esp_index;
                }
            } else if sum > score {
                score = sum;
                best_q_start = current_q_start;
                best_s_start = current_s_start;
                best_q_end = query;
                best_s_end = subject;
                best_start_esp_index = current_start_esp_index;
                best_end_esp_index = index;
                best_end_esp_num = op_index;
            }
        }
    }
    score /= factor;
    if best_start_esp_index < size && best_end_esp_index < size {
        let mut qp = best_q_start;
        let mut sp = best_s_start;
        let mut ext = 0;
        while qp > 0 && sp > 0 {
            qp -= 1;
            sp -= 1;
            if q[qp as usize] == s[sp as usize] && q[qp as usize] < 4 {
                ext += 1;
            } else {
                break;
            }
        }
        best_q_start -= ext as isize;
        best_s_start -= ext as isize;
        esp.ops[best_start_esp_index].1 += ext;
        if best_end_esp_index == best_start_esp_index {
            best_end_esp_num += ext;
        }
        score += ext * params.reward;
        let mut qp = best_q_end;
        let mut sp = best_s_end;
        let mut ext = 0;
        while (qp as i32) < qlen && (sp as i32) < slen && q[qp as usize] < 4 {
            let eq = q[qp as usize] == s[sp as usize];
            qp += 1;
            sp += 1;
            if eq {
                ext += 1;
            } else {
                break;
            }
        }
        best_q_end += ext as isize;
        best_s_end += ext as isize;
        esp.ops[best_end_esp_index].1 += ext;
        best_end_esp_num += ext;
        score += ext * params.reward;
    }
    // s_UpdateReevaluatedHSP
    h.score = score;
    if h.score >= cutoff_score {
        h.query.offset = best_q_start as i32;
        h.query.end = h.query.offset + (best_q_end - best_q_start) as i32;
        h.subject.offset = best_s_start as i32;
        h.subject.end = h.subject.offset + (best_s_end - best_s_start) as i32;
        let last = esp.ops.len() - 1;
        if best_end_esp_index != last || best_start_esp_index > 0 {
            esp.ops = esp.ops[best_start_esp_index..=best_end_esp_index].to_vec();
        }
        let last = esp.ops.len() - 1;
        esp.ops[last].1 = best_end_esp_num;
        false
    } else {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hsp(score: i32, qo: i32, qe: i32, so: i32, se: i32) -> Hsp {
        Hsp::new(qo, qe, so, se, 0, 0, 0, 1, 1, score, None)
    }

    #[test]
    fn score_order_is_total_on_distinct_hsps() {
        let mut v = [
            hsp(10, 0, 10, 0, 10),
            hsp(12, 5, 17, 5, 17),
            hsp(10, 0, 10, 2, 12),
            hsp(10, 0, 11, 0, 10),
        ];
        v.sort_by(score_compare_hsps);
        assert_eq!(v[0].score, 12);
        assert_eq!((v[1].subject.offset, v[1].query.end), (0, 11)); // longer query end first
        assert_eq!(v[3].subject.offset, 2);
    }

    #[test]
    fn evalue_epsilon_ties_fall_through_to_score() {
        let mut a = hsp(50, 0, 50, 0, 50);
        let mut b = hsp(60, 0, 60, 0, 60);
        a.evalue = 1e-190;
        b.evalue = 1e-200;
        assert_eq!(evalue_compare_hsps(&a, &b), Ordering::Greater);
        let mut la = HspList::default();
        la.hsps.push(a);
        la.oid = 3;
        la.best_evalue = la.best_evalue();
        let mut lb = HspList::default();
        lb.hsps.push(b);
        lb.oid = 1;
        lb.best_evalue = lb.best_evalue();
        assert_eq!(evalue_compare_hsp_lists(&la, &lb), Ordering::Greater);
        // same evalue and score: the higher oid sorts first
        lb.hsps[0].score = 50;
        assert_eq!(evalue_compare_hsp_lists(&la, &lb), Ordering::Less);
    }

    #[test]
    fn hit_list_heap_keeps_best_and_higher_oid_on_ties() {
        let mut hl = HitList::new(2);
        for (oid, e) in [(0, 1e-5), (1, 1e-3)] {
            let mut l = HspList::default();
            let mut h = hsp(30, 0, 30, 0, 30);
            h.evalue = e;
            l.hsps.push(h);
            l.oid = oid;
            hl.update(l);
        }
        // a third list ties with the worst (oid 1): higher oid replaces it
        let mut l = HspList::default();
        let mut h = hsp(30, 0, 30, 0, 30);
        h.evalue = 1e-3;
        l.hsps.push(h);
        l.oid = 7;
        hl.update(l);
        assert!(hl.heapified);
        let mut oids: Vec<i32> = hl.hsplists.iter().map(|l| l.oid).collect();
        oids.sort();
        assert_eq!(oids, vec![0, 7]);
        // a worse one is discarded
        let mut l = HspList::default();
        let mut h = hsp(20, 0, 20, 0, 20);
        h.evalue = 1e-2;
        l.hsps.push(h);
        l.oid = 9;
        hl.update(l);
        let mut oids: Vec<i32> = hl.hsplists.iter().map(|l| l.oid).collect();
        oids.sort();
        assert_eq!(oids, vec![0, 7]);
    }

    #[test]
    fn purge_common_endpoints_removes_duplicates() {
        let mut l = HspList::default();
        l.hsps.push(hsp(20, 0, 20, 0, 20));
        l.hsps.push(hsp(15, 0, 15, 0, 15)); // same start
        l.hsps.push(hsp(10, 5, 20, 5, 20)); // same end
        let n = l.purge_common_endpoints(true);
        assert_eq!(n, 1);
        assert_eq!(l.hsps.len(), 1);
        assert_eq!(l.hsps[0].score, 20);
    }

    /// blastn's cut (not purged) HSPs sit after the live range in reverse order of moving:
    /// the C writes each to `hsp_array[hsp_count]` after decrementing the count.
    #[test]
    fn purge_common_endpoints_cut_order() {
        let mut l = HspList::default();
        let long = |score: i32, qe: i32| {
            let mut h = hsp(score, 0, qe, 0, qe);
            h.gap_info = Some(EditScript {
                ops: vec![(EditOp::Sub, qe)],
            });
            h
        };
        l.hsps.push(long(20, 20)); // kept
        l.hsps.push(long(15, 40)); // moved first (higher score sorts first)
        l.hsps.push(long(10, 60)); // moved second
        let n = l.purge_common_endpoints(false);
        assert_eq!(n, 1);
        assert_eq!(l.hsps.len(), 3);
        assert_eq!(l.hsps[0].score, 20);
        assert_eq!(l.hsps[1].score, 10, "the last moved comes first");
        assert_eq!(l.hsps[2].score, 15);
        assert_eq!((l.hsps[1].query.offset, l.hsps[1].query.end), (20, 60));
    }

    #[test]
    fn length_and_gaps() {
        let mut h = hsp(10, 0, 10, 0, 12);
        h.gap_info = Some(EditScript {
            ops: vec![(EditOp::Sub, 4), (EditOp::Del, 2), (EditOp::Sub, 6)],
        });
        assert_eq!(calc_length_and_gaps(&h), (12, 2, 1));
        assert_eq!(prelim_hitlist_size(500, true), 550);
        assert_eq!(prelim_hitlist_size(1, true), 10);
    }
}
