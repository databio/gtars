//! Traceback stage: final alignments, re-evaluation, e-values and the hit-list pruning.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/core/blast_traceback.c` (`Blast_TracebackFromHSPList`
//! for the greedy and the dynamic-programming paths, `Blast_HSPUpdateWithTraceback`, `s_HSPListPostTracebackUpdate`,
//! `s_FilterBlastResults`, `s_BlastPruneExtraHits`, `BLAST_ComputeTraceback_MT` in one
//! thread, `s_SThreadLocalDataArraySetGapXDropoffFinal`), `blast_hspstream.c::
//! BlastHSPStreamClose` (HSP lists grouped by OID) and `blast_gapalign.c`
//! (`BlastGetOffsetsForGappedAlignment`, `BlastGetStartForGappedAlignmentNucl`,
//! `AdjustSubjectRange`). blastn uses the `Nucl` start finder for both tasks; the generic
//! `BlastGetStartForGappedAlignment` is reached only by protein searches.

use super::gapalign::{MAX_SUBJECT_OFFSET, MAX_TOTAL_GAPS, adjust_subject_range};
use super::hsp::{
    Hsp, HspList, HspResults, reevaluate_with_ambiguities_gapped, test_identity_and_length,
};
use super::itree::IntervalTree;
use super::prelim::{GapAligner, SearchContext};
use super::subject::StoreSubjects;
use crate::error::Error;
use std::collections::HashMap;
use std::ops::Range;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, mpsc};

const HSP_MAX_WINDOW: i32 = 11;

/// `BlastGetOffsetsForGappedAlignment` (nucleotide matrix): `Some((q, s))` or `None`.
fn offsets_for_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    matrix: &[[i32; 16]; 16],
    hsp: &Hsp,
) -> Option<(i32, i32)> {
    let q_length = hsp.query.end - hsp.query.offset;
    let s_length = hsp.subject.end - hsp.subject.offset;
    let q_start = hsp.query.offset;
    let s_start = hsp.subject.offset;
    if q_length <= HSP_MAX_WINDOW {
        return Some((q_start + q_length / 2, s_start + q_length / 2));
    }
    let m = |qi: i32, si: i32| matrix[query[qi as usize] as usize][subject[si as usize] as usize];
    let mut hsp_end = q_start + HSP_MAX_WINDOW;
    let mut score = 0;
    let mut qi = q_start;
    let mut si = s_start;
    for _ in q_start..hsp_end {
        score += m(qi, si);
        qi += 1;
        si += 1;
    }
    let mut max_score = score;
    let mut max_offset = hsp_end - 1;
    hsp_end = q_start + q_length.min(s_length);
    let mut index1 = q_start + HSP_MAX_WINDOW;
    while index1 < hsp_end {
        score -= m(qi - HSP_MAX_WINDOW, si - HSP_MAX_WINDOW);
        score += m(qi, si);
        if score > max_score {
            max_score = score;
            max_offset = index1;
        }
        qi += 1;
        si += 1;
        index1 += 1;
    }
    if max_score > 0 {
        return Some((max_offset, (max_offset - q_start) + s_start));
    }
    let mut score = 0;
    let mut qi = q_start + q_length - HSP_MAX_WINDOW;
    let mut si = s_start + s_length - HSP_MAX_WINDOW;
    for _ in (hsp.query.end - HSP_MAX_WINDOW)..hsp.query.end {
        score += m(qi, si);
        qi += 1;
        si += 1;
    }
    if score > 0 {
        return Some((
            hsp.query.end - HSP_MAX_WINDOW / 2,
            hsp.subject.end - HSP_MAX_WINDOW / 2,
        ));
    }
    None
}

/// `BlastGetStartForGappedAlignmentNucl`: moves `gapped_start` into a run of identities.
fn start_for_gapped_alignment_nucl(query: &[u8], subject: &[u8], hsp: &mut Hsp) {
    let mut hsp_max_ident_run = 10;
    let offset = (hsp.subject.gapped_start - hsp.subject.offset)
        .min(hsp.query.gapped_start - hsp.query.offset);
    let mut q_start = hsp.query.gapped_start;
    let mut s_start = hsp.subject.gapped_start;
    let mut score: i32 = -1;
    let mut q = q_start;
    let mut s = s_start;
    let q_len = hsp.query.end;
    while q < q_len && query[q as usize] == subject[s as usize] {
        q += 1;
        s += 1;
        score += 1;
        if score > hsp_max_ident_run {
            return;
        }
    }
    let mut q = q_start;
    let mut s = s_start;
    // the C reads subject[-1], the sentinel before the sequence, which never matches
    while q >= 0 && s >= 0 && query[q as usize] == subject[s as usize] {
        q -= 1;
        s -= 1;
        score += 1;
        if score > hsp_max_ident_run {
            return;
        }
    }
    hsp_max_ident_run = (hsp_max_ident_run as f64 * 1.5) as i32;
    q_start = hsp.query.gapped_start - offset;
    s_start = hsp.subject.gapped_start - offset;
    let q_len = (hsp.subject.end - s_start).min(hsp.query.end - q_start);
    let mut q = q_start;
    let mut s = s_start;
    let mut max_score = 0;
    let mut max_offset = q_start;
    let mut score = 0;
    let mut m = false;
    let mut prev_match = false;
    let mut index = q_start;
    while index < q_start + q_len {
        m = query[q as usize] == subject[s as usize];
        q += 1;
        s += 1;
        if m != prev_match {
            prev_match = m;
            if m {
                score = 1;
            } else if score > max_score {
                max_score = score;
                max_offset = index - score / 2;
            }
        } else if m {
            score += 1;
            if score > hsp_max_ident_run {
                max_offset = index - hsp_max_ident_run / 2;
                hsp.query.gapped_start = max_offset;
                hsp.subject.gapped_start = max_offset + s_start - q_start;
                return;
            }
        }
        index += 1;
    }
    if m && score > max_score {
        max_score = score;
        max_offset = index - score / 2;
    }
    if max_score > 0 {
        hsp.query.gapped_start = max_offset;
        hsp.subject.gapped_start = max_offset + s_start - q_start;
    }
}

/// The `sequence_nomask` letters of context `c`.
fn query_nomask(qb: &super::query::QueryBlock, c: usize) -> &[u8] {
    match &qb.nomask {
        Some(nm) => {
            let o = qb.contexts[c].query_offset as usize + 1;
            &nm[o..o + qb.contexts[c].query_length as usize]
        }
        None => qb.context_seq(c),
    }
}

/// The letters of a long subject that one HSP list's traceback reads: around each HSP, its
/// subject span widened by the query length plus `MAX_TOTAL_GAPS` (what `AdjustSubjectRange`
/// lets the final alignment reach on either side) and a little slack, clipped to the subject;
/// sorted, overlaps merged. The start-finding scans stay inside an HSP's span, the alignment
/// inside its adjusted range and the identity counts inside the aligned span, so the
/// traceback reads no other letter (blastn's SeqDB likewise fetches only HSP ranges of a long
/// subject).
pub fn hsp_windows(ctx: &SearchContext, hsp_list: &HspList, slen: i32) -> Vec<Range<i32>> {
    let mut w: Vec<Range<i32>> = hsp_list
        .hsps
        .iter()
        .map(|h| {
            let qlen = ctx.qb.contexts[h.context as usize].query_length;
            let margin = qlen.saturating_add(MAX_TOTAL_GAPS).saturating_add(16);
            let lo = h.subject.offset.min(h.subject.gapped_start);
            let hi = h.subject.end.max(h.subject.gapped_start);
            lo.saturating_sub(margin).max(0)..hi.saturating_add(margin).min(slen)
        })
        .collect();
    w.sort_by_key(|r| r.start);
    let mut merged: Vec<Range<i32>> = Vec::with_capacity(w.len());
    for r in w {
        match merged.last_mut() {
            Some(m) if r.start <= m.end => m.end = m.end.max(r.end),
            _ => merged.push(r),
        }
    }
    merged
}

/// `Blast_TracebackFromHSPList`. `subject` is the whole subject in blastna (`slen` letters),
/// followed by a sentinel.
pub fn traceback_from_hsp_list(
    ctx: &SearchContext,
    hsp_list: &mut HspList,
    subject: &[u8],
    slen: i32,
    gap: &mut GapAligner,
) {
    if hsp_list.hsps.is_empty() {
        return;
    }
    let qb = ctx.qb;
    let min_diag_sep = ctx.options.min_diag_separation;
    let matrix = &ctx.sbp.matrix;
    let percent_identity = ctx.options.percent_identity;
    let greedy_traceback = matches!(gap, GapAligner::Greedy(_));
    let mut tree = IntervalTree::new(0, qb.length + 1, 0, slen + 1);
    let mut kept: Vec<Hsp> = Vec::with_capacity(hsp_list.hsps.len());
    let initial = std::mem::take(&mut hsp_list.hsps);
    for mut hsp in initial {
        let ctx_i = hsp.context as usize;
        let query = qb.context_seq(ctx_i);
        let query_length = qb.contexts[ctx_i].query_length;
        if tree.contains_hsp(&hsp, &kept, qb, min_diag_sep) {
            continue;
        }
        let (q_start, mut s_start);
        if hsp.query.gapped_start == 0 && hsp.subject.gapped_start == 0 {
            match offsets_for_gapped_alignment(query, subject, matrix, &hsp) {
                Some((q, s)) => {
                    q_start = q;
                    s_start = s;
                    hsp.query.gapped_start = q;
                    hsp.subject.gapped_start = s;
                }
                None => continue,
            }
        } else {
            start_for_gapped_alignment_nucl(query, subject, &mut hsp);
            q_start = hsp.query.gapped_start;
            s_start = hsp.subject.gapped_start;
        }
        // AdjustSubjectRange: long subjects are aligned inside a window around the start
        let (adj_s_start, adjusted_s_length, start_shift) =
            adjust_subject_range(s_start, slen, q_start, query_length);
        s_start = adj_s_start;
        hsp.subject.gapped_start = s_start;
        let adjusted_subject = &subject[start_shift as usize..];
        let (score, qs, qe, ss, se, esp) = match gap {
            GapAligner::Greedy(greedy) => {
                let r = greedy.align(
                    query,
                    adjusted_subject,
                    query_length,
                    adjusted_s_length,
                    q_start,
                    s_start,
                    false,
                    true,
                );
                (
                    r.score,
                    r.query_start,
                    r.query_stop,
                    r.subject_start,
                    r.subject_stop,
                    r.edit_script,
                )
            }
            GapAligner::Dp(dp) => {
                let r = dp.gapped_alignment_with_traceback(
                    qb.context_seq_with_sentinel(ctx_i),
                    adjusted_subject,
                    q_start,
                    s_start,
                    query_length,
                    adjusted_s_length,
                    matrix,
                );
                (
                    r.score,
                    r.query_start,
                    r.query_stop,
                    r.subject_start,
                    r.subject_stop,
                    r.edit_script,
                )
            }
        };
        // Blast_HSPUpdateWithTraceback
        hsp.score = score;
        hsp.query.offset = qs;
        hsp.subject.offset = ss;
        hsp.query.end = qe;
        hsp.subject.end = se;
        if esp.is_some() {
            hsp.gap_info = esp;
        }
        // Blast_HSPAdjustSubjectOffset
        hsp.subject.offset += start_shift;
        hsp.subject.end += start_shift;
        hsp.subject.gapped_start += start_shift;
        if !greedy_traceback {
            // the DP traceback is exact: identities and the -perc_identity test here
            let delete = test_identity_and_length(
                &mut hsp,
                query_nomask(qb, ctx_i),
                subject,
                percent_identity,
            );
            if delete {
                continue;
            }
        }
        kept.push(hsp);
        let idx = kept.len() - 1;
        tree.add_hsp(idx, &kept, qb);
    }
    hsp_list.hsps = kept;
    // common endpoints: cut, not purged (blastn). The cut HSPs sit after the live range;
    // only they are re-evaluated by the DP path, while the greedy path (which ignored
    // ambiguities) re-evaluates everything.
    let mut extra_start = hsp_list.purge_common_endpoints(false);
    if greedy_traceback {
        extra_start = 0;
    }
    let mut survivors: Vec<Hsp> = Vec::with_capacity(hsp_list.hsps.len());
    for (i, mut hsp) in std::mem::take(&mut hsp_list.hsps).into_iter().enumerate() {
        if i < extra_start {
            survivors.push(hsp);
            continue;
        }
        let ctx_i = hsp.context as usize;
        let query = qb.context_seq(ctx_i);
        let qlen = qb.contexts[ctx_i].query_length;
        let mut delete = reevaluate_with_ambiguities_gapped(
            &mut hsp, query, qlen, subject, slen, ctx.params, ctx.sbp,
        );
        if !delete {
            delete = test_identity_and_length(
                &mut hsp,
                query_nomask(qb, ctx_i),
                subject,
                percent_identity,
            );
        }
        if !delete {
            survivors.push(hsp);
        }
    }
    hsp_list.hsps = survivors;
    hsp_list.purge_common_endpoints(true);
    hsp_list.sort_by_score();
    tree.reset();
    let mut kept: Vec<Hsp> = Vec::with_capacity(hsp_list.hsps.len());
    for hsp in std::mem::take(&mut hsp_list.hsps) {
        if tree.contains_hsp(&hsp, &kept, qb, min_diag_sep) {
            continue;
        }
        kept.push(hsp);
        let idx = kept.len() - 1;
        tree.add_hsp(idx, &kept, qb);
    }
    hsp_list.hsps = kept;
    // s_HSPListPostTracebackUpdate
    hsp_list.get_evalues(&ctx.sbp.kbp_gap_std, &ctx.eff_searchsp, ctx.sbp.round_down);
    hsp_list.reap_by_evalue(ctx.options.expect_value);
    hsp_list.get_bit_scores(&ctx.sbp.kbp_gap_std);
}

/// Subjects at least this long are decoded per HSP list, only around its HSPs
/// (`hsp_windows`); shorter ones are decoded whole, once per subject, and shared by their
/// lists (`MAX_SUBJECT_OFFSET`, below which `AdjustSubjectRange` does nothing either).
pub const WINDOW_MIN_SUBJECT: i32 = MAX_SUBJECT_OFFSET;

/// One subject's blastna letters, decoded once by the first thread that needs them and shared
/// by every list of that subject; dropped when its last list is done.
struct SubjectSlot {
    letters: Mutex<Option<Arc<Vec<u8>>>>,
    remaining: AtomicUsize,
}

/// What the traceback threads share.
struct TracebackShared<'a> {
    ctx: &'a SearchContext<'a>,
    subjects: &'a StoreSubjects,
    /// The HSP lists in OID order (within an OID, query order); a worker takes a list out,
    /// runs the traceback and puts it back in the same slot, so results are placed by index.
    items: Vec<Mutex<Option<HspList>>>,
    slots: HashMap<i32, SubjectSlot>,
    next: AtomicUsize,
    cancel: AtomicBool,
    max_subject_length: i32,
}

/// Messages from a traceback worker to the main thread.
enum TbMsg {
    /// A subject's last list is done.
    SubjectDone,
    Error(Error),
}

fn traceback_worker(sh: &TracebackShared, tx: &mpsc::Sender<TbMsg>) {
    let mut gap = GapAligner::new(sh.ctx, sh.max_subject_length);
    // s_SThreadLocalDataArraySetGapXDropoffFinal
    gap.set_gap_x_dropoff(sh.ctx.params.gap_x_dropoff_final);
    let mut letters: Vec<u8> = Vec::new();
    loop {
        if sh.cancel.load(Ordering::Relaxed) {
            return;
        }
        let i = sh.next.fetch_add(1, Ordering::Relaxed);
        if i >= sh.items.len() {
            return;
        }
        let mut list = sh.items[i]
            .lock()
            .unwrap()
            .take()
            .expect("each traceback item is taken once");
        let oid = list.oid;
        let slot = &sh.slots[&oid];
        let slen = sh.subjects.seq_len(oid);
        if slen >= WINDOW_MIN_SUBJECT {
            let windows = hsp_windows(sh.ctx, &list, slen);
            match sh.subjects.blastna_sparse(oid, &windows, &mut letters) {
                Ok(subject) => traceback_from_hsp_list(sh.ctx, &mut list, &subject, slen, &mut gap),
                Err(e) => {
                    let _ = tx.send(TbMsg::Error(e));
                    return;
                }
            }
        } else {
            let subject = {
                let mut g = slot.letters.lock().unwrap();
                match &*g {
                    Some(s) => s.clone(),
                    None => {
                        let mut v = Vec::new();
                        if let Err(e) = sh.subjects.blastna_window(oid, 0..slen, &mut v) {
                            let _ = tx.send(TbMsg::Error(e));
                            return;
                        }
                        let s = Arc::new(v);
                        *g = Some(s.clone());
                        s
                    }
                }
            };
            traceback_from_hsp_list(sh.ctx, &mut list, &subject, slen, &mut gap);
        }
        *sh.items[i].lock().unwrap() = Some(list);
        if slot.remaining.fetch_sub(1, Ordering::AcqRel) == 1 {
            *slot.letters.lock().unwrap() = None;
            if tx.send(TbMsg::SubjectDone).is_err() {
                return;
            }
        }
    }
}

/// `BlastHSPStreamClose` + `BLAST_ComputeTraceback_MT`: consume the preliminary results, run
/// the traceback of every (query, hit list entry) pair on `threads` threads, build the final
/// results in OID order (then query order), filter and prune to `hitlist_size`. The
/// per-subject work is independent, and `Blast_HSPResultsInsertHSPList` is order-independent,
/// so thread timing cannot change the output. `progress` is called once per subject; `false`
/// cancels.
pub fn compute_traceback(
    ctx: &SearchContext,
    prelim: HspResults,
    subjects: &StoreSubjects,
    threads: usize,
    max_subject_length: i32,
    progress: &mut dyn FnMut(u32, u32) -> bool,
) -> Result<HspResults, Error> {
    // gather every HSP list with its query index, grouped by OID ascending
    let mut lists: Vec<HspList> = Vec::new();
    for (qi, hl) in prelim.hitlists.into_iter().enumerate() {
        if let Some(hl) = hl {
            for mut l in hl.hsplists {
                l.query_index = qi as i32;
                lists.push(l);
            }
        }
    }
    lists.sort_by_key(|l| l.oid);
    super::search::debug_shuffle(&mut lists);
    let mut slots: HashMap<i32, SubjectSlot> = HashMap::new();
    for l in &lists {
        slots
            .entry(l.oid)
            .or_insert_with(|| SubjectSlot {
                letters: Mutex::new(None),
                remaining: AtomicUsize::new(0),
            })
            .remaining
            .fetch_add(1, Ordering::Relaxed);
    }
    let total = slots.len() as u32;
    let hitlist_size = ctx.options.hitlist_size;
    let n_items = lists.len();
    let sh = TracebackShared {
        ctx,
        subjects,
        items: lists.into_iter().map(|l| Mutex::new(Some(l))).collect(),
        slots,
        next: AtomicUsize::new(0),
        cancel: AtomicBool::new(false),
        max_subject_length,
    };
    let workers = threads.clamp(1, n_items.max(1));
    let mut first_error: Option<Error> = None;
    let mut cancelled = false;
    if n_items > 0 {
        let (tx, rx) = mpsc::channel::<TbMsg>();
        std::thread::scope(|s| {
            for t in 0..workers {
                let tx = tx.clone();
                let sh = &sh;
                std::thread::Builder::new()
                    .name(format!("rgblastn-tb-{t}"))
                    .stack_size(super::search::WORKER_STACK)
                    .spawn_scoped(s, move || traceback_worker(sh, &tx))
                    .expect("spawn traceback worker");
            }
            drop(tx);
            let mut done = 0u32;
            for msg in rx {
                match msg {
                    TbMsg::SubjectDone => {
                        done += 1;
                        if !cancelled && !first_error.is_some() && !progress(done, total) {
                            cancelled = true;
                            sh.cancel.store(true, Ordering::Relaxed);
                        }
                    }
                    TbMsg::Error(e) => {
                        if first_error.is_none() {
                            first_error = Some(e);
                        }
                        sh.cancel.store(true, Ordering::Relaxed);
                    }
                }
            }
        });
    }
    if let Some(e) = first_error {
        return Err(e);
    }
    if cancelled {
        return Err(Error::Cancelled);
    }
    let mut results = HspResults::new(ctx.qb.num_queries());
    for slot in sh.items {
        let l = slot
            .into_inner()
            .unwrap()
            .expect("every traceback item was put back");
        if !l.hsps.is_empty() {
            results.insert_hsp_list(l, hitlist_size);
        }
    }
    // s_FilterBlastResults: -max_hsps
    let max_hsps = ctx.options.max_hsps_per_subject;
    if max_hsps > 0 {
        for hl in results.hitlists.iter_mut().flatten() {
            for l in hl.hsplists.iter_mut() {
                l.trim_by_max_hsps(max_hsps);
            }
        }
    }
    results.sort_by_evalue();
    // s_BlastPruneExtraHits
    for hl in results.hitlists.iter_mut().flatten() {
        hl.hsplists.truncate(hitlist_size as usize);
    }
    Ok(results)
}
