//! Preliminary stage: one subject in, its HSP list out.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/core/blast_gapalign.c::BLAST_GetGappedScore` (the greedy
//! and the dynamic-programming branches), `blast_engine.c` (`s_BackupSubject`, `s_GetNextSubjectChunk`,
//! `s_BlastSearchEngineOneContext`, `s_BlastSearchEngineCore`, the subject loop of
//! `BLAST_PreliminarySearchEngine`) and `hspfilter_collector.c::s_BlastHSPCollectorRun`.

use super::gapalign::DpAligner;
use super::greedy::GreedyAligner;
use super::hsp::{Hsp, HspList, HspResults, merge_chunk_lists};
use super::itree::IntervalTree;
use super::lookup::Lookup;
use super::options::{BlastOptions, DBSEQ_CHUNK_OVERLAP, MAX_DBSEQ_LEN, Parameters, PrelimGapExt};
use super::query::QueryBlock;
use super::stat::ScoreBlk;
use super::ungapped::{
    SubjectBlk, WordFinderContext, WordFinderScratch, bsearch_context, na_word_finder,
};

/// Everything the preliminary search reads and never writes.
pub struct SearchContext<'a> {
    pub qb: &'a QueryBlock,
    pub lookup: &'a Lookup,
    pub sbp: &'a ScoreBlk,
    pub params: &'a Parameters,
    pub options: &'a BlastOptions,
    /// `eff_searchsp` per context.
    pub eff_searchsp: Vec<i64>,
}

impl SearchContext<'_> {
    pub fn word_finder_ctx(&self) -> WordFinderContext<'_> {
        WordFinderContext {
            qb: self.qb,
            lookup: self.lookup,
            matrix: &self.sbp.matrix,
            params: self.params,
        }
    }
}

/// `BlastGapAlignStruct`: the greedy aligner (megablast) or the dynamic-programming one
/// (blastn), chosen by `ePrelimGapExt` as `BLAST_GapAlignStructNew` does.
pub enum GapAligner {
    Greedy(GreedyAligner),
    Dp(DpAligner),
}

impl GapAligner {
    /// `BLAST_GapAlignStructNew`: the aligner `ePrelimGapExt` names, sized for the longest
    /// subject (one per thread; the preliminary stage and the traceback stage each build
    /// their own).
    pub fn new(ctx: &SearchContext, max_subject_length: i32) -> GapAligner {
        let p = ctx.params;
        match ctx.options.prelim_gap_ext {
            PrelimGapExt::GreedyScoreOnly => GapAligner::Greedy(GreedyAligner::new(
                p.reward,
                p.penalty,
                p.gap_x_dropoff,
                p.gap_x_dropoff_final,
                max_subject_length,
            )),
            PrelimGapExt::DynProgScoreOnly => {
                GapAligner::Dp(DpAligner::new(p.gap_open, p.gap_extend, p.gap_x_dropoff))
            }
        }
    }

    /// `gap_align->gap_x_dropoff = ...` (the traceback stage sets the final x-drop).
    pub fn set_gap_x_dropoff(&mut self, x: i32) {
        match self {
            GapAligner::Greedy(g) => g.gap_x_dropoff = x,
            GapAligner::Dp(d) => d.gap_x_dropoff = x,
        }
    }
}

/// Per-thread scratch for the preliminary stage.
pub struct Scratch {
    pub wf: WordFinderScratch,
    pub gap: GapAligner,
    tree_hsps: Vec<Hsp>,
    /// `BlastUngappedStats::good_init_extends`: the sum of `init_hitlist->total` over the
    /// word finder calls (`Blast_UngappedStatsUpdate` in `BlastNaWordFinder`), an `Int4`
    /// that `CLocalBlast::GetNumExtensions` hands to the query batch mixer.
    pub good_init_extends: i32,
}

impl Scratch {
    pub fn new(ctx: &SearchContext, max_subject_length: i32) -> Scratch {
        Scratch {
            wf: WordFinderScratch::new(ctx.qb, ctx.lookup, ctx.params),
            gap: GapAligner::new(ctx, max_subject_length),
            tree_hsps: Vec::new(),
            good_init_extends: 0,
        }
    }
}

/// `BLAST_GetGappedScore` (score-only) on the sorted init hit list.
fn get_gapped_score(ctx: &SearchContext, subject: &SubjectBlk, scratch: &mut Scratch) -> HspList {
    let mut hsp_list = HspList::default();
    let qb = ctx.qb;
    let mut tree = IntervalTree::new(0, qb.length + 1, 0, subject.length + 1);
    scratch.tree_hsps.clear();
    let min_diag_sep = ctx.options.min_diag_separation;
    let init_hits = std::mem::take(&mut scratch.wf.init_hitlist);
    for init in &init_hits {
        let context = bsearch_context(qb, init.q_off);
        let cinfo = qb.contexts[context];
        let query_start = cinfo.query_offset;
        let u = init.ungapped;
        // s_AdjustInitialHSPOffsets: offsets relative to the context
        let q_start = u.q_start - query_start;
        let q_end = q_start + u.length;
        let s_start = u.s_start;
        let s_end = s_start + u.length;
        let tmp_hsp = Hsp::new(
            q_start,
            q_end,
            s_start,
            s_end,
            0,
            0,
            context as i32,
            cinfo.frame,
            1,
            u.score,
            None,
        );
        if tree.contains_hsp(&tmp_hsp, &scratch.tree_hsps, qb, min_diag_sep) {
            continue;
        }
        let cutoff = ctx.params.hit_cutoffs[context].cutoff_score;
        let query = qb.context_seq(context);
        // (score, box, gapped start) from the aligner the task uses
        let (score, qs, qe, ss, se, seed_q, seed_s) = match &mut scratch.gap {
            GapAligner::Greedy(greedy) => {
                let q_off = q_start + u.length / 2;
                let s_off = s_start + u.length / 2;
                let r = greedy.align(
                    query,
                    subject.seq,
                    cinfo.query_length,
                    subject.length,
                    q_off,
                    s_off,
                    true,
                    false,
                );
                (
                    r.score,
                    r.query_start,
                    r.query_stop,
                    r.subject_start,
                    r.subject_stop,
                    r.seed_q,
                    r.seed_s,
                )
            }
            GapAligner::Dp(dp) => {
                // start inside the first byte of an 8-letter exact seed when the ungapped
                // alignment is long enough; s_BlastDynProgNtGappedAlignment rounds up to the
                // next 4-base subject boundary anyway
                let mut q_off = init.q_off - query_start;
                let mut s_off = init.s_off;
                if s_end >= s_off + 8 {
                    s_off += 3;
                    q_off += 3;
                }
                let r = dp.dyn_prog_nt_gapped_alignment(
                    &qb.buf,
                    query_start as usize + 1,
                    cinfo.query_length,
                    subject.seq,
                    subject.length,
                    q_off,
                    s_off,
                    u.score,
                    &ctx.sbp.matrix,
                );
                (
                    r.score,
                    r.query_start,
                    r.query_stop,
                    r.subject_start,
                    r.subject_stop,
                    q_off,
                    s_off,
                )
            }
        };
        if score >= cutoff {
            let new_hsp = Hsp::new(
                qs,
                qe,
                ss,
                se,
                seed_q,
                seed_s,
                context as i32,
                cinfo.frame,
                1,
                score,
                None,
            );
            scratch.tree_hsps.push(new_hsp.clone());
            let idx = scratch.tree_hsps.len() - 1;
            hsp_list.hsps.push(new_hsp);
            tree.add_hsp(idx, &scratch.tree_hsps, qb);
        }
    }
    scratch.wf.init_hitlist = init_hits;
    hsp_list
}

/// `s_BackupSubject` / `s_GetNextSubjectChunk`: the (offset, length) of every chunk of a
/// subject of `length` bases, in order. A subject up to `MAX_DBSEQ_LEN` is one chunk; longer
/// ones are cut every `MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP` bases with the overlap repeated.
/// Both constants are multiples of 4, so the C's `residual` is always 0 and every chunk
/// starts on a packed byte.
pub fn subject_chunks(length: i32) -> Vec<(i32, i32)> {
    let mut chunks = Vec::new();
    let mut next = 0i32;
    let full_right = length;
    while next < full_right {
        let residual = next % 4;
        let offset = next - residual;
        let chunk_len;
        // (the C adds Int4s; a subject within 5 Mbp of 2^31 would overflow there, so compare
        // in 64 bits)
        if (offset as i64 + MAX_DBSEQ_LEN as i64) < full_right as i64 {
            chunk_len = MAX_DBSEQ_LEN;
            next = offset + MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP;
        } else {
            chunk_len = full_right - offset;
            next = full_right;
        }
        chunks.push((offset, chunk_len));
    }
    chunks
}

/// One chunk of `s_BlastSearchEngineCore`: the word finder and the gapped scoring on the
/// chunk at `offset` (`chunk_len` bases) of a subject, whose packed bytes start `seq`, then
/// `Blast_HSPListPurgeHSPsWithCommonEndpoints`, the sort by score and
/// `Blast_HSPListAdjustOffsets`. Returns the chunk's list in subject coordinates, or `None`
/// when nothing scored (the caller's `continue`). Adds the saved initial hits to
/// `scratch.good_init_extends`.
pub fn search_chunk(
    ctx: &SearchContext,
    seq: &[u8],
    offset: i32,
    chunk_len: i32,
    oid: i32,
    scratch: &mut Scratch,
) -> Option<HspList> {
    let wctx = ctx.word_finder_ctx();
    let ranges = [(0, chunk_len)];
    let blk = SubjectBlk {
        seq,
        length: chunk_len,
        seq_ranges: &ranges,
    };
    na_word_finder(&wctx, &blk, &mut scratch.wf);
    scratch.good_init_extends = scratch
        .good_init_extends
        .wrapping_add(scratch.wf.init_hitlist.len() as i32);
    if scratch.wf.init_hitlist.is_empty() {
        return None;
    }
    let mut hsp_list = get_gapped_score(ctx, &blk, scratch);
    hsp_list.purge_common_endpoints(true);
    hsp_list.sort_by_score();
    if hsp_list.hsps.is_empty() {
        return None;
    }
    hsp_list.oid = oid;
    hsp_list.adjust_offsets(offset);
    Some(hsp_list)
}

/// `Blast_HSPListsMerge` of one chunk's list (from `search_chunk`, chunk at `offset`) into the
/// subject's combined list. Chunks must be merged in order.
pub fn merge_chunk(
    ctx: &SearchContext,
    hsp_list: HspList,
    combined: &mut Option<HspList>,
    offset: i32,
) {
    let overlap = if offset == 0 { 0 } else { DBSEQ_CHUNK_OVERLAP };
    merge_chunk_lists(
        hsp_list,
        combined,
        offset,
        overlap,
        ctx.options.gapped_calculation,
    );
}

/// The end of `s_BlastSearchEngineCore` after the chunk loop: preliminary e-values and
/// reaping. `None` when nothing survives.
pub fn finish_subject(ctx: &SearchContext, combined: Option<HspList>) -> Option<HspList> {
    let mut out = combined?;
    if out.hsps.is_empty() {
        return None;
    }
    out.get_evalues(&ctx.sbp.kbp_gap_std, &ctx.eff_searchsp, ctx.sbp.round_down);
    out.reap_by_evalue(ctx.params.prelim_evalue);
    if out.hsps.is_empty() {
        return None;
    }
    Some(out)
}

/// `s_BlastHSPCollectorRun`: split by query and update the per-query hit lists.
pub fn collect(results: &mut HspResults, hsp_list: HspList, prelim_hitlist_size: i32) {
    let nq = results.hitlists.len();
    if nq > 1 {
        let mut per_query: Vec<Option<HspList>> = vec![None; nq];
        let oid = hsp_list.oid;
        for hsp in hsp_list.hsps {
            let qi = (hsp.context / 2) as usize;
            let l = per_query[qi].get_or_insert_with(|| HspList {
                oid,
                query_index: qi as i32,
                ..Default::default()
            });
            l.hsps.push(hsp);
        }
        for (qi, l) in per_query.into_iter().enumerate() {
            if let Some(l) = l {
                let hl = results.hitlists[qi]
                    .get_or_insert_with(|| super::hsp::HitList::new(prelim_hitlist_size));
                hl.update(l);
            }
        }
    } else if !hsp_list.hsps.is_empty() {
        let mut l = hsp_list;
        l.query_index = 0;
        let hl = results.hitlists[0]
            .get_or_insert_with(|| super::hsp::HitList::new(prelim_hitlist_size));
        hl.update(l);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chunk_geometry_matches_the_c_loop() {
        let step = MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP;
        for len in [
            1,
            4,
            4999,
            MAX_DBSEQ_LEN - 1,
            MAX_DBSEQ_LEN,
            MAX_DBSEQ_LEN + 1,
            step,
            step + 1,
            2 * step,
            2 * step + DBSEQ_CHUNK_OVERLAP,
            2 * step + DBSEQ_CHUNK_OVERLAP + 1,
            248_956_422,
            i32::MAX - 1,
        ] {
            let chunks = subject_chunks(len);
            // every chunk starts on a packed byte, the first at 0, each next one `step` on
            for (k, (offset, clen)) in chunks.iter().enumerate() {
                assert_eq!(*offset as i64, k as i64 * step as i64, "length {len}");
                assert_eq!(offset % 4, 0);
                assert!(*clen > 0 && *clen <= MAX_DBSEQ_LEN);
            }
            let (last_off, last_len) = *chunks.last().unwrap();
            assert_eq!(last_off + last_len, len, "length {len}");
        }
        assert_eq!(subject_chunks(MAX_DBSEQ_LEN).len(), 1);
        assert_eq!(subject_chunks(MAX_DBSEQ_LEN + 1).len(), 2);
    }
}
