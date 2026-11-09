//! Interval tree over HSPs, used to skip gapped extensions of seeds already inside a
//! higher-scoring alignment.
//!
//! Source: BLAST+ 2.17.0 `algo/blast/core/blast_itree.c` (`Blast_IntervalTreeInit`,
//! `Blast_IntervalTreeReset`, `s_IntervalNodeInit`, `s_IntervalRootNodeInit`,
//! `s_GetQueryStrandOffset`, `s_HSPsHaveCommonEndpoint`, `s_MidpointTreeHasHSPEndpoint`,
//! `s_IntervalTreeHasHSPEndpoint`, `BlastIntervalTreeAddHSP` (`eQueryAndSubject`),
//! `s_HSPIsContained`, `s_MidpointTreeContainsHSP`, `BlastIntervalTreeContainsHSP`) with
//! `CONTAINED_IN_HSP` (blast_hits_priv.h) and `MB_HSP_CLOSE` (blast_gapalign_priv.h).
//!
//! HSPs are referenced by index into the caller's list (the C stores pointers).

use super::hsp::Hsp;
use super::query::QueryBlock;

#[derive(Debug, Clone, Copy, Default)]
struct Node {
    leftptr: i32,
    midptr: i32,
    rightptr: i32,
    leftend: i32,
    rightend: i32,
    /// Index of the HSP in a leaf, -1 for internal nodes.
    hsp: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Dir {
    Left,
    Right,
    Neither,
}

pub struct IntervalTree {
    nodes: Vec<Node>,
    s_min: i32,
    s_max: i32,
}

/// `s_GetQueryStrandOffset`.
pub fn query_strand_offset(qb: &QueryBlock, context: i32) -> i32 {
    let mut c = context as usize;
    while c != 0 {
        let frame = qb.contexts[c].frame;
        if frame == 0 || frame.signum() != qb.contexts[c - 1].frame.signum() {
            break;
        }
        c -= 1;
    }
    qb.contexts[c].query_offset
}

#[inline]
fn contained_in_hsp(a: i32, b: i32, c: i32, d: i32, e: i32, f: i32) -> bool {
    (a <= c && b >= c) && (d <= f && e >= f)
}

#[inline]
fn mb_hsp_close(q1: i32, s1: i32, q2: i32, s2: i32, c: i32) -> bool {
    ((q1 - s1) - (q2 - s2)).abs() < c
}

/// `s_HSPIsContained`.
fn hsp_is_contained(
    in_hsp: &Hsp,
    in_q_start: i32,
    tree_hsp: &Hsp,
    tree_q_start: i32,
    min_diag_separation: i32,
) -> bool {
    if in_q_start != tree_q_start {
        return false;
    }
    if in_hsp.score <= tree_hsp.score
        && in_hsp.subject.frame.signum() == tree_hsp.subject.frame.signum()
        && contained_in_hsp(
            tree_hsp.query.offset,
            tree_hsp.query.end,
            in_hsp.query.offset,
            tree_hsp.subject.offset,
            tree_hsp.subject.end,
            in_hsp.subject.offset,
        )
        && contained_in_hsp(
            tree_hsp.query.offset,
            tree_hsp.query.end,
            in_hsp.query.end,
            tree_hsp.subject.offset,
            tree_hsp.subject.end,
            in_hsp.subject.end,
        )
    {
        if min_diag_separation == 0 {
            return true;
        }
        if mb_hsp_close(
            tree_hsp.query.offset,
            tree_hsp.subject.offset,
            in_hsp.query.offset,
            in_hsp.subject.offset,
            min_diag_separation,
        ) || mb_hsp_close(
            tree_hsp.query.end,
            tree_hsp.subject.end,
            in_hsp.query.end,
            in_hsp.subject.end,
            min_diag_separation,
        ) {
            return true;
        }
    }
    false
}

/// `s_HSPsHaveCommonEndpoint`: Some(true) = the tree HSP is better (or equal), Some(false) =
/// the input is better, None = no common endpoint.
fn have_common_endpoint(
    in_hsp: &Hsp,
    in_q_start: i32,
    tree_hsp: &Hsp,
    tree_q_start: i32,
    left: bool,
) -> Option<bool> {
    if in_q_start != tree_q_start {
        return None;
    }
    if in_hsp.subject.frame.signum() != tree_hsp.subject.frame.signum() {
        return None;
    }
    let m = if left {
        in_hsp.query.offset == tree_hsp.query.offset
            && in_hsp.subject.offset == tree_hsp.subject.offset
    } else {
        in_hsp.query.end == tree_hsp.query.end && in_hsp.subject.end == tree_hsp.subject.end
    };
    if !m {
        return None;
    }
    if in_hsp.score > tree_hsp.score {
        return Some(false);
    }
    if in_hsp.score < tree_hsp.score {
        return Some(true);
    }
    let in_q = in_hsp.query.end - in_hsp.query.offset;
    let tr_q = tree_hsp.query.end - tree_hsp.query.offset;
    if in_q > tr_q {
        return Some(true);
    }
    if in_q < tr_q {
        return Some(false);
    }
    let in_s = in_hsp.subject.end - in_hsp.subject.offset;
    let tr_s = tree_hsp.subject.end - tree_hsp.subject.offset;
    if in_s > tr_s {
        return Some(true);
    }
    if in_s < tr_s {
        return Some(false);
    }
    Some(true)
}

impl IntervalTree {
    /// `Blast_IntervalTreeInit`.
    pub fn new(q_start: i32, q_end: i32, s_start: i32, s_end: i32) -> IntervalTree {
        let mut t = IntervalTree {
            nodes: Vec::with_capacity(100),
            s_min: s_start,
            s_max: s_end,
        };
        t.root_node_init(q_start, q_end);
        t
    }

    /// `Blast_IntervalTreeReset`.
    pub fn reset(&mut self) {
        self.nodes.truncate(1);
        let r = &mut self.nodes[0];
        r.leftptr = 0;
        r.midptr = 0;
        r.rightptr = 0;
        r.hsp = -1;
    }

    fn node_init(&mut self, parent: i32, dir: Dir) -> i32 {
        let idx = self.nodes.len() as i32;
        self.nodes.push(Node {
            hsp: -1,
            ..Default::default()
        });
        if dir == Dir::Neither {
            return idx;
        }
        let p = self.nodes[parent as usize];
        let midpt = ((p.leftend as i64 + p.rightend as i64) / 2) as i32;
        let n = &mut self.nodes[idx as usize];
        if dir == Dir::Left {
            n.leftend = p.leftend;
            n.rightend = midpt;
        } else {
            n.leftend = midpt + 1;
            n.rightend = p.rightend;
        }
        idx
    }

    fn root_node_init(&mut self, region_start: i32, region_end: i32) -> i32 {
        let idx = self.node_init(0, Dir::Neither);
        let n = &mut self.nodes[idx as usize];
        n.leftend = region_start;
        n.rightend = region_end;
        idx
    }

    #[inline]
    fn middle(&self, idx: i32) -> i64 {
        let n = &self.nodes[idx as usize];
        (n.leftend as i64 + n.rightend as i64) / 2
    }

    /// `s_MidpointTreeHasHSPEndpoint`.
    fn midpoint_tree_has_endpoint(
        &mut self,
        root_index: i32,
        in_hsp: &Hsp,
        in_q_start: i32,
        left: bool,
        hsps: &[Hsp],
    ) -> bool {
        let mut root = root_index;
        let target_offset = if left {
            in_hsp.subject.offset
        } else {
            in_hsp.subject.end
        };
        loop {
            // midpoint list
            let mut tmp_index = self.nodes[root as usize].midptr;
            let mut list_node = root;
            while tmp_index != 0 {
                let next = self.nodes[tmp_index as usize];
                let best = have_common_endpoint(
                    in_hsp,
                    in_q_start,
                    &hsps[next.hsp as usize],
                    next.leftptr,
                    left,
                );
                let after = next.midptr;
                match best {
                    Some(true) => return true,
                    Some(false) => self.nodes[list_node as usize].midptr = after,
                    None => {}
                }
                // list_node advances to the examined node even if it was unlinked (C does this)
                list_node = tmp_index;
                tmp_index = after;
            }
            let midpt = self.middle(root);
            let mut tmp_index = 0;
            if (target_offset as i64) < midpt {
                tmp_index = self.nodes[root as usize].leftptr;
            } else if (target_offset as i64) > midpt {
                tmp_index = self.nodes[root as usize].rightptr;
            }
            if tmp_index == 0 {
                return false;
            }
            let next = self.nodes[tmp_index as usize];
            if next.hsp >= 0 {
                let best = have_common_endpoint(
                    in_hsp,
                    in_q_start,
                    &hsps[next.hsp as usize],
                    next.leftptr,
                    left,
                );
                match best {
                    Some(true) => return true,
                    Some(false) => {
                        if (target_offset as i64) < midpt {
                            self.nodes[root as usize].leftptr = 0;
                        } else if (target_offset as i64) > midpt {
                            self.nodes[root as usize].rightptr = 0;
                        }
                        return false;
                    }
                    None => {}
                }
                break;
            }
            root = tmp_index;
        }
        false
    }

    /// `s_IntervalTreeHasHSPEndpoint`.
    fn has_endpoint(&mut self, in_hsp: &Hsp, in_q_start: i32, left: bool, hsps: &[Hsp]) -> bool {
        let mut root = 0i32;
        let target_offset = if left {
            in_q_start + in_hsp.query.offset
        } else {
            in_q_start + in_hsp.query.end
        };
        loop {
            let mid = self.nodes[root as usize].midptr;
            if mid != 0 && self.midpoint_tree_has_endpoint(mid, in_hsp, in_q_start, left, hsps) {
                return true;
            }
            let midpt = self.middle(root);
            let mut tmp_index = 0;
            if (target_offset as i64) < midpt {
                tmp_index = self.nodes[root as usize].leftptr;
            } else if (target_offset as i64) > midpt {
                tmp_index = self.nodes[root as usize].rightptr;
            }
            if tmp_index == 0 {
                return false;
            }
            let next = self.nodes[tmp_index as usize];
            if next.hsp >= 0 {
                let best = have_common_endpoint(
                    in_hsp,
                    in_q_start,
                    &hsps[next.hsp as usize],
                    next.leftptr,
                    left,
                );
                match best {
                    Some(true) => return true,
                    Some(false) => {
                        if (target_offset as i64) < midpt {
                            self.nodes[root as usize].leftptr = 0;
                        } else if (target_offset as i64) > midpt {
                            self.nodes[root as usize].rightptr = 0;
                        }
                        return false;
                    }
                    None => {}
                }
                break;
            }
            root = tmp_index;
        }
        false
    }

    /// `BlastIntervalTreeAddHSP` with `eQueryAndSubject`. `hsp_index` indexes `hsps`.
    pub fn add_hsp(&mut self, hsp_index: usize, hsps: &[Hsp], qb: &QueryBlock) {
        let hsp = &hsps[hsp_index];
        let query_start = query_strand_offset(qb, hsp.context);
        let mut region_start = query_start + hsp.query.offset;
        let mut region_end = query_start + hsp.query.end;
        if self.has_endpoint(hsp, query_start, true, hsps) {
            return;
        }
        if self.has_endpoint(hsp, query_start, false, hsps) {
            return;
        }
        let mut index_subject_range = false;
        let mut root_index = 0i32;
        let new_index = self.node_init(0, Dir::Neither);
        {
            let n = &mut self.nodes[new_index as usize];
            n.leftptr = query_start;
            n.midptr = 0;
            n.hsp = hsp_index as i32;
        }
        loop {
            let middle = self.middle(root_index);
            let (old_index, which_half);
            if (region_end as i64) < middle {
                if self.nodes[root_index as usize].leftptr == 0 {
                    self.nodes[root_index as usize].leftptr = new_index;
                    return;
                }
                old_index = self.nodes[root_index as usize].leftptr;
                if self.nodes[old_index as usize].hsp < 0 {
                    root_index = old_index;
                    continue;
                }
                which_half = Dir::Left;
            } else if (region_start as i64) > middle {
                if self.nodes[root_index as usize].rightptr == 0 {
                    self.nodes[root_index as usize].rightptr = new_index;
                    return;
                }
                old_index = self.nodes[root_index as usize].rightptr;
                if self.nodes[old_index as usize].hsp < 0 {
                    root_index = old_index;
                    continue;
                }
                which_half = Dir::Right;
            } else {
                if index_subject_range {
                    let m = self.nodes[root_index as usize].midptr;
                    self.nodes[new_index as usize].midptr = m;
                    self.nodes[root_index as usize].midptr = new_index;
                    return;
                }
                index_subject_range = true;
                if self.nodes[root_index as usize].midptr == 0 {
                    let (smin, smax) = (self.s_min, self.s_max);
                    let mid_index = self.root_node_init(smin, smax);
                    self.nodes[root_index as usize].midptr = mid_index;
                }
                root_index = self.nodes[root_index as usize].midptr;
                region_start = hsp.subject.offset;
                region_end = hsp.subject.end;
                continue;
            }
            // split the leaf
            let mid_index = self.node_init(root_index, which_half);
            let old_hsp = &hsps[self.nodes[old_index as usize].hsp as usize];
            if which_half == Dir::Left {
                self.nodes[root_index as usize].leftptr = mid_index;
            } else {
                self.nodes[root_index as usize].rightptr = mid_index;
            }
            let (old_region_start, old_region_end) = if index_subject_range {
                (old_hsp.subject.offset, old_hsp.subject.end)
            } else {
                let lp = self.nodes[old_index as usize].leftptr;
                (lp + old_hsp.query.offset, lp + old_hsp.query.end)
            };
            root_index = mid_index;
            let middle = self.middle(root_index);
            if (old_region_end as i64) < middle {
                self.nodes[mid_index as usize].leftptr = old_index;
            } else if (old_region_start as i64) > middle {
                self.nodes[mid_index as usize].rightptr = old_index;
            } else if index_subject_range {
                self.nodes[mid_index as usize].midptr = old_index;
            } else {
                let (smin, smax) = (self.s_min, self.s_max);
                let mid_index2 = self.root_node_init(smin, smax);
                let (ors, ore) = (old_hsp.subject.offset, old_hsp.subject.end);
                self.nodes[mid_index as usize].midptr = mid_index2;
                let middle2 = self.middle(mid_index2);
                if (ore as i64) < middle2 {
                    self.nodes[mid_index2 as usize].leftptr = old_index;
                } else if (ors as i64) > middle2 {
                    self.nodes[mid_index2 as usize].rightptr = old_index;
                } else {
                    self.nodes[mid_index2 as usize].midptr = old_index;
                }
            }
        }
    }

    /// `s_MidpointTreeContainsHSP`.
    fn midpoint_tree_contains(
        &self,
        root_index: i32,
        in_hsp: &Hsp,
        in_q_start: i32,
        min_diag_separation: i32,
        hsps: &[Hsp],
    ) -> bool {
        let mut node = root_index;
        let region_start = in_hsp.subject.offset;
        let region_end = in_hsp.subject.end;
        while self.nodes[node as usize].hsp < 0 {
            let mut tmp_index = self.nodes[node as usize].midptr;
            while tmp_index != 0 {
                let tn = self.nodes[tmp_index as usize];
                if hsp_is_contained(
                    in_hsp,
                    in_q_start,
                    &hsps[tn.hsp as usize],
                    tn.leftptr,
                    min_diag_separation,
                ) {
                    return true;
                }
                tmp_index = tn.midptr;
            }
            let middle = self.middle(node);
            let mut tmp_index = 0;
            if (region_end as i64) < middle {
                tmp_index = self.nodes[node as usize].leftptr;
            } else if (region_start as i64) > middle {
                tmp_index = self.nodes[node as usize].rightptr;
            }
            if tmp_index == 0 {
                return false;
            }
            node = tmp_index;
        }
        let n = self.nodes[node as usize];
        hsp_is_contained(
            in_hsp,
            in_q_start,
            &hsps[n.hsp as usize],
            n.leftptr,
            min_diag_separation,
        )
    }

    /// `BlastIntervalTreeContainsHSP`.
    pub fn contains_hsp(
        &self,
        hsp: &Hsp,
        hsps: &[Hsp],
        qb: &QueryBlock,
        min_diag_separation: i32,
    ) -> bool {
        let query_start = query_strand_offset(qb, hsp.context);
        let region_start = query_start + hsp.query.offset;
        let region_end = query_start + hsp.query.end;
        let mut node = 0i32;
        while self.nodes[node as usize].hsp < 0 {
            let mid = self.nodes[node as usize].midptr;
            if mid > 0
                && self.midpoint_tree_contains(mid, hsp, query_start, min_diag_separation, hsps)
            {
                return true;
            }
            let middle = self.middle(node);
            let mut tmp_index = 0;
            if (region_end as i64) < middle {
                tmp_index = self.nodes[node as usize].leftptr;
            } else if (region_start as i64) > middle {
                tmp_index = self.nodes[node as usize].rightptr;
            }
            if tmp_index == 0 {
                return false;
            }
            node = tmp_index;
        }
        let n = self.nodes[node as usize];
        hsp_is_contained(
            hsp,
            query_start,
            &hsps[n.hsp as usize],
            n.leftptr,
            min_diag_separation,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::options::{BlastOptions, Strand};
    use crate::engine::query::QuerySet;

    #[test]
    fn contained_and_endpoint_rules() {
        let qs = QuerySet::from_fasta_bytes(&{
            let mut v = b">a\n".to_vec();
            v.extend(std::iter::repeat_n(b'A', 1000));
            v
        })
        .unwrap();
        let mut o = BlastOptions::megablast();
        o.dust = None;
        o.strand = Strand::Plus;
        let qb = crate::engine::query::QueryBlock::build(&qs, &o);
        let big = Hsp::new(100, 500, 1000, 1400, 0, 0, 0, 1, 1, 300, None);
        let inside = Hsp::new(200, 300, 1100, 1200, 0, 0, 0, 1, 1, 50, None);
        let off_diag = Hsp::new(200, 300, 1150, 1250, 0, 0, 0, 1, 1, 50, None);
        let hsps = vec![big, inside, off_diag];
        let mut t = IntervalTree::new(0, qb.length + 1, 0, 5000);
        t.add_hsp(0, &hsps, &qb);
        assert!(t.contains_hsp(&hsps[1], &hsps, &qb, 6));
        // same box but 50 diagonals away: not "close" with min_diag_separation 6
        assert!(!t.contains_hsp(&hsps[2], &hsps, &qb, 6));
        assert!(t.contains_hsp(&hsps[2], &hsps, &qb, 0));
        // same start, lower score: add refuses it (common endpoint with a better HSP)
        let same_start = Hsp::new(100, 200, 1000, 1100, 0, 0, 0, 1, 1, 60, None);
        let hsps2 = vec![hsps[0].clone(), same_start];
        let n_before = t.nodes.len();
        t.add_hsp(1, &hsps2, &qb);
        assert_eq!(t.nodes.len(), n_before);
    }
}
