//! Dropping the HSPs a report cannot use, as they stream in.
//!
//! A search can return tens of millions of HSPs (many per subject), but the report keeps one
//! hit per genome and origin's stage 1 one score per group. Every HSP of one subject for one
//! query fans out to the same genomes, so of each such run of HSPs only these can matter:
//! - the first one passing the filters that is best by the report's order (`rank::order`);
//!   among equals the earliest wins, as in `build_report`;
//! - the first one passing the filters with the highest score (stage 1's candidate cut);
//! - the first HSP, when none passes (the subject still counts toward `truncated`).
//!
//! The kept HSPs stay in their input order, so ties between subjects resolve as before and
//! the report is the same as on every HSP.

use super::parse::Hit;
use super::rank::{RankOptions, order};
use std::cmp::Ordering;

fn score(h: &Hit) -> f64 {
    h.score.unwrap_or(h.bits)
}

fn same_subject(a: &Hit, b: &Hit) -> bool {
    a.query == b.query && a.subject == b.subject && a.sallseqid == b.sallseqid
}

/// One subject's run of HSPs: the first, and the best and top passing ones with their
/// positions in the run.
struct Run {
    first: Hit,
    best: Option<(usize, Hit)>,
    top: Option<(usize, Hit)>,
    n: usize,
}

/// Keeps the HSPs a report under `opts` can use (see the module docs).
pub struct HitReducer {
    opts: RankOptions,
    run: Option<Run>,
    out: Vec<Hit>,
}

impl HitReducer {
    pub fn new(opts: &RankOptions) -> Self {
        HitReducer {
            opts: opts.clone(),
            run: None,
            out: Vec::new(),
        }
    }

    pub fn push(&mut self, h: Hit) {
        match &mut self.run {
            Some(r) if same_subject(&r.first, &h) => {
                let i = r.n;
                r.n += 1;
                Self::consider(&self.opts, r, i, h);
            }
            _ => {
                self.flush();
                let mut r = Run {
                    first: h.clone(),
                    best: None,
                    top: None,
                    n: 1,
                };
                Self::consider(&self.opts, &mut r, 0, h);
                self.run = Some(r);
            }
        }
    }

    fn consider(opts: &RankOptions, r: &mut Run, i: usize, h: Hit) {
        if !opts.passes(&h) {
            return;
        }
        let better = r
            .best
            .as_ref()
            .is_none_or(|(_, b)| order(b, &h) == Ordering::Greater);
        let higher = r.top.as_ref().is_none_or(|(_, t)| score(&h) > score(t));
        match (better, higher) {
            (true, true) => {
                r.top = Some((i, h.clone()));
                r.best = Some((i, h));
            }
            (true, false) => r.best = Some((i, h)),
            (false, true) => r.top = Some((i, h)),
            (false, false) => {}
        }
    }

    fn flush(&mut self) {
        let Some(r) = self.run.take() else { return };
        match (r.best, r.top) {
            (None, _) => self.out.push(r.first),
            (Some((bi, b)), Some((ti, t))) if ti != bi => {
                if ti < bi {
                    self.out.extend([t, b]);
                } else {
                    self.out.extend([b, t]);
                }
            }
            (Some((_, b)), _) => self.out.push(b),
        }
    }

    pub fn finish(mut self) -> Vec<Hit> {
        self.flush();
        self.out
    }
}

/// `hits` reduced to what a report under `opts` can use.
pub fn reduce_hits(hits: impl IntoIterator<Item = Hit>, opts: &RankOptions) -> Vec<Hit> {
    let mut r = HitReducer::new(opts);
    for h in hits {
        r.push(h);
    }
    r.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(q: &str, s: &str, bits: f64, score: f64, pident: f64, sstart: u64) -> Hit {
        Hit {
            query: q.into(),
            subject: Some(s.into()),
            pident,
            pident_s: format!("{pident:.3}"),
            length: 100,
            qstart: 1,
            qend: 100,
            sstart,
            send: sstart + 99,
            evalue: 1e-10,
            evalue_s: "1e-10".into(),
            bits,
            bits_s: format!("{bits}"),
            qlen: 100,
            score: Some(score),
            sallseqid: vec![format!("SQ.{s}")],
            salltitles: None,
            staxids: Vec::new(),
        }
    }

    #[test]
    fn keeps_best_top_and_order() {
        let hs = vec![
            hit("q", "a", 50.0, 30.0, 90.0, 1),
            hit("q", "a", 80.0, 40.0, 90.0, 5),
            hit("q", "b", 80.0, 40.0, 90.0, 1),
            hit("q", "a", 80.0, 40.0, 95.0, 9),
            hit("q", "a", 70.0, 45.0, 90.0, 3),
        ];
        let got = reduce_hits(hs.clone(), &RankOptions::default());
        // a: the bits 80 hit is best and top; b alone; a again (not consecutive, a new run):
        // best by order (bits 80) and top by score (45) are different hits, both kept
        assert_eq!(
            got,
            vec![hs[1].clone(), hs[2].clone(), hs[3].clone(), hs[4].clone()]
        );
    }

    #[test]
    fn filtered_subject_keeps_first() {
        let hs = vec![
            hit("q", "a", 50.0, 30.0, 80.0, 1),
            hit("q", "a", 60.0, 35.0, 85.0, 2),
            hit("q", "a", 90.0, 50.0, 99.0, 3),
        ];
        let opts = RankOptions {
            min_pident: Some(99.5),
            ..Default::default()
        };
        assert_eq!(reduce_hits(hs.clone(), &opts), vec![hs[0].clone()]);
        let opts = RankOptions {
            min_pident: Some(84.0),
            ..Default::default()
        };
        assert_eq!(reduce_hits(hs.clone(), &opts), vec![hs[2].clone()]);
    }

    #[test]
    fn top_score_kept_when_order_differs() {
        // higher raw score but lower bits (not possible from one engine, but kept anyway)
        let hs = vec![
            hit("q", "a", 90.0, 40.0, 90.0, 1),
            hit("q", "a", 80.0, 60.0, 90.0, 2),
        ];
        assert_eq!(reduce_hits(hs.clone(), &RankOptions::default()), hs);
    }
}
