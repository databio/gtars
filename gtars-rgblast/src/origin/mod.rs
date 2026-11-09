//! `rgblast origin`: species first, then fan-out.
//!
//! Stage 1 searches one representative per group (taxid) of the store selection, whose
//! deflines still list every selected genome holding each sequence. Each query's candidate
//! groups (stage 1 genome score within `--fanout-percent` of its best) are then searched in
//! stage 2: the candidates' sequences that stage 1 did not cover. Every stage searches under
//! the full selection's statistics (rgblastn's `STATS_NSEQ`/`STATS_TOTLEN` override), so
//! E-values equal a full search. Staged mode needs rgblastn; blastn on a database the user
//! built runs one exhaustive search.
//!
//! When staged and exhaustive results can differ (the contract, also in `--help`):
//! 1. the best sequence is in no candidate group's representative (no stage 1 hit at all gives
//!    `none` + `reps-only`);
//! 2. the true group's representative scores below the fan-out cut;
//! 3. `max_target_seqs` truncation sees different subject pools (flag `truncated`).
//!
//! E-values are not a difference.

pub mod queries;
pub mod reps;

use crate::error::Error;
use crate::plan::{DbTotals, ResolvedSelection, SubjectPlan};
use crate::progress::ProgressEvent;
use crate::report::Hit;
use crate::report::fanout::{DbMeta, fanout};
use crate::report::parse::read_hits_file;
use crate::report::rank::{QueryInput, RankOptions, Report, ReportMeta, StageHits, build_report};
use crate::report::taxonomy::Taxonomy;
use crate::run::outfmt::ORIGIN_OUTFMT;
use crate::run::{BlastRunner, EngineOutput, RunSpec, SearchTarget};
use crate::store::OpenedStore;
use queries::{QueryRec, write_queries};
use reps::{GroupKey, RepColl, RepsTable};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// `--exhaustive` / `--staged` / neither (decided by the cost prediction).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum ModeChoice {
    #[default]
    Auto,
    Staged,
    Exhaustive,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OriginOptions {
    pub mode: ModeChoice,
    pub fanout_percent: f64,
    /// `--rep` collections (digests).
    pub reps: Vec<String>,
    pub rep_namespace: String,
    pub rank: RankOptions,
}

impl Default for OriginOptions {
    fn default() -> Self {
        OriginOptions {
            mode: ModeChoice::Auto,
            fanout_percent: 20.0,
            reps: Vec::new(),
            rep_namespace: "reference".into(),
            rank: RankOptions::default(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StageMeta {
    pub stage: String,
    pub queries: usize,
    pub subjects: usize,
    pub subjects_bases: u64,
    /// The database searched (blastn).
    pub db: Option<PathBuf>,
    pub exit_code: i32,
    pub args: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct QueryCandidates {
    pub query: String,
    pub groups: Vec<String>,
    /// Collections a stage 2 hit fans out to (digests).
    pub scope: Vec<String>,
}

/// `origin.meta.json`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OriginMeta {
    pub mode: String,
    pub stats: DbTotals,
    /// The representatives (store searches).
    pub reps: Option<RepsTable>,
    pub stages: Vec<StageMeta>,
    pub candidates: Vec<QueryCandidates>,
    /// The stderr summary lines.
    pub notes: Vec<String>,
}

/// What origin searches.
pub enum Source<'a> {
    /// rgblastn on a store selection: staged or exhaustive.
    Store {
        store: &'a mut OpenedStore,
        selection: &'a ResolvedSelection,
        /// The selection's plan (deflines, totals; no sequence reads).
        plan: &'a SubjectPlan,
    },
    /// blastn on a database the user built: exhaustive.
    Db(&'a Path),
}

/// What an origin run needs.
pub struct OriginRun<'a> {
    pub runner: &'a dyn BlastRunner,
    pub source: Source<'a>,
    /// The collections and members searched.
    pub meta: DbMeta,
    /// The full selection's (or database's) statistics.
    pub stats: DbTotals,
    pub queries: &'a [QueryRec],
    /// Pass-through arguments (no `-outfmt`; origin adds its own).
    pub args: Vec<String>,
    pub threads_added: bool,
    pub max_target_seqs: usize,
    pub rundir: &'a Path,
}

fn group_key(
    meta: &DbMeta,
    coll: usize,
    tax: &mut Option<&mut Taxonomy>,
    rank: &Option<String>,
) -> GroupKey {
    match meta.collections[coll].taxid.filter(|t| *t > 0) {
        None => GroupKey::Unknown(meta.collections[coll].collection.clone()),
        Some(t) => {
            let t = t as u32;
            match tax {
                Some(tx) => {
                    let t = tx.resolve(t);
                    GroupKey::Taxid(
                        rank.as_ref()
                            .and_then(|r| tx.ancestor_at(t, r))
                            .unwrap_or(t),
                    )
                }
                None => GroupKey::Taxid(t),
            }
        }
    }
}

/// The representative input: every plan collection with its sequences as (OID, length).
pub fn rep_input(
    plan: &SubjectPlan,
    meta: &DbMeta,
    rep_namespace: &str,
    mut tax: Option<&mut Taxonomy>,
    group_rank: &Option<String>,
) -> Result<Vec<RepColl>, Error> {
    let mut seqs: Vec<Vec<(u32, u64)>> = vec![Vec::new(); plan.collections.len()];
    for (oid, s) in plan.subjects.iter().enumerate() {
        for o in &s.occurrences {
            let v = &mut seqs[o.collection as usize];
            if v.last().map(|x| x.0) != Some(oid as u32) {
                v.push((oid as u32, s.length as u64));
            }
        }
    }
    let keys: Vec<GroupKey> = (0..plan.collections.len())
        .map(|i| group_key(meta, i, &mut tax, group_rank))
        .collect();
    if let Some(tx) = tax.as_deref_mut() {
        let want: HashSet<u32> = keys
            .iter()
            .filter_map(|k| match k {
                GroupKey::Taxid(t) => Some(*t),
                _ => None,
            })
            .collect();
        tx.load_names(&want)?;
    }
    Ok(plan
        .collections
        .iter()
        .enumerate()
        .map(|(i, c)| RepColl {
            digest: c.digest.clone(),
            display: c.display.clone(),
            group_name: match &keys[i] {
                GroupKey::Taxid(t) => tax
                    .as_deref()
                    .and_then(|tx| tx.name(*t).map(str::to_string))
                    .or_else(|| c.species.clone())
                    .unwrap_or_else(|| t.to_string()),
                GroupKey::Unknown(_) => c.display.clone(),
            },
            group: keys[i].clone(),
            n_fraction: c.n_fraction,
            in_namespace: c.aliases.iter().any(|(ns, _)| ns == rep_namespace),
            seqs: std::mem::take(&mut seqs[i]),
        })
        .collect())
}

fn qindex(id: &str) -> Option<usize> {
    id.strip_prefix('Q')?.parse::<usize>().ok()?.checked_sub(1)
}

fn gbp(b: u64) -> String {
    format!("{:.2}", b as f64 / 1e9)
}

fn plural(n: usize, one: &str, many: &str) -> String {
    format!("{n} {}", if n == 1 { one } else { many })
}

/// Run one stage on the source (`subjects` limits a store search to these digests); returns
/// its hits split by query index.
#[allow(clippy::too_many_arguments)]
fn run_stage(
    r: &mut OriginRun,
    stage: &str,
    subjects: Option<Arc<Vec<String>>>,
    stats: Option<DbTotals>,
    only: Option<&HashSet<usize>>,
    subjects_bases: u64,
    metas: &mut Vec<StageMeta>,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<Vec<Vec<Hit>>, Error> {
    let qfile = r.rundir.join(format!("{stage}.queries.fa"));
    write_queries(r.queries, &qfile, only)?;
    let out = r.rundir.join(format!("{stage}.tsv"));
    let nq = only.map(|o| o.len()).unwrap_or(r.queries.len());
    let n_subjects = subjects
        .as_ref()
        .map(|d| d.len())
        .unwrap_or(r.stats.num_seqs as usize);
    if !progress(&ProgressEvent::Stage {
        stage: stage.into(),
        queries: nq,
        subjects_gbp: subjects_bases as f64 / 1e9,
    }) {
        return Err(Error::Cancelled);
    }
    let mut args = r.args.clone();
    args.push("-outfmt".into());
    args.push(ORIGIN_OUTFMT.into());
    let (target, store) = match &mut r.source {
        Source::Store {
            store, selection, ..
        } => (
            SearchTarget::Store {
                selection: Box::new((*selection).clone()),
                subjects,
            },
            Some(&mut **store),
        ),
        Source::Db(db) => (SearchTarget::Db(db.to_path_buf()), None),
    };
    let spec = RunSpec {
        stage: stage.into(),
        target,
        args,
        stats,
        query: Some(qfile),
        output: EngineOutput::File(out.clone()),
        threads_added: r.threads_added,
    };
    let outcome = r.runner.run(&spec, store, progress)?;
    if outcome.exit_code != 0 {
        return Err(Error::Engine(format!(
            "{stage}: {} exited with status {}: {}",
            r.runner.describe(),
            outcome.exit_code,
            outcome.engine_stderr.join(" / ")
        )));
    }
    let parsed = read_hits_file(&out, &ORIGIN_OUTFMT[2..])?;
    let mut hits = vec![Vec::new(); r.queries.len()];
    for h in parsed.hits {
        if let Some(i) = qindex(&h.query).filter(|i| *i < hits.len()) {
            hits[i].push(h);
        }
    }
    metas.push(StageMeta {
        stage: stage.into(),
        queries: nq,
        subjects: n_subjects,
        subjects_bases,
        db: outcome.db.clone(),
        exit_code: outcome.exit_code,
        args: outcome.args.clone(),
    });
    Ok(hits)
}

/// Run origin: decide the mode, run the stages, build the report.
pub fn run(
    mut r: OriginRun,
    opts: &OriginOptions,
    mut tax: Option<&mut Taxonomy>,
    mut info: ReportMeta,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<(Report, OriginMeta), Error> {
    let meta = std::mem::take(&mut r.meta);
    let n_coll = meta.collections.len();
    let nq = r.queries.len();
    let stats = r.stats;
    let full_bases = stats.total_length;
    let plan = match &r.source {
        Source::Store { plan, .. } => Some(*plan),
        Source::Db(_) => None,
    };
    // the representatives (store searches only)
    let reps = match plan {
        Some(plan) => {
            let colls = rep_input(
                plan,
                &meta,
                &opts.rep_namespace,
                tax.as_deref_mut(),
                &opts.rank.group_rank,
            )?;
            let table = reps::choose(&colls, &opts.reps, plan.totals.total_length)?;
            Some((plan, colls, table))
        }
        None => None,
    };
    let mut notes = Vec::new();
    let staged = match (&reps, opts.mode) {
        (None, _) | (_, ModeChoice::Exhaustive) => false,
        (Some(_), ModeChoice::Staged) => true,
        (Some((_, _, t)), ModeChoice::Auto) => {
            let s = &t.summary;
            if !s.staged_genome {
                notes.push(format!(
                    "origin: a staged search would scan up to {:.0}% of the selection ({} groups for {} collections); running one full search",
                    s.worst_genome_bases as f64 / full_bases.max(1) as f64 * 100.0,
                    s.groups,
                    s.collections
                ));
            }
            s.staged_genome
        }
    };
    let mut metas: Vec<StageMeta> = Vec::new();
    let mut inputs: Vec<QueryInput> = r
        .queries
        .iter()
        .map(|q| QueryInput {
            name: q.name.clone(),
            qlen: Some(q.seq.len() as u64),
            max_target_seqs: r.max_target_seqs,
            ..Default::default()
        })
        .collect();
    let mut candidates = Vec::new();
    match reps.as_ref().filter(|_| staged) {
        None => {
            let hits = run_stage(
                &mut r, "full", None, None, None, full_bases, &mut metas, progress,
            )?;
            for (i, h) in hits.into_iter().enumerate() {
                inputs[i].runs.push(StageHits {
                    stage: "full".into(),
                    hits: h,
                    scope: None,
                });
                inputs[i].searched = (n_coll, n_coll);
            }
            notes.insert(
                0,
                format!(
                    "origin: exhaustive; 1 search of {} collections, {} Gbp, {}",
                    n_coll,
                    gbp(full_bases),
                    plural(nq, "query", "queries")
                ),
            );
            info.mode = "exhaustive".into();
        }
        Some((plan, colls, table)) => {
            info.mode = "staged".into();
            let s = &table.summary;
            let rep_oids = table.rep_seqs(colls);
            let mut rep_digests: Vec<String> = rep_oids
                .iter()
                .map(|o| plan.subjects[*o as usize].digest.clone())
                .collect();
            rep_digests.sort();
            let reps_n = table.groups.len();
            let mut st1 = run_stage(
                &mut r,
                "stage1",
                Some(Arc::new(rep_digests)),
                Some(stats),
                None,
                s.reps_bases,
                &mut metas,
                progress,
            )?;
            let rep_colls: HashSet<usize> = table.groups.iter().map(|g| g.rep_index).collect();
            let group_of: Vec<&GroupKey> = colls.iter().map(|c| &c.group).collect();
            let no_ids = HashMap::new();
            let ro = &opts.rank;
            let mut cand_colls: Vec<Option<HashSet<usize>>> = vec![None; nq];
            let mut cand_names: Vec<Vec<String>> = vec![Vec::new(); nq];
            for (i, hits) in st1.iter().enumerate() {
                let mut best: HashMap<&GroupKey, f64> = HashMap::new();
                for h in hits {
                    if ro.min_pident.is_some_and(|m| h.pident < m)
                        || ro.min_qcov.is_some_and(|m| h.qcov() < m)
                        || ro.max_evalue.is_some_and(|m| h.evalue > m)
                    {
                        continue;
                    }
                    let sc = h.score.unwrap_or(h.bits);
                    for p in fanout(&meta, h, None, &no_ids)?.pairs {
                        let e = best.entry(group_of[p.coll]).or_insert(f64::NEG_INFINITY);
                        *e = e.max(sc);
                    }
                }
                let Some(top) = best.values().copied().reduce(f64::max) else {
                    continue;
                };
                let cut = top * (1.0 - opts.fanout_percent / 100.0);
                let groups: BTreeSet<&GroupKey> = best
                    .iter()
                    .filter(|(_, v)| **v >= cut)
                    .map(|(k, _)| *k)
                    .collect();
                cand_names[i] = table
                    .groups
                    .iter()
                    .filter(|g| g.key.as_ref().is_some_and(|k| groups.contains(k)))
                    .map(|g| g.taxon.clone())
                    .collect();
                cand_colls[i] = Some(
                    (0..n_coll)
                        .filter(|c| groups.contains(group_of[*c]))
                        .collect(),
                );
            }
            let with_hits = cand_colls.iter().filter(|c| c.is_some()).count();
            let mut line = format!(
                "origin: staged; stage 1: {} representatives, {} Gbp, {}, {with_hits} with hits",
                reps_n,
                gbp(s.reps_bases),
                plural(nq, "query", "queries"),
            );
            // R = the candidate collections' sequences minus every representative sequence
            let per_query: Vec<HashSet<u32>> = cand_colls
                .iter()
                .map(|c| {
                    c.iter()
                        .flatten()
                        .flat_map(|ci| colls[*ci].seqs.iter().map(|x| x.0))
                        .filter(|o| !rep_oids.contains(o))
                        .collect()
                })
                .collect();
            let r_all: BTreeSet<u32> = per_query.iter().flatten().copied().collect();
            let only: HashSet<usize> = (0..nq).filter(|i| !per_query[*i].is_empty()).collect();
            let mut st2: Option<Vec<Vec<Hit>>> = None;
            if !r_all.is_empty() {
                let bases: u64 = r_all
                    .iter()
                    .map(|o| plan.subjects[*o as usize].length as u64)
                    .sum();
                let mut digests: Vec<String> = r_all
                    .iter()
                    .map(|o| plan.subjects[*o as usize].digest.clone())
                    .collect();
                digests.sort();
                let names: BTreeSet<&String> =
                    only.iter().flat_map(|i| cand_names[*i].iter()).collect();
                line += &format!(
                    "; stage 2: {} remainder {} Gbp for {}",
                    names.into_iter().cloned().collect::<Vec<_>>().join(", "),
                    gbp(bases),
                    plural(only.len(), "query", "queries")
                );
                st2 = Some(run_stage(
                    &mut r,
                    "stage2",
                    Some(Arc::new(digests)),
                    Some(stats),
                    Some(&only),
                    bases,
                    &mut metas,
                    progress,
                )?);
            }
            let no_hit = cand_colls.iter().filter(|c| c.is_none()).count();
            if no_hit > 0 {
                line += &format!(
                    "; {} with no stage 1 hit",
                    plural(no_hit, "query", "queries")
                );
            }
            notes.insert(0, line);
            if no_hit > 0 {
                notes.push(format!(
                    "{} no hit among representatives; rerun with --exhaustive",
                    plural(no_hit, "query had", "queries had")
                ));
            }
            for (i, inp) in inputs.iter_mut().enumerate() {
                inp.runs.push(StageHits {
                    stage: "reps".into(),
                    hits: std::mem::take(&mut st1[i]),
                    scope: None,
                });
                if let Some(h) = st2.as_mut() {
                    inp.runs.push(StageHits {
                        stage: "fanout".into(),
                        hits: std::mem::take(&mut h[i]),
                        scope: cand_colls[i].clone(),
                    });
                }
                if cand_colls[i].is_none() {
                    inp.flags.push("reps-only".into());
                    inp.searched = (reps_n, n_coll);
                } else {
                    let mut searched = rep_colls.clone();
                    searched.extend(cand_colls[i].iter().flatten().copied());
                    inp.searched = (searched.len(), n_coll);
                }
                candidates.push(QueryCandidates {
                    query: r.queries[i].name.clone(),
                    groups: cand_names[i].clone(),
                    scope: cand_colls[i]
                        .iter()
                        .flatten()
                        .map(|c| meta.collections[*c].collection.clone())
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect(),
                });
            }
        }
    }
    let report = build_report(&meta, inputs, &opts.rank, tax, info)?;
    Ok((
        report,
        OriginMeta {
            mode: if staged { "staged" } else { "exhaustive" }.into(),
            stats,
            reps: reps.map(|(_, _, t)| t),
            stages: metas,
            candidates,
            notes,
        },
    ))
}
