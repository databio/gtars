//! Best hit per genome, ties, the taxon call and its status. The rules are fixed and total:
//! - best HSP per genome: bits desc, evalue asc, pident desc, length desc, sstart asc;
//! - top score T: raw `score` when every hit has it, else bits;
//! - tied genomes: score >= T x (1 - P/100), P = `--tie-percent`;
//! - status: `taxon` (one group), `lca` (several, taxonomy loaded), `multi` (several, no
//!   taxonomy), `none` (no HSP survives the filters);
//! - top genome: the first tied genome by pident desc, length desc, display bytewise.

use super::fanout::{DbMeta, fanout};
use super::parse::Hit;
use super::taxonomy::Taxonomy;
use crate::error::Error;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct RankOptions {
    pub tie_percent: f64,
    pub min_pident: Option<f64>,
    pub min_qcov: Option<f64>,
    pub max_evalue: Option<f64>,
    /// Fold taxids to their ancestor at this rank (needs a taxonomy).
    pub group_rank: Option<String>,
}

/// The hits of one engine run for one query, and the run's fan-out scope.
#[derive(Debug, Clone, Default)]
pub struct StageHits {
    /// "full", "reps", "fanout".
    pub stage: String,
    pub hits: Vec<Hit>,
    /// Collection indices an `SQ.`-only hit fans out to; None = all.
    pub scope: Option<HashSet<usize>>,
}

#[derive(Debug, Clone, Default)]
pub struct QueryInput {
    pub name: String,
    pub qlen: Option<u64>,
    pub runs: Vec<StageHits>,
    /// (genomes fully searched, selection size).
    pub searched: (usize, usize),
    pub flags: Vec<String>,
    pub max_target_seqs: usize,
}

/// A number (`9606`) or text (`unknown:asm1`, `9606,10090`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum TaxId {
    Num(u32),
    Text(String),
}

impl std::fmt::Display for TaxId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TaxId::Num(n) => write!(f, "{n}"),
            TaxId::Text(s) => f.write_str(s),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Call {
    pub taxid: TaxId,
    pub name: String,
    pub rank: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GenomeHit {
    pub genome: String,
    pub collection: String,
    pub seq_name: String,
    pub alt_names: Vec<String>,
    pub seq_digest: String,
    pub sstart: u64,
    pub send: u64,
    pub strand: String,
    pub pident: String,
    pub length: u64,
    pub qcov: f64,
    pub evalue: String,
    pub bits: String,
    pub tied: bool,
    pub stage: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Tied {
    pub genomes: usize,
    pub taxa: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Searched {
    pub genomes: usize,
    pub of: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TaxonGroup {
    pub taxid: TaxId,
    pub name: String,
    pub best_bits: String,
    pub best_pident: String,
    pub tied: usize,
    pub hit: usize,
    pub genomes: Vec<GenomeHit>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct QueryResult {
    pub query: String,
    pub qlen: Option<u64>,
    pub status: String,
    pub call: Option<Call>,
    pub top: Option<GenomeHit>,
    pub tied: Tied,
    pub searched: Searched,
    pub flags: Vec<String>,
    pub taxa: Vec<TaxonGroup>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ReportMeta {
    pub engine: Option<String>,
    pub db: Option<String>,
    pub mode: String,
    pub taxonomy: Option<String>,
    pub tie_percent: f64,
    /// For example unknown taxids met in the taxonomy.
    pub warnings: Vec<String>,
}

/// The `--format json` document.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Report {
    pub rgblast_report: u32,
    pub meta: ReportMeta,
    pub queries: Vec<QueryResult>,
}

/// The TSV summary's columns, in order (`SummaryRow`'s fields). Serde maps lose field order, so
/// the bindings use this to order the columns of the rows they return.
pub const SUMMARY_COLUMNS: [&str; 15] = [
    "query",
    "status",
    "taxid",
    "taxon",
    "rank",
    "top_genome",
    "top_seq",
    "pident",
    "qcov",
    "bits",
    "evalue",
    "tied_genomes",
    "tied_taxa",
    "searched",
    "flags",
];

/// One line of the TSV summary (all fields as printed).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SummaryRow {
    pub query: String,
    pub status: String,
    pub taxid: String,
    pub taxon: String,
    pub rank: String,
    pub top_genome: String,
    pub top_seq: String,
    pub pident: String,
    pub qcov: String,
    pub bits: String,
    pub evalue: String,
    pub tied_genomes: usize,
    pub tied_taxa: usize,
    pub searched: String,
    pub flags: String,
}

impl Report {
    /// The step-16 TSV rows.
    pub fn summary_rows(&self) -> Vec<SummaryRow> {
        self.queries
            .iter()
            .map(|q| {
                let dash = || "-".to_string();
                let top = q.top.as_ref();
                SummaryRow {
                    query: q.query.clone(),
                    status: q.status.clone(),
                    taxid: q
                        .call
                        .as_ref()
                        .map(|c| c.taxid.to_string())
                        .unwrap_or_else(dash),
                    taxon: q.call.as_ref().map(|c| c.name.clone()).unwrap_or_else(dash),
                    rank: q
                        .call
                        .as_ref()
                        .and_then(|c| c.rank.clone())
                        .unwrap_or_else(dash),
                    top_genome: top.map(|t| t.genome.clone()).unwrap_or_else(dash),
                    top_seq: top.map(|t| t.seq_name.clone()).unwrap_or_else(dash),
                    pident: top.map(|t| t.pident.clone()).unwrap_or_else(dash),
                    qcov: top.map(|t| format!("{:.1}", t.qcov)).unwrap_or_else(dash),
                    bits: top.map(|t| t.bits.clone()).unwrap_or_else(dash),
                    evalue: top.map(|t| t.evalue.clone()).unwrap_or_else(dash),
                    tied_genomes: q.tied.genomes,
                    tied_taxa: q.tied.taxa,
                    searched: format!("{}/{}", q.searched.genomes, q.searched.of),
                    flags: if q.flags.is_empty() {
                        "-".into()
                    } else {
                        q.flags.join(",")
                    },
                }
            })
            .collect()
    }
}

#[derive(Debug, Clone)]
struct Cand {
    coll: usize,
    name: String,
    alt: Vec<String>,
    digest: String,
    hit: Hit,
    stage: String,
    score: f64,
}

/// `a` before `b`: bits desc, evalue asc, pident desc, length desc, sstart asc.
fn order(a: &Hit, b: &Hit) -> Ordering {
    b.bits
        .total_cmp(&a.bits)
        .then(a.evalue.total_cmp(&b.evalue))
        .then(b.pident.total_cmp(&a.pident))
        .then(b.length.cmp(&a.length))
        .then(a.sstart.cmp(&b.sstart))
}

/// A genome's best hit, its group and whether it ties.
type Scored = (Cand, Group, bool);

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Group {
    Tax(u32),
    Unknown(usize),
}

fn round1(x: f64) -> f64 {
    (x * 10.0).round() / 10.0
}

/// Build the report: fan every hit out, keep the best per genome, call each query.
pub fn build_report(
    meta: &DbMeta,
    queries: Vec<QueryInput>,
    opts: &RankOptions,
    mut tax: Option<&mut Taxonomy>,
    info: ReportMeta,
) -> Result<Report, Error> {
    // ids that need ids.tsv, looked up in one pass
    let wanted: HashSet<String> = queries
        .iter()
        .flat_map(|q| q.runs.iter())
        .flat_map(|r| r.hits.iter())
        .filter_map(|h| meta.unresolved(h))
        .collect();
    let ids = meta.lookup_ids(&wanted)?;
    let group_of = |coll: usize, tax: &mut Option<&mut Taxonomy>| -> Group {
        match meta.collections[coll].taxid.filter(|t| *t > 0) {
            None => Group::Unknown(coll),
            Some(t) => {
                let t = t as u32;
                match tax {
                    Some(tx) => {
                        let t = tx.resolve(t);
                        match &opts.group_rank {
                            Some(r) => Group::Tax(tx.ancestor_at(t, r).unwrap_or(t)),
                            None => Group::Tax(t),
                        }
                    }
                    None => Group::Tax(t),
                }
            }
        }
    };
    struct Computed {
        q: QueryInput,
        qlen: Option<u64>,
        genomes: Vec<Scored>,
        status: String,
        call: Option<(Vec<Group>, Option<u32>)>,
        top: Option<usize>,
        tied_taxa: usize,
        truncated: bool,
    }
    let mut computed = Vec::new();
    let mut want_names: HashSet<u32> = HashSet::new();
    for q in queries {
        let mut best: HashMap<usize, Cand> = HashMap::new();
        let mut truncated = false;
        let mut qlen = q.qlen;
        for run in &q.runs {
            let mut subjects: HashSet<&str> = HashSet::new();
            for h in &run.hits {
                subjects.insert(h.subject.as_deref().unwrap_or(""));
                qlen = qlen.or(Some(h.qlen));
                if opts.min_pident.is_some_and(|m| h.pident < m)
                    || opts.min_qcov.is_some_and(|m| h.qcov() < m)
                    || opts.max_evalue.is_some_and(|m| h.evalue > m)
                {
                    continue;
                }
                let f = fanout(meta, h, run.scope.as_ref(), &ids)?;
                for p in f.pairs {
                    let c = Cand {
                        coll: p.coll,
                        name: p.name,
                        alt: p.alt_names,
                        digest: f.digest.clone(),
                        hit: h.clone(),
                        stage: run.stage.clone(),
                        score: 0.0,
                    };
                    match best.get(&p.coll) {
                        Some(b) if order(&b.hit, &c.hit) != Ordering::Greater => {}
                        _ => {
                            best.insert(p.coll, c);
                        }
                    }
                }
            }
            if q.max_target_seqs > 0 && subjects.len() >= q.max_target_seqs {
                truncated = true;
            }
        }
        let mut cands: Vec<Cand> = best.into_values().collect();
        let raw = !cands.is_empty() && cands.iter().all(|c| c.hit.score.is_some());
        for c in &mut cands {
            c.score = if raw {
                c.hit.score.unwrap()
            } else {
                c.hit.bits
            };
        }
        cands.sort_by(|a, b| {
            order(&a.hit, &b.hit).then_with(|| {
                meta.collections[a.coll]
                    .display
                    .cmp(&meta.collections[b.coll].display)
            })
        });
        let t = cands
            .iter()
            .map(|c| c.score)
            .fold(f64::NEG_INFINITY, f64::max);
        let thr = t * (1.0 - opts.tie_percent / 100.0);
        let genomes: Vec<Scored> = cands
            .into_iter()
            .map(|c| {
                let g = group_of(c.coll, &mut tax);
                let tied = c.score >= thr;
                (c, g, tied)
            })
            .collect();
        let mut tied_groups: BTreeMap<Group, usize> = BTreeMap::new();
        for (_, g, tied) in &genomes {
            if *tied {
                *tied_groups.entry(g.clone()).or_default() += 1;
            }
            if let Group::Tax(t) = g {
                want_names.insert(*t);
            }
        }
        let (status, call) = if genomes.is_empty() {
            ("none".to_string(), None)
        } else if tied_groups.len() == 1 {
            let g = tied_groups.keys().next().unwrap().clone();
            ("taxon".to_string(), Some((vec![g], None)))
        } else {
            let all_known = tied_groups.keys().all(|g| matches!(g, Group::Tax(_)));
            let mut gs: Vec<(Group, usize)> = tied_groups.clone().into_iter().collect();
            gs.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
            let groups: Vec<Group> = gs.into_iter().map(|(g, _)| g).collect();
            match (&tax, all_known) {
                (Some(tx), true) => {
                    let ts: Vec<u32> = groups
                        .iter()
                        .map(|g| match g {
                            Group::Tax(t) => *t,
                            Group::Unknown(_) => 1,
                        })
                        .collect();
                    let l = tx.lca(&ts);
                    want_names.insert(l);
                    ("lca".to_string(), Some((groups, Some(l))))
                }
                _ => ("multi".to_string(), Some((groups, None))),
            }
        };
        let top = genomes
            .iter()
            .enumerate()
            .filter(|(_, g)| g.2)
            .min_by(|(_, a), (_, b)| {
                b.0.hit
                    .pident
                    .total_cmp(&a.0.hit.pident)
                    .then(b.0.hit.length.cmp(&a.0.hit.length))
                    .then(
                        meta.collections[a.0.coll]
                            .display
                            .as_bytes()
                            .cmp(meta.collections[b.0.coll].display.as_bytes()),
                    )
            })
            .map(|(i, _)| i);
        computed.push(Computed {
            q,
            qlen,
            tied_taxa: tied_groups.len(),
            genomes,
            status,
            call,
            top,
            truncated,
        });
    }
    if let Some(tx) = tax.as_deref_mut() {
        tx.load_names(&want_names)?;
    }
    let tax = tax.as_deref();
    let species_of = |t: u32| -> Option<String> {
        meta.collections
            .iter()
            .find(|c| c.taxid == Some(t as i32))
            .and_then(|c| c.species.clone())
    };
    let name_of = |g: &Group| -> String {
        match g {
            Group::Unknown(c) => meta.collections[*c].display.clone(),
            Group::Tax(t) => tax
                .and_then(|tx| tx.name(*t).map(str::to_string))
                .or_else(|| species_of(*t))
                .unwrap_or_else(|| t.to_string()),
        }
    };
    let id_of = |g: &Group| -> TaxId {
        match g {
            Group::Unknown(c) => TaxId::Text(format!("unknown:{}", meta.collections[*c].display)),
            Group::Tax(t) => TaxId::Num(*t),
        }
    };
    let gh = |c: &Cand, tied: bool| -> GenomeHit {
        let cm = &meta.collections[c.coll];
        GenomeHit {
            genome: cm.display.clone(),
            collection: cm.collection.clone(),
            seq_name: c.name.clone(),
            alt_names: c.alt.clone(),
            seq_digest: format!("SQ.{}", c.digest),
            sstart: c.hit.sstart,
            send: c.hit.send,
            strand: if c.hit.sstart > c.hit.send {
                "minus"
            } else {
                "plus"
            }
            .into(),
            pident: c.hit.pident_s.clone(),
            length: c.hit.length,
            qcov: round1(c.hit.qcov()),
            evalue: c.hit.evalue_s.clone(),
            bits: c.hit.bits_s.clone(),
            tied,
            stage: c.stage.clone(),
        }
    };
    let mut out = Vec::new();
    for c in computed {
        let call = c
            .call
            .as_ref()
            .map(|(groups, lca)| match (c.status.as_str(), lca) {
                ("lca", Some(l)) => Call {
                    taxid: TaxId::Num(*l),
                    name: tax
                        .and_then(|tx| tx.name(*l).map(str::to_string))
                        .unwrap_or_else(|| l.to_string()),
                    rank: tax.map(|tx| tx.rank(*l).to_string()),
                },
                ("multi", _) => Call {
                    taxid: TaxId::Text(
                        groups
                            .iter()
                            .map(|g| id_of(g).to_string())
                            .collect::<Vec<_>>()
                            .join(","),
                    ),
                    name: groups.iter().map(&name_of).collect::<Vec<_>>().join(","),
                    rank: None,
                },
                _ => {
                    let g = &groups[0];
                    Call {
                        taxid: id_of(g),
                        name: name_of(g),
                        rank: match (g, tax) {
                            (Group::Tax(t), Some(tx)) => Some(tx.rank(*t).to_string()),
                            _ => None,
                        },
                    }
                }
            });
        // per-group detail, groups by best score then name, genomes by score then display
        let mut groups: Vec<(Group, Vec<&Scored>)> = Vec::new();
        for g in &c.genomes {
            match groups.iter_mut().find(|(k, _)| *k == g.1) {
                Some((_, v)) => v.push(g),
                None => groups.push((g.1.clone(), vec![g])),
            }
        }
        for (_, v) in groups.iter_mut() {
            v.sort_by(|a, b| {
                b.0.score.total_cmp(&a.0.score).then_with(|| {
                    meta.collections[a.0.coll]
                        .display
                        .cmp(&meta.collections[b.0.coll].display)
                })
            });
        }
        groups.sort_by(|a, b| {
            b.1[0]
                .0
                .score
                .total_cmp(&a.1[0].0.score)
                .then_with(|| name_of(&a.0).cmp(&name_of(&b.0)))
        });
        let taxa = groups
            .iter()
            .map(|(g, v)| TaxonGroup {
                taxid: id_of(g),
                name: name_of(g),
                best_bits: v[0].0.hit.bits_s.clone(),
                best_pident: v[0].0.hit.pident_s.clone(),
                tied: v.iter().filter(|x| x.2).count(),
                hit: v.len(),
                genomes: v.iter().map(|x| gh(&x.0, x.2)).collect(),
            })
            .collect();
        let mut flags = c.q.flags.clone();
        if c.truncated {
            flags.insert(0, "truncated".into());
        }
        out.push(QueryResult {
            query: c.q.name.clone(),
            qlen: c.qlen,
            status: c.status.clone(),
            call,
            top: c.top.map(|i| gh(&c.genomes[i].0, true)),
            tied: Tied {
                genomes: c.genomes.iter().filter(|g| g.2).count(),
                taxa: c.tied_taxa,
            },
            searched: Searched {
                genomes: c.q.searched.0,
                of: c.q.searched.1,
            },
            flags,
            taxa,
        });
    }
    let mut info = info;
    if let Some(tx) = tax {
        info.warnings.extend(tx.warnings.iter().cloned());
    }
    Ok(Report {
        rgblast_report: 1,
        meta: info,
        queries: out,
    })
}
