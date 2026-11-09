//! The representative rule: one assembly per group (taxid, or `--group-rank` ancestor), chosen
//! by the first rule that applies:
//! 1. `override`: `--rep COLLECTION`;
//! 2. `namespace`: the only candidate with an alias in `--rep-namespace` (default `reference`);
//!    several such candidates go to rule 3;
//! 3. `most-shared`: the largest sum over its sequences of length x (candidates holding it - 1),
//!    ties to the larger total length, then the smaller digest;
//! 4. `largest`: rule 3 when every score is 0; `only` for a group of one.
//!
//! Candidates are collections with N fraction below 0.25 (unknown counts as below); when no
//! collection of a group qualifies, all are candidates.

use crate::error::Error;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};

pub const MAX_CANDIDATE_N: f64 = 0.25;
/// Staged at genome level when reps + the largest remainder is at most this share of a full search.
pub const GENOME_THRESHOLD: f64 = 0.8;

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum GroupKey {
    Taxid(u32),
    /// A collection with no taxid (its digest): its own group.
    Unknown(String),
}

impl std::fmt::Display for GroupKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GroupKey::Taxid(t) => write!(f, "{t}"),
            GroupKey::Unknown(d) => write!(f, "unknown:{d}"),
        }
    }
}

/// One selected collection, with sequences as (index, length); indices stand for digests.
#[derive(Debug, Clone, PartialEq)]
pub struct RepColl {
    pub digest: String,
    pub display: String,
    pub group: GroupKey,
    pub group_name: String,
    pub n_fraction: Option<f64>,
    /// Has an alias in the representative namespace.
    pub in_namespace: bool,
    pub seqs: Vec<(u32, u64)>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RepsRow {
    pub group: String,
    pub taxon: String,
    pub collections: usize,
    pub rep: String,
    pub rep_collection: String,
    pub rule: String,
    pub rep_bases: u64,
    /// Unique bases of the group.
    pub group_bases: u64,
    pub not_in_rep_seqs: u64,
    pub not_in_rep_bases: u64,
    /// Indices (into the input) of the group's collections, and of the representative.
    #[serde(skip)]
    pub members: Vec<usize>,
    #[serde(skip)]
    pub rep_index: usize,
    #[serde(skip)]
    pub key: Option<GroupKey>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RepsSummary {
    pub representatives: usize,
    pub collections: usize,
    pub groups: usize,
    /// Unique bases of all representatives.
    pub reps_bases: u64,
    /// The full selection's total length.
    pub full_bases: u64,
    /// reps + the largest group remainder (group minus every representative sequence).
    pub worst_genome_bases: u64,
    pub staged_genome: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RepsTable {
    pub groups: Vec<RepsRow>,
    pub summary: RepsSummary,
}

fn gbp(b: u64) -> String {
    format!("{:.2}", b as f64 / 1e9)
}

impl RepsTable {
    pub const HEADER: &'static str = "#group\ttaxon\tcollections\trep\trep_collection\trule\trep_gbp\tgroup_gbp\tnot_in_rep_seqs\tnot_in_rep_gbp";

    /// The `rgblast reps` table and its two summary lines.
    pub fn render(&self) -> String {
        let mut s = format!("{}\n", Self::HEADER);
        for r in &self.groups {
            s += &format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                r.group,
                r.taxon,
                r.collections,
                r.rep,
                r.rep_collection,
                r.rule,
                gbp(r.rep_bases),
                gbp(r.group_bases),
                r.not_in_rep_seqs,
                gbp(r.not_in_rep_bases)
            );
        }
        let m = &self.summary;
        let full = m.full_bases.max(1) as f64;
        s += &format!(
            "# representatives: {} of {} collections, {} of {} Gbp ({:.1}x smaller)\n",
            m.representatives,
            m.collections,
            gbp(m.reps_bases),
            gbp(m.full_bases),
            full / m.reps_bases.max(1) as f64
        );
        s += &format!(
            "# worst-case staged query (genome level): {} Gbp = {:.2} of a full search -> staged mode {}\n",
            gbp(m.worst_genome_bases),
            m.worst_genome_bases as f64 / full,
            if m.staged_genome { "on" } else { "off" }
        );
        s
    }

    /// Sequence indices of every representative.
    pub fn rep_seqs(&self, colls: &[RepColl]) -> HashSet<u32> {
        self.groups
            .iter()
            .flat_map(|g| colls[g.rep_index].seqs.iter().map(|s| s.0))
            .collect()
    }
}

fn unique_bases(seqs: impl Iterator<Item = (u32, u64)>) -> (HashMap<u32, u64>, u64) {
    let m: HashMap<u32, u64> = seqs.collect();
    let total = m.values().sum();
    (m, total)
}

/// Choose the representatives. `overrides` are collection digests (`--rep`), `full_bases` the
/// full selection's total length (0 = the union of the collections' sequences).
pub fn choose(
    colls: &[RepColl],
    overrides: &[String],
    full_bases: u64,
) -> Result<RepsTable, Error> {
    let mut groups: BTreeMap<GroupKey, Vec<usize>> = BTreeMap::new();
    for (i, c) in colls.iter().enumerate() {
        groups.entry(c.group.clone()).or_default().push(i);
    }
    for o in overrides {
        if !colls.iter().any(|c| &c.digest == o) {
            return Err(Error::InvalidOptions(format!(
                "--rep {o} is not in the selection"
            )));
        }
    }
    let mut rows = Vec::new();
    for (key, members) in &groups {
        let picked: Vec<usize> = members
            .iter()
            .copied()
            .filter(|&i| overrides.contains(&colls[i].digest))
            .collect();
        if picked.len() > 1 {
            return Err(Error::InvalidOptions(format!(
                "--rep names two collections of group {} ({}): {} and {}",
                key,
                colls[picked[0]].group_name,
                colls[picked[0]].display,
                colls[picked[1]].display
            )));
        }
        let (rep, rule) = if let Some(&p) = picked.first() {
            (p, "override")
        } else if members.len() == 1 {
            (members[0], "only")
        } else {
            let mut cands: Vec<usize> = members
                .iter()
                .copied()
                .filter(|&i| colls[i].n_fraction.is_none_or(|f| f < MAX_CANDIDATE_N))
                .collect();
            if cands.is_empty() {
                cands = members.clone();
            }
            let ns: Vec<usize> = cands
                .iter()
                .copied()
                .filter(|&i| colls[i].in_namespace)
                .collect();
            if ns.len() == 1 {
                (ns[0], "namespace")
            } else {
                if ns.len() > 1 {
                    cands = ns;
                }
                most_shared(colls, &cands)
            }
        };
        let (group_seqs, group_bases) =
            unique_bases(members.iter().flat_map(|&i| colls[i].seqs.iter().copied()));
        let (rep_seqs, rep_bases) = unique_bases(colls[rep].seqs.iter().copied());
        let missing: Vec<u64> = group_seqs
            .iter()
            .filter(|(s, _)| !rep_seqs.contains_key(s))
            .map(|(_, l)| *l)
            .collect();
        rows.push(RepsRow {
            group: key.to_string(),
            taxon: colls[rep].group_name.clone(),
            collections: members.len(),
            rep: colls[rep].display.clone(),
            rep_collection: colls[rep].digest.clone(),
            rule: rule.into(),
            rep_bases,
            group_bases,
            not_in_rep_seqs: missing.len() as u64,
            not_in_rep_bases: missing.iter().sum(),
            members: members.clone(),
            rep_index: rep,
            key: Some(key.clone()),
        });
    }
    let (_, reps_bases) = unique_bases(
        rows.iter()
            .flat_map(|r| colls[r.rep_index].seqs.iter().copied()),
    );
    let all_reps: HashSet<u32> = rows
        .iter()
        .flat_map(|r| colls[r.rep_index].seqs.iter().map(|s| s.0))
        .collect();
    let (_, union) = unique_bases(colls.iter().flat_map(|c| c.seqs.iter().copied()));
    let full = if full_bases > 0 { full_bases } else { union };
    let max_rem = rows
        .iter()
        .map(|r| {
            let (m, _) = unique_bases(
                r.members
                    .iter()
                    .flat_map(|&i| colls[i].seqs.iter().copied()),
            );
            m.iter()
                .filter(|(s, _)| !all_reps.contains(s))
                .map(|(_, l)| *l)
                .sum::<u64>()
        })
        .max()
        .unwrap_or(0);
    let worst = reps_bases + max_rem;
    let summary = RepsSummary {
        representatives: rows.len(),
        collections: colls.len(),
        groups: rows.len(),
        reps_bases,
        full_bases: full,
        worst_genome_bases: worst,
        staged_genome: (worst as f64) <= GENOME_THRESHOLD * full as f64,
    };
    Ok(RepsTable {
        groups: rows,
        summary,
    })
}

/// Rule 3 (and 4) among `cands`.
fn most_shared(colls: &[RepColl], cands: &[usize]) -> (usize, &'static str) {
    let mut share: HashMap<u32, u64> = HashMap::new();
    let sets: Vec<HashMap<u32, u64>> = cands
        .iter()
        .map(|&i| colls[i].seqs.iter().copied().collect())
        .collect();
    for s in &sets {
        for k in s.keys() {
            *share.entry(*k).or_default() += 1;
        }
    }
    let scored: Vec<(usize, u64, u64)> = cands
        .iter()
        .zip(&sets)
        .map(|(&i, s)| {
            let score = s.iter().map(|(k, l)| l * (share[k] - 1)).sum();
            let total = s.values().sum();
            (i, score, total)
        })
        .collect();
    let best = scored
        .iter()
        .min_by(|a, b| {
            b.1.cmp(&a.1).then(b.2.cmp(&a.2)).then(
                colls[a.0]
                    .digest
                    .as_bytes()
                    .cmp(colls[b.0].digest.as_bytes()),
            )
        })
        .unwrap();
    let rule = if scored.iter().all(|s| s.1 == 0) {
        "largest"
    } else {
        "most-shared"
    };
    (best.0, rule)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(d: &str, g: u32, seqs: &[(u32, u64)]) -> RepColl {
        RepColl {
            digest: d.into(),
            display: format!("asm-{d}"),
            group: GroupKey::Taxid(g),
            group_name: format!("taxon{g}"),
            n_fraction: Some(0.01),
            in_namespace: false,
            seqs: seqs.to_vec(),
        }
    }

    fn rule_of(t: &RepsTable, g: &str) -> (String, String) {
        let r = t.groups.iter().find(|r| r.group == g).unwrap();
        (r.rep_collection.clone(), r.rule.clone())
    }

    #[test]
    fn rules_and_tie_breaks() {
        // most-shared: b shares 1 and 2 with a and c
        let colls = vec![
            c("a", 1, &[(1, 100), (9, 50)]),
            c("b", 1, &[(1, 100), (2, 80)]),
            c("c", 1, &[(2, 80), (3, 500)]),
        ];
        let t = choose(&colls, &[], 0).unwrap();
        assert_eq!(rule_of(&t, "1"), ("b".into(), "most-shared".into()));
        // override
        let t = choose(&colls, &["c".into()], 0).unwrap();
        assert_eq!(rule_of(&t, "1"), ("c".into(), "override".into()));
        assert!(choose(&colls, &["a".into(), "c".into()], 0).is_err());
        assert!(choose(&colls, &["zz".into()], 0).is_err());
        // namespace: one candidate wins outright; several go to most-shared
        let mut ns = colls.clone();
        ns[0].in_namespace = true;
        assert_eq!(
            rule_of(&choose(&ns, &[], 0).unwrap(), "1"),
            ("a".into(), "namespace".into())
        );
        ns[2].in_namespace = true;
        // a and c: no shared sequence between them -> largest (c)
        assert_eq!(
            rule_of(&choose(&ns, &[], 0).unwrap(), "1"),
            ("c".into(), "largest".into())
        );
        // largest: no sharing; ties on size go to the smaller digest
        let colls = vec![c("y", 2, &[(5, 10)]), c("x", 2, &[(6, 10)])];
        assert_eq!(
            rule_of(&choose(&colls, &[], 0).unwrap(), "2"),
            ("x".into(), "largest".into())
        );
        // equal score, size decides
        let colls = vec![
            c("p", 3, &[(1, 10), (2, 10)]),
            c("q", 3, &[(1, 10), (2, 10), (7, 5)]),
        ];
        assert_eq!(
            rule_of(&choose(&colls, &[], 0).unwrap(), "3"),
            ("q".into(), "most-shared".into())
        );
        // only; and a collection with no taxid is its own group
        let mut u = c("u", 0, &[(8, 10)]);
        u.group = GroupKey::Unknown("u".into());
        let t = choose(&[c("o", 4, &[(1, 1)]), u], &[], 0).unwrap();
        assert_eq!(rule_of(&t, "4"), ("o".into(), "only".into()));
        assert_eq!(rule_of(&t, "unknown:u"), ("u".into(), "only".into()));
    }

    #[test]
    fn n_fraction_candidates() {
        let mut colls = vec![
            c("m", 1, &[(1, 100), (2, 100), (3, 100)]),
            c("n", 1, &[(1, 100), (2, 100)]),
            c("k", 1, &[(1, 100)]),
        ];
        colls[0].n_fraction = Some(0.5); // hard-masked: not a candidate
        let t = choose(&colls, &[], 0).unwrap();
        assert_eq!(rule_of(&t, "1").0, "n");
        // all masked: everyone is a candidate again
        for x in &mut colls {
            x.n_fraction = Some(0.5);
        }
        assert_eq!(rule_of(&choose(&colls, &[], 0).unwrap(), "1").0, "m");
    }

    #[test]
    fn vgp_shaped_table_is_exhaustive() {
        // 585 groups for 605 collections, no shared sequences
        let mut colls = Vec::new();
        let mut s = 0u32;
        for g in 0..585u32 {
            let n = if g < 20 { 2 } else { 1 };
            for k in 0..n {
                s += 1;
                colls.push(c(&format!("{g}-{k}"), g + 1, &[(s, 1_000_000 + s as u64)]));
            }
        }
        let t = choose(&colls, &[], 0).unwrap();
        assert_eq!(t.summary.representatives, 585);
        assert!(!t.summary.staged_genome);
    }
}
