//! From a hit's `sallseqid` (or its subject digest) to the (collection, name) pairs it stands for.

use super::parse::Hit;
use crate::blastdb::sidecars::{self, CollectionMeta};
use crate::error::Error;
use crate::plan::SubjectPlan;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Member {
    pub coll: usize,
    pub digest: String,
    pub name: String,
}

/// The collections and members of a database (its sidecars) or of a store selection (its plan).
#[derive(Debug, Clone, Default)]
pub struct DbMeta {
    pub collections: Vec<CollectionMeta>,
    pub members: Vec<Member>,
    by_digest: HashMap<String, usize>,
    by_tag: HashMap<String, usize>,
    by_idname: HashMap<(usize, String), usize>,
    by_seq: HashMap<String, Vec<usize>>,
    /// `<db>.ids.tsv`, the last resort for ids that are neither `SQ.` nor `gnl|<db_tag>|`.
    pub ids_path: Option<PathBuf>,
}

impl DbMeta {
    fn index(collections: Vec<CollectionMeta>, rows: Vec<(usize, String, String, String)>) -> Self {
        let mut m = DbMeta {
            by_digest: collections
                .iter()
                .enumerate()
                .map(|(i, c)| (c.collection.clone(), i))
                .collect(),
            by_tag: collections
                .iter()
                .enumerate()
                .map(|(i, c)| (c.db_tag.clone(), i))
                .collect(),
            collections,
            ..Default::default()
        };
        for (coll, digest, name, id_name) in rows {
            let i = m.members.len();
            m.by_idname.entry((coll, id_name)).or_insert(i);
            m.by_seq.entry(digest.clone()).or_default().push(i);
            m.members.push(Member { coll, digest, name });
        }
        m
    }

    /// From `<base>.collections.tsv` and `<base>.members.tsv`.
    pub fn from_sidecars(base: &Path) -> Result<Self, Error> {
        let [ids, colls, members, ..] = sidecars::paths(base);
        let collections = sidecars::read_collections(&colls)?;
        let by: HashMap<String, usize> = collections
            .iter()
            .enumerate()
            .map(|(i, c)| (c.collection.clone(), i))
            .collect();
        let mut rows = Vec::new();
        for r in sidecars::read_members(&members)? {
            let Some(&c) = by.get(&r.collection) else {
                return Err(Error::InvalidOptions(format!(
                    "{}: collection {} is not in collections.tsv",
                    members.display(),
                    r.collection
                )));
            };
            rows.push((c, r.digest, r.name, r.id_name));
        }
        let mut m = Self::index(collections, rows);
        m.ids_path = ids.exists().then_some(ids);
        Ok(m)
    }

    /// From a selection's plan (the same tables, built from the store directly).
    pub fn from_plan(plan: &SubjectPlan) -> Self {
        let collections = plan
            .collections
            .iter()
            .map(|c| CollectionMeta {
                collection: c.digest.clone(),
                db_tag: c.db_tag.clone(),
                display: c.display.clone(),
                aliases: c.aliases.clone(),
                taxid: c.taxid,
                species: c.species.clone(),
                assembly: c.assembly.clone(),
                n_sequences: c.n_sequences,
                total_length: c.total_length,
                n_fraction: c.n_fraction,
            })
            .collect();
        let mut rows = Vec::new();
        for s in &plan.subjects {
            for o in &s.occurrences {
                rows.push((
                    o.collection as usize,
                    s.digest.clone(),
                    o.name.clone(),
                    o.id_name.clone(),
                ));
            }
        }
        Self::index(collections, rows)
    }

    pub fn coll_index(&self, digest: &str) -> Option<usize> {
        self.by_digest.get(digest).copied()
    }

    fn gnl(&self, id: &str) -> Option<usize> {
        let rest = id.strip_prefix("gnl|")?;
        let (tag, name) = rest.split_once('|')?;
        let c = *self.by_tag.get(tag)?;
        self.by_idname.get(&(c, name.to_string())).copied()
    }

    /// Digests for ids that are neither `SQ.` nor resolvable `gnl|` ids, from `ids.tsv` in one
    /// streaming pass that keeps only `wanted`.
    pub fn lookup_ids(&self, wanted: &HashSet<String>) -> Result<HashMap<String, String>, Error> {
        let mut out = HashMap::new();
        let Some(p) = &self.ids_path else {
            return Ok(out);
        };
        if wanted.is_empty() {
            return Ok(out);
        }
        let f = std::fs::File::open(p).map_err(|e| Error::io(p, e))?;
        for line in BufReader::new(f).lines() {
            let line = line.map_err(|e| Error::io(p, e))?;
            if line.starts_with('#') {
                continue;
            }
            let mut it = line.split('\t');
            if let (Some(id), Some(_), Some(d)) = (it.next(), it.next(), it.next())
                && wanted.contains(id)
            {
                out.insert(id.to_string(), d.to_string());
            }
        }
        Ok(out)
    }

    /// The first id of a hit that does not resolve without `ids.tsv`, if any.
    pub fn unresolved(&self, hit: &Hit) -> Option<String> {
        if hit.sallseqid.iter().any(|i| i.starts_with("SQ.")) {
            return None;
        }
        let first = hit.sallseqid.first()?;
        if self.gnl(first).is_some() {
            return None;
        }
        Some(first.clone())
    }
}

/// One (collection, name) a hit stands for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pair {
    pub coll: usize,
    pub name: String,
    /// Further names of the same sequence in the same collection.
    pub alt_names: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fanned {
    pub digest: String,
    pub pairs: Vec<Pair>,
}

fn add(pairs: &mut Vec<Pair>, coll: usize, name: &str) {
    if let Some(p) = pairs.iter_mut().find(|p| p.coll == coll) {
        if p.name != name && !p.alt_names.iter().any(|a| a == name) {
            p.alt_names.push(name.to_string());
        }
    } else {
        pairs.push(Pair {
            coll,
            name: name.to_string(),
            alt_names: Vec::new(),
        });
    }
}

/// Fan a hit out to the pairs within `scope` (collection indices; None = all). With `gnl|`
/// entries, they are exactly the pairs BLAST kept (after any view filter). With only
/// `SQ.<digest>` (an `SQ.` seqid list), every member of the digest.
pub fn fanout(
    meta: &DbMeta,
    hit: &Hit,
    scope: Option<&HashSet<usize>>,
    ids: &HashMap<String, String>,
) -> Result<Fanned, Error> {
    let unresolved = |id: &str| {
        Error::InvalidOptions(format!(
            "hits: subject id '{id}' is not in the database's sidecars (wrong --db or --store?)"
        ))
    };
    let gnls: Vec<&String> = hit
        .sallseqid
        .iter()
        .filter(|i| i.starts_with("gnl|"))
        .collect();
    let digest = if let Some(sq) = hit.sallseqid.iter().find_map(|i| i.strip_prefix("SQ.")) {
        sq.to_string()
    } else if let Some(g) = gnls.first() {
        match meta.gnl(g) {
            Some(i) => meta.members[i].digest.clone(),
            None => ids.get(*g).cloned().ok_or_else(|| unresolved(g))?,
        }
    } else {
        let first = hit.sallseqid.first().map(|s| s.as_str()).unwrap_or("");
        ids.get(first).cloned().ok_or_else(|| unresolved(first))?
    };
    let mut pairs = Vec::new();
    if !gnls.is_empty() {
        for g in gnls {
            let i = meta.gnl(g).ok_or_else(|| unresolved(g))?;
            let m = &meta.members[i];
            if scope.is_none_or(|s| s.contains(&m.coll)) {
                add(&mut pairs, m.coll, &m.name);
            }
        }
    } else {
        for &i in meta.by_seq.get(&digest).into_iter().flatten() {
            let m = &meta.members[i];
            if scope.is_none_or(|s| s.contains(&m.coll)) {
                add(&mut pairs, m.coll, &m.name);
            }
        }
    }
    Ok(Fanned { digest, pairs })
}
