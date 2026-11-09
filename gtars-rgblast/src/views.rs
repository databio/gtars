//! `rgmakeblastdb view`: subset views of a built database, with no rebuild. Each view is a
//! `.bsl` (the member rows' `<collection digest>:<id_name>` ids, never `SQ.`, so hits show the
//! collection's own ids) plus a `.nal` naming the database, with the exact `NSEQ` and `LENGTH`
//! of the subset, so `blastdbcmd -info`, the blastn footer and every E-value equal a fresh
//! build of the same selection.
//!
//! [`DbCatalog`] serves selections from the database's own sidecars and manifest: no store and
//! no network are needed.

use crate::blastdb::alias::{self, NalSpec};
use crate::blastdb::output::{db_dir, with_ext};
pub use crate::blastdb::sidecars::ViewKind;
use crate::blastdb::sidecars::{self, CollectionMeta, Manifest, MemberRow, OidRow};
use crate::error::Error;
use crate::plan::CollectionInfo;
use crate::progress::ProgressEvent;
use crate::select::resolve::{CatalogRow, Names};
use crate::select::{Catalog, ResolveOptions, Selection};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Component, Path, PathBuf};

/// [`Catalog`] over a built database's `collections.tsv`, `members.tsv` and `oids.tsv`.
pub struct DbCatalog {
    pub base: PathBuf,
    pub manifest: Manifest,
    pub collections: Vec<CollectionMeta>,
    pub members: Vec<MemberRow>,
    pub oids: Vec<OidRow>,
    by_digest: HashMap<String, usize>,
    members_of: HashMap<String, Vec<usize>>,
    oid_of: HashMap<String, u32>,
    namespaces: Vec<String>,
}

impl DbCatalog {
    /// Open the sidecars of the database at `base`; one without a manifest is
    /// `NotFound`.
    pub fn open(base: &Path) -> Result<Self, Error> {
        let manifest = sidecars::read_manifest(base)?;
        let collections = sidecars::read_collections(&with_ext(base, "collections.tsv"))?;
        let members_tsv = with_ext(base, "members.tsv");
        let members = sidecars::read_members(&members_tsv)?;
        let oids = sidecars::read_oids(&with_ext(base, "oids.tsv"))?;
        let by_digest = collections
            .iter()
            .enumerate()
            .map(|(i, c)| (c.collection.clone(), i))
            .collect();
        let mut members_of: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, m) in members.iter().enumerate() {
            members_of.entry(m.collection.clone()).or_default().push(i);
        }
        let oid_of = oids.iter().map(|r| (r.digest.clone(), r.oid)).collect();
        let mut namespaces: Vec<String> = Vec::new();
        for c in &collections {
            for (ns, _) in &c.aliases {
                if !namespaces.contains(ns) {
                    namespaces.push(ns.clone());
                }
            }
        }
        Ok(DbCatalog {
            base: base.to_path_buf(),
            manifest,
            collections,
            members,
            oids,
            by_digest,
            members_of,
            oid_of,
            namespaces,
        })
    }

    pub fn name(&self) -> String {
        crate::blastdb::writer::file_name(&self.base)
    }

    pub fn meta(&self, digest: &str) -> Option<&CollectionMeta> {
        self.by_digest.get(digest).map(|&i| &self.collections[i])
    }

    /// Member rows of a collection, in OID order.
    pub fn members_of(&self, digest: &str) -> Vec<&MemberRow> {
        self.members_of
            .get(digest)
            .map(|v| v.iter().map(|&i| &self.members[i]).collect())
            .unwrap_or_default()
    }
}

impl Names for DbCatalog {
    fn label(&self) -> String {
        format!("database {}", self.name())
    }

    fn has_collection(&self, digest: &str) -> bool {
        self.by_digest.contains_key(digest)
    }

    fn namespaces(&self) -> Vec<String> {
        self.namespaces.clone()
    }

    fn lookup_alias(&self, ns: &str, alias: &str) -> Option<String> {
        self.collections
            .iter()
            .find(|c| c.aliases.iter().any(|(n, a)| n == ns && a == alias))
            .map(|c| c.collection.clone())
    }
}

impl Catalog for DbCatalog {
    fn collections(&self) -> Vec<String> {
        self.collections
            .iter()
            .map(|c| c.collection.clone())
            .collect()
    }

    fn aliases(&self, digest: &str) -> Vec<(String, String)> {
        self.meta(digest)
            .map(|c| c.aliases.clone())
            .unwrap_or_default()
    }

    fn info(&mut self, digest: &str) -> Result<(CollectionInfo, Option<String>), Error> {
        let c = self.meta(digest);
        Ok((
            CollectionInfo {
                taxid: c.and_then(|c| c.taxid),
                species: c.and_then(|c| c.species.clone()),
                assembly: c.and_then(|c| c.assembly.clone()),
                n_fraction: c.and_then(|c| c.n_fraction),
            },
            None,
        ))
    }

    fn rows(&mut self, digest: &str) -> Result<Vec<CatalogRow>, Error> {
        Ok(self
            .members_of(digest)
            .into_iter()
            .map(|m| CatalogRow {
                name: m.name.clone(),
                length: self.oids.get(m.oid as usize).map(|r| r.length).unwrap_or(0),
                digest: m.digest.clone(),
                alphabet: None,
            })
            .collect())
    }

    fn n_fractions(
        &mut self,
        digests: &[String],
        _scan: bool,
        _progress: &mut dyn FnMut(&ProgressEvent) -> bool,
    ) -> Result<HashMap<String, Option<f64>>, Error> {
        Ok(digests
            .iter()
            .map(|d| (d.clone(), self.meta(d).and_then(|c| c.n_fraction)))
            .collect())
    }

    fn has_sequence(&mut self, digest: &str) -> Result<bool, Error> {
        Ok(self.oid_of.contains_key(digest))
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ViewRequest {
    /// Database base path.
    pub db: PathBuf,
    pub kind: ViewKind,
    /// The view name (named views only).
    pub name: Option<String>,
    /// Which collections (required for named views; limits the others).
    pub selection: Option<Selection>,
    /// Default `<base>.views/`.
    pub out_dir: Option<PathBuf>,
    /// Title prefix (per-collection) or title (named); default the database title.
    pub title: Option<String>,
    /// Also write `STATS_NSEQ`/`STATS_TOTLEN` of the whole database, for E-values on its scale.
    pub full_db_stats: bool,
}

impl ViewRequest {
    pub fn validate(&self) -> Result<(), Error> {
        let bad = |m: &str| Err(Error::InvalidOptions(m.to_string()));
        if self.db.as_os_str().is_empty() {
            return bad("--db BASE is required");
        }
        match (self.kind, &self.name, &self.selection) {
            (ViewKind::Named, None, _) => return bad("a named view needs --name NAME"),
            (ViewKind::Named, Some(_), None) => {
                return bad(
                    "--name requires a selection (--collection, --collections-file, --namespace, --taxon or --all)",
                );
            }
            (ViewKind::PerCollection, Some(_), _) => {
                return bad("exactly one of --per-collection and --name is required");
            }
            _ => {}
        }
        if let Some(n) = &self.name
            && safe_name(n) != *n
        {
            return Err(Error::InvalidOptions(format!(
                "view name '{n}' may hold only letters, digits, '.', '_' and '-', and may not start with '.' or '-'"
            )));
        }
        if let Some(s) = &self.selection {
            s.validate()?;
        }
        Ok(())
    }
}

/// One view written.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ViewInfo {
    pub name: String,
    pub kind: ViewKind,
    pub collections: Vec<String>,
    pub taxid: Option<i32>,
    pub nseq: u64,
    pub length: u64,
    pub nal: PathBuf,
    pub bsl: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViewSummary {
    pub out_dir: PathBuf,
    pub index: PathBuf,
    pub views: Vec<ViewInfo>,
}

/// A file name from a tag: `[A-Za-z0-9._-]` kept, anything else `_`, a leading `-` or `.`
/// becomes `_`.
pub fn safe_name(s: &str) -> String {
    let mut v: String = s
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-') {
                c
            } else {
                '_'
            }
        })
        .collect();
    if v.starts_with('-') || v.starts_with('.') {
        v.replace_range(0..1, "_");
    }
    if v.is_empty() {
        v.push('_');
    }
    v
}

/// `to` relative to directory `from` (both made absolute first).
fn relative(from: &Path, to: &Path) -> PathBuf {
    let abs = |p: &Path| {
        std::fs::canonicalize(p).unwrap_or_else(|_| {
            std::env::current_dir()
                .map(|c| c.join(p))
                .unwrap_or_else(|_| p.to_path_buf())
        })
    };
    let from = abs(from);
    let to_dir = abs(&db_dir(to));
    let to = to_dir.join(to.file_name().unwrap_or_default());
    let f: Vec<Component> = from.components().collect();
    let t: Vec<Component> = to.components().collect();
    let common = f.iter().zip(&t).take_while(|(a, b)| a == b).count();
    let mut r = PathBuf::new();
    for _ in common..f.len() {
        r.push("..");
    }
    for c in &t[common..] {
        r.push(c.as_os_str());
    }
    r
}

struct Planned {
    name: String,
    kind: ViewKind,
    title: String,
    collections: Vec<String>,
    taxid: Option<i32>,
}

/// Write the views of a request (`<base>.views/` unless `out_dir` is given).
pub fn view(
    req: &ViewRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<ViewSummary, Error> {
    req.validate()?;
    let base = req.db.as_path();
    let out_dir = req.out_dir.clone().unwrap_or_else(|| {
        db_dir(base).join(format!(
            "{}.{}",
            crate::blastdb::writer::file_name(base),
            crate::blastdb::output::VIEWS_DIR
        ))
    });
    let mut cat = DbCatalog::open(base)?;
    let dbname = cat.name();
    // which collections (and which rows)
    let (chosen, sequences): (Vec<String>, Option<HashSet<String>>) = match &req.selection {
        None => (cat.collections(), None),
        Some(sel) => {
            let r = crate::select::resolve(&mut cat, sel, &ResolveOptions::default(), progress)
                .map_err(|e| match e {
                    Error::UnknownCollection { key, .. } => Error::InvalidOptions(format!(
                        "collection {key} is not in database {dbname}; to add it, rebuild the database with rgmakeblastdb build --force"
                    )),
                    e => e,
                })?;
            (r.collections, r.sequences.map(|v| v.into_iter().collect()))
        }
    };
    let prefix = match (&req.kind, &req.title) {
        (ViewKind::Named, _) | (_, None) => cat.manifest.title.clone(),
        (_, Some(t)) => t.clone(),
    };
    let mut chosen = chosen;
    chosen.sort_by(|a, b| {
        let (ma, mb) = (cat.meta(a), cat.meta(b));
        let da = ma.map(|m| m.display.clone()).unwrap_or_default();
        let db = mb.map(|m| m.display.clone()).unwrap_or_default();
        (da, a).cmp(&(db, b))
    });
    let mut planned: Vec<Planned> = Vec::new();
    match req.kind {
        ViewKind::PerCollection => {
            let mut used: HashSet<String> = HashSet::new();
            for d in &chosen {
                let m = cat.meta(d).unwrap();
                let mut name = safe_name(&m.db_tag);
                if !used.insert(name.to_lowercase()) {
                    name = format!(
                        "{name}.{}",
                        safe_name(&d.chars().take(8).collect::<String>())
                    );
                    used.insert(name.to_lowercase());
                }
                let sp = m
                    .species
                    .as_ref()
                    .map(|s| format!("{s} "))
                    .unwrap_or_default();
                planned.push(Planned {
                    name,
                    kind: req.kind,
                    title: format!("{prefix}: {sp}{}", m.display),
                    collections: vec![d.clone()],
                    taxid: m.taxid,
                });
            }
        }
        ViewKind::Named => {
            let taxa: BTreeSet<Option<i32>> = chosen
                .iter()
                .map(|d| cat.meta(d).and_then(|m| m.taxid))
                .collect();
            let name = req.name.clone().unwrap_or_default();
            planned.push(Planned {
                kind: req.kind,
                title: req
                    .title
                    .clone()
                    .unwrap_or_else(|| format!("{prefix}: {name}")),
                name,
                collections: chosen.clone(),
                taxid: if taxa.len() == 1 {
                    taxa.into_iter().next().flatten()
                } else {
                    None
                },
            });
        }
    }
    std::fs::create_dir_all(&out_dir).map_err(|e| Error::io(&out_dir, e))?;
    let dblist = vec![relative(&out_dir, base).to_string_lossy().into_owned()];
    let date = cat.manifest.date.clone();
    let full = req.full_db_stats.then_some((
        cat.manifest.totals.num_seqs as u64,
        cat.manifest.totals.total_length,
    ));
    let mut infos = Vec::new();
    for p in planned {
        let mut ids = Vec::new();
        let mut oids: BTreeSet<u32> = BTreeSet::new();
        for d in &p.collections {
            for m in cat.members_of(d) {
                if sequences.as_ref().is_some_and(|s| !s.contains(&m.digest)) {
                    continue;
                }
                ids.push(format!("{d}:{}", m.id_name));
                oids.insert(m.oid);
            }
        }
        if oids.is_empty() {
            if p.kind == ViewKind::Named {
                return Err(Error::empty_selection(
                    0,
                    &format!("view {} has no rows in database {dbname}", p.name),
                ));
            }
            continue;
        }
        let length: u64 = oids
            .iter()
            .map(|&o| cat.oids.get(o as usize).map(|r| r.length).unwrap_or(0))
            .sum();
        let bsl = out_dir.join(format!("{}.bsl", p.name));
        let nal = out_dir.join(format!("{}.nal", p.name));
        let bsl_name = format!("{}.bsl", p.name);
        let bytes = alias::bsl_bytes(&ids, &p.name, &date);
        alias::write_atomic(&bsl, &bytes)?;
        let text = alias::nal_text(&NalSpec {
            title: &p.title,
            dblist: &dblist,
            seqidlist: Some(&bsl_name),
            nseq: Some(oids.len() as u64),
            length: Some(length),
            full_stats: full,
            date: &date,
        });
        alias::write_atomic(&nal, text.as_bytes())?;
        infos.push(ViewInfo {
            name: p.name,
            kind: p.kind,
            collections: p.collections,
            taxid: p.taxid,
            nseq: oids.len() as u64,
            length,
            nal,
            bsl,
        });
    }
    let index = write_index(&out_dir, req.kind, &infos)?;
    Ok(ViewSummary {
        out_dir,
        index,
        views: infos,
    })
}

const INDEX_HEADER: &str = "#view\tkind\tcollections\ttaxid\tnseq\tlength\tnal";
/// Rewrite `views.tsv`: rows of other views are kept; this batch replaces earlier rows of the
/// same name (and, for per-collection views, every earlier per-collection row).
fn write_index(out_dir: &Path, kind: ViewKind, infos: &[ViewInfo]) -> Result<PathBuf, Error> {
    let p = out_dir.join("views.tsv");
    let new_names: HashSet<&str> = infos.iter().map(|v| v.name.as_str()).collect();
    let mut rows: Vec<String> = std::fs::read_to_string(&p)
        .unwrap_or_default()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .filter(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f.len() == 7
                && !new_names.contains(f[0])
                && !(kind == ViewKind::PerCollection && f[1] == kind.as_str())
        })
        .map(String::from)
        .collect();
    for v in infos {
        rows.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            v.name,
            v.kind.as_str(),
            v.collections.join(","),
            v.taxid.map(|t| t.to_string()).unwrap_or_default(),
            v.nseq,
            v.length,
            v.nal.file_name().unwrap_or_default().to_string_lossy()
        ));
    }
    rows.sort();
    let mut t = format!("{INDEX_HEADER}\n");
    for r in rows {
        t += &r;
        t.push('\n');
    }
    alias::write_atomic(&p, t.as_bytes())?;
    Ok(p)
}

/// The rows of a `views.tsv` (name, kind, collections, taxid, nseq, length, nal file).
pub fn read_index(path: &Path) -> Result<Vec<Vec<String>>, Error> {
    let t = std::fs::read_to_string(path).map_err(|e| Error::io(path, e))?;
    let mut out = Vec::new();
    for l in t.lines() {
        if l.starts_with('#') || l.is_empty() {
            continue;
        }
        out.push(l.split('\t').map(String::from).collect());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn safe_names() {
        assert_eq!(safe_name("hg38-p14-ucsc"), "hg38-p14-ucsc");
        assert_eq!(safe_name("-e70JAQq4"), "_e70JAQq4");
        assert_eq!(safe_name(".x"), "_x");
        assert_eq!(safe_name("GRCh38 p14/x:y"), "GRCh38_p14_x_y");
        assert_eq!(safe_name(""), "_");
    }

    #[test]
    fn relative_paths() {
        let d = tempfile::tempdir().unwrap();
        let base = d.path().join("db");
        let v = d.path().join("db.views");
        std::fs::create_dir_all(&v).unwrap();
        assert_eq!(relative(&v, &base), PathBuf::from("../db"));
        let deep = d.path().join("a/b");
        std::fs::create_dir_all(&deep).unwrap();
        assert_eq!(relative(&deep, &base), PathBuf::from("../../db"));
        assert_eq!(relative(d.path(), &base), PathBuf::from("db"));
    }
}
