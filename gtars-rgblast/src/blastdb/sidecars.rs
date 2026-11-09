//! Sidecar tables written next to the database (not listed in `.njs`), and readers for the
//! ones later tools consume (`collections.tsv`, `members.tsv`, `oids.tsv`, the manifest).

use super::output::{PartialOutput, with_ext};
use crate::convert::BaseCounts;
use crate::error::Error;
use crate::plan::taxonomy::TaxidMapEntry;
use crate::plan::{DbTotals, Excluded, OnOversize, SubjectPlan};
use crate::select::Selection;
use serde::{Deserialize, Serialize};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Changes only when the `.nsq` encoding changes (recorded in the manifest).
pub const FORMAT_VERSION: u32 = 1;

pub const IDS_HELP: &str = "# ids are looked up with blastdbcmd -entry '<id>' or -entry_batch FILE; for collection digests starting with '-' use the gnl|<digest>|<name> form";

/// One row of `<db>.collections.tsv`. `taxid`, `species`, `assembly` and `n_fraction` are
/// None (an empty field) when unknown.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CollectionMeta {
    pub collection: String,
    pub db_tag: String,
    pub display: String,
    /// (namespace, alias), from the comma-joined `ns:alias` column.
    pub aliases: Vec<(String, String)>,
    pub taxid: Option<i32>,
    pub species: Option<String>,
    pub assembly: Option<String>,
    pub n_sequences: u32,
    pub total_length: u64,
    pub n_fraction: Option<f64>,
}

/// One row of `<db>.members.tsv`: a (collection, name) in the database.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MemberRow {
    pub oid: u32,
    pub digest: String,
    pub collection: String,
    /// Original name.
    pub name: String,
    /// Sanitized tag used in ids.
    pub id_name: String,
}

/// Paths of the five sidecars of `base`, in a fixed order: ids, collections, members, sanitized, skipped.
pub fn paths(base: &Path) -> [PathBuf; 5] {
    super::output::SIDECARS.map(|s| with_ext(base, s))
}

fn create(out: &PartialOutput, p: &Path) -> Result<BufWriter<std::fs::File>, Error> {
    out.record(p);
    Ok(BufWriter::new(
        std::fs::File::create(p).map_err(|e| Error::io(p, e))?,
    ))
}

fn opt<T: ToString>(v: &Option<T>) -> String {
    v.as_ref().map(|x| x.to_string()).unwrap_or_default()
}

/// Write the five sidecars (header-only when empty). `n_fractions` (one per plan collection)
/// overrides each collection's `n_fraction` when given. Returns the paths written.
pub fn write_all(
    base: &Path,
    plan: &SubjectPlan,
    n_fractions: Option<&[Option<f64>]>,
    out: &PartialOutput,
) -> Result<Vec<PathBuf>, Error> {
    let [ids, colls, members, sanitized, skipped] = paths(base);
    let io = |p: &Path| {
        let p = p.to_path_buf();
        move |e| Error::io(p, e)
    };
    write_ids(&ids, plan, out)?;
    let mut w = create(out, &colls)?;
    (|| -> std::io::Result<()> {
        writeln!(w, "#collection\tdb_tag\tdisplay\taliases\ttaxid\tspecies\tassembly\tn_sequences\ttotal_length\tn_fraction")?;
        for (ci, c) in plan.collections.iter().enumerate() {
            let aliases: Vec<String> = c.aliases.iter().map(|(ns, a)| format!("{ns}:{a}")).collect();
            let tsv = |s: &Option<String>| s.as_ref().map(|x| x.replace(['\t', '\n'], " "));
            let m = CollectionMeta {
                collection: c.digest.clone(),
                db_tag: c.db_tag.clone(),
                display: c.display.clone(),
                aliases: Vec::new(),
                taxid: c.taxid,
                species: tsv(&c.species),
                assembly: tsv(&c.assembly),
                n_sequences: c.n_sequences,
                total_length: c.total_length,
                n_fraction: n_fractions
                    .and_then(|f| f.get(ci).copied())
                    .unwrap_or(c.n_fraction),
            };
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                m.collection,
                m.db_tag,
                m.display,
                aliases.join(","),
                opt(&m.taxid),
                opt(&m.species),
                opt(&m.assembly),
                m.n_sequences,
                m.total_length,
                opt(&m.n_fraction)
            )?;
        }
        w.flush()
    })()
    .map_err(io(&colls))?;

    {
        let mut w = create(out, &members)?;
        (|| -> std::io::Result<()> {
            writeln!(w, "#oid\tdigest\tcollection\tname\tid_name")?;
            for (oid, s) in plan.subjects.iter().enumerate() {
                for o in &s.occurrences {
                    let c = &plan.collections[o.collection as usize];
                    writeln!(
                        w,
                        "{oid}\t{}\t{}\t{}\t{}",
                        s.digest, c.digest, o.name, o.id_name
                    )?;
                }
            }
            w.flush()
        })()
        .map_err(io(&members))?;
    }

    let mut w = create(out, &sanitized)?;
    (|| -> std::io::Result<()> {
        writeln!(w, "#kind\tcollection\toriginal\tsanitized")?;
        for s in &plan.sanitized {
            writeln!(
                w,
                "{}\t{}\t{}\t{}",
                s.kind.as_str(),
                opt(&s.collection),
                s.original,
                s.sanitized
            )?;
        }
        w.flush()
    })()
    .map_err(io(&sanitized))?;

    let mut w = create(out, &skipped)?;
    (|| -> std::io::Result<()> {
        writeln!(w, "#sequence_digest\tlength\talphabet\treason\tcollections")?;
        for s in &plan.skipped {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}",
                s.digest,
                s.length,
                s.alphabet,
                s.reason.as_str(),
                s.n_collections
            )?;
        }
        w.flush()
    })()
    .map_err(io(&skipped))?;

    Ok(vec![ids, colls, members, sanitized, skipped])
}

fn write_ids(ids: &Path, plan: &SubjectPlan, out: &PartialOutput) -> Result<(), Error> {
    let mut w = create(out, ids)?;
    (|| -> std::io::Result<()> {
        writeln!(w, "#id\toid\tsequence_digest")?;
        writeln!(w, "{IDS_HELP}")?;
        for (oid, s) in plan.subjects.iter().enumerate() {
            for f in plan.lookup_forms(oid as u32) {
                if f.cli_safe {
                    writeln!(w, "{}\t{oid}\t{}", f.id, s.digest)?;
                }
            }
        }
        w.flush()
    })()
    .map_err(|e| Error::io(ids, e))
}

fn rows(path: &Path, ncol: usize) -> Result<Vec<Vec<String>>, Error> {
    let text = std::fs::read_to_string(path).map_err(|e| Error::io(path, e))?;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<String> = line.split('\t').map(String::from).collect();
        if f.len() != ncol {
            return Err(Error::Store(anyhow::anyhow!(
                "{} line {}: {} columns, expected {ncol}",
                path.display(),
                i + 1,
                f.len()
            )));
        }
        out.push(f);
    }
    Ok(out)
}

fn num<T: std::str::FromStr>(path: &Path, s: &str) -> Result<T, Error> {
    s.parse()
        .map_err(|_| Error::Store(anyhow::anyhow!("{}: '{s}' is not a number", path.display())))
}

fn opt_num<T: std::str::FromStr>(path: &Path, s: &str) -> Result<Option<T>, Error> {
    if s.is_empty() {
        Ok(None)
    } else {
        num(path, s).map(Some)
    }
}

pub fn read_collections(path: &Path) -> Result<Vec<CollectionMeta>, Error> {
    rows(path, 10)?
        .into_iter()
        .map(|f| {
            let aliases = f[3]
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    let (ns, a) = s.split_once(':').unwrap_or(("", s));
                    (ns.to_string(), a.to_string())
                })
                .collect();
            Ok(CollectionMeta {
                collection: f[0].clone(),
                db_tag: f[1].clone(),
                display: f[2].clone(),
                aliases,
                taxid: opt_num(path, &f[4])?,
                species: (!f[5].is_empty()).then(|| f[5].clone()),
                assembly: (!f[6].is_empty()).then(|| f[6].clone()),
                n_sequences: num(path, &f[7])?,
                total_length: num(path, &f[8])?,
                n_fraction: opt_num(path, &f[9])?,
            })
        })
        .collect()
}

pub fn read_members(path: &Path) -> Result<Vec<MemberRow>, Error> {
    rows(path, 5)?
        .into_iter()
        .map(|f| {
            Ok(MemberRow {
                oid: num(path, &f[0])?,
                digest: f[1].clone(),
                collection: f[2].clone(),
                name: f[3].clone(),
                id_name: f[4].clone(),
            })
        })
        .collect()
}

/// `(id, oid, sequence digest)` rows of `<db>.ids.tsv`.
pub fn read_ids(path: &Path) -> Result<Vec<(String, u32, String)>, Error> {
    rows(path, 3)?
        .into_iter()
        .map(|f| Ok((f[0].clone(), num(path, &f[1])?, f[2].clone())))
        .collect()
}

/// One row of `<db>.oids.tsv`: the record of OID order and exact lengths (for view `LENGTH`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OidRow {
    pub oid: u32,
    pub digest: String,
    pub length: u64,
    pub n: u64,
    pub ambiguous: u64,
}

/// Write `<db>.oids.tsv` (`#oid digest length n ambiguous`, OID order) to a temp name, then
/// rename it into place.
pub fn write_oids(base: &Path, rows: &[OidRow], out: &PartialOutput) -> Result<PathBuf, Error> {
    let p = with_ext(base, "oids.tsv");
    let mut s = String::with_capacity(64 * rows.len() + 40);
    s += "#oid\tdigest\tlength\tn\tambiguous\n";
    for r in rows {
        s += &format!(
            "{}\t{}\t{}\t{}\t{}\n",
            r.oid, r.digest, r.length, r.n, r.ambiguous
        );
    }
    out.record(&p);
    super::alias::write_atomic(&p, s.as_bytes())?;
    Ok(p)
}

/// The `oids.tsv` rows of a plan, with per-OID base counts.
pub fn oid_rows(plan: &SubjectPlan, counts: &[BaseCounts]) -> Vec<OidRow> {
    plan.subjects
        .iter()
        .enumerate()
        .map(|(oid, s)| {
            let c = counts.get(oid).copied().unwrap_or_default();
            OidRow {
                oid: oid as u32,
                digest: s.digest.clone(),
                length: s.length as u64,
                n: c.n,
                ambiguous: c.ambiguous,
            }
        })
        .collect()
}

pub fn read_oids(path: &Path) -> Result<Vec<OidRow>, Error> {
    let v: Vec<OidRow> = rows(path, 5)?
        .into_iter()
        .map(|f| {
            Ok(OidRow {
                oid: num(path, &f[0])?,
                digest: f[1].clone(),
                length: num(path, &f[2])?,
                n: num(path, &f[3])?,
                ambiguous: num(path, &f[4])?,
            })
        })
        .collect::<Result<_, Error>>()?;
    for (i, r) in v.iter().enumerate() {
        if r.oid as usize != i {
            return Err(Error::Store(anyhow::anyhow!(
                "{}: OID {} on row {i}; rows must be in OID order",
                path.display(),
                r.oid
            )));
        }
    }
    Ok(v)
}

/// The kinds of view `rgmakeblastdb view` writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ViewKind {
    PerCollection,
    Named,
}

impl ViewKind {
    pub fn as_str(&self) -> &'static str {
        match self {
            ViewKind::PerCollection => "per-collection",
            ViewKind::Named => "named",
        }
    }
}

/// The store a database was built from: where, and the `rgstore.json` content digests.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ManifestStore {
    pub location: String,
    /// Whichever of the four `rgstore.json` digests the store has.
    #[serde(default)]
    pub digests: std::collections::BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ManifestSelection {
    pub given: Selection,
    pub canonical: serde_json::Value,
    pub digest: String,
    pub excluded: Vec<Excluded>,
}

/// The build options a database was written with.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ManifestOptions {
    pub max_file_size: u64,
    pub max_subject_len: u64,
    pub on_oversize: OnOversize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ManifestVolume {
    /// Volume stem (file name without extension), e.g. `db`, `db.00`.
    pub stem: String,
    pub first_oid: u32,
    pub n_oids: u32,
    pub letters: u64,
}

/// `<db>.manifest.json`: how a database was built, and what `view` needs to know about it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Manifest {
    pub format_version: u32,
    pub subject_plan_version: u32,
    pub tool: String,
    pub title: String,
    /// makeblastdb date text.
    pub date: String,
    pub store: ManifestStore,
    pub selection: ManifestSelection,
    pub taxid_map: Vec<TaxidMapEntry>,
    pub options: ManifestOptions,
    pub volumes: Vec<ManifestVolume>,
    /// The volume, LMDB and OID files, by name.
    pub files: Vec<String>,
    pub totals: DbTotals,
}

pub fn manifest_path(base: &Path) -> PathBuf {
    with_ext(base, "manifest.json")
}

/// Write the manifest to a temp name, then rename it into place.
pub fn write_manifest(
    base: &Path,
    m: &Manifest,
    out: Option<&PartialOutput>,
) -> Result<PathBuf, Error> {
    let p = manifest_path(base);
    let mut t = serde_json::to_string_pretty(m).expect("manifest serializes");
    t.push('\n');
    if let Some(o) = out {
        o.record(&p);
    }
    super::alias::write_atomic(&p, t.as_bytes())?;
    Ok(p)
}

/// Read the manifest; a database without one is `NotFound`.
pub fn read_manifest(base: &Path) -> Result<Manifest, Error> {
    let p = manifest_path(base);
    let t = match std::fs::read_to_string(&p) {
        Ok(t) => t,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            return Err(Error::NotFound(format!(
                "{} is not a database written by rgmakeblastdb with a manifest (no {}); rebuild with rgmakeblastdb build --force",
                base.display(),
                p.file_name().unwrap_or_default().to_string_lossy()
            )));
        }
        Err(e) => return Err(Error::io(&p, e)),
    };
    serde_json::from_str(&t)
        .map_err(|e| Error::Store(anyhow::anyhow!("{}: bad manifest: {e}", p.display())))
}
