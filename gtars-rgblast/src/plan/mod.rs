//! The canonical `SubjectPlan`: which sequences become OIDs, in what order, with which
//! deflines. Both `rgmakeblastdb` and the `rgblastn` engine use it, so they agree on OID
//! order, ids and database totals (blastn's hit list keeps the higher OID on ties, and its
//! statistics use the database length and sequence count).

mod build;
pub mod defline;
pub mod ids;
pub mod order;
pub mod taxonomy;
pub mod titles;

pub use defline::{Defline, LookupForm, SeqId, numeric_tag};
use gtars_refget::digest::AlphabetType;
use gtars_refget::store::StorageMode;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Bump whenever OID order, ids, titles, taxids or sidecar columns change (recorded in the
/// manifest).
pub const SUBJECT_PLAN_VERSION: u32 = 2;

/// A selection resolved against a store (or a database's sidecars): the collections, the
/// subject filter, and what each collection contributes to deflines (taxid, species, assembly,
/// N fraction). Made by `select::resolve`; `SubjectPlan::build` reads it.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ResolvedSelection {
    /// Collection digests, sorted bytewise (`SubjectPlan::build` applies its own order).
    pub collections: Vec<String>,
    /// True only for `--all`: sequences in no collection are subjects too.
    pub include_orphans: bool,
    /// `--sequences-file`: subject digests (sorted, no `SQ.`); None = every sequence.
    pub sequences: Option<Vec<String>>,
    /// Per selected collection: taxid, species, assembly, N fraction.
    pub info: BTreeMap<String, CollectionInfo>,
    /// Collections an N filter removed (attribute and all modes).
    pub excluded: Vec<Excluded>,
    /// The selection as given.
    pub given: crate::select::Selection,
    /// Canonical JSON of the resolved selection (stored in the manifest).
    pub canonical: serde_json::Value,
    /// sha512t24u of the canonical JSON.
    pub digest: String,
    /// `--taxid-map` entries, resolved to digests (every line, also those outside the selection).
    pub taxid_map: Vec<crate::plan::taxonomy::TaxidMapEntry>,
    /// For example "3 collections have no FHR taxon: ...".
    pub warnings: Vec<String>,
}

impl ResolvedSelection {
    /// A bare selection of collection digests with no metadata (every taxid unknown).
    pub fn of(collections: Vec<String>, include_orphans: bool) -> Self {
        let mut collections = collections;
        collections.sort();
        collections.dedup();
        ResolvedSelection {
            collections,
            include_orphans,
            ..Default::default()
        }
    }

    /// One line for the preflight: `67 of 79 collections (12 excluded: N ≥ 0.25, e.g.
    /// hg38-p14-masked-ucsc 54.5%); taxa: Homo sapiens 45, Mus musculus 22`.
    pub fn summary(&self, store_collections: usize) -> String {
        let mut s = format!(
            "{} of {} collections",
            crate::progress::commas(self.collections.len() as u64),
            crate::progress::commas(store_collections as u64)
        );
        if !self.excluded.is_empty() {
            let max = self.given.max_n_fraction.unwrap_or(0.0);
            let worst = self
                .excluded
                .iter()
                .max_by(|a, b| {
                    a.n_fraction
                        .total_cmp(&b.n_fraction)
                        .then(b.display.cmp(&a.display))
                })
                .unwrap();
            s += &format!(
                " ({} excluded: N \u{2265} {}, e.g. {} {:.1}%)",
                self.excluded.len(),
                crate::error::fmt_fraction(max),
                worst.display,
                worst.n_fraction * 100.0
            );
        }
        let mut taxa: BTreeMap<(i32, String), usize> = BTreeMap::new();
        let mut unknown = 0usize;
        for d in &self.collections {
            match self.info.get(d) {
                Some(CollectionInfo {
                    taxid: Some(t),
                    species,
                    ..
                }) => {
                    *taxa
                        .entry((*t, species.clone().unwrap_or_else(|| t.to_string())))
                        .or_default() += 1
                }
                _ => unknown += 1,
            }
        }
        if !taxa.is_empty() || unknown > 0 {
            let mut v: Vec<String> = taxa
                .iter()
                .map(|((_, name), n)| format!("{name} {n}"))
                .collect();
            if unknown > 0 {
                v.push(format!("unknown {unknown}"));
            }
            s += &format!("; taxa: {}", v.join(", "));
        }
        s
    }
}

/// What a collection contributes beyond its rows: copied into `PlannedCollection` and
/// written to the `collections.tsv` columns. Each field is None when unknown.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct CollectionInfo {
    pub taxid: Option<i32>,
    pub species: Option<String>,
    pub assembly: Option<String>,
    pub n_fraction: Option<f64>,
}

/// A collection an N filter left out.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Excluded {
    pub digest: String,
    pub display: String,
    pub n_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SubjectPlan {
    pub store_mode: StorageMode,
    /// Selected collections, in canonical order.
    pub collections: Vec<PlannedCollection>,
    /// OID = index.
    pub subjects: Vec<Subject>,
    pub totals: DbTotals,
    pub skipped: Vec<SkippedSequence>,
    pub sanitized: Vec<SanitizedId>,
    /// Sanitized collection alias strings left out of the ids because more than one
    /// collection in the store uses them (sorted, upper-case duplicates removed).
    pub ambiguous_aliases: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlannedCollection {
    pub digest: String,
    /// (namespace, alias), canonical order.
    pub aliases: Vec<(String, String)>,
    /// Unique per store (case-insensitively) and independent of the selection.
    pub db_tag: String,
    /// First alias in alias order (unsanitized), else the digest.
    pub display: String,
    /// NCBI taxid from FHR (or `--taxid-map`); None when unknown (deflines then carry 0).
    pub taxid: Option<i32>,
    /// FHR `taxon.name`, trimmed.
    pub species: Option<String>,
    /// Public assembly accession: FHR `accessionID.name`, else the first refseq, else insdc alias.
    pub assembly: Option<String>,
    /// Σ N / Σ length over the collection's nucleotide rows, when known.
    pub n_fraction: Option<f64>,
    /// All rows, including skipped sequences.
    pub n_names: u32,
    /// Rows present in the database.
    pub n_sequences: u32,
    /// Rows left out as oversize (the build does not read them, so it can't fill `n_fraction`).
    pub n_oversize: u32,
    /// Summed length of those rows.
    pub total_length: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Subject {
    /// sha512t24u, no "SQ.".
    pub digest: String,
    pub md5: String,
    pub alphabet: AlphabetType,
    /// Always < 2^31, so it fits BLAST's Int4.
    pub length: u32,
    /// Encoded: ceil(len*bits/8); Raw: len; Zstd: None.
    pub stored_len: Option<u64>,
    /// [0] = SQ.<digest> defline, then one per occurrence.
    pub deflines: Vec<Defline>,
    /// Aligned with deflines[1..].
    pub occurrences: Vec<Occurrence>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Occurrence {
    /// Index into SubjectPlan.collections.
    pub collection: u32,
    /// Original name.
    pub name: String,
    /// Sanitized tag used in ids.
    pub id_name: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub struct DbTotals {
    pub num_seqs: u32,
    pub total_length: u64,
    pub max_len: u32,
    pub min_len: u32,
}

impl DbTotals {
    pub fn avg_len(&self) -> u32 {
        (self.total_length / self.num_seqs.max(1) as u64) as u32
    }

    pub fn of(subjects: &[Subject]) -> DbTotals {
        DbTotals {
            num_seqs: subjects.len() as u32,
            total_length: subjects.iter().map(|s| s.length as u64).sum(),
            max_len: subjects.iter().map(|s| s.length).max().unwrap_or(0),
            min_len: subjects.iter().map(|s| s.length).min().unwrap_or(0),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default, clap::ValueEnum)]
pub enum OnOversize {
    #[default]
    Skip,
    Fail,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlanOptions {
    pub on_oversize: OnOversize,
    /// Longest sequence kept; default `i32::MAX` (BLAST lengths are Int4).
    pub max_subject_len: u64,
}

impl Default for PlanOptions {
    fn default() -> Self {
        PlanOptions {
            on_oversize: OnOversize::Skip,
            max_subject_len: i32::MAX as u64,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SkippedSequence {
    pub digest: String,
    pub length: u64,
    pub alphabet: AlphabetType,
    pub reason: SkipReason,
    pub n_collections: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SkipReason {
    Oversize,
    NonNucleotide,
    Empty,
}

impl SkipReason {
    pub fn as_str(&self) -> &'static str {
        match self {
            SkipReason::Oversize => "oversize",
            SkipReason::NonNucleotide => "non_nucleotide",
            SkipReason::Empty => "empty",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SanitizedId {
    pub kind: SanitizedKind,
    /// Collection digest, for names and collection aliases.
    pub collection: Option<String>,
    pub original: String,
    pub sanitized: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SanitizedKind {
    Name,
    CollectionAlias,
    SequenceAlias,
}

impl SanitizedKind {
    pub fn as_str(&self) -> &'static str {
        match self {
            SanitizedKind::Name => "name",
            SanitizedKind::CollectionAlias => "collection_alias",
            SanitizedKind::SequenceAlias => "sequence_alias",
        }
    }
}

impl SubjectPlan {
    /// Build the plan from store metadata only (indexes and aliases; never sequence bytes).
    /// The selected collections and the sequence index must be loaded.
    pub fn build(
        store: &gtars_refget::store::ReadonlyRefgetStore,
        sel: &ResolvedSelection,
        opts: &PlanOptions,
    ) -> Result<SubjectPlan, crate::Error> {
        build::build(store, sel, opts)
    }

    /// Every id string blastdbcmd accepts for this OID, in defline order:
    /// Local -> the string; General -> "gnl|db|tag" then "db:tag".
    pub fn lookup_forms(&self, oid: u32) -> Vec<LookupForm> {
        let mut out = Vec::new();
        for d in &self.subjects[oid as usize].deflines {
            for id in &d.ids {
                id.lookup_forms(&mut out);
            }
        }
        out
    }

    /// Keep the subjects whose digest passes `keep`, in their original relative order (OIDs are
    /// renumbered densely). Each kept subject keeps its full deflines and occurrences, so every
    /// selected collection that holds it is still listed. Recompute `totals` and each collection's
    /// `n_sequences` and `total_length`. `collections`, `skipped` and `sanitized` are kept.
    pub fn retain(&self, keep: impl Fn(&str) -> bool) -> SubjectPlan {
        let subjects: Vec<Subject> = self
            .subjects
            .iter()
            .filter(|s| keep(&s.digest))
            .cloned()
            .collect();
        let mut collections = self.collections.clone();
        for c in &mut collections {
            c.n_sequences = 0;
            c.total_length = 0;
        }
        for s in &subjects {
            for o in &s.occurrences {
                let c = &mut collections[o.collection as usize];
                c.n_sequences += 1;
                c.total_length += s.length as u64;
            }
        }
        SubjectPlan {
            store_mode: self.store_mode,
            totals: DbTotals::of(&subjects),
            collections,
            subjects,
            skipped: self.skipped.clone(),
            sanitized: self.sanitized.clone(),
            ambiguous_aliases: self.ambiguous_aliases.clone(),
        }
    }

    pub fn n_deflines(&self) -> u64 {
        self.subjects.iter().map(|s| s.deflines.len() as u64).sum()
    }

    /// Number of lookup forms over all OIDs (the size of the id index users can query).
    pub fn n_lookup_ids(&self) -> u64 {
        self.subjects
            .iter()
            .flat_map(|s| s.deflines.iter())
            .flat_map(|d| d.ids.iter())
            .map(|id| match id {
                SeqId::Local(_) => 1,
                SeqId::General { .. } => 2,
            })
            .sum()
    }

    pub fn n_names(&self) -> u64 {
        self.collections.iter().map(|c| c.n_names as u64).sum()
    }

    pub fn n_orphans(&self) -> u64 {
        self.subjects
            .iter()
            .filter(|s| s.occurrences.is_empty())
            .count() as u64
    }

    /// Selected collections whose digest starts with '-' (their `db:tag` ids look like flags).
    pub fn n_dash_leading_collections(&self) -> u64 {
        self.collections
            .iter()
            .filter(|c| c.digest.starts_with('-'))
            .count() as u64
    }

    pub fn n_skipped(&self, reason: SkipReason) -> u64 {
        self.skipped.iter().filter(|s| s.reason == reason).count() as u64
    }
}
