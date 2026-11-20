//! Resolving a [`Selection`] against a [`Catalog`]: a store ([`StoreCatalog`]) or a built
//! database's sidecars (`views::DbCatalog`).

use super::nstats::{self, NstatsOptions};
use super::{CollectionSet, Selection};
use crate::error::Error;
use crate::plan::taxonomy::{self, TaxidMapEntry};
use crate::plan::{CollectionInfo, Excluded, ResolvedSelection, order};
use crate::progress::ProgressEvent;
use crate::store::OpenedStore;
use gtars_refget::digest::AlphabetType;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::PathBuf;

/// One row of a collection, for resolution.
#[derive(Debug, Clone, PartialEq)]
pub struct CatalogRow {
    pub name: String,
    pub length: u64,
    pub digest: String,
    /// None when the catalog does not know it (a built database holds nucleotides only).
    pub alphabet: Option<AlphabetType>,
}

/// What resolving a selector needs: collection digests and alias lookups.
pub trait Names {
    /// Where the collections come from, for messages.
    fn label(&self) -> String;
    fn has_collection(&self, digest: &str) -> bool;
    /// Collection alias namespaces, in `rgstore.json` order.
    fn namespaces(&self) -> Vec<String>;
    /// The collection an alias names.
    fn lookup_alias(&self, ns: &str, alias: &str) -> Option<String>;
}

/// What resolution needs to know about the collections it chooses from.
pub trait Catalog: Names {
    /// Every collection digest.
    fn collections(&self) -> Vec<String>;
    /// A collection's aliases, in canonical (ranked) order.
    fn aliases(&self, digest: &str) -> Vec<(String, String)>;
    /// Taxid, species and assembly from FHR (before `--taxid-map`), plus a warning when
    /// `taxon.uri` can't be parsed.
    fn info(&mut self, digest: &str) -> Result<(CollectionInfo, Option<String>), Error>;
    /// A collection's rows.
    fn rows(&mut self, digest: &str) -> Result<Vec<CatalogRow>, Error>;
    /// N fraction of each collection; `scan` = false uses only counts already known (cache,
    /// `--nstats`), leaving the others None.
    fn n_fractions(
        &mut self,
        digests: &[String],
        scan: bool,
        progress: &mut dyn FnMut(&ProgressEvent) -> bool,
    ) -> Result<HashMap<String, Option<f64>>, Error>;
    /// Whether a sequence digest is anywhere in the catalog (in a collection or an orphan).
    fn has_sequence(&mut self, digest: &str) -> Result<bool, Error>;
}

/// How to resolve: taxid overrides and where N counts come from.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct ResolveOptions {
    /// `--taxid-map FILE`.
    pub taxid_map_file: Option<PathBuf>,
    pub nstats: NstatsOptions,
    /// Workers for an N scan; 0 = all cores.
    pub threads: usize,
}

/// Resolve a selector: a collection digest; `NS:ALIAS` when the part before the first `:` is a
/// collection alias namespace; otherwise a bare alias unique across namespaces.
pub fn resolve_selector<C: Names + ?Sized>(cat: &C, key: &str) -> Result<String, Error> {
    let key = key.trim();
    if cat.has_collection(key) {
        return Ok(key.to_string());
    }
    let namespaces = cat.namespaces();
    if let Some((ns, alias)) = key.split_once(':')
        && namespaces.iter().any(|n| n == ns)
    {
        return cat
            .lookup_alias(ns, alias)
            .ok_or_else(|| Error::UnknownCollection {
                key: key.into(),
                store: cat.label(),
                namespaces: namespaces.clone(),
            });
    }
    let mut matches: Vec<(String, String)> = Vec::new();
    for ns in &namespaces {
        if let Some(d) = cat.lookup_alias(ns, key) {
            matches.push((ns.clone(), d));
        }
    }
    let mut digests: Vec<&String> = matches.iter().map(|(_, d)| d).collect();
    digests.sort();
    digests.dedup();
    match digests.len() {
        1 => Ok(digests[0].clone()),
        0 => match key.split_once(':') {
            Some((ns, _)) if !ns.is_empty() => Err(Error::NotFound(
                crate::error::unknown_namespace(ns, &namespaces),
            )),
            _ => Err(Error::UnknownCollection {
                key: key.into(),
                store: cat.label(),
                namespaces,
            }),
        },
        _ => Err(Error::InvalidOptions(crate::error::ambiguous_alias(
            key, &matches,
        ))),
    }
}

/// `display` of a collection: its first ranked alias, else its digest.
pub fn display(cat: &dyn Catalog, digest: &str) -> String {
    cat.aliases(digest)
        .first()
        .map(|(_, a)| a.clone())
        .unwrap_or_else(|| digest.to_string())
}

fn norm_species(s: &str) -> String {
    s.trim().replace('_', " ").to_lowercase()
}

/// Resolve `--taxid-map` lines to entries (an unknown selector is an error).
pub fn resolve_taxid_map<C: Names + ?Sized>(
    cat: &C,
    opts: &ResolveOptions,
) -> Result<Vec<TaxidMapEntry>, Error> {
    let Some(p) = &opts.taxid_map_file else {
        return Ok(Vec::new());
    };
    let mut out: Vec<TaxidMapEntry> = Vec::new();
    for l in taxonomy::read_taxid_map(p)? {
        let d = resolve_selector(cat, &l.selector).map_err(|e| {
            Error::InvalidOptions(format!("--taxid-map {} line {}: {e}", p.display(), l.line))
        })?;
        out.retain(|e| e.collection != d);
        out.push(TaxidMapEntry {
            collection: d,
            taxid: l.taxid,
            species: l.species,
        });
    }
    out.sort_by(|a, b| a.collection.cmp(&b.collection));
    Ok(out)
}

/// A collection's info after `--taxid-map`: the map's taxid; its species, else the FHR species
/// when FHR has the same taxid.
fn apply_map(mut i: CollectionInfo, m: Option<&TaxidMapEntry>) -> CollectionInfo {
    if let Some(m) = m {
        let species = m.species.clone().or_else(|| {
            (i.taxid == Some(m.taxid))
                .then(|| i.species.clone())
                .flatten()
        });
        i.taxid = Some(m.taxid);
        i.species = species;
    }
    i
}

fn sha256_hex(b: &[u8]) -> String {
    use sha2::{Digest, Sha256};
    Sha256::digest(b)
        .iter()
        .map(|x| format!("{x:02x}"))
        .collect()
}

/// Resolve a selection (steps: explicit/attributes/all, the N filter, `--sequences-file`, the
/// empty check, collection info).
pub fn resolve(
    cat: &mut dyn Catalog,
    sel: &Selection,
    opts: &ResolveOptions,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<ResolvedSelection, Error> {
    sel.validate()?;
    let map = resolve_taxid_map(&*cat, opts)?;
    let map_of: HashMap<&str, &TaxidMapEntry> =
        map.iter().map(|e| (e.collection.as_str(), e)).collect();
    let mut infos: HashMap<String, CollectionInfo> = HashMap::new();
    let mut warnings_by: BTreeMap<String, String> = BTreeMap::new();
    let mut info_of = |cat: &mut dyn Catalog, d: &str| -> Result<CollectionInfo, Error> {
        if let Some(i) = infos.get(d) {
            return Ok(i.clone());
        }
        let (i, w) = cat.info(d)?;
        if let Some(w) = w {
            warnings_by.insert(d.to_string(), w);
        }
        let i = apply_map(i, map_of.get(d).copied());
        infos.insert(d.to_string(), i.clone());
        Ok(i)
    };

    // 1-3: candidates in a stable order (explicit: as given, first occurrence of a digest)
    let mode: Value;
    let candidates: Vec<String> = match &sel.set {
        CollectionSet::Explicit(selectors) => {
            let mut seen = HashSet::new();
            let mut v = Vec::new();
            for s in selectors {
                let d = resolve_selector(&*cat, s)?;
                if seen.insert(d.clone()) {
                    v.push(d);
                }
            }
            let mut sorted = v.clone();
            sorted.sort();
            mode = json!({ "explicit": sorted });
            v
        }
        CollectionSet::Attributes { namespaces, taxa } => {
            let have = cat.namespaces();
            for ns in namespaces {
                if !have.iter().any(|h| h == ns) {
                    return Err(Error::NotFound(crate::error::unknown_namespace(ns, &have)));
                }
            }
            let mut all = cat.collections();
            all.sort();
            if !namespaces.is_empty() {
                all.retain(|d| cat.aliases(d).iter().any(|(n, _)| namespaces.contains(n)));
            }
            let mut want_taxa: BTreeSet<i32> = BTreeSet::new();
            if !taxa.is_empty() {
                // every collection's taxon (after --taxid-map), to match and to list
                let mut by_taxid: BTreeMap<i32, (BTreeMap<String, usize>, usize)> = BTreeMap::new();
                for d in cat.collections() {
                    let i = info_of(cat, &d)?;
                    if let Some(t) = i.taxid {
                        let e = by_taxid.entry(t).or_default();
                        e.1 += 1;
                        if let Some(sp) = &i.species {
                            *e.0.entry(sp.clone()).or_default() += 1;
                        }
                    }
                }
                let available: Vec<(i32, String, usize)> = by_taxid
                    .iter()
                    .map(|(t, (names, n))| {
                        let name = names
                            .iter()
                            .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
                            .map(|(s, _)| s.clone())
                            .unwrap_or_default();
                        (*t, name, *n)
                    })
                    .collect();
                for t in taxa {
                    let t = t.trim();
                    if !t.is_empty() && t.bytes().all(|b| b.is_ascii_digit()) {
                        match t.parse::<i32>().ok().filter(|x| by_taxid.contains_key(x)) {
                            Some(x) => {
                                want_taxa.insert(x);
                            }
                            None => {
                                return Err(Error::NotFound(crate::error::unknown_taxon(
                                    t, &available,
                                )));
                            }
                        }
                    } else {
                        let n = norm_species(t);
                        let hits: Vec<i32> = by_taxid
                            .iter()
                            .filter(|(_, (names, _))| names.keys().any(|s| norm_species(s) == n))
                            .map(|(t, _)| *t)
                            .collect();
                        if hits.is_empty() {
                            return Err(Error::NotFound(crate::error::unknown_taxon(
                                t, &available,
                            )));
                        }
                        want_taxa.extend(hits);
                    }
                }
                let mut keep = Vec::new();
                for d in all {
                    if info_of(cat, &d)?
                        .taxid
                        .is_some_and(|t| want_taxa.contains(&t))
                    {
                        keep.push(d);
                    }
                }
                all = keep;
            }
            let mut ns_sorted = namespaces.clone();
            ns_sorted.sort();
            ns_sorted.dedup();
            mode = json!({ "attributes": { "namespaces": ns_sorted, "taxa": want_taxa.iter().collect::<Vec<_>>() } });
            all
        }
        CollectionSet::All => {
            mode = json!("all");
            let mut all = cat.collections();
            all.sort();
            all
        }
    };
    let all_set = matches!(sel.set, CollectionSet::All);

    // 4: N fractions (scanned only under an N filter)
    let fractions = cat.n_fractions(&candidates, sel.max_n_fraction.is_some(), progress)?;
    let mut selected = Vec::new();
    let mut excluded = Vec::new();
    for d in &candidates {
        let f = fractions.get(d).copied().flatten();
        if let (Some(max), Some(f)) = (sel.max_n_fraction, f)
            && f >= max
        {
            let disp = display(&*cat, d);
            if matches!(sel.set, CollectionSet::Explicit(_)) {
                let ns = cat
                    .aliases(d)
                    .first()
                    .map(|(n, a)| format!("{n}:{a}"))
                    .unwrap_or_else(|| d.clone());
                return Err(Error::InvalidOptions(format!(
                    "{ns} is {:.1}% N, at or above --max-n-fraction {}; remove it or drop the N filter",
                    f * 100.0,
                    crate::error::fmt_fraction(max)
                )));
            }
            excluded.push(Excluded {
                digest: d.clone(),
                display: disp,
                n_fraction: f,
            });
            continue;
        }
        selected.push(d.clone());
    }
    excluded.sort_by(|a, b| {
        b.n_fraction
            .total_cmp(&a.n_fraction)
            .then(a.display.cmp(&b.display))
    });

    // 6: empty
    if selected.is_empty() {
        let reason = match sel.max_n_fraction {
            Some(max) if !excluded.is_empty() => {
                format!("N fraction \u{2265} {}", crate::error::fmt_fraction(max))
            }
            _ => "no collection matches".to_string(),
        };
        return Err(Error::empty_selection(excluded.len(), &reason));
    }

    // 5: --sequences-file
    let sequences = match &sel.sequences {
        None => None,
        Some(seqs) => {
            let mut want: Vec<String> = seqs.iter().map(|s| super::strip_sq(s)).collect();
            want.sort();
            want.dedup();
            let mut held: HashSet<String> = HashSet::new();
            for d in &selected {
                for r in cat.rows(d)? {
                    held.insert(r.digest);
                }
            }
            for w in &want {
                if !(held.contains(w) || all_set && excluded.is_empty() && cat.has_sequence(w)?) {
                    return Err(Error::InvalidOptions(format!(
                        "--sequences-file: SQ.{w} is in no selected collection"
                    )));
                }
            }
            Some(want)
        }
    };

    // 7: info
    let mut info = BTreeMap::new();
    let mut no_taxon = Vec::new();
    let mut warnings = Vec::new();
    let mut sorted = selected.clone();
    sorted.sort();
    for d in &sorted {
        let mut i = info_of(cat, d)?;
        i.n_fraction = fractions.get(d).copied().flatten();
        if i.taxid.is_none() {
            no_taxon.push(display(&*cat, d));
        }
        info.insert(d.clone(), i);
    }
    for d in &sorted {
        if let Some(w) = warnings_by.get(d) {
            warnings.push(w.clone());
        }
    }
    if !no_taxon.is_empty() {
        no_taxon.sort();
        let n = no_taxon.len();
        let shown: Vec<&str> = no_taxon.iter().take(5).map(|s| s.as_str()).collect();
        warnings.push(format!(
            "{n} collection{} {} no FHR taxon (taxid 0): {}{}",
            if n == 1 { "" } else { "s" },
            if n == 1 { "has" } else { "have" },
            shown.join(", "),
            if n > 5 { ", ..." } else { "" }
        ));
    }
    // Orphans come back with --all, but not when collections were excluded: the plan cannot tell
    // an excluded collection's sequences from true orphans, so they would return through it.
    let include_orphans = all_set && excluded.is_empty();
    let canonical = json!({
        "collections": sorted,
        "include_orphans": include_orphans,
        "max_n_fraction": sel.max_n_fraction,
        "mode": mode,
        "sequences": sequences.as_ref().map(|s| sha256_hex(s.join("\n").as_bytes())),
    });
    let digest = gtars_refget::digest::sha512t24u(
        gtars_refget::digest::canonicalize_json(&canonical).as_bytes(),
    );
    Ok(ResolvedSelection {
        collections: sorted,
        include_orphans,
        sequences,
        info,
        excluded,
        given: sel.clone(),
        canonical,
        digest,
        taxid_map: map,
        warnings,
    })
}

/// [`Catalog`] over an opened store: aliases ranked once, FHR read (local sidecars) or pulled
/// (remote) per collection, rows from the collection indexes, N counts through `nstats`.
pub struct StoreCatalog<'a> {
    pub store: &'a mut OpenedStore,
    aliases: order::ReverseAliases,
    nstats: NstatsOptions,
    threads: usize,
}

impl<'a> StoreCatalog<'a> {
    pub fn new(store: &'a mut OpenedStore, opts: &ResolveOptions) -> Self {
        let aliases = order::collection_aliases(&store.store);
        StoreCatalog {
            store,
            aliases,
            nstats: opts.nstats.clone(),
            threads: crate::api::resolve_threads(opts.threads),
        }
    }
}

impl Names for StoreCatalog<'_> {
    fn label(&self) -> String {
        self.store.location_string()
    }

    fn has_collection(&self, digest: &str) -> bool {
        Names::has_collection(&*self.store, digest)
    }

    fn namespaces(&self) -> Vec<String> {
        Names::namespaces(&*self.store)
    }

    fn lookup_alias(&self, ns: &str, alias: &str) -> Option<String> {
        Names::lookup_alias(&*self.store, ns, alias)
    }
}

impl Catalog for StoreCatalog<'_> {
    fn collections(&self) -> Vec<String> {
        self.store
            .store
            .list_collections(0, usize::MAX, &[])
            .map(|p| p.results.into_iter().map(|m| m.digest).collect())
            .unwrap_or_default()
    }

    fn aliases(&self, digest: &str) -> Vec<(String, String)> {
        self.aliases.get(digest).cloned().unwrap_or_default()
    }

    fn info(&mut self, digest: &str) -> Result<(CollectionInfo, Option<String>), Error> {
        let fhr = self.store.fhr(digest)?;
        let aliases = self.aliases(digest);
        let disp = aliases
            .first()
            .map(|(_, a)| a.clone())
            .unwrap_or_else(|| digest.to_string());
        Ok(taxonomy::info_from_fhr(&disp, fhr.as_ref(), &aliases))
    }

    fn rows(&mut self, digest: &str) -> Result<Vec<CatalogRow>, Error> {
        let rows = self.store.collection_rows(digest)?;
        Ok(rows
            .into_iter()
            .map(|r| {
                let alphabet = self
                    .store
                    .store
                    .get_sequence_metadata(&r.digest)
                    .map(|m| m.alphabet);
                CatalogRow {
                    name: r.name,
                    length: r.length,
                    digest: r.digest,
                    alphabet,
                }
            })
            .collect())
    }

    fn n_fractions(
        &mut self,
        digests: &[String],
        scan: bool,
        progress: &mut dyn FnMut(&ProgressEvent) -> bool,
    ) -> Result<HashMap<String, Option<f64>>, Error> {
        let mut rows_of = Vec::with_capacity(digests.len());
        let mut seqs: Vec<String> = Vec::new();
        for d in digests {
            let rows: Vec<CatalogRow> = self
                .rows(d)?
                .into_iter()
                .filter(|r| {
                    r.alphabet
                        .and_then(crate::convert::nucleotide_bits)
                        .is_some()
                })
                .collect();
            seqs.extend(rows.iter().map(|r| r.digest.clone()));
            rows_of.push((d.clone(), rows));
        }
        seqs.sort();
        seqs.dedup();
        let opts = if scan {
            self.nstats.clone()
        } else {
            NstatsOptions {
                write_cache: false,
                ..self.nstats.clone()
            }
        };
        let counts = if scan {
            nstats::ensure(self.store, &seqs, &opts, self.threads, progress)?
        } else {
            known_counts(self.store, &seqs, &opts)?
        };
        Ok(rows_of
            .into_iter()
            .map(|(d, rows)| {
                let mut n = 0u64;
                let mut len = 0u64;
                for r in &rows {
                    match counts.get(&r.digest) {
                        Some(c) => {
                            n += c.n;
                            len += c.length;
                        }
                        None => return (d, None),
                    }
                }
                (d, (len > 0).then(|| n as f64 / len as f64))
            })
            .collect())
    }

    fn has_sequence(&mut self, digest: &str) -> Result<bool, Error> {
        self.store.store.load_sequence_index()?;
        Ok(self.store.store.get_sequence_metadata(digest).is_some())
    }
}

/// N counts known without reading sequences: dna2bit (zero), the cache, and `--nstats`.
fn known_counts(
    store: &OpenedStore,
    digests: &[String],
    opts: &NstatsOptions,
) -> Result<HashMap<String, nstats::NCount>, Error> {
    let mut out = HashMap::new();
    let mut need: Vec<(String, u64)> = Vec::new();
    for d in digests {
        let Some(m) = store.store.get_sequence_metadata(d) else {
            continue;
        };
        let len = m.length as u64;
        if m.alphabet == AlphabetType::Dna2bit {
            out.insert(
                d.clone(),
                nstats::NCount {
                    length: len,
                    n: 0,
                    ambiguous: 0,
                },
            );
        } else if nstats::needs_read(m.alphabet) {
            need.push((d.clone(), len));
        }
    }
    if need.is_empty() {
        return Ok(out);
    }
    if let Some(c) = opts.cache() {
        let ds: Vec<&str> = need.iter().map(|x| x.0.as_str()).collect();
        let got = c.read(&ds);
        need.retain(|(d, len)| {
            match got.get(d).and_then(|v| v.iter().find(|c| c.length == *len)) {
                Some(x) => {
                    out.insert(d.clone(), *x);
                    false
                }
                None => true,
            }
        });
    }
    Ok(out)
}
