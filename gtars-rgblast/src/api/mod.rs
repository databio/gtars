//! The binding surface: `preflight`, `build`, `verify`, `view`, `nstats`. Plain data
//! in and out, no printing, progress through a callback that is always called on the calling
//! thread. Returning `false` from the callback cancels: the call returns `Error::Cancelled` and
//! removes partial output.

use crate::blastdb::preflight::{self, PreflightSummary};
use crate::blastdb::sidecars::{
    FORMAT_VERSION, Manifest, ManifestOptions, ManifestSelection, ManifestStore,
};
use crate::blastdb::verify::{VerifyReport, verify_db};
use crate::blastdb::{self, BlastDate, BuildSummary, WriteOptions};
use crate::error::Error;
use crate::plan::{
    DbTotals, OnOversize, PlanOptions, ResolvedSelection, SUBJECT_PLAN_VERSION, SubjectPlan,
};
use crate::progress::ProgressEvent;
use crate::select::nstats::NstatsOptions;
use crate::select::{Catalog, ResolveOptions, Selection, StoreCatalog};
use crate::store::{OpenedStore, RemoteOptions, StoreFiles};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::path::PathBuf;
use std::time::Instant;

mod origin;
pub mod rgblastn;

pub use crate::views::{ViewRequest, ViewSummary, view};
pub use origin::{
    OriginOutcome, OriginRequest, ReportRequest, RepsRequest, SearchOutcome, SearchRequest,
    StoreMeta, origin, origin_with, report, reps, search,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BuildRequest {
    /// Local path or http(s) URL.
    pub store: String,
    /// The cache root a remote store downloads into (None = gtars-refget's shared default).
    pub cache_dir: Option<PathBuf>,
    pub selection: Selection,
    pub options: BuildOptions,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct BuildOptions {
    pub out: PathBuf,
    /// None = "<store label>: <N> collections, <M> sequences".
    pub title: Option<String>,
    /// Volume size limit, as makeblastdb `-max_file_sz`.
    pub max_file_size: u64,
    /// Conversion workers (local and remote stores); 0 = all cores.
    pub threads: usize,
    /// How a remote store is downloaded.
    pub remote: RemoteOptions,
    /// Converted bytes held in flight (`--budget-mb`).
    pub budget_bytes: u64,
    pub on_oversize: OnOversize,
    pub force: bool,
    pub check_free_space: bool,
    /// makeblastdb date text; None = SOURCE_DATE_EPOCH, else now.
    pub date: Option<String>,
    pub progress_interval_ms: u64,
    /// Longest sequence kept (default `i32::MAX`, BLAST's limit); lower only for testing.
    pub max_subject_len: u64,
    /// `--taxid-map FILE`: per-collection taxids that override FHR.
    pub taxid_map: Option<PathBuf>,
    /// Where N counts are cached and read from.
    pub nstats: NstatsOptions,
}

impl Default for BuildOptions {
    fn default() -> Self {
        BuildOptions {
            out: PathBuf::new(),
            title: None,
            max_file_size: 3_000_000_000,
            threads: 0,
            remote: RemoteOptions::default(),
            budget_bytes: 2 << 30,
            on_oversize: OnOversize::Skip,
            force: false,
            check_free_space: true,
            date: None,
            progress_interval_ms: 2000,
            max_subject_len: i32::MAX as u64,
            taxid_map: None,
            nstats: NstatsOptions::default(),
        }
    }
}

impl BuildOptions {
    /// Local worker count with 0 turned into the number of cores.
    pub fn resolved_threads(&self) -> usize {
        resolve_threads(self.threads)
    }

    pub fn plan_options(&self) -> PlanOptions {
        PlanOptions {
            on_oversize: self.on_oversize,
            max_subject_len: self.max_subject_len,
        }
    }

    pub fn resolve_options(&self) -> ResolveOptions {
        ResolveOptions {
            taxid_map_file: self.taxid_map.clone(),
            nstats: self.nstats.clone(),
            threads: self.threads,
        }
    }

    fn validate(&self) -> Result<(), Error> {
        if self.out.as_os_str().is_empty() {
            return Err(Error::InvalidOptions(
                "an output base path (--out) is required".into(),
            ));
        }
        if self.max_file_size == 0 {
            return Err(Error::InvalidOptions(
                "--max-file-size must be greater than 0".into(),
            ));
        }
        if self.budget_bytes == 0 {
            return Err(Error::InvalidOptions(
                "--budget-mb must be at least 1".into(),
            ));
        }
        Ok(())
    }
}

pub fn resolve_threads(t: usize) -> usize {
    if t > 0 {
        t
    } else {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4)
    }
}

/// An opened store and the plan for a request.
pub struct Planned {
    pub store: OpenedStore,
    pub plan: SubjectPlan,
    pub write: WriteOptions,
    pub resolved: ResolvedSelection,
    pub plan_s: f64,
}

/// Open the store, emitting `Opened`. `cache_dir` is the cache root (None = the shared
/// default).
pub(crate) fn open_store_with(
    spec: &str,
    cache_dir: Option<&std::path::Path>,
    remote: &RemoteOptions,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<OpenedStore, Error> {
    let store = OpenedStore::open_with(spec, cache_dir, remote)?;
    let opened = ProgressEvent::Opened {
        store: store.location_string(),
        mode: format!("{:?}", store.mode()),
        remote: store.is_remote(),
    };
    if !progress(&opened) {
        return Err(Error::Cancelled);
    }
    Ok(store)
}

/// Resolve a selection against an opened store; returns it and the store's collection count.
pub fn resolve_in_store(
    store: &mut OpenedStore,
    sel: &Selection,
    opts: &ResolveOptions,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<(ResolvedSelection, usize), Error> {
    let mut cat = StoreCatalog::new(store, opts);
    let n = cat.collections().len();
    let r = crate::select::resolve(&mut cat, sel, opts, progress)?;
    Ok((r, n))
}

/// Build the plan for a resolved selection (loading what it needs), applying
/// `--sequences-file` with `SubjectPlan::retain`.
pub fn plan_for_selection(
    store: &mut OpenedStore,
    sel: &ResolvedSelection,
    opts: &PlanOptions,
) -> Result<SubjectPlan, Error> {
    store.load_for(sel)?;
    let plan = SubjectPlan::build(&store.store, sel, opts)?;
    Ok(match &sel.sequences {
        Some(seqs) => {
            let keep: HashSet<&str> = seqs.iter().map(|s| s.as_str()).collect();
            plan.retain(|d| keep.contains(d))
        }
        None => plan,
    })
}

/// The manifest of a build (the writer fills volumes, files and totals).
pub(crate) fn manifest_template(
    store: &OpenedStore,
    sel: &ResolvedSelection,
    title: &str,
    o: &ManifestOptions,
) -> Manifest {
    Manifest {
        format_version: FORMAT_VERSION,
        subject_plan_version: SUBJECT_PLAN_VERSION,
        tool: format!("rgmakeblastdb {}", env!("CARGO_PKG_VERSION")),
        title: title.to_string(),
        date: String::new(),
        store: ManifestStore {
            location: store.location_string(),
            digests: store.digests(),
        },
        selection: ManifestSelection {
            given: sel.given.clone(),
            canonical: sel.canonical.clone(),
            digest: sel.digest.clone(),
            excluded: sel.excluded.clone(),
        },
        taxid_map: sel.taxid_map.clone(),
        options: o.clone(),
        volumes: Vec::new(),
        files: Vec::new(),
        totals: DbTotals::default(),
    }
}

/// Open the store, resolve and load the selection, build the plan, resolve write options.
pub fn plan_request(
    req: &BuildRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<Planned, Error> {
    let t0 = Instant::now();
    req.options.validate()?;
    req.selection.validate()?;
    let o = &req.options;
    let mut store = open_store_with(&req.store, req.cache_dir.as_deref(), &o.remote, progress)?;
    let (sel, n_store) =
        resolve_in_store(&mut store, &req.selection, &o.resolve_options(), progress)?;
    let plan = plan_for_selection(&mut store, &sel, &o.plan_options())?;
    let write = write_options(&store, &sel, &plan, n_store, o)?;
    Ok(Planned {
        store,
        plan,
        write,
        resolved: sel,
        plan_s: t0.elapsed().as_secs_f64(),
    })
}

/// The writer's options and the manifest template for a plan.
pub(crate) fn write_options(
    store: &OpenedStore,
    sel: &ResolvedSelection,
    plan: &SubjectPlan,
    n_store: usize,
    o: &BuildOptions,
) -> Result<WriteOptions, Error> {
    let date = BlastDate::resolve(o.date.as_deref())?;
    let title = o.title.clone().unwrap_or_else(|| {
        format!(
            "{}: {} collections, {} sequences",
            store.label(),
            plan.collections.len(),
            plan.subjects.len()
        )
    });
    let mopts = ManifestOptions {
        max_file_size: o.max_file_size,
        max_subject_len: o.max_subject_len,
        on_oversize: o.on_oversize,
    };
    let manifest = manifest_template(store, sel, &title, &mopts);
    Ok(WriteOptions {
        out: o.out.clone(),
        title,
        max_file_size: o.max_file_size,
        budget_bytes: o.budget_bytes,
        threads: o.resolved_threads(),
        batch_bytes: crate::pipeline::DEFAULT_BATCH_BYTES,
        force: o.force,
        check_free_space: o.check_free_space,
        date,
        progress_interval_ms: o.progress_interval_ms,
        store: store.location_string(),
        store_label: store.label(),
        remote: store.is_remote(),
        selection: sel.summary(n_store),
        nstats: o.nstats.clone(),
        manifest: Some(manifest),
    })
}

/// Everything `build` would do before writing: open, plan, size. Writes nothing.
pub fn preflight(req: &BuildRequest) -> Result<PreflightSummary, Error> {
    let p = plan_request(req, &mut |_| true)?;
    let (_, mut s) = blastdb::plan_preflight(&p.plan, &p.write);
    if p.write.check_free_space {
        let dir = s.out.parent().map(PathBuf::from).unwrap_or_default();
        s.free_bytes = preflight::free_space(&dir);
    }
    Ok(s)
}

/// Open the store, resolve the selection, build the `SubjectPlan`, preflight and write.
pub fn build(
    req: &BuildRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<BuildSummary, Error> {
    let t0 = Instant::now();
    let p = plan_request(req, progress)?;
    write_planned(&p, t0, progress)
}

/// Pull a remote store's subjects into the cache (emitting `Pulling`); `None` for a local
/// store.
pub(crate) fn pull_plan(
    store: &OpenedStore,
    subjects: &[crate::plan::Subject],
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<Option<gtars_refget::store::FetchCounters>, Error> {
    if !store.is_remote() {
        return Ok(None);
    }
    let ev = ProgressEvent::Pulling {
        sequences: subjects.len() as u64,
        bases: subjects.iter().map(|s| s.length as u64).sum(),
    };
    if !progress(&ev) {
        return Err(Error::Cancelled);
    }
    let digests: Vec<String> = subjects.iter().map(|s| s.digest.clone()).collect();
    store.pull(&digests)
}

/// Write a planned database (a remote store's subjects pulled into the cache first); emits
/// `Finished`.
pub(crate) fn write_planned(
    p: &Planned,
    t0: Instant,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<BuildSummary, Error> {
    let notes: Vec<String> = p.resolved.warnings.clone();
    for n in &notes {
        if !progress(&ProgressEvent::Warning(n.clone())) {
            return Err(Error::Cancelled);
        }
    }
    let network = pull_plan(&p.store, &p.plan.subjects, progress)?;
    let fetch = StoreFiles::new(&p.store);
    let mut s = blastdb::write_blastdb(&p.plan, &fetch, &p.write, progress)?;
    s.network = network.as_ref().map(blastdb::NetworkReport::from);
    s.warnings.splice(0..0, notes);
    s.timings.plan_s = p.plan_s;
    s.elapsed_s = t0.elapsed().as_secs_f64();
    progress(&ProgressEvent::Finished(Box::new(s.clone())));
    Ok(s)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VerifyRequest {
    pub db: PathBuf,
    pub store: String,
    /// The cache root a remote store downloads into (None = the shared default).
    pub cache_dir: Option<PathBuf>,
    /// Workers; 0 = all cores.
    pub threads: usize,
    /// How a remote store is downloaded.
    #[serde(default)]
    pub fetch: RemoteOptions,
}

/// Check a database against its store: every id in `<db>.ids.tsv` resolves to its OID, every
/// sequence has the store's md5, and the totals agree.
pub fn verify(
    req: &VerifyRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<VerifyReport, Error> {
    let mut store = open_store_with(&req.store, req.cache_dir.as_deref(), &req.fetch, progress)?;
    store.store.load_sequence_index()?;
    let threads = resolve_threads(req.threads);
    let fetch = StoreFiles::new(&store);
    verify_db(&req.db, &store.store, &fetch, threads, progress)
}

/// `rgmakeblastdb nstats`: fill the N-count cache for a selection and report each collection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NstatsRequest {
    pub store: String,
    /// The cache root a remote store downloads into (None = the shared default).
    pub cache_dir: Option<PathBuf>,
    /// How a remote store is downloaded.
    #[serde(default)]
    pub fetch: RemoteOptions,
    pub selection: Selection,
    pub nstats: NstatsOptions,
    /// `--taxid-map`, for the species column.
    pub taxid_map: Option<PathBuf>,
    /// Workers; 0 = all cores.
    pub threads: usize,
}

/// One collection of an nstats report.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NstatsRow {
    pub digest: String,
    pub display: String,
    pub species: Option<String>,
    /// Nucleotide rows.
    pub rows: u64,
    pub bases: u64,
    pub n: u64,
    pub n_fraction: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NstatsSummary {
    /// Sorted by N fraction, descending.
    pub collections: Vec<NstatsRow>,
    /// Distinct nucleotide sequences counted.
    pub sequences: u64,
    pub excluded: Vec<crate::plan::Excluded>,
}

pub fn nstats(
    req: &NstatsRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<NstatsSummary, Error> {
    req.selection.validate()?;
    let mut store = open_store_with(&req.store, req.cache_dir.as_deref(), &req.fetch, progress)?;
    store.store.load_sequence_index()?;
    let ro = ResolveOptions {
        taxid_map_file: req.taxid_map.clone(),
        nstats: req.nstats.clone(),
        threads: req.threads,
    };
    let (sel, _) = resolve_in_store(&mut store, &req.selection, &ro, progress)?;
    let aliases = crate::plan::order::collection_aliases(&store.store);
    let mut rows_of = Vec::new();
    let mut all: Vec<String> = Vec::new();
    for d in &sel.collections {
        let rows: Vec<crate::store::CollectionRow> = store
            .collection_rows(d)?
            .into_iter()
            .filter(|r| {
                store
                    .store
                    .get_sequence_metadata(&r.digest)
                    .and_then(|m| crate::convert::nucleotide_bits(m.alphabet))
                    .is_some()
            })
            .collect();
        all.extend(rows.iter().map(|r| r.digest.clone()));
        rows_of.push((d.clone(), rows));
    }
    all.sort();
    all.dedup();
    let counts = crate::select::nstats::ensure(
        &store,
        &all,
        &req.nstats,
        resolve_threads(req.threads),
        progress,
    )?;
    let mut out = Vec::new();
    for (d, rows) in rows_of {
        let (mut n, mut bases, mut known) = (0u64, 0u64, true);
        for r in &rows {
            match counts.get(&r.digest) {
                Some(c) => {
                    n += c.n;
                    bases += c.length;
                }
                None => {
                    known = false;
                    bases += r.length;
                }
            }
        }
        out.push(NstatsRow {
            display: aliases
                .get(&d)
                .and_then(|v| v.first())
                .map(|(_, a)| a.clone())
                .unwrap_or_else(|| d.clone()),
            species: sel.info.get(&d).and_then(|i| i.species.clone()),
            rows: rows.len() as u64,
            bases,
            n,
            n_fraction: (known && bases > 0).then(|| n as f64 / bases as f64),
            digest: d,
        });
    }
    out.sort_by(|a, b| {
        b.n_fraction
            .unwrap_or(-1.0)
            .total_cmp(&a.n_fraction.unwrap_or(-1.0))
            .then(a.display.cmp(&b.display))
    });
    Ok(NstatsSummary {
        collections: out,
        sequences: all.len() as u64,
        excluded: sel.excluded,
    })
}
