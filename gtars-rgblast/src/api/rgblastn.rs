//! `api::rgblastn`: the one in-process call the `rgblastn` binary, `run::RgblastnRunner` and
//! the bindings make. Opens the store, resolves the selection, builds the `SubjectPlan`, maps
//! the subjects (`store::access`), runs `engine::search::search` and fills the aligned
//! sequence strings when the output asks for them. Prints nothing;
//! progress and blastn's warnings go through the callback; the caller writes the tabular
//! output with `format::tabular::write`.

use crate::engine::query::QuerySet;
use crate::engine::search::{SearchParams, SearchResults, search};
use crate::engine::subject::StoreSubjects;
use crate::error::Error;
use crate::format::tabular::TabularSpec;
use crate::plan::{PlanOptions, ResolvedSelection, SubjectPlan};
use crate::progress::ProgressEvent;
use crate::select::{ResolveOptions, Selection};
use crate::store::{OpenedStore, RemoteOptions};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum QuerySource {
    Path(PathBuf),
    Bytes(Vec<u8>),
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RgblastnRequest {
    pub store: String,
    /// `--cache-dir`: the cache root a remote store's metadata and prefetched sequences are
    /// kept under (None: gtars-refget's default root). Unused for a local store.
    pub cache_dir: Option<PathBuf>,
    pub selection: Selection,
    pub query: QuerySource,
    pub params: SearchParams,
    /// blastn's `-outfmt` string (parsed by `format::tabular::TabularSpec::parse_outfmt`).
    pub outfmt: String,
    /// How a remote store is downloaded: the subjects are pulled into the cache and the engine
    /// maps the cached files.
    #[serde(default)]
    pub remote: RemoteOptions,
}

impl RgblastnRequest {
    pub fn new(store: impl Into<String>, selection: Selection, query: QuerySource) -> Self {
        RgblastnRequest {
            store: store.into(),
            cache_dir: None,
            selection,
            query,
            params: SearchParams::default(),
            outfmt: "6".into(),
            remote: RemoteOptions::default(),
        }
    }
}

pub struct RgblastnOutcome {
    pub results: SearchResults,
    /// The subjects searched (after `--sequences-file`), with the deflines the output shows.
    pub plan: Arc<SubjectPlan>,
    pub subjects: u32,
    pub total_length: u64,
    pub store_label: String,
    /// The parsed `-outfmt`.
    pub spec: TabularSpec,
}

/// Open the store, resolve the selection and run.
pub fn rgblastn(
    req: &RgblastnRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<RgblastnOutcome, Error> {
    req.selection.validate()?;
    let _ = TabularSpec::parse_outfmt(&req.outfmt)?;
    let mut store =
        super::open_store_with(&req.store, req.cache_dir.as_deref(), &req.remote, progress)?;
    let (sel, _n) = super::resolve_in_store(
        &mut store,
        &req.selection,
        &ResolveOptions::default(),
        progress,
    )?;
    run_resolved(&mut store, &sel, req, progress)
}

/// Run on an opened store with a selection already resolved (what `run::RgblastnRunner`
/// has). `req.store` and `req.selection` are not read.
pub fn run_resolved(
    store: &mut OpenedStore,
    sel: &ResolvedSelection,
    req: &RgblastnRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<RgblastnOutcome, Error> {
    let spec = TabularSpec::parse_outfmt(&req.outfmt)?;
    let label = store.label();
    let plan = super::plan_for_selection(store, sel, &PlanOptions::default())?;
    if plan.subjects.is_empty() {
        return Err(Error::empty_selection(
            0,
            "the selection holds no searchable sequences",
        ));
    }
    let plan = Arc::new(plan);
    let access = crate::store::access::open(store, &plan)?;
    let subjects = StoreSubjects::new(plan.clone(), access, req.params.stats);
    let queries = match &req.query {
        QuerySource::Path(p) => QuerySet::from_fasta_path(p)?,
        QuerySource::Bytes(b) => QuerySet::from_fasta_bytes(b)?,
    };
    let totals = plan.totals;
    let mut results = search(&req.params, &queries, &subjects, progress)?;
    if spec.needs_sequences() {
        crate::format::alignment::fill(&mut results, &queries, &subjects)?;
    }
    Ok(RgblastnOutcome {
        results,
        plan,
        subjects: totals.num_seqs,
        total_length: totals.total_length,
        store_label: label,
        spec,
    })
}

/// Render an outcome's results with `format::tabular::write`.
pub fn write_tabular(out: &RgblastnOutcome, w: &mut dyn std::io::Write) -> Result<(), Error> {
    crate::format::tabular::write(&out.results, &out.plan, &out.spec, w)
}
