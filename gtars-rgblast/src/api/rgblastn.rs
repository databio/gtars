//! `api::rgblastn`: the one in-process call the `rgblastn` binary, `run::RgblastnRunner` and
//! the bindings make. Opens the store, resolves the selection, builds the `SubjectPlan`, maps
//! the subjects (`store::access`), runs `engine::search::search` and fills the aligned
//! sequence strings when the output asks for them. Prints nothing;
//! progress and blastn's warnings go through the callback; the caller writes the tabular
//! output with `format::tabular::write`.

use crate::engine::query::QuerySet;
use crate::engine::search::{QueryResult, SearchParams, SearchResults, search_each};
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
    let plan = super::plan_for_selection(store, sel, &PlanOptions::default())?;
    run_planned(store, Arc::new(plan), req, progress)
}

/// Run on an opened store with the subjects' plan already built (`plan_for_selection` with
/// the default `PlanOptions`; origin builds it once for all its stages).
pub fn run_planned(
    store: &mut OpenedStore,
    plan: Arc<SubjectPlan>,
    req: &RgblastnRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<RgblastnOutcome, Error> {
    let mut all = Vec::new();
    let mut out = run_planned_each(store, plan, req, progress, &mut |_, _, q| {
        all.push(q);
        Ok(())
    })?;
    out.results.queries = all;
    Ok(out)
}

/// A query's finished result, with the plan and output spec of its search.
pub type QuerySink<'a> =
    dyn FnMut(&SubjectPlan, &TabularSpec, QueryResult) -> Result<(), Error> + 'a;

/// `run_resolved`, handing each query's result to `sink` as it is done
/// (`engine::search::search_each`); the outcome's `results.queries` is empty.
pub fn run_resolved_each(
    store: &mut OpenedStore,
    sel: &ResolvedSelection,
    req: &RgblastnRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
    sink: &mut QuerySink,
) -> Result<RgblastnOutcome, Error> {
    let plan = super::plan_for_selection(store, sel, &PlanOptions::default())?;
    run_planned_each(store, Arc::new(plan), req, progress, sink)
}

/// `run_planned`, handing each query's result to `sink` as it is done
/// (`engine::search::search_each`); the outcome's `results.queries` is empty.
pub fn run_planned_each(
    store: &mut OpenedStore,
    plan: Arc<SubjectPlan>,
    req: &RgblastnRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
    sink: &mut QuerySink,
) -> Result<RgblastnOutcome, Error> {
    let spec = TabularSpec::parse_outfmt(&req.outfmt)?;
    let label = store.label();
    if plan.subjects.is_empty() {
        return Err(Error::empty_selection(
            0,
            "the selection holds no searchable sequences",
        ));
    }
    let access = crate::store::access::open(store, &plan, progress)?;
    let subjects = StoreSubjects::new(plan.clone(), access, req.params.stats);
    let queries = match &req.query {
        QuerySource::Path(p) => QuerySet::from_fasta_path(p)?,
        QuerySource::Bytes(b) => QuerySet::from_fasta_bytes(b)?,
    };
    let totals = plan.totals;
    let mut params = req.params.clone();
    params.drop_edit_scripts |= !spec.needs_edit_scripts();
    let fill = spec.needs_sequences();
    let mut window = Vec::new();
    let results = search_each(&params, &queries, &subjects, progress, &mut |qi, mut q| {
        if fill {
            let plus = &queries.queries[qi].seq;
            crate::format::alignment::fill_query(&mut q, plus, &subjects, &mut window)?;
        }
        sink(&plan, &spec, q)
    })?;
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
