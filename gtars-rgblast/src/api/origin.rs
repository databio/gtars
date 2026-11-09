//! The `rgblast` library calls: `search`, `report`, `origin` and `reps`. They print nothing;
//! the engine's stderr comes back as progress events.

use crate::api::{BuildOptions, open_store_with, plan_for_selection, resolve_in_store};
use crate::blastdb::sidecars;
use crate::error::Error;
use crate::origin::queries::{read_queries, write_queries};
use crate::origin::reps::RepsTable;
use crate::origin::{self, ModeChoice, OriginMeta, OriginOptions, OriginRun, Source};
use crate::plan::ResolvedSelection;
use crate::progress::ProgressEvent;
use crate::report::fanout::DbMeta;
use crate::report::parse::{fasta_names, read_hits_file};
use crate::report::rank::{QueryInput, RankOptions, Report, ReportMeta, StageHits, build_report};
use crate::report::taxonomy::{Taxonomy, taxdump_dir};
use crate::run::outfmt;
use crate::run::passthrough::{self, check_origin, check_search, taxonomy4blast_present};
use crate::run::{EngineChoice, EngineInfo, EngineOutput, RunSpec, SearchTarget, resolve_engine};
use crate::select::nstats::NstatsOptions;
use crate::select::{ResolveOptions, Selection};
use crate::store::{OpenedStore, RemoteOptions};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

fn note(progress: &mut dyn FnMut(&ProgressEvent) -> bool, msg: String) -> Result<(), Error> {
    if progress(&ProgressEvent::Warning(msg)) {
        Ok(())
    } else {
        Err(Error::Cancelled)
    }
}

/// A run directory for stage queries and outputs (removed on drop).
fn rundir() -> Result<tempfile::TempDir, Error> {
    tempfile::Builder::new()
        .prefix("rgblast-run-")
        .tempdir()
        .map_err(|e| Error::io(std::env::temp_dir(), e))
}

/// Open the store (remote stores download into the cache root `cache_dir`) and resolve the
/// selection.
fn open_and_resolve(
    store: &str,
    fetch: &RemoteOptions,
    selection: &Selection,
    cache_dir: Option<&Path>,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<(OpenedStore, ResolvedSelection), Error> {
    selection.validate()?;
    let mut st = open_store_with(store, cache_dir, fetch, progress)?;
    let opts = ResolveOptions {
        nstats: NstatsOptions {
            cache_dir: cache_dir.map(Path::to_path_buf),
            ..Default::default()
        },
        ..Default::default()
    };
    let (sel, _) = resolve_in_store(&mut st, selection, &opts, progress)?;
    for w in &sel.warnings {
        note(progress, w.clone())?;
    }
    Ok((st, sel))
}

// ---------------------------------------------------------------------------------------------
// search

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SearchRequest {
    /// The store rgblastn searches (empty with `EngineChoice::Blastn`, which searches its `db`).
    pub store: String,
    pub selection: Selection,
    pub engine: EngineChoice,
    /// The cache root: a remote store's downloads and N counts (None = the shared default).
    pub cache_dir: Option<PathBuf>,
    /// How a remote store is downloaded.
    #[serde(default)]
    pub fetch: RemoteOptions,
    /// Append the columns `rgblast report` needs to `-outfmt 6/7` (default true).
    pub extra_fields: bool,
    /// Pass-through engine arguments.
    pub args: Vec<String>,
    pub output: EngineOutput,
    /// Where to write the run metadata (`--meta`, or `<-out>.meta.json`).
    pub meta: Option<PathBuf>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StoreMeta {
    pub location: String,
    #[serde(flatten)]
    pub digests: BTreeMap<String, String>,
    /// The selected collections (digests).
    pub collections: Vec<String>,
}

/// The `--meta` JSON (plus captured output and the engine's last stderr lines, which are
/// not written to it).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SearchOutcome {
    pub tool: String,
    pub version: String,
    pub engine: EngineInfo,
    /// The store and selection searched (rgblastn).
    pub store: Option<StoreMeta>,
    /// The database searched (blastn).
    pub db: Option<PathBuf>,
    pub max_target_seqs: usize,
    pub args: Vec<String>,
    pub exit_code: i32,
    #[serde(skip)]
    pub output: Option<String>,
    #[serde(skip)]
    pub engine_stderr: Vec<String>,
}

/// Search the store's selection with rgblastn, or a database with blastn.
pub fn search(
    req: &SearchRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<SearchOutcome, Error> {
    let mut args = req.args.clone();
    let taxonomy4blast = match &req.engine {
        EngineChoice::Blastn { db, .. } => taxonomy4blast_present(db.parent()),
        // rgblastn refuses the taxid flags itself
        EngineChoice::Rgblastn => true,
    };
    check_search(&args, taxonomy4blast)?;
    let runner = resolve_engine(&req.engine)?;
    let (mut opened, target, store) = match &req.engine {
        EngineChoice::Blastn { db, .. } => (None, SearchTarget::Db(db.clone()), None),
        EngineChoice::Rgblastn => {
            let (st, sel) = open_and_resolve(
                &req.store,
                &req.fetch,
                &req.selection,
                req.cache_dir.as_deref(),
                progress,
            )?;
            let meta = StoreMeta {
                location: st.location_string(),
                digests: st.digests(),
                collections: sel.collections.clone(),
            };
            let target = SearchTarget::Store {
                selection: Box::new(sel),
                subjects: None,
            };
            (Some(st), target, Some(meta))
        }
    };
    if req.extra_fields {
        let o = outfmt::rewrite(&mut args);
        if let Some(n) = o.note() {
            note(progress, n)?;
        }
    }
    let threads_added = passthrough::ensure_threads(&mut args).is_some();
    let spec = RunSpec {
        stage: "search".into(),
        target,
        args: args.clone(),
        stats: None,
        query: None,
        output: req.output.clone(),
        threads_added,
    };
    let o = runner.run(&spec, opened.as_mut(), progress)?;
    let outcome = SearchOutcome {
        tool: "rgblast search".into(),
        version: env!("CARGO_PKG_VERSION").into(),
        engine: runner.info(),
        store,
        db: o.db,
        max_target_seqs: passthrough::max_target_seqs(&args),
        args,
        exit_code: o.exit_code,
        output: o.output,
        engine_stderr: o.engine_stderr,
    };
    if let Some(m) = &req.meta {
        let text = serde_json::to_string_pretty(&outcome).unwrap() + "\n";
        std::fs::write(m, text).map_err(|e| Error::io(m, e))?;
    }
    Ok(outcome)
}

// ---------------------------------------------------------------------------------------------
// report

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ReportRequest {
    /// Hit files (`-` = stdin).
    pub hits: Vec<PathBuf>,
    /// A database's sidecars (`--db`); default the meta's database, else its store.
    pub db: Option<PathBuf>,
    /// `--meta`, else `<first hits file>.meta.json` when it exists.
    pub meta: Option<PathBuf>,
    /// The query FASTA, for the query list and order.
    pub query: Option<PathBuf>,
    /// `-outfmt 6` columns (blastn syntax).
    pub fields: String,
    pub taxdump: Option<PathBuf>,
    pub rank: RankOptions,
    /// `--max-target-seqs` (else from the meta, else 500).
    pub max_target_seqs: Option<usize>,
}

fn load_taxonomy(dir: Option<&Path>) -> Result<Option<Taxonomy>, Error> {
    match taxdump_dir(dir) {
        Some(d) => Ok(Some(Taxonomy::load(&d)?)),
        None => Ok(None),
    }
}

/// The collections and members of the store selection a search ran on.
fn store_meta(m: &StoreMeta) -> Result<DbMeta, Error> {
    let sel = Selection::explicit(&m.collections);
    let (mut st, rs) = open_and_resolve(
        &m.location,
        &RemoteOptions::default(),
        &sel,
        None,
        &mut |_| true,
    )?;
    let plan = plan_for_selection(&mut st, &rs, &BuildOptions::default().plan_options())?;
    Ok(DbMeta::from_plan(&plan))
}

/// Turn hit files into the report.
pub fn report(req: &ReportRequest) -> Result<Report, Error> {
    if req.hits.is_empty() {
        return Err(Error::InvalidOptions(
            "report needs at least one hits file".into(),
        ));
    }
    let mut hits = Vec::new();
    let mut listed: Vec<String> = Vec::new();
    for h in &req.hits {
        let r = read_hits_file(h, &req.fields)?;
        hits.extend(r.hits);
        listed.extend(r.queries);
    }
    let meta_path = req.meta.clone().or_else(|| {
        let first = req.hits[0].display().to_string();
        let p = PathBuf::from(format!("{first}.meta.json"));
        (first != "-" && p.exists()).then_some(p)
    });
    let run_meta: Option<SearchOutcome> = match &meta_path {
        Some(p) => {
            let t = std::fs::read_to_string(p).map_err(|e| Error::io(p, e))?;
            Some(serde_json::from_str(&t).map_err(|e| {
                Error::InvalidOptions(format!(
                    "{}: not an 'rgblast search' meta file: {e}",
                    p.display()
                ))
            })?)
        }
        None => None,
    };
    let run_db = run_meta.as_ref().and_then(|m| m.db.as_ref());
    let run_store = run_meta.as_ref().and_then(|m| m.store.as_ref());
    let (meta, db_label) = if let Some(db) = req.db.as_ref().or(run_db) {
        (DbMeta::from_sidecars(db)?, db.display().to_string())
    } else if let Some(s) = run_store {
        (store_meta(s)?, s.location.clone())
    } else {
        return Err(Error::InvalidOptions(
            "report needs the database's collections: give --db BASE (or keep the <HITS>.meta.json that 'rgblast search -- -out HITS' writes)".into(),
        ));
    };
    // a store search covers its selection; a database search, the whole database
    let scope: Option<HashSet<usize>> = run_store.filter(|_| req.db.is_none()).map(|s| {
        s.collections
            .iter()
            .filter_map(|d| meta.coll_index(d))
            .collect()
    });
    let n = scope
        .as_ref()
        .map(|s| s.len())
        .unwrap_or(meta.collections.len());
    let mts = req
        .max_target_seqs
        .or(run_meta.as_ref().map(|m| m.max_target_seqs))
        .unwrap_or(500);
    // query order: --query, else "# Query:" lines, else first appearance
    let mut order: Vec<String> = match &req.query {
        Some(q) => fasta_names(q)?,
        None => listed,
    };
    let known = !order.is_empty();
    let mut by_query: HashMap<String, Vec<crate::report::Hit>> = HashMap::new();
    for h in hits {
        if !known && !by_query.contains_key(&h.query) {
            order.push(h.query.clone());
        }
        by_query.entry(h.query.clone()).or_default().push(h);
    }
    let mut seen = HashSet::new();
    order.retain(|q| seen.insert(q.clone()));
    let inputs: Vec<QueryInput> = order
        .iter()
        .map(|q| QueryInput {
            name: q.clone(),
            qlen: None,
            runs: vec![StageHits {
                stage: "full".into(),
                hits: by_query.remove(q).unwrap_or_default(),
                scope: scope.clone(),
            }],
            searched: (n, n),
            flags: Vec::new(),
            max_target_seqs: mts,
        })
        .collect();
    let mut tax = load_taxonomy(req.taxdump.as_deref())?;
    let info = ReportMeta {
        engine: run_meta
            .as_ref()
            .map(|m| format!("{} {}", m.engine.name, m.engine.version)),
        db: Some(db_label),
        mode: "exhaustive".into(),
        taxonomy: tax.as_ref().map(|t| t.describe()),
        tie_percent: req.rank.tie_percent,
        warnings: Vec::new(),
    };
    build_report(&meta, inputs, &req.rank, tax.as_mut(), info)
}

// ---------------------------------------------------------------------------------------------
// origin

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OriginRequest {
    /// The store rgblastn searches (empty with `EngineChoice::Blastn`, which searches its `db`).
    pub store: String,
    pub selection: Selection,
    /// rgblastn (staged or exhaustive), or blastn on a database (exhaustive only).
    pub engine: EngineChoice,
    /// The cache root: a remote store's downloads and N counts (None = the shared default).
    pub cache_dir: Option<PathBuf>,
    /// How a remote store is downloaded.
    #[serde(default)]
    pub fetch: RemoteOptions,
    /// The query FASTA (`-` = stdin).
    pub query: PathBuf,
    /// Pass-through engine arguments (no `-query`).
    pub args: Vec<String>,
    /// `--rep` selectors (digest or `ns:alias`).
    pub reps: Vec<String>,
    pub options: OriginOptions,
    pub taxdump: Option<PathBuf>,
    pub keep_hits: Option<PathBuf>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OriginOutcome {
    pub report: Report,
    pub meta: OriginMeta,
}

/// Species first, then fan-out (or one exhaustive search), then the report.
pub fn origin(
    req: &OriginRequest,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<OriginOutcome, Error> {
    let runner = resolve_engine(&req.engine)?;
    origin_with(req, runner.as_ref(), progress)
}

/// `origin` with a given runner (tests replay recorded hits through a fake one).
pub fn origin_with(
    req: &OriginRequest,
    runner: &dyn crate::run::BlastRunner,
    progress: &mut dyn FnMut(&ProgressEvent) -> bool,
) -> Result<OriginOutcome, Error> {
    check_origin(&req.args)?;
    let blastn_db = match &req.engine {
        EngineChoice::Blastn { db, .. } => {
            if req.options.mode == ModeChoice::Staged || !req.reps.is_empty() {
                return Err(Error::InvalidOptions(
                    "staged origin needs the default engine, rgblastn; with --engine blastn, origin runs one search of --db".into(),
                ));
            }
            Some(db.clone())
        }
        EngineChoice::Rgblastn => None,
    };
    let queries = read_queries(&req.query)?;
    let mut opts = req.options.clone();
    let mut args = req.args.clone();
    let threads_added = passthrough::ensure_threads(&mut args).is_some();
    let mts = passthrough::max_target_seqs(&args);
    let mut tax = load_taxonomy(req.taxdump.as_deref())?;
    let dir = rundir()?;
    write_queries(&queries, &dir.path().join("queries.fa"), None)?;
    let info = ReportMeta {
        engine: Some(runner.describe()),
        db: None,
        mode: String::new(),
        taxonomy: tax.as_ref().map(|t| t.describe()),
        tie_percent: opts.rank.tie_percent,
        warnings: Vec::new(),
    };
    let (report, meta) = match blastn_db {
        Some(db) => {
            opts.mode = ModeChoice::Exhaustive;
            let run = OriginRun {
                runner,
                meta: DbMeta::from_sidecars(&db)?,
                stats: sidecars::read_manifest(&db)?.totals,
                source: Source::Db(&db),
                queries: &queries,
                args,
                threads_added,
                max_target_seqs: mts,
                rundir: dir.path(),
            };
            origin::run(run, &opts, tax.as_mut(), info, progress)?
        }
        None => {
            let (mut st, sel) = open_and_resolve(
                &req.store,
                &req.fetch,
                &req.selection,
                req.cache_dir.as_deref(),
                progress,
            )?;
            for r in &req.reps {
                opts.reps.push(st.resolve_collection(r)?);
            }
            let plan = plan_for_selection(&mut st, &sel, &BuildOptions::default().plan_options())?;
            let run = OriginRun {
                runner,
                meta: DbMeta::from_plan(&plan),
                stats: plan.totals,
                source: Source::Store {
                    store: &mut st,
                    selection: &sel,
                    plan: &plan,
                },
                queries: &queries,
                args,
                threads_added,
                max_target_seqs: mts,
                rundir: dir.path(),
            };
            origin::run(run, &opts, tax.as_mut(), info, progress)?
        }
    };
    if let Some(keep) = &req.keep_hits {
        std::fs::create_dir_all(keep).map_err(|e| Error::io(keep, e))?;
        for e in std::fs::read_dir(dir.path()).map_err(|e| Error::io(dir.path(), e))? {
            let e = e.map_err(|e| Error::io(dir.path(), e))?;
            let to = keep.join(e.file_name());
            std::fs::copy(e.path(), &to).map_err(|err| Error::io(&to, err))?;
        }
        let m = keep.join("origin.meta.json");
        std::fs::write(&m, serde_json::to_string_pretty(&meta).unwrap() + "\n")
            .map_err(|e| Error::io(&m, e))?;
    }
    Ok(OriginOutcome { report, meta })
}

// ---------------------------------------------------------------------------------------------
// reps

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RepsRequest {
    pub store: String,
    pub selection: Selection,
    pub taxdump: Option<PathBuf>,
    pub group_rank: Option<String>,
    pub reps: Vec<String>,
    pub rep_namespace: String,
    /// The cache root: a remote store's downloads and N counts (None = the shared default).
    pub cache_dir: Option<PathBuf>,
    /// How a remote store is downloaded.
    #[serde(default)]
    pub fetch: RemoteOptions,
}

/// The representative table and the staged-cost prediction (no search).
pub fn reps(req: &RepsRequest) -> Result<RepsTable, Error> {
    let (mut st, sel) = open_and_resolve(
        &req.store,
        &req.fetch,
        &req.selection,
        req.cache_dir.as_deref(),
        &mut |_| true,
    )?;
    let overrides = req
        .reps
        .iter()
        .map(|r| st.resolve_collection(r))
        .collect::<Result<Vec<_>, _>>()?;
    let plan = plan_for_selection(&mut st, &sel, &BuildOptions::default().plan_options())?;
    let meta = DbMeta::from_plan(&plan);
    let mut tax = load_taxonomy(req.taxdump.as_deref())?;
    let colls = origin::rep_input(
        &plan,
        &meta,
        &req.rep_namespace,
        tax.as_mut(),
        &req.group_rank,
    )?;
    origin::reps::choose(&colls, &overrides, plan.totals.total_length)
}
