//! Running a search engine: rgblastn in process on a store selection (the default), or stock
//! blastn on a database the user built (`--engine blastn --db BASE`). `search` and `origin`
//! describe a run as a [`RunSpec`]; a [`BlastRunner`] turns it into hits.

pub mod blastn;
pub mod outfmt;
pub mod passthrough;
pub mod rgblastn;

use crate::error::Error;
use crate::plan::{DbTotals, ResolvedSelection};
use crate::progress::ProgressEvent;
use crate::store::OpenedStore;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::sync::Arc;

pub use blastn::BlastnRunner;
pub use rgblastn::RgblastnRunner;

/// What a run searches.
#[derive(Debug, Clone)]
pub enum SearchTarget {
    /// rgblastn: the selection on the store; `subjects` limits it to these digests (no `SQ.`).
    Store {
        selection: Box<ResolvedSelection>,
        subjects: Option<Arc<Vec<String>>>,
    },
    /// blastn: a database the user built, searched as given.
    Db(PathBuf),
}

/// Where the engine's tabular output goes.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub enum EngineOutput {
    /// The caller's stdout (the CLI without `-out`).
    #[default]
    Inherit,
    File(PathBuf),
    /// Returned in the outcome (the bindings' default).
    Capture,
}

#[derive(Debug, Clone)]
pub struct RunSpec {
    /// "search", "stage1", "stage2", "full": names the run in progress events.
    pub stage: String,
    pub target: SearchTarget,
    /// Pass-through arguments, `-outfmt` already final.
    pub args: Vec<String>,
    /// Statistics to search under (rgblastn's `STATS_NSEQ`/`STATS_TOTLEN` override), None =
    /// the target's own.
    pub stats: Option<DbTotals>,
    /// `-query` set by origin (its rewritten query file).
    pub query: Option<PathBuf>,
    pub output: EngineOutput,
    /// `-num_threads` in `args` was added, not given (for the summary line).
    pub threads_added: bool,
}

#[derive(Debug, Clone, Default)]
pub struct RunOutcome {
    pub exit_code: i32,
    /// The database searched (blastn).
    pub db: Option<PathBuf>,
    /// Captured output (`EngineOutput::Capture`).
    pub output: Option<String>,
    /// The last 50 stderr lines.
    pub engine_stderr: Vec<String>,
    /// The engine's final argument list (after `blastn`).
    pub args: Vec<String>,
}

/// The engine, for summaries and `--meta`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EngineInfo {
    pub name: String,
    pub version: String,
    /// "/usr/bin/blastn", "in process".
    pub via: String,
}

pub trait BlastRunner {
    /// "blastn 2.17.0+ (/usr/bin/blastn)"
    fn describe(&self) -> String;
    fn info(&self) -> EngineInfo;
    /// Run one search. `store` is the opened store of a `SearchTarget::Store` run.
    fn run(
        &self,
        spec: &RunSpec,
        store: Option<&mut OpenedStore>,
        progress: &mut dyn FnMut(&ProgressEvent) -> bool,
    ) -> Result<RunOutcome, Error>;
}

/// `--engine`, `--db`, `--blastn`.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub enum EngineChoice {
    /// rgblastn, in process, on the store and selection.
    #[default]
    Rgblastn,
    /// Stock blastn on `db`, a database the user built, searched as given. `blastn` is
    /// `--blastn PATH` (else `$GTARS_RGBLAST_BLASTN`, else blastn on PATH).
    Blastn {
        db: PathBuf,
        blastn: Option<PathBuf>,
    },
}

/// The runner for the engine flags; a missing blastn is `Error::EngineNotFound`.
pub fn resolve_engine(c: &EngineChoice) -> Result<Box<dyn BlastRunner>, Error> {
    match c {
        EngineChoice::Rgblastn => Ok(Box::new(RgblastnRunner)),
        EngineChoice::Blastn { blastn, .. } => {
            Ok(Box::new(BlastnRunner::resolve(blastn.as_deref())?))
        }
    }
}
