//! clap argument structs, in the library so every front end (the binaries, `gtars blast`, and
//! the Python and R bindings through `run_json`) parses the same flags. The `TryFrom`
//! conversions check what clap cannot (values that depend on each other or on the store).

mod blast;
mod json;
mod nstats;
mod origin;
mod report;
mod rgblast;
pub mod rgblastn;
mod search;
mod selection;
mod view;

pub use blast::{BlastCommand, db_printer, main_with_verbs};
pub use json::{main_argv, run_json};
pub use nstats::NstatsArgs;
pub use origin::{ORIGIN_AFTER_HELP, OriginArgs, RepsArgs};
pub use report::{FormatArg, RankArgs, ReportArgs, ReportOutputArgs};
pub use rgblast::{
    install_sigint, printer, render_report, run_origin, run_report, run_reps, run_search,
};
pub use rgblastn::{RgblastnArgs, rgblastn_main};
pub use search::{BLASTN_SEARCHES_DB, EngineArg, EngineArgs, SearchArgs};
pub use selection::SelectionArgs;
pub use view::ViewArgs;

use crate::api::{BuildOptions, BuildRequest, VerifyRequest};
use crate::error::Error;
use crate::plan::OnOversize;
use crate::select::Selection;
use crate::select::nstats::NstatsOptions;
use crate::store::{RemoteOptions, is_remote_spec};
use std::path::PathBuf;

/// The N-count options of build and nstats: the cache under `cache_dir` (`--cache-dir`).
pub(crate) fn nstats_options(cache_dir: Option<&PathBuf>) -> NstatsOptions {
    NstatsOptions {
        cache_dir: cache_dir.cloned(),
        write_cache: true,
    }
}

/// The cache root and how a remote store is downloaded: the flags every gtars tool that reads
/// remote stores shares (gtars CLI, refget CLI, rgmakeblastdb, rgblast).
#[derive(clap::Args, Debug, Clone, PartialEq, Default)]
pub struct FetchFlags {
    /// Cache root: sequences and metadata pulled from remote stores, plus N counts in DIR/blast
    /// (default: $REFGET_CACHE_DIR, else $XDG_CACHE_HOME/refget, else ~/.cache/refget)
    #[arg(long, value_name = "DIR")]
    pub cache_dir: Option<PathBuf>,
    /// Remote stores: keep pulled sequences in the cache, so rerunning the same command
    /// downloads only what is missing (always on here: a remote store is pulled into the cache
    /// and read from there)
    #[arg(long)]
    pub resume: bool,
    /// Remote stores: parallel HTTP requests (default 32)
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u64).range(1..))]
    pub fetch_jobs: Option<u64>,
}

impl FetchFlags {
    /// The download flags apply to remote stores only. `--cache-dir` applies to every store:
    /// N counts live under it too.
    pub fn check(&self, store: &str) -> Result<(), Error> {
        if !is_remote_spec(store) {
            let flag = if self.resume {
                Some("--resume")
            } else if self.fetch_jobs.is_some() {
                Some("--fetch-jobs")
            } else {
                None
            };
            if let Some(flag) = flag {
                return Err(invalid(format!(
                    "{flag} only applies to remote stores (--store http(s)://…)"
                )));
            }
        }
        Ok(())
    }

    /// The download settings, with defaults for flags not given.
    pub fn remote(&self) -> RemoteOptions {
        RemoteOptions {
            fetch_jobs: self
                .fetch_jobs
                .map(|n| n as usize)
                .unwrap_or(RemoteOptions::default().fetch_jobs),
        }
    }
}

/// `rgmakeblastdb build`
#[derive(clap::Args, Debug, Clone, PartialEq)]
pub struct BuildArgs {
    /// Local store directory, or an http(s) store URL
    #[arg(long, value_name = "PATH|URL")]
    pub store: String,
    #[command(flatten)]
    pub selection: SelectionArgs,
    /// Per-collection taxids that override FHR: lines of <collection><TAB><taxid>[<TAB><species>].
    /// Unlike makeblastdb's per-id -taxid_map, it tells collections with the same names apart
    #[arg(long, value_name = "FILE")]
    pub taxid_map: Option<PathBuf>,
    /// Database base path (makeblastdb -out)
    #[arg(long, value_name = "BASE")]
    pub out: PathBuf,
    /// Database title (default: "<store>: N collections, M sequences")
    #[arg(long, value_name = "TEXT")]
    pub title: Option<String>,
    /// Volume size limit (as makeblastdb -max_file_sz)
    #[arg(long, value_name = "BYTES", default_value_t = 3_000_000_000, value_parser = clap::value_parser!(u64).range(1..))]
    pub max_file_size: u64,
    /// Conversion threads (default: all cores)
    #[arg(long, value_name = "N")]
    pub threads: Option<usize>,
    #[command(flatten)]
    pub fetch: FetchFlags,
    /// Converted bytes held in memory, in MiB
    #[arg(long, value_name = "N", default_value_t = 2048, value_parser = clap::value_parser!(u64).range(1..))]
    pub budget_mb: u64,
    /// Sequences of 2^31 bp or more (listed in BASE.skipped.tsv)
    #[arg(long, value_enum, default_value_t = OnOversize::Skip)]
    pub on_oversize: OnOversize,
    /// Print the preflight summary and stop
    #[arg(long)]
    pub dry_run: bool,
    /// Replace an existing database at BASE
    #[arg(long)]
    pub force: bool,
    /// Skip the free-space check
    #[arg(long)]
    pub no_space_check: bool,
    /// No progress lines
    #[arg(short, long, conflicts_with = "verbose")]
    pub quiet: bool,
    /// Timing breakdown, network report, file list
    #[arg(short, long)]
    pub verbose: bool,
}

/// `rgmakeblastdb verify`
#[derive(clap::Args, Debug, Clone, PartialEq, Default)]
pub struct VerifyArgs {
    /// Database base path
    #[arg(long, value_name = "BASE")]
    pub db: PathBuf,
    /// The store the database was built from (local path or http(s) URL)
    #[arg(long, value_name = "PATH|URL")]
    pub store: String,
    #[command(flatten)]
    pub fetch: FetchFlags,
    /// Worker threads (default: all cores)
    #[arg(long, value_name = "N")]
    pub threads: Option<usize>,
}

pub(crate) fn invalid(msg: impl Into<String>) -> Error {
    Error::InvalidOptions(msg.into())
}

impl TryFrom<&BuildArgs> for BuildRequest {
    type Error = Error;

    fn try_from(a: &BuildArgs) -> Result<Self, Error> {
        let selection = Selection::try_from(&a.selection)?;
        a.fetch.check(&a.store)?;
        let d = BuildOptions::default();
        Ok(BuildRequest {
            store: a.store.clone(),
            cache_dir: a.fetch.cache_dir.clone(),
            selection,
            options: BuildOptions {
                out: a.out.clone(),
                title: a.title.clone(),
                max_file_size: a.max_file_size,
                threads: a.threads.unwrap_or(0),
                remote: a.fetch.remote(),
                budget_bytes: a.budget_mb << 20,
                on_oversize: a.on_oversize,
                force: a.force,
                check_free_space: !a.no_space_check,
                taxid_map: a.taxid_map.clone(),
                nstats: nstats_options(a.fetch.cache_dir.as_ref()),
                ..d
            },
        })
    }
}

impl TryFrom<&VerifyArgs> for VerifyRequest {
    type Error = Error;

    fn try_from(a: &VerifyArgs) -> Result<Self, Error> {
        a.fetch.check(&a.store)?;
        Ok(VerifyRequest {
            db: a.db.clone(),
            store: a.store.clone(),
            cache_dir: a.fetch.cache_dir.clone(),
            threads: a.threads.unwrap_or(0),
            fetch: a.fetch.remote(),
        })
    }
}
