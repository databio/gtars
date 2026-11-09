//! `rgblast search`, and the engine flags shared with `origin`.

use super::{FetchFlags, SelectionArgs, invalid};
use crate::api::SearchRequest;
use crate::error::Error;
use crate::run::passthrough::value;
use crate::run::{EngineChoice, EngineOutput};
use crate::select::Selection;
use std::path::PathBuf;

#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EngineArg {
    #[default]
    Rgblastn,
    Blastn,
}

/// What `--engine blastn` says about selection flags.
pub const BLASTN_SEARCHES_DB: &str = "with --engine blastn, select collections when building the database, or build a view with `rgmakeblastdb view`, and pass it as --db";

/// `--engine`, `--db`, `--blastn`.
#[derive(clap::Args, Debug, Clone, PartialEq, Default)]
pub struct EngineArgs {
    /// Search engine: rgblastn, in process on the store and selection (megablast, blastn and
    /// blastn-short; see rgblastn -help), or stock blastn on a database you built (--db)
    #[arg(long, value_enum, default_value_t = EngineArg::Rgblastn)]
    pub engine: EngineArg,
    /// --engine blastn: the database to search, as built by 'rgmakeblastdb build' or 'view'
    #[arg(long, value_name = "BASE")]
    pub db: Option<PathBuf>,
    /// --engine blastn: the blastn binary (default: $GTARS_RGBLAST_BLASTN, else blastn on PATH)
    #[arg(long, value_name = "PATH")]
    pub blastn: Option<PathBuf>,
}

impl EngineArgs {
    /// The engine, the store and the selection: rgblastn needs `--store` and a selection;
    /// blastn needs `--db` and takes no store or selection (`--all` alone is accepted).
    pub fn target(
        &self,
        store: Option<&str>,
        sel: &SelectionArgs,
    ) -> Result<(EngineChoice, String, Selection), Error> {
        match self.engine {
            EngineArg::Rgblastn => {
                if self.db.is_some() || self.blastn.is_some() {
                    return Err(invalid(
                        "--db and --blastn go with --engine blastn; rgblastn searches --store",
                    ));
                }
                let store = match store {
                    Some(s) if !s.is_empty() => s.to_string(),
                    _ => {
                        return Err(invalid(
                            "--store PATH|URL is required (or --engine blastn --db BASE)",
                        ));
                    }
                };
                Ok((EngineChoice::Rgblastn, store, Selection::try_from(sel)?))
            }
            EngineArg::Blastn => {
                let Some(db) = self.db.clone() else {
                    return Err(invalid(
                        "--engine blastn needs --db BASE, a database built with 'rgmakeblastdb build' or 'view'",
                    ));
                };
                let all_only = SelectionArgs {
                    all: true,
                    ..Default::default()
                };
                if store.is_some() || !(sel.is_empty() || *sel == all_only) {
                    return Err(invalid(BLASTN_SEARCHES_DB));
                }
                let engine = EngineChoice::Blastn {
                    db,
                    blastn: self.blastn.clone(),
                };
                Ok((engine, String::new(), Selection::all()))
            }
        }
    }
}

pub(crate) fn check_store(store: &str) -> Result<(), Error> {
    if store.is_empty() {
        return Err(invalid("--store PATH|URL is required"));
    }
    Ok(())
}

/// `rgblast search`
#[derive(clap::Args, Debug, Clone, PartialEq, Default)]
pub struct SearchArgs {
    /// Local store directory, or an http(s) store URL (rgblastn)
    #[arg(long, value_name = "PATH|URL")]
    pub store: Option<String>,
    #[command(flatten)]
    pub selection: SelectionArgs,
    #[command(flatten)]
    pub engine: EngineArgs,
    #[command(flatten)]
    pub fetch: FetchFlags,
    /// Leave -outfmt 6/7 as given (do not append qlen score sallseqid salltitles staxids)
    #[arg(long)]
    pub no_extra_fields: bool,
    /// Write the run metadata (JSON) here; with -out FILE it goes to FILE.meta.json by default
    #[arg(long, value_name = "FILE")]
    pub meta: Option<PathBuf>,
    /// Only warnings and the engine's own messages on stderr
    #[arg(long)]
    pub quiet: bool,
    /// blastn arguments (after '--', or from the first single-dash flag)
    #[arg(
        trailing_var_arg = true,
        allow_hyphen_values = true,
        value_name = "BLASTN_ARGS"
    )]
    pub blastn_args: Vec<String>,
}

impl TryFrom<&SearchArgs> for SearchRequest {
    type Error = Error;

    fn try_from(a: &SearchArgs) -> Result<Self, Error> {
        let (engine, store, selection) = a.engine.target(a.store.as_deref(), &a.selection)?;
        a.fetch.check(&store)?;
        let meta = a.meta.clone().or_else(|| {
            value(&a.blastn_args, "-out")
                .filter(|o| *o != "-")
                .map(|o| PathBuf::from(format!("{o}.meta.json")))
        });
        Ok(SearchRequest {
            store,
            selection,
            engine,
            cache_dir: a.fetch.cache_dir.clone(),
            fetch: a.fetch.remote(),
            extra_fields: !a.no_extra_fields,
            args: a.blastn_args.clone(),
            output: EngineOutput::Inherit,
            meta,
        })
    }
}
