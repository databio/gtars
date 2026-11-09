//! `rgblast origin` and `rgblast reps`.

use super::report::{RankArgs, ReportOutputArgs, percent};
use super::search::{EngineArgs, check_store};
use super::{FetchFlags, SelectionArgs, invalid};
use crate::api::{OriginRequest, RepsRequest};
use crate::error::Error;
use crate::origin::{ModeChoice, OriginOptions};
use crate::run::EngineChoice;
use crate::run::passthrough::{check_origin, value};
use crate::select::Selection;
use std::path::PathBuf;

pub const ORIGIN_AFTER_HELP: &str = "\
Staged mode (the default when predicted to pay off) searches one representative assembly per
species, then the rest of each candidate species. Every stage uses the full selection's
statistics, so E-values equal a full search. Staged and --exhaustive results can differ when:
  1. the best-matching sequence is in no candidate species' representative (a decoy, alt or
     patch contig, or a region only another build has): with no stage 1 hit the query is
     'none' with flag reps-only (rerun it with --exhaustive);
  2. the true species' representative scores below the --fanout-percent cut;
  3. -max_target_seqs truncates the hit list differently (flag truncated).
The 'searched' column and the reps-only flag show each case per query. Staged mode needs the
default engine, rgblastn; with --engine blastn, origin runs one search of --db.";

/// `rgblast origin`
#[derive(clap::Args, Debug, Clone, PartialEq)]
pub struct OriginArgs {
    /// Local store directory, or an http(s) store URL (rgblastn)
    #[arg(long, value_name = "PATH|URL")]
    pub store: Option<String>,
    #[command(flatten)]
    pub selection: SelectionArgs,
    /// Query FASTA ('-' = stdin); or pass -query after '--'
    #[arg(long, value_name = "FILE")]
    pub query: Option<PathBuf>,
    /// One search of the whole selection
    #[arg(long)]
    pub exhaustive: bool,
    /// Staged even when the cost prediction says it does not pay off
    #[arg(long)]
    pub staged: bool,
    /// Fan out to species whose best stage 1 score is within F% of the query's best
    #[arg(long, value_name = "F", default_value_t = 20.0)]
    pub fanout_percent: f64,
    /// The representative of its species (digest or NS:ALIAS; repeatable)
    #[arg(long = "rep", value_name = "COLLECTION")]
    pub rep: Vec<String>,
    /// Alias namespace that marks preferred representatives
    #[arg(long, value_name = "NS", default_value = "reference")]
    pub rep_namespace: String,
    #[command(flatten)]
    pub ranking: RankArgs,
    #[command(flatten)]
    pub out: ReportOutputArgs,
    /// Keep the stage hits, queries and origin.meta.json in DIR
    #[arg(long, value_name = "DIR")]
    pub keep_hits: Option<PathBuf>,
    #[command(flatten)]
    pub engine: EngineArgs,
    #[command(flatten)]
    pub fetch: FetchFlags,
    /// Only warnings and the summary on stderr
    #[arg(long)]
    pub quiet: bool,
    /// blastn arguments for every stage (after '--', or from the first single-dash flag)
    #[arg(
        trailing_var_arg = true,
        allow_hyphen_values = true,
        value_name = "BLASTN_ARGS"
    )]
    pub blastn_args: Vec<String>,
}

impl TryFrom<&OriginArgs> for OriginRequest {
    type Error = Error;

    fn try_from(a: &OriginArgs) -> Result<Self, Error> {
        let (engine, store, selection) = a.engine.target(a.store.as_deref(), &a.selection)?;
        a.fetch.check(&store)?;
        if matches!(engine, EngineChoice::Blastn { .. }) && (a.staged || !a.rep.is_empty()) {
            return Err(invalid(
                "staged origin (--staged, --rep) needs the default engine, rgblastn; with --engine blastn, origin runs one search of --db",
            ));
        }
        if a.exhaustive && a.staged {
            return Err(invalid(
                "the argument '--exhaustive' cannot be used with '--staged'",
            ));
        }
        percent("--fanout-percent", a.fanout_percent)?;
        let mut args = a.blastn_args.clone();
        let pass_query = value(&args, "-query").map(PathBuf::from);
        let query = match (&a.query, pass_query) {
            (Some(_), Some(_)) => {
                return Err(invalid(
                    "give the query once: --query FILE or -query FILE, not both",
                ));
            }
            (None, None) => {
                return Err(invalid(
                    "origin needs a query: --query FILE (or -query FILE)",
                ));
            }
            (Some(q), None) => q.clone(),
            (None, Some(q)) => {
                let i = args.iter().position(|x| x == "-query").unwrap();
                args.drain(i..i + 2);
                q
            }
        };
        check_origin(&args)?;
        let blastn = matches!(engine, EngineChoice::Blastn { .. });
        Ok(OriginRequest {
            store,
            selection,
            engine,
            cache_dir: a.fetch.cache_dir.clone(),
            fetch: a.fetch.remote(),
            query,
            args,
            reps: a.rep.clone(),
            options: OriginOptions {
                mode: if a.exhaustive || blastn {
                    ModeChoice::Exhaustive
                } else if a.staged {
                    ModeChoice::Staged
                } else {
                    ModeChoice::Auto
                },
                fanout_percent: a.fanout_percent,
                reps: Vec::new(),
                rep_namespace: a.rep_namespace.clone(),
                rank: a.ranking.options()?,
            },
            taxdump: a.ranking.taxdump.clone(),
            keep_hits: a.keep_hits.clone(),
        })
    }
}

/// `rgblast reps`
#[derive(clap::Args, Debug, Clone, PartialEq)]
pub struct RepsArgs {
    /// Local store directory, or an http(s) store URL
    #[arg(long, value_name = "PATH|URL")]
    pub store: String,
    #[command(flatten)]
    pub selection: SelectionArgs,
    /// NCBI taxdump directory (for --group-rank and taxon names)
    #[arg(long, value_name = "DIR")]
    pub taxdump: Option<PathBuf>,
    /// Group by the ancestor at this rank (needs a taxdump)
    #[arg(long, value_name = "RANK")]
    pub group_rank: Option<String>,
    /// The representative of its species (digest or NS:ALIAS; repeatable)
    #[arg(long = "rep", value_name = "COLLECTION")]
    pub rep: Vec<String>,
    /// Alias namespace that marks preferred representatives
    #[arg(long, value_name = "NS", default_value = "reference")]
    pub rep_namespace: String,
    #[command(flatten)]
    pub fetch: FetchFlags,
}

impl TryFrom<&RepsArgs> for RepsRequest {
    type Error = Error;

    fn try_from(a: &RepsArgs) -> Result<Self, Error> {
        check_store(&a.store)?;
        a.fetch.check(&a.store)?;
        let ranking = RankArgs {
            taxdump: a.taxdump.clone(),
            group_rank: a.group_rank.clone(),
            ..Default::default()
        };
        ranking.options()?;
        Ok(RepsRequest {
            store: a.store.clone(),
            selection: Selection::try_from(&a.selection)?,
            taxdump: a.taxdump.clone(),
            group_rank: a.group_rank.clone(),
            reps: a.rep.clone(),
            rep_namespace: a.rep_namespace.clone(),
            cache_dir: a.fetch.cache_dir.clone(),
            fetch: a.fetch.remote(),
        })
    }
}
