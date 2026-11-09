//! `rgblast report`, and the ranking and report-output flags shared with `origin`.

use super::invalid;
use crate::api::ReportRequest;
use crate::error::Error;
use crate::report::parse::DEFAULT_FIELDS;
use crate::report::rank::RankOptions;
use crate::report::taxonomy::taxdump_dir;
use std::path::PathBuf;

/// `--taxdump`, `--group-rank`, `--tie-percent` and the HSP filters.
#[derive(clap::Args, Debug, Clone, PartialEq, Default)]
pub struct RankArgs {
    /// NCBI taxdump directory (nodes.dmp, names.dmp, merged.dmp) for LCA calls and ranks
    /// (default: $GTARS_RGBLAST_TAXDUMP)
    #[arg(long, value_name = "DIR")]
    pub taxdump: Option<PathBuf>,
    /// Group genomes by their ancestor at this rank (e.g. species folds subspecies; needs a taxdump)
    #[arg(long, value_name = "RANK")]
    pub group_rank: Option<String>,
    /// Genomes within P% of the top score tie (0 = exact ties only)
    #[arg(long, value_name = "P", default_value_t = 0.0)]
    pub tie_percent: f64,
    /// Drop HSPs below this percent identity
    #[arg(long, value_name = "X")]
    pub min_pident: Option<f64>,
    /// Drop HSPs covering less than this percent of the query
    #[arg(long, value_name = "X")]
    pub min_qcov: Option<f64>,
    /// Drop HSPs above this E-value
    #[arg(long, value_name = "E")]
    pub max_evalue: Option<f64>,
}

pub(crate) fn percent(flag: &str, p: f64) -> Result<(), Error> {
    if (0.0..100.0).contains(&p) {
        Ok(())
    } else {
        Err(invalid(format!("{flag} must be at least 0 and below 100")))
    }
}

impl RankArgs {
    pub fn options(&self) -> Result<RankOptions, Error> {
        percent("--tie-percent", self.tie_percent)?;
        if self.group_rank.is_some() && taxdump_dir(self.taxdump.as_deref()).is_none() {
            return Err(invalid(
                "--group-rank needs a taxonomy: give --taxdump DIR or set GTARS_RGBLAST_TAXDUMP",
            ));
        }
        Ok(RankOptions {
            tie_percent: self.tie_percent,
            min_pident: self.min_pident,
            min_qcov: self.min_qcov,
            max_evalue: self.max_evalue,
            group_rank: self.group_rank.clone(),
        })
    }
}

#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FormatArg {
    #[default]
    Tsv,
    Json,
}

/// Where the report goes.
#[derive(clap::Args, Debug, Clone, PartialEq, Default)]
pub struct ReportOutputArgs {
    /// Summary format (one line per query for tsv)
    #[arg(long, value_enum, default_value_t = FormatArg::Tsv)]
    pub format: FormatArg,
    /// Write one row per (query, genome), grouped by taxon
    #[arg(long, value_name = "FILE")]
    pub detail: Option<PathBuf>,
    /// Keep the best N genomes per taxon in --detail
    #[arg(long, value_name = "N")]
    pub detail_top: Option<usize>,
    /// Write the summary here instead of stdout
    #[arg(long = "output", value_name = "FILE")]
    pub output: Option<PathBuf>,
}

/// `rgblast report`
#[derive(clap::Args, Debug, Clone, PartialEq)]
pub struct ReportArgs {
    /// Tabular hits (-outfmt 6 or 7); '-' is stdin
    #[arg(value_name = "HITS", required = true)]
    pub hits: Vec<PathBuf>,
    /// The database the hits came from (reads its sidecars; default: the one HITS.meta.json
    /// names, else the store and collections it names)
    #[arg(long, value_name = "BASE")]
    pub db: Option<PathBuf>,
    /// The run metadata of 'rgblast search' (default: HITS.meta.json when it exists)
    #[arg(long, value_name = "FILE")]
    pub meta: Option<PathBuf>,
    /// The query FASTA: lists queries with no hits, in input order
    #[arg(long, value_name = "FILE")]
    pub query: Option<PathBuf>,
    /// Columns of -outfmt 6 input (blastn syntax; outfmt 7 names its own)
    #[arg(long, value_name = "LIST", default_value = DEFAULT_FIELDS)]
    pub fields: String,
    #[command(flatten)]
    pub ranking: RankArgs,
    /// The -max_target_seqs of the run (for the 'truncated' flag; default from the meta, else 500)
    #[arg(long, value_name = "N")]
    pub max_target_seqs: Option<usize>,
    #[command(flatten)]
    pub out: ReportOutputArgs,
}

impl TryFrom<&ReportArgs> for ReportRequest {
    type Error = Error;

    fn try_from(a: &ReportArgs) -> Result<Self, Error> {
        if a.hits.is_empty() {
            return Err(invalid(
                "the following required arguments were not provided: <HITS>...",
            ));
        }
        Ok(ReportRequest {
            hits: a.hits.clone(),
            db: a.db.clone(),
            meta: a.meta.clone(),
            query: a.query.clone(),
            fields: a.fields.clone(),
            taxdump: a.ranking.taxdump.clone(),
            rank: a.ranking.options()?,
            max_target_seqs: a.max_target_seqs,
        })
    }
}
