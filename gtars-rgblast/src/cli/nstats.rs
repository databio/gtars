//! `NstatsArgs`: `rgmakeblastdb nstats`.

use super::{FetchFlags, SelectionArgs, nstats_options};
use crate::api::NstatsRequest;
use crate::error::Error;
use crate::select::Selection;

/// `rgmakeblastdb nstats`
#[derive(clap::Args, Debug, Clone, PartialEq, Default)]
pub struct NstatsArgs {
    /// Local store directory, or an http(s) store URL
    #[arg(long, value_name = "PATH|URL")]
    pub store: String,
    #[command(flatten)]
    pub selection: SelectionArgs,
    /// Per-collection taxids that override FHR (for the species column)
    #[arg(long, value_name = "FILE")]
    pub taxid_map: Option<std::path::PathBuf>,
    /// Counting threads (default: all cores)
    #[arg(long, value_name = "N")]
    pub threads: Option<usize>,
    #[command(flatten)]
    pub fetch: FetchFlags,
}

impl TryFrom<&NstatsArgs> for NstatsRequest {
    type Error = Error;

    fn try_from(a: &NstatsArgs) -> Result<Self, Error> {
        let selection = Selection::try_from(&a.selection)?;
        a.fetch.check(&a.store)?;
        Ok(NstatsRequest {
            store: a.store.clone(),
            cache_dir: a.fetch.cache_dir.clone(),
            fetch: a.fetch.remote(),
            selection,
            nstats: nstats_options(a.fetch.cache_dir.as_ref()),
            taxid_map: a.taxid_map.clone(),
            threads: a.threads.unwrap_or(0),
        })
    }
}
