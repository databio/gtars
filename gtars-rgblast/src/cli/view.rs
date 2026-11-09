//! `ViewArgs`: `rgmakeblastdb view`.

use super::{SelectionArgs, invalid};
use crate::blastdb::sidecars::ViewKind;
use crate::error::Error;
use crate::select::Selection;
use crate::views::ViewRequest;
use std::path::PathBuf;

/// `rgmakeblastdb view`
#[derive(clap::Args, Debug, Clone, PartialEq, Default)]
#[command(group(clap::ArgGroup::new("view_kind").required(true).args(["per_collection", "name"])))]
pub struct ViewArgs {
    /// Database base path (built by rgmakeblastdb)
    #[arg(long, value_name = "BASE")]
    pub db: PathBuf,
    /// One view per collection (selection flags limit which)
    #[arg(long)]
    pub per_collection: bool,
    /// One view of the selection, with this name (needs selection flags)
    #[arg(long, value_name = "NAME")]
    pub name: Option<String>,
    /// Where to write the views (default: BASE.views/)
    #[arg(long, value_name = "DIR")]
    pub out_dir: Option<PathBuf>,
    /// Title prefix (per-collection) or title (--name); default the database title
    #[arg(long, value_name = "TEXT")]
    pub title: Option<String>,
    /// Also write the whole database's sequence count and length as statistics (E-values on
    /// the scale of the whole database; NSEQ and LENGTH stay the view's own)
    #[arg(long)]
    pub full_db_stats: bool,
    #[command(flatten)]
    pub selection: SelectionArgs,
}

impl TryFrom<&ViewArgs> for ViewRequest {
    type Error = Error;

    fn try_from(a: &ViewArgs) -> Result<Self, Error> {
        let selection = if a.selection.is_empty() {
            if a.name.is_some() {
                return Err(invalid(
                    "--name requires a selection (--collection, --collections-file, --namespace, --taxon or --all)",
                ));
            }
            None
        } else {
            Some(Selection::try_from(&a.selection)?)
        };
        let kind = if a.per_collection {
            ViewKind::PerCollection
        } else {
            ViewKind::Named
        };
        let r = ViewRequest {
            db: a.db.clone(),
            kind,
            name: a.name.clone(),
            selection,
            out_dir: a.out_dir.clone(),
            title: a.title.clone(),
            full_db_stats: a.full_db_stats,
        };
        r.validate()?;
        Ok(r)
    }
}
