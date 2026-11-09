//! `SelectionArgs`: the selection flags shared by `rgmakeblastdb build|view|nstats`,
//! `rgblast` and `rgblastn`.

use crate::error::Error;
use crate::plan::ResolvedSelection;
use crate::progress::ProgressEvent;
use crate::select::{
    CollectionSet, EXCLUDE_MASKED_FRACTION, Flags, NO_MODE, ResolveOptions, Selection,
    read_list_file, strip_sq,
};
use crate::store::OpenedStore;
use std::path::PathBuf;

fn parse_fraction(s: &str) -> Result<f64, String> {
    let f: f64 = s.parse().map_err(|_| format!("'{s}' is not a number"))?;
    if f > 0.0 && f <= 1.0 {
        Ok(f)
    } else {
        Err("must be greater than 0 and at most 1".into())
    }
}

/// Which collections, and optionally which sequences. See the README's "Choosing
/// collections" for the grammar and the conflict table.
#[derive(clap::Args, Debug, Clone, PartialEq, Default)]
pub struct SelectionArgs {
    /// A collection: its digest, NS:ALIAS, or an alias unique across namespaces (repeatable;
    /// a value may start with '-')
    #[arg(
        long,
        value_name = "SEL",
        allow_hyphen_values = true,
        conflicts_with_all = ["namespace", "taxon", "all"]
    )]
    pub collection: Vec<String>,
    /// A file of collection selectors, one per line ('-' = stdin; '#' comments; first field)
    #[arg(long, value_name = "FILE", conflicts_with_all = ["namespace", "taxon", "all"])]
    pub collections_file: Vec<PathBuf>,
    /// Collections with an alias in this namespace (repeatable, OR'd; AND'd with --taxon)
    #[arg(long, value_name = "NS", conflicts_with = "all")]
    pub namespace: Vec<String>,
    /// Collections of this taxon: an NCBI taxid, or a species name from FHR, matched ignoring
    /// case with '_' as a space (repeatable, OR'd)
    #[arg(long, value_name = "TAXID|NAME", conflicts_with = "all")]
    pub taxon: Vec<String>,
    /// Every collection, plus sequences in no collection
    #[arg(long)]
    pub all: bool,
    /// Leave out collections whose N fraction is F or more (0 < F <= 1); listing one
    /// explicitly is an error
    #[arg(long, value_name = "F", value_parser = parse_fraction, conflicts_with = "exclude_masked")]
    pub max_n_fraction: Option<f64>,
    /// Same as --max-n-fraction 0.25 (leaves out hard-masked assemblies)
    #[arg(long)]
    pub exclude_masked: bool,
    /// Only these sequences (digests, 'SQ.' optional, one per line; '-' = stdin), with every
    /// selected collection's defline
    #[arg(long, value_name = "FILE")]
    pub sequences_file: Option<PathBuf>,
}

impl SelectionArgs {
    /// True when no selection flag is given (view then uses its default).
    pub fn is_empty(&self) -> bool {
        *self == SelectionArgs::default()
    }

    pub fn flags(&self) -> Flags {
        Flags {
            collection: !self.collection.is_empty(),
            collections_file: !self.collections_file.is_empty(),
            namespace: !self.namespace.is_empty(),
            taxon: !self.taxon.is_empty(),
            all: self.all,
            max_n_fraction: self.max_n_fraction,
            exclude_masked: self.exclude_masked,
        }
    }

    /// Every rule violation, as clap words them (empty when the flags are valid).
    pub fn violations(&self) -> Vec<String> {
        let f = self.flags();
        let mut v = f.violations();
        if !f.any_mode() {
            v.push(NO_MODE.to_string());
        }
        v
    }

    /// Convert and resolve against a store (default options).
    pub fn resolve(&self, store: &mut OpenedStore) -> Result<ResolvedSelection, Error> {
        self.resolve_with(store, &ResolveOptions::default(), &mut |_| true)
    }

    /// Convert and resolve against a store.
    pub fn resolve_with(
        &self,
        store: &mut OpenedStore,
        opts: &ResolveOptions,
        progress: &mut dyn FnMut(&ProgressEvent) -> bool,
    ) -> Result<ResolvedSelection, Error> {
        let sel = Selection::try_from(self)?;
        Ok(crate::api::resolve_in_store(store, &sel, opts, progress)?.0)
    }
}

/// The authority for the selection rules: reads the files, strips `SQ.`, re-checks every rule
/// clap enforces and reports all violations in one `InvalidOptions`.
impl TryFrom<&SelectionArgs> for Selection {
    type Error = Error;

    fn try_from(a: &SelectionArgs) -> Result<Self, Error> {
        let v = a.violations();
        if !v.is_empty() {
            return Err(Error::InvalidOptions(v.join("; ")));
        }
        let set = if a.all {
            CollectionSet::All
        } else if !a.collection.is_empty() || !a.collections_file.is_empty() {
            let mut sel = a.collection.clone();
            for f in &a.collections_file {
                let entries = read_list_file(f)?;
                if entries.is_empty() {
                    return Err(Error::InvalidOptions(format!(
                        "--collections-file {} lists no collections",
                        f.display()
                    )));
                }
                sel.extend(entries);
            }
            CollectionSet::Explicit(sel)
        } else {
            CollectionSet::Attributes {
                namespaces: a.namespace.clone(),
                taxa: a.taxon.clone(),
            }
        };
        let sequences = match &a.sequences_file {
            Some(f) => {
                let v: Vec<String> = read_list_file(f)?.iter().map(|s| strip_sq(s)).collect();
                if v.is_empty() {
                    return Err(Error::InvalidOptions(format!(
                        "--sequences-file {} lists no sequences",
                        f.display()
                    )));
                }
                Some(v)
            }
            None => None,
        };
        let s = Selection {
            set,
            max_n_fraction: if a.exclude_masked {
                Some(EXCLUDE_MASKED_FRACTION)
            } else {
                a.max_n_fraction
            },
            sequences,
        };
        s.validate()?;
        Ok(s)
    }
}
