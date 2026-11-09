//! Which collections, and optionally which sequences, go into a database or a search.
//!
//! [`Selection`] is plain serde data: the CLI's `SelectionArgs` converts into it, and bindings
//! may build it from keyword arguments. [`resolve`] turns it into a
//! [`ResolvedSelection`](crate::plan::ResolvedSelection) against a store (or a built
//! database's sidecars). This module imports only `plan`, `store` and `convert`, never
//! `blastdb`, so the rgblastn engine can use it.

pub mod nstats;
pub mod resolve;

use crate::error::Error;
pub use crate::plan::{CollectionInfo, Excluded, ResolvedSelection};
pub use resolve::{Catalog, ResolveOptions, StoreCatalog, resolve};
use serde::{Deserialize, Serialize};
use std::io::Read;
use std::path::Path;

/// `--exclude-masked` means exactly this `--max-n-fraction`.
pub const EXCLUDE_MASKED_FRACTION: f64 = 0.25;

/// The collection part of a selection; the modes exclude each other.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum CollectionSet {
    /// Selectors from `--collection` / `--collections-file` (digest, `NS:ALIAS` or bare alias).
    Explicit(Vec<String>),
    /// `--namespace` (values OR'd) AND `--taxon` (TAXID or species name, values OR'd).
    Attributes {
        namespaces: Vec<String>,
        taxa: Vec<String>,
    },
    /// `--all`: every collection, plus sequences in no collection.
    All,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Selection {
    pub set: CollectionSet,
    /// `--max-n-fraction` / `--exclude-masked`: collections with N fraction ≥ this are left out
    /// (attribute and all modes) or refused (explicit mode).
    pub max_n_fraction: Option<f64>,
    /// `--sequences-file`: subject digests (no `SQ.`).
    pub sequences: Option<Vec<String>>,
}

impl Default for Selection {
    fn default() -> Self {
        Selection {
            set: CollectionSet::All,
            max_n_fraction: None,
            sequences: None,
        }
    }
}

impl Selection {
    pub fn all() -> Self {
        Selection::default()
    }

    pub fn explicit<S: AsRef<str>>(selectors: &[S]) -> Self {
        Selection {
            set: CollectionSet::Explicit(
                selectors.iter().map(|s| s.as_ref().to_string()).collect(),
            ),
            ..Default::default()
        }
    }

    pub fn attributes<S: AsRef<str>>(namespaces: &[S], taxa: &[S]) -> Self {
        Selection {
            set: CollectionSet::Attributes {
                namespaces: namespaces.iter().map(|s| s.as_ref().to_string()).collect(),
                taxa: taxa.iter().map(|s| s.as_ref().to_string()).collect(),
            },
            ..Default::default()
        }
    }

    pub fn with_max_n_fraction(mut self, f: f64) -> Self {
        self.max_n_fraction = Some(f);
        self
    }

    pub fn with_sequences<S: AsRef<str>>(mut self, digests: &[S]) -> Self {
        self.sequences = Some(digests.iter().map(|s| strip_sq(s.as_ref())).collect());
        self
    }

    /// The rule checks `TryFrom<&SelectionArgs>` runs, for a `Selection` built directly.
    pub fn validate(&self) -> Result<(), Error> {
        let mut v = Vec::new();
        let flags = match &self.set {
            CollectionSet::Explicit(sel) => {
                if sel.is_empty() {
                    v.push("--collection: the list of collections is empty".to_string());
                }
                if sel.iter().any(|s| s.trim().is_empty()) {
                    v.push("--collection: a collection selector is empty".to_string());
                }
                Flags {
                    collection: true,
                    ..Default::default()
                }
            }
            CollectionSet::Attributes { namespaces, taxa } => {
                if namespaces.iter().any(|s| s.trim().is_empty()) {
                    v.push("--namespace: a namespace is empty".to_string());
                }
                if taxa.iter().any(|s| s.trim().is_empty()) {
                    v.push("--taxon: a taxon is empty".to_string());
                }
                Flags {
                    namespace: !namespaces.is_empty(),
                    taxon: !taxa.is_empty(),
                    ..Default::default()
                }
            }
            CollectionSet::All => Flags {
                all: true,
                ..Default::default()
            },
        };
        let flags = Flags {
            max_n_fraction: self.max_n_fraction,
            ..flags
        };
        v.extend(flags.violations());
        if let Some(seqs) = &self.sequences {
            if seqs
                .iter()
                .any(|s| s.is_empty() || s.contains(char::is_whitespace))
            {
                v.push("--sequences-file: a sequence digest is empty or holds whitespace".into());
            }
            if seqs.iter().any(|s| s.starts_with("SQ.")) {
                v.push("--sequences-file: digests are stored without 'SQ.'".into());
            }
        }
        if v.is_empty() {
            Ok(())
        } else {
            Err(Error::InvalidOptions(v.join("; ")))
        }
    }
}

/// Which selection flags were given, for the rule table (shared by `Selection::validate` and
/// the CLI conversions, so both report the same wording).
#[derive(Debug, Clone, Default)]
pub struct Flags {
    /// `--collection` or `--collections-file`.
    pub collection: bool,
    pub collections_file: bool,
    pub namespace: bool,
    pub taxon: bool,
    pub all: bool,
    pub max_n_fraction: Option<f64>,
    pub exclude_masked: bool,
}

/// clap's spelling of each flag in its messages.
pub mod flag {
    pub const COLLECTION: &str = "--collection <SEL>";
    pub const COLLECTIONS_FILE: &str = "--collections-file <FILE>";
    pub const NAMESPACE: &str = "--namespace <NS>";
    pub const TAXON: &str = "--taxon <TAXID|NAME>";
    pub const ALL: &str = "--all";
    pub const MAX_N_FRACTION: &str = "--max-n-fraction <F>";
    pub const EXCLUDE_MASKED: &str = "--exclude-masked";
}

/// clap's wording for a conflict.
pub fn conflict(a: &str, b: &str) -> String {
    format!("the argument '{a}' cannot be used with '{b}'")
}

/// The message when no mode is given.
pub const NO_MODE: &str =
    "choose collections with --collection, --collections-file, --namespace, --taxon or --all";

impl Flags {
    pub fn any_mode(&self) -> bool {
        self.collection || self.collections_file || self.namespace || self.taxon || self.all
    }

    /// Every rule violation, in a fixed order, with clap's wording.
    pub fn violations(&self) -> Vec<String> {
        use flag::*;
        let mut v = Vec::new();
        let explicit: Vec<&str> = [
            (self.collection, COLLECTION),
            (self.collections_file, COLLECTIONS_FILE),
        ]
        .iter()
        .filter(|x| x.0)
        .map(|x| x.1)
        .collect();
        let attrs: Vec<&str> = [(self.namespace, NAMESPACE), (self.taxon, TAXON)]
            .iter()
            .filter(|x| x.0)
            .map(|x| x.1)
            .collect();
        for e in &explicit {
            for a in &attrs {
                v.push(conflict(e, a));
            }
            if self.all {
                v.push(conflict(e, ALL));
            }
        }
        if self.all {
            for a in &attrs {
                v.push(conflict(a, ALL));
            }
        }
        if self.max_n_fraction.is_some() && self.exclude_masked {
            v.push(conflict(MAX_N_FRACTION, EXCLUDE_MASKED));
        }
        if let Some(f) = self.max_n_fraction
            && !(f > 0.0 && f <= 1.0)
        {
            v.push(format!(
                "invalid value '{f}' for '{MAX_N_FRACTION}': must be greater than 0 and at most 1"
            ));
        }
        v
    }
}

/// Strip an optional `SQ.` prefix.
pub fn strip_sq(s: &str) -> String {
    s.strip_prefix("SQ.").unwrap_or(s).to_string()
}

/// Read a list file (`-` = stdin): one entry per line, the first whitespace-separated field;
/// `#` comments and blank lines are skipped.
pub fn read_list_file(path: &Path) -> Result<Vec<String>, Error> {
    let text = if path == Path::new("-") {
        let mut s = String::new();
        std::io::stdin()
            .read_to_string(&mut s)
            .map_err(|e| Error::io(path, e))?;
        s
    } else {
        std::fs::read_to_string(path).map_err(|e| Error::io(path, e))?
    };
    Ok(text
        .lines()
        .filter_map(|l| {
            let t = l.trim();
            if t.is_empty() || t.starts_with('#') {
                None
            } else {
                t.split_whitespace().next().map(String::from)
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_direct_selections() {
        assert!(Selection::all().validate().is_ok());
        assert!(Selection::explicit(&["a"]).validate().is_ok());
        assert!(Selection::explicit::<&str>(&[]).validate().is_err());
        assert!(Selection::attributes(&["refseq"], &[]).validate().is_ok());
        assert!(Selection::all().with_max_n_fraction(1.0).validate().is_ok());
        for bad in [0.0, -0.1, 1.5, f64::NAN] {
            let e = Selection::all()
                .with_max_n_fraction(bad)
                .validate()
                .unwrap_err();
            assert!(e.to_string().contains("--max-n-fraction"), "{e}");
        }
        let s = Selection::all().with_sequences(&["SQ.abc", "def"]);
        assert_eq!(
            s.sequences,
            Some(vec!["abc".to_string(), "def".to_string()])
        );
        assert!(s.validate().is_ok());
    }

    #[test]
    fn list_files() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("l.txt");
        std::fs::write(
            &p,
            "# header\n\n a  extra\nname:b\t#x\n  # indented comment\nc\n",
        )
        .unwrap();
        assert_eq!(read_list_file(&p).unwrap(), vec!["a", "name:b", "c"]);
    }

    #[test]
    fn flag_rules() {
        let f = Flags {
            all: true,
            collection: true,
            namespace: true,
            max_n_fraction: Some(0.3),
            exclude_masked: true,
            ..Default::default()
        };
        assert_eq!(
            f.violations(),
            vec![
                conflict(flag::COLLECTION, flag::NAMESPACE),
                conflict(flag::COLLECTION, flag::ALL),
                conflict(flag::NAMESPACE, flag::ALL),
                conflict(flag::MAX_N_FRACTION, flag::EXCLUDE_MASKED),
            ]
        );
    }
}
