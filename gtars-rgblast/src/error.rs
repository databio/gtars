//! Typed errors. Every public function returns [`Error`]; bindings map
//! [`Error::kind`] to their native exception types.

use serde::{Deserialize, Serialize};
use std::fmt;
use std::path::PathBuf;

/// What an id was made from, for collision messages.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum IdOrigin {
    /// A sequence name in a collection (`display` is the collection's first alias, else its digest).
    Name {
        collection: String,
        display: String,
        name: String,
    },
    /// A sequence digest (`SQ.<digest>`), or one of its sequence aliases.
    Sequence {
        digest: String,
        alias: Option<String>,
    },
}

impl fmt::Display for IdOrigin {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IdOrigin::Name { display, name, .. } => {
                write!(f, "name '{name}' in collection {display}")
            }
            IdOrigin::Sequence {
                digest,
                alias: None,
            } => write!(f, "sequence SQ.{digest}"),
            IdOrigin::Sequence {
                digest,
                alias: Some(a),
            } => {
                write!(f, "sequence alias '{a}' of SQ.{digest}")
            }
        }
    }
}

/// Fieldless error category, for bindings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ErrorKind {
    NotFound,
    InvalidInput,
    Conflict,
    Io,
    Network,
    Data,
    Cancelled,
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A store, database, collection namespace or taxon that is not there.
    #[error("{0}")]
    NotFound(String),

    #[error("{}", unknown_collection(key, store, namespaces))]
    UnknownCollection {
        key: String,
        store: String,
        namespaces: Vec<String>,
    },

    #[error("{0}")]
    InvalidOptions(String),

    #[error("{}", id_collision(id, first, second, *sanitized))]
    IdCollision {
        id: String,
        first: IdOrigin,
        second: IdOrigin,
        sanitized: bool,
    },

    #[error("{}", oversize(*count, examples))]
    Oversize {
        count: u64,
        examples: Vec<(String, u64)>,
    },

    #[error("{} already exists (a database is already at this path)", path.display())]
    OutputExists { path: PathBuf },

    #[error(
        "not enough free space in {}: need about {}, {} available",
        dir.display(), crate::progress::gb(*need), crate::progress::gb(*available)
    )]
    InsufficientSpace {
        dir: PathBuf,
        need: u64,
        available: u64,
    },

    #[error("sequence SQ.{digest}: {msg}")]
    Sequence { digest: String, msg: String },

    #[error("{}: {source}", path.display())]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("{what}: {msg}")]
    Fetch { what: String, msg: String },

    #[error("{}: LMDB error: {source}", path.display())]
    Lmdb {
        path: PathBuf,
        source: crate::blastdb::lmdb::LmdbError,
    },

    #[error("{0:#}")]
    Store(anyhow::Error),

    /// No usable engine (blastn not found, or `--engine rgblastn` not available yet).
    #[error("{hint}")]
    EngineNotFound { hint: String },

    /// The engine could not be run or is unsuitable (for example older than 2.10).
    #[error("{0}")]
    Engine(String),

    /// A query FASTA problem; `pos` is the 0-based letter position within query `id`.
    #[error("query {id}: {msg}")]
    Query { id: String, pos: usize, msg: String },

    #[error("cancelled")]
    Cancelled,
}

impl Error {
    pub fn kind(&self) -> ErrorKind {
        use Error::*;
        match self {
            NotFound(_) | UnknownCollection { .. } | EngineNotFound { .. } => ErrorKind::NotFound,
            InvalidOptions(_) => ErrorKind::InvalidInput,
            IdCollision { .. } | OutputExists { .. } => ErrorKind::Conflict,
            InsufficientSpace { .. } | Io { .. } | Lmdb { .. } => ErrorKind::Io,
            Fetch { .. } => ErrorKind::Network,
            Oversize { .. } | Sequence { .. } | Store(_) | Engine(_) | Query { .. } => {
                ErrorKind::Data
            }
            Cancelled => ErrorKind::Cancelled,
        }
    }

    /// A blastn flag rgblastn does not implement.
    pub fn unsupported(flag: &str) -> Self {
        Error::InvalidOptions(format!(
            "rgblastn: {flag} is not supported by rgblastn; build a database with rgmakeblastdb and use blastn"
        ))
    }

    /// A remote store file that is not there (or not readable).
    pub fn remote_not_found(url: &str, status: u16) -> Self {
        Error::NotFound(format!("{url}: not found or access denied (HTTP {status})"))
    }

    /// The selection matched nothing.
    pub fn empty_selection(excluded: usize, reason: &str) -> Self {
        Error::InvalidOptions(empty_selection(excluded, reason))
    }

    /// `Error::Io` with the path that failed.
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        Error::Io {
            path: path.into(),
            source,
        }
    }

    pub fn lmdb(path: impl Into<PathBuf>, source: crate::blastdb::lmdb::LmdbError) -> Self {
        Error::Lmdb {
            path: path.into(),
            source,
        }
    }

    /// A data error about one sequence; `digest` may be filled in later with [`Error::with_digest`].
    pub fn sequence(digest: &str, msg: impl Into<String>) -> Self {
        Error::Sequence {
            digest: digest.to_string(),
            msg: msg.into(),
        }
    }

    /// Fill in the digest of a `Sequence` error raised by code that did not know it.
    pub fn with_digest(self, d: &str) -> Self {
        match self {
            Error::Sequence { digest, msg } if digest.is_empty() => Error::Sequence {
                digest: d.to_string(),
                msg,
            },
            e => e,
        }
    }

    /// A one-line hint the CLI prints under the error, where one helps.
    pub fn hint(&self) -> Option<String> {
        match self {
            Error::OutputExists { .. } => Some("use --force to replace it".into()),
            Error::InsufficientSpace { .. } => Some(
                "free some space, choose another --out, or skip the check with --no-space-check"
                    .into(),
            ),
            Error::Oversize { .. } => Some(
                "use --on-oversize skip to leave them out (they are listed in BASE.skipped.tsv)"
                    .into(),
            ),
            _ => None,
        }
    }
}

impl From<anyhow::Error> for Error {
    fn from(e: anyhow::Error) -> Self {
        Error::Store(e)
    }
}

fn unknown_collection(key: &str, store: &str, namespaces: &[String]) -> String {
    let ns = if namespaces.is_empty() {
        "the store has no collection aliases".to_string()
    } else {
        format!("collection alias namespaces: {}", namespaces.join(", "))
    };
    format!("unknown collection '{key}' in {store}: not a collection digest or alias ({ns})")
}

pub(crate) fn unknown_namespace(ns: &str, namespaces: &[String]) -> String {
    let list = if namespaces.is_empty() {
        "the store has no collection aliases".to_string()
    } else {
        format!("collection alias namespaces: {}", namespaces.join(", "))
    };
    format!("unknown collection alias namespace '{ns}' ({list})")
}

pub(crate) fn unknown_taxon(taxon: &str, available: &[(i32, String, usize)]) -> String {
    if available.is_empty() {
        return format!("unknown taxon '{taxon}': no collection in this store has an FHR taxon");
    }
    let list: Vec<String> = available
        .iter()
        .map(|(t, name, n)| {
            format!(
                "{t} {name} ({n} collection{})",
                if *n == 1 { "" } else { "s" }
            )
        })
        .collect();
    format!(
        "unknown taxon '{taxon}'; taxa in this store: {}",
        list.join(", ")
    )
}

fn empty_selection(excluded: usize, reason: &str) -> String {
    if excluded > 0 {
        format!(
            "the selection is empty: {excluded} collection{} excluded ({reason})",
            if excluded == 1 { "" } else { "s" }
        )
    } else {
        format!("the selection is empty ({reason})")
    }
}

/// A fraction as given on the command line: `0.25`, `1`, `0.125`.
pub fn fmt_fraction(f: f64) -> String {
    let s = format!("{f}");
    s
}

pub(crate) fn ambiguous_alias(alias: &str, matches: &[(String, String)]) -> String {
    let m: Vec<String> = matches
        .iter()
        .map(|(ns, d)| format!("{ns}:{alias} -> {d}"))
        .collect();
    format!(
        "collection alias '{alias}' names more than one collection ({}); give the digest",
        m.join(", ")
    )
}

fn id_collision(id: &str, first: &IdOrigin, second: &IdOrigin, sanitized: bool) -> String {
    let what = match (first, second) {
        (
            IdOrigin::Name {
                collection: c1,
                display,
                name: n1,
            },
            IdOrigin::Name {
                collection: c2,
                name: n2,
                ..
            },
        ) if c1 == c2 => format!("names '{n1}' and '{n2}' in collection {display}"),
        (
            IdOrigin::Sequence {
                alias: Some(a1), ..
            },
            IdOrigin::Sequence {
                alias: Some(a2), ..
            },
        ) => {
            format!("sequence aliases '{a1}' and '{a2}'")
        }
        _ => format!("{first} and {second}"),
    };
    let why = if sanitized {
        "after replacing characters BLAST can't store"
    } else {
        "because BLAST compares ids ignoring case"
    };
    format!("{what} both become id '{id}' {why}; rename one in the store")
}

fn oversize(count: u64, examples: &[(String, u64)]) -> String {
    let ex: Vec<String> = examples
        .iter()
        .map(|(d, l)| format!("SQ.{d} ({l} bp)"))
        .collect();
    format!(
        "{count} sequences are 2^31 bp or longer, which BLAST can't store: {}{}",
        ex.join(", "),
        if (count as usize) > examples.len() {
            ", ..."
        } else {
            ""
        }
    )
}
