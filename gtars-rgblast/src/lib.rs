//! BLAST nucleotide databases built straight from a RefgetStore.
//!
//! - [`api`]: the binding surface (`preflight`, `build`, `verify`, `view`, `nstats`); no
//!   printing, progress through a callback, typed errors, serializable options and summaries.
//! - [`select`]: which collections and sequences (`Selection`, resolution, N counts).
//! - [`plan`], [`store`], [`convert`]: the canonical `SubjectPlan` (OID order, ids, taxids,
//!   titles, totals), reading stores through gtars-refget, and sequence conversion. The
//!   rgblastn engine uses only these and `select`; they never import [`blastdb`].
//! - [`blastdb`]: the version 5 database writer, its reader, alias files and `verify`.
//! - [`views`]: subset views without a rebuild.

// Errors are the cold path; `Error` carries enough context to explain itself.
#![allow(clippy::result_large_err)]

pub mod api;
pub mod blastdb;
pub mod cache_root;
pub mod cli;
pub mod convert;
pub mod engine;
pub mod error;
pub mod format;
pub mod origin;
pub mod pipeline;
pub mod plan;
pub mod progress;
pub mod report;
pub mod run;
pub mod select;
pub mod store;
pub mod views;

#[cfg(feature = "reference-fasta")]
pub mod reference_fasta;

pub use api::{
    BuildOptions, BuildRequest, NstatsRequest, VerifyRequest, build, nstats, preflight, verify,
    view,
};
pub use blastdb::preflight::PreflightSummary;
pub use blastdb::verify::VerifyReport;
pub use blastdb::{BuildSummary, NetworkReport, TaxonCount, Timings};
pub use error::{Error, ErrorKind, IdOrigin};
pub use plan::{DbTotals, OnOversize, ResolvedSelection, SUBJECT_PLAN_VERSION, SubjectPlan};
pub use progress::{ProgressEvent, describe, event_json};
pub use select::{CollectionSet, Selection};
pub use store::OpenedStore;
