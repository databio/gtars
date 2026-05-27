//! Errors produced by HGVS parsing and the VRS bridge.

use thiserror::Error;

use super::ast::ReferenceType;

#[derive(Debug, Error)]
pub enum HgvsError {
    #[error("HGVS parse error at position {pos}: {msg} (in: {input:?})")]
    Parse {
        input: String,
        pos: usize,
        msg: String,
    },

    #[error("Unsupported reference type {0:?} for VRS conversion")]
    UnsupportedReferenceType(ReferenceType),

    #[error("Unsupported edit type: {0}")]
    UnsupportedEdit(String),

    #[error("Accession {accession} not found in collection {collection}")]
    AccessionNotFound {
        accession: String,
        collection: String,
    },

    #[error("Transcript not found: {0}")]
    TranscriptNotFound(String),

    #[error("No MANE Select transcript found for gene: {0}")]
    NoManeTranscript(String),

    #[error("Transcript provider required for {0} but none was supplied")]
    NoTranscriptProvider(String),

    #[error("Coordinate mapping failed for {accession}: {detail}")]
    MappingError { accession: String, detail: String },

    #[error("Invalid chromosome accession format: {0}")]
    InvalidChromAccession(String),

    #[error("Refget error: {0}")]
    RefgetError(String),
}
