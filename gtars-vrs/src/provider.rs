//! Transcript provider trait for HGVS coordinate resolution.
//!
//! This trait abstracts over different transcript annotation sources
//! (e.g., `gtars-reftx`, UTA, etc.). HGVS parsers use it to resolve
//! transcript-relative coordinates (c./n.) to genomic SequenceLocations.

use thiserror::Error;

use crate::models::SequenceLocation;

/// Errors from transcript operations.
#[derive(Debug, Error)]
pub enum ProviderError {
    #[error("Transcript not found: {0}")]
    TranscriptNotFound(String),

    #[error("No MANE Select transcript for gene: {0}")]
    NoManeTranscript(String),

    #[error("Invalid coordinate: {0}")]
    InvalidCoordinate(String),

    #[error("Mapping error: {0}")]
    MappingError(String),

    #[error("Non-coding transcript")]
    NonCodingTranscript,
}

/// Trait for looking up transcript annotations.
///
/// Implemented by `gtars-reftx::ReadonlyTxStore` when the `vrs` feature is
/// enabled on the `gtars-reftx` crate.
pub trait TranscriptProvider {
    /// Map a coding coordinate (c.) to genomic SequenceLocation.
    ///
    /// c. coordinates are 1-based, relative to CDS start.
    fn c_to_genomic(&self, accession: &str, c_pos: i64) -> Result<SequenceLocation, ProviderError>;

    /// Map a non-coding coordinate (n.) to genomic SequenceLocation.
    ///
    /// n. coordinates are 1-based, relative to transcript start.
    fn n_to_genomic(&self, accession: &str, n_pos: u64) -> Result<SequenceLocation, ProviderError>;

    /// Get the chromosome refget accession for a transcript.
    fn get_chrom_accession(&self, accession: &str) -> Result<String, ProviderError>;

    /// Get transcript strand (+1 or -1).
    fn get_strand(&self, accession: &str) -> Result<i8, ProviderError>;

    /// Map a coding coordinate (c.) with optional intronic offset and an
    /// `is_cds_end` flag (for `c.*N` 3' UTR coordinates).
    ///
    /// Default impl falls back to `c_to_genomic` when `offset == 0` and
    /// `is_cds_end == false`; otherwise it returns
    /// `ProviderError::InvalidCoordinate` so implementations that don't yet
    /// support intronic offsets surface a clear error.
    fn c_to_genomic_full(
        &self,
        accession: &str,
        c_pos: i64,
        offset: i64,
        is_cds_end: bool,
    ) -> Result<SequenceLocation, ProviderError> {
        if offset == 0 && !is_cds_end {
            self.c_to_genomic(accession, c_pos)
        } else {
            Err(ProviderError::InvalidCoordinate(format!(
                "c. position with offset {} or CDS-end flag {} not supported by this provider",
                offset, is_cds_end
            )))
        }
    }

    /// Map an n. coordinate with optional intronic offset.
    fn n_to_genomic_full(
        &self,
        accession: &str,
        n_pos: i64,
        offset: i64,
    ) -> Result<SequenceLocation, ProviderError> {
        if offset == 0 && n_pos > 0 {
            self.n_to_genomic(accession, n_pos as u64)
        } else {
            Err(ProviderError::InvalidCoordinate(format!(
                "n. position with offset {} or non-positive base {} not supported by this provider",
                offset, n_pos
            )))
        }
    }

    /// Resolve a gene symbol to its MANE Select transcript accession.
    ///
    /// Default impl returns `None` (no MANE awareness).
    fn gene_to_mane_accession(&self, _gene: &str) -> Option<String> {
        None
    }
}

/// No-op provider that rejects transcript lookups. Use when only `g.` /
/// `m.` variants are expected.
pub struct NoTranscriptProvider;

impl TranscriptProvider for NoTranscriptProvider {
    fn c_to_genomic(&self, accession: &str, _: i64) -> Result<SequenceLocation, ProviderError> {
        Err(ProviderError::TranscriptNotFound(accession.to_string()))
    }

    fn n_to_genomic(&self, accession: &str, _: u64) -> Result<SequenceLocation, ProviderError> {
        Err(ProviderError::TranscriptNotFound(accession.to_string()))
    }

    fn get_chrom_accession(&self, accession: &str) -> Result<String, ProviderError> {
        Err(ProviderError::TranscriptNotFound(accession.to_string()))
    }

    fn get_strand(&self, accession: &str) -> Result<i8, ProviderError> {
        Err(ProviderError::TranscriptNotFound(accession.to_string()))
    }
}
