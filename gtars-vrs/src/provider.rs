//! Transcript provider trait and concrete implementations for HGVS coordinate resolution.
//!
//! This module defines the `TranscriptProvider` trait, which abstracts over
//! different transcript annotation sources. HGVS parsers use it to resolve
//! transcript-relative coordinates (c./n.) to genomic `SequenceLocation`s.
//!
//! ## Concrete provider
//!
//! When the `transcripts` feature is enabled on `gtars-vrs`, this module also
//! provides `TxProvider`, a concrete `Arc<ReadonlyTxStore>`-backed provider
//! implemented over the `gtars-refget` transcript store. This impl is available
//! on **both native and WASM targets** — it is gated on the transcripts CORE
//! feature, NOT on `filesystem`, so `c.`/`n.` resolution is possible in the
//! browser once a real provider is wired in.
//!
//! On native targets the same `TxProvider` wraps a mmap-backed store opened via
//! `gtars_refget::transcripts::TxStore::open(...).into_readonly()` (the
//! opener types are native-only, behind `gtars-refget`'s `filesystem` gate). On
//! WASM, `TxProvider` wraps an in-memory store constructed via the WASM-safe
//! in-memory constructor.

use thiserror::Error;

use crate::models::SequenceLocation;

// Concrete provider impl — available on all targets (WASM-safe core).
#[cfg(feature = "transcripts")]
use std::sync::Arc;
#[cfg(feature = "transcripts")]
use gtars_refget::transcripts::{CoordinateMapper, ReadonlyTxStore, Strand};

#[cfg(feature = "transcripts")]
use crate::models::SequenceReference;

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
/// Implemented in this crate (when the `transcripts` feature is enabled) by
/// `ReadonlyTxStore` from `gtars-refget`. That impl is available on both native
/// and WASM targets — see the module-level doc for details.
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

// ============================================================================
// Concrete provider over the gtars-refget transcript store.
// Gated on the `transcripts` CORE feature — WASM-safe (no mmap, no std::fs).
// ============================================================================

/// `TranscriptProvider` impl for the readonly transcript store from
/// `gtars-refget`. Delegates each method to a `CoordinateMapper` over the
/// store.
///
/// This impl is available on ALL targets (native and WASM) when the
/// `transcripts` feature is enabled. It contains no `memmap2`, no `std::fs`,
/// and no file-open calls — those live behind the store's own `filesystem`
/// gate.
///
/// On native, the store is typically opened via
/// `gtars_refget::transcripts::TxStore::open(...).into_readonly()` (the
/// opener types are native-only, behind `filesystem`). On WASM, the same
/// `ReadonlyTxStore` is constructed from in-memory bytes.
#[cfg(feature = "transcripts")]
impl TranscriptProvider for ReadonlyTxStore {
    fn c_to_genomic(&self, accession: &str, c_pos: i64) -> Result<SequenceLocation, ProviderError> {
        let mapper = CoordinateMapper::new(self);
        let result = mapper
            .c_to_g(accession, c_pos)
            .map_err(|e| ProviderError::MappingError(e.to_string()))?;

        Ok(SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: format!("SQ.{}", base64_url::encode(&result.chrom_digest)),
            },
            start: result.position,
            end: result.position + 1,
        })
    }

    fn n_to_genomic(&self, accession: &str, n_pos: u64) -> Result<SequenceLocation, ProviderError> {
        let mapper = CoordinateMapper::new(self);
        let result = mapper
            .n_to_g(accession, n_pos)
            .map_err(|e| ProviderError::MappingError(e.to_string()))?;

        Ok(SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: format!("SQ.{}", base64_url::encode(&result.chrom_digest)),
            },
            start: result.position,
            end: result.position + 1,
        })
    }

    fn get_chrom_accession(&self, accession: &str) -> Result<String, ProviderError> {
        let tx = self
            .lookup(accession)
            .ok_or_else(|| ProviderError::TranscriptNotFound(accession.to_string()))?;
        Ok(format!("SQ.{}", base64_url::encode(&tx.chrom_digest)))
    }

    fn get_strand(&self, accession: &str) -> Result<i8, ProviderError> {
        let tx = self
            .lookup(accession)
            .ok_or_else(|| ProviderError::TranscriptNotFound(accession.to_string()))?;
        Ok(match tx.strand {
            Strand::Forward => 1,
            Strand::Reverse => -1,
        })
    }

    fn c_to_genomic_full(
        &self,
        accession: &str,
        c_pos: i64,
        offset: i64,
        is_cds_end: bool,
    ) -> Result<SequenceLocation, ProviderError> {
        let mapper = CoordinateMapper::new(self);
        let result = mapper
            .c_to_g_full(accession, c_pos, offset, is_cds_end)
            .map_err(|e| ProviderError::MappingError(e.to_string()))?;
        Ok(SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: format!("SQ.{}", base64_url::encode(&result.chrom_digest)),
            },
            start: result.position,
            end: result.position + 1,
        })
    }

    fn n_to_genomic_full(
        &self,
        accession: &str,
        n_pos: i64,
        offset: i64,
    ) -> Result<SequenceLocation, ProviderError> {
        let mapper = CoordinateMapper::new(self);
        let result = mapper
            .n_to_g_full(accession, n_pos, offset)
            .map_err(|e| ProviderError::MappingError(e.to_string()))?;
        Ok(SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: format!("SQ.{}", base64_url::encode(&result.chrom_digest)),
            },
            start: result.position,
            end: result.position + 1,
        })
    }

    fn gene_to_mane_accession(&self, gene: &str) -> Option<String> {
        self.lookup_mane(gene).map(|tx| tx.accession)
    }
}

/// Thread-safe transcript provider wrapping `Arc<ReadonlyTxStore>`.
///
/// For concurrent HGVS parsing, share a single provider across workers:
///
/// ```ignore
/// // Native path (opener types are native-only, behind gtars-refget `filesystem` gate):
/// use gtars_refget::transcripts::TxStore;
/// use std::sync::Arc;
/// use gtars_vrs::TxProvider;
///
/// let store = TxStore::open("transcripts.reftx")?.into_readonly();
/// let provider = TxProvider::new(Arc::new(store));
///
/// // Clone is cheap (just bumps an Arc).
/// for hgvs in variants {
///     let loc = provider.c_to_genomic(&hgvs.accession, hgvs.c_pos)?;
/// }
///
/// // WASM path: same TxProvider wraps an in-memory store built via the
/// // WASM-safe in-memory constructor from gtars-refget.
/// ```
#[cfg(feature = "transcripts")]
#[derive(Clone)]
pub struct TxProvider {
    store: Arc<ReadonlyTxStore>,
}

#[cfg(feature = "transcripts")]
impl TxProvider {
    pub fn new(store: Arc<ReadonlyTxStore>) -> Self {
        Self { store }
    }

    /// Access the underlying store.
    pub fn store(&self) -> &Arc<ReadonlyTxStore> {
        &self.store
    }
}

#[cfg(feature = "transcripts")]
impl TranscriptProvider for TxProvider {
    fn c_to_genomic(&self, accession: &str, c_pos: i64) -> Result<SequenceLocation, ProviderError> {
        self.store.c_to_genomic(accession, c_pos)
    }

    fn n_to_genomic(&self, accession: &str, n_pos: u64) -> Result<SequenceLocation, ProviderError> {
        self.store.n_to_genomic(accession, n_pos)
    }

    fn get_chrom_accession(&self, accession: &str) -> Result<String, ProviderError> {
        self.store.get_chrom_accession(accession)
    }

    fn get_strand(&self, accession: &str) -> Result<i8, ProviderError> {
        self.store.get_strand(accession)
    }

    fn c_to_genomic_full(
        &self,
        accession: &str,
        c_pos: i64,
        offset: i64,
        is_cds_end: bool,
    ) -> Result<SequenceLocation, ProviderError> {
        self.store.c_to_genomic_full(accession, c_pos, offset, is_cds_end)
    }

    fn n_to_genomic_full(
        &self,
        accession: &str,
        n_pos: i64,
        offset: i64,
    ) -> Result<SequenceLocation, ProviderError> {
        self.store.n_to_genomic_full(accession, n_pos, offset)
    }

    fn gene_to_mane_accession(&self, gene: &str) -> Option<String> {
        self.store.gene_to_mane_accession(gene)
    }
}
