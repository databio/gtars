//! TranscriptProvider implementation for ReadonlyTxStore.

use std::sync::Arc;

use gtars_vrs::models::{SequenceLocation, SequenceReference};
use gtars_vrs::provider::{ProviderError, TranscriptProvider};

use crate::mapper::CoordinateMapper;
use crate::models::Strand;
use crate::store::ReadonlyTxStore;

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

/// Thread-safe transcript provider using `Arc<ReadonlyTxStore>`.
///
/// For concurrent HGVS parsing, share a single provider across workers:
///
/// ```ignore
/// let store = TxStore::open("transcripts.reftx")?.into_readonly();
/// let provider = ReftxProvider::new(Arc::new(store));
///
/// // Clone is cheap (just bumps an Arc).
/// for hgvs in variants {
///     let loc = provider.c_to_genomic(&hgvs.accession, hgvs.c_pos)?;
/// }
/// ```
#[derive(Clone)]
pub struct ReftxProvider {
    store: Arc<ReadonlyTxStore>,
}

impl ReftxProvider {
    pub fn new(store: Arc<ReadonlyTxStore>) -> Self {
        Self { store }
    }

    /// Access the underlying store.
    pub fn store(&self) -> &Arc<ReadonlyTxStore> {
        &self.store
    }
}

impl TranscriptProvider for ReftxProvider {
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
