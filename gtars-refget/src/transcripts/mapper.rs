//! Coordinate mapping between transcript and genomic coordinates.

use thiserror::Error;

use crate::transcripts::models::{Strand, Transcript};
use crate::transcripts::store::ReadonlyTxStore;

/// Result of coordinate mapping.
#[derive(Debug, Clone)]
pub struct MappingResult {
    /// Genomic position (0-based).
    pub position: u64,
    /// Chromosome refget digest (24 bytes, truncated SHA-512).
    pub chrom_digest: [u8; 24],
}

/// Errors that can occur during coordinate mapping.
#[derive(Debug, Error)]
pub enum MappingError {
    #[error("Transcript not found: {0}")]
    TranscriptNotFound(String),

    #[error("No MANE Select transcript for gene: {0}")]
    NoManeTranscript(String),

    #[error("Position {0} is outside transcript")]
    OutsideTranscript(i64),

    #[error("Position {0} is outside CDS")]
    OutsideCds(i64),

    #[error("Position falls in intron (offset {0} from exon boundary)")]
    IntronicPosition(i64),

    #[error("Intronic offset {offset} at transcript position {pos} is invalid (not at exon boundary)")]
    InvalidIntronicOffset { pos: i64, offset: i64 },

    #[error("5' UTR position c.{0} extends beyond transcript start")]
    FivePrimeUtrOverflow(i64),

    #[error("3' UTR position c.*{0} extends beyond transcript end")]
    ThreePrimeUtrOverflow(i64),

    #[error("Non-coding transcript has no CDS")]
    NonCodingTranscript,
}

// ============================================================================
// CoordinateMapper: Simple API (allocates per call)
// ============================================================================

/// Coordinate mapper for one-off conversions.
pub struct CoordinateMapper<'a> {
    store: &'a ReadonlyTxStore,
}

impl<'a> CoordinateMapper<'a> {
    /// Create a coordinate mapper backed by a read-only transcript store.
    pub fn new(store: &'a ReadonlyTxStore) -> Self {
        Self { store }
    }

    /// Map a coding coordinate (c.) to genomic coordinate (g.).
    pub fn c_to_g(&self, accession: &str, c_pos: i64) -> Result<MappingResult, MappingError> {
        let tx = self
            .store
            .lookup(accession)
            .ok_or_else(|| MappingError::TranscriptNotFound(accession.to_string()))?;

        let cds_start = match (tx.cds_start, tx.cds_end) {
            (Some(s), Some(_)) => s,
            _ => return Err(MappingError::NonCodingTranscript),
        };

        let mut offsets = Vec::with_capacity(tx.exons.len());
        build_exon_offsets_into(&tx, &mut offsets);

        let cds_tx_start = genomic_to_transcript_fast(&tx, cds_start as u64, &offsets)
            .ok_or(MappingError::OutsideCds(c_pos))?;

        let tx_pos = if c_pos > 0 {
            cds_tx_start + (c_pos as u64) - 1
        } else {
            cds_tx_start
                .checked_sub((-c_pos) as u64)
                .ok_or(MappingError::OutsideCds(c_pos))?
        };

        let g_pos = transcript_to_genomic_fast(&tx, tx_pos, &offsets)?;

        Ok(MappingResult {
            position: g_pos,
            chrom_digest: tx.chrom_digest,
        })
    }

    /// Map a non-coding coordinate (n.) to genomic coordinate.
    pub fn n_to_g(&self, accession: &str, n_pos: u64) -> Result<MappingResult, MappingError> {
        if n_pos == 0 {
            return Err(MappingError::OutsideTranscript(0));
        }
        let tx = self
            .store
            .lookup(accession)
            .ok_or_else(|| MappingError::TranscriptNotFound(accession.to_string()))?;

        let mut offsets = Vec::with_capacity(tx.exons.len());
        build_exon_offsets_into(&tx, &mut offsets);
        let g_pos = transcript_to_genomic_fast(&tx, n_pos - 1, &offsets)?;

        Ok(MappingResult {
            position: g_pos,
            chrom_digest: tx.chrom_digest,
        })
    }

    /// Map a full HGVS c. coordinate to genomic.
    ///
    /// - `c_pos`: 1-based CDS position (negative for 5' UTR, e.g. `c.-14`)
    /// - `offset`: intronic offset (`c.93+5` → offset `+5`; `c.94-3` → offset `-3`)
    /// - `is_cds_end`: true for `c.*N` 3' UTR positions
    pub fn c_to_g_full(
        &self,
        accession: &str,
        c_pos: i64,
        offset: i64,
        is_cds_end: bool,
    ) -> Result<MappingResult, MappingError> {
        let tx = self
            .store
            .lookup(accession)
            .ok_or_else(|| MappingError::TranscriptNotFound(accession.to_string()))?;
        let mut offsets = Vec::with_capacity(tx.exons.len());
        build_exon_offsets_into(&tx, &mut offsets);
        c_to_g_inner(&tx, &offsets, c_pos, offset, is_cds_end)
    }

    /// Map a full HGVS n. coordinate (with optional intronic offset) to genomic.
    pub fn n_to_g_full(
        &self,
        accession: &str,
        n_pos: i64,
        offset: i64,
    ) -> Result<MappingResult, MappingError> {
        let tx = self
            .store
            .lookup(accession)
            .ok_or_else(|| MappingError::TranscriptNotFound(accession.to_string()))?;
        let mut offsets = Vec::with_capacity(tx.exons.len());
        build_exon_offsets_into(&tx, &mut offsets);
        n_to_g_inner(&tx, &offsets, n_pos, offset)
    }

    /// Map a 0-based genomic position to a 0-based offset on the mature
    /// (spliced) mRNA for `accession`. Returns `None` if the position is
    /// intronic or otherwise outside every exon (i.e. has no coordinate on
    /// the mature mRNA).
    ///
    /// This is the inverse of the forward c./n. -> genomic projection: it lets
    /// a caller take a genomic interbase coordinate (e.g. produced by
    /// [`Self::c_to_g_full`]) and re-anchor it onto the derived mature mRNA.
    /// `Some(off)` is an exonic offset; `None` signals an intronic / out-of-
    /// exon position which has no mature-mRNA coordinate.
    pub fn g_to_transcript_offset(
        &self,
        accession: &str,
        g_pos: u64,
    ) -> Result<Option<u64>, MappingError> {
        let tx = self
            .store
            .lookup(accession)
            .ok_or_else(|| MappingError::TranscriptNotFound(accession.to_string()))?;
        let mut offsets = Vec::with_capacity(tx.exons.len());
        build_exon_offsets_into(&tx, &mut offsets);
        Ok(genomic_to_transcript_fast(&tx, g_pos, &offsets))
    }

    /// Map a c. coordinate by gene symbol via MANE Select.
    ///
    /// Returns `(accession_used, result)` so the caller can record which
    /// transcript MANE resolved to.
    pub fn c_to_g_by_gene(
        &self,
        gene: &str,
        c_pos: i64,
        offset: i64,
        is_cds_end: bool,
    ) -> Result<(String, MappingResult), MappingError> {
        let tx = self
            .store
            .lookup_mane(gene)
            .ok_or_else(|| MappingError::NoManeTranscript(gene.to_string()))?;
        let mut offsets = Vec::with_capacity(tx.exons.len());
        build_exon_offsets_into(&tx, &mut offsets);
        let acc = tx.accession.clone();
        let result = c_to_g_inner(&tx, &offsets, c_pos, offset, is_cds_end)?;
        Ok((acc, result))
    }
}

// ============================================================================
// CoordinateMapperWriter: Zero-allocation API (reusable buffers)
// ============================================================================

/// Reusable coordinate mapper that avoids per-call allocations.
pub struct CoordinateMapperWriter<'a> {
    store: &'a ReadonlyTxStore,
    exon_offsets: Vec<ExonOffset>,
    /// Reusable string buffer for formatting positions. Kept WASM-safe (no
    /// `itoa`, which is a native-only builder dependency).
    fmt_buf: String,
}

/// Pre-computed exon offset for coordinate mapping.
#[derive(Clone, Copy)]
struct ExonOffset {
    tx_start: u64,
    tx_end: u64,
    g_start: u32,
    g_end: u32,
}

impl<'a> CoordinateMapperWriter<'a> {
    /// Create a reusable coordinate mapper that keeps scratch buffers between calls.
    pub fn new(store: &'a ReadonlyTxStore) -> Self {
        Self {
            store,
            exon_offsets: Vec::with_capacity(64),
            fmt_buf: String::with_capacity(24),
        }
    }

    pub fn c_to_g(&mut self, accession: &str, c_pos: i64) -> Result<MappingResult, MappingError> {
        let tx = self
            .store
            .lookup(accession)
            .ok_or_else(|| MappingError::TranscriptNotFound(accession.to_string()))?;

        let cds_start = match (tx.cds_start, tx.cds_end) {
            (Some(s), Some(_)) => s,
            _ => return Err(MappingError::NonCodingTranscript),
        };

        self.exon_offsets.clear();
        build_exon_offsets_into(&tx, &mut self.exon_offsets);

        let cds_tx_start = genomic_to_transcript_fast(&tx, cds_start as u64, &self.exon_offsets)
            .ok_or(MappingError::OutsideCds(c_pos))?;

        let tx_pos = if c_pos > 0 {
            cds_tx_start + (c_pos as u64) - 1
        } else {
            cds_tx_start
                .checked_sub((-c_pos) as u64)
                .ok_or(MappingError::OutsideCds(c_pos))?
        };

        let g_pos = transcript_to_genomic_fast(&tx, tx_pos, &self.exon_offsets)?;

        Ok(MappingResult {
            position: g_pos,
            chrom_digest: tx.chrom_digest,
        })
    }

    pub fn n_to_g(&mut self, accession: &str, n_pos: u64) -> Result<MappingResult, MappingError> {
        if n_pos == 0 {
            return Err(MappingError::OutsideTranscript(0));
        }
        let tx = self
            .store
            .lookup(accession)
            .ok_or_else(|| MappingError::TranscriptNotFound(accession.to_string()))?;

        self.exon_offsets.clear();
        build_exon_offsets_into(&tx, &mut self.exon_offsets);

        let g_pos = transcript_to_genomic_fast(&tx, n_pos - 1, &self.exon_offsets)?;

        Ok(MappingResult {
            position: g_pos,
            chrom_digest: tx.chrom_digest,
        })
    }

    #[inline]
    pub fn format_position(&mut self, pos: u64) -> &str {
        use std::fmt::Write;
        self.fmt_buf.clear();
        let _ = write!(self.fmt_buf, "{}", pos);
        &self.fmt_buf
    }
}

// ============================================================================
// Helpers
// ============================================================================

fn build_exon_offsets_into(tx: &Transcript, out: &mut Vec<ExonOffset>) {
    let mut tx_pos = 0u64;

    match tx.strand {
        Strand::Forward => {
            for exon in &tx.exons {
                let len = exon.len() as u64;
                out.push(ExonOffset {
                    tx_start: tx_pos,
                    tx_end: tx_pos + len,
                    g_start: exon.start,
                    g_end: exon.end,
                });
                tx_pos += len;
            }
        }
        Strand::Reverse => {
            for exon in tx.exons.iter().rev() {
                let len = exon.len() as u64;
                out.push(ExonOffset {
                    tx_start: tx_pos,
                    tx_end: tx_pos + len,
                    g_start: exon.start,
                    g_end: exon.end,
                });
                tx_pos += len;
            }
        }
    }
}

fn transcript_to_genomic_fast(
    tx: &Transcript,
    tx_pos: u64,
    offsets: &[ExonOffset],
) -> Result<u64, MappingError> {
    for eo in offsets {
        if tx_pos >= eo.tx_start && tx_pos < eo.tx_end {
            let offset = tx_pos - eo.tx_start;
            let g_pos = match tx.strand {
                Strand::Forward => eo.g_start as u64 + offset,
                Strand::Reverse => eo.g_end as u64 - 1 - offset,
            };
            return Ok(g_pos);
        }
    }
    Err(MappingError::OutsideTranscript(tx_pos as i64))
}

/// Find CDS start/end in transcript coordinates.
fn cds_tx_bounds(tx: &Transcript, offsets: &[ExonOffset]) -> Option<(u64, u64)> {
    let cds_g_start = tx.cds_start?;
    let cds_g_end = tx.cds_end?;

    // cds_g_end is exclusive in genomic; map the last included base (cds_g_end - 1).
    let last_cds_g = cds_g_end.checked_sub(1)?;

    let tx_a = genomic_to_transcript_fast(tx, cds_g_start as u64, offsets)?;
    let tx_b = genomic_to_transcript_fast(tx, last_cds_g as u64, offsets)?;
    let (lo, hi) = if tx_a <= tx_b { (tx_a, tx_b) } else { (tx_b, tx_a) };
    // tx-coordinate CDS is [lo, hi+1) (half-open, like genomic).
    Some((lo, hi + 1))
}

/// Verify that `tx_pos` sits on an exon boundary in transcript space, where
/// `offset_positive=true` requires it to be the last base of an exon and
/// `offset_positive=false` requires it to be the first base of an exon that
/// is not the very first exon.
fn is_exon_boundary(tx_pos: u64, offsets: &[ExonOffset], offset_positive: bool) -> bool {
    for (i, eo) in offsets.iter().enumerate() {
        if offset_positive && tx_pos + 1 == eo.tx_end && i + 1 < offsets.len() {
            return true;
        }
        if !offset_positive && tx_pos == eo.tx_start && i > 0 {
            return true;
        }
    }
    false
}

/// Core c. → g. mapper.
fn c_to_g_inner(
    tx: &Transcript,
    offsets: &[ExonOffset],
    c_pos: i64,
    offset: i64,
    is_cds_end: bool,
) -> Result<MappingResult, MappingError> {
    if !tx.is_coding() {
        return Err(MappingError::NonCodingTranscript);
    }
    let (cds_tx_start, cds_tx_end) =
        cds_tx_bounds(tx, offsets).ok_or(MappingError::NonCodingTranscript)?;
    let tx_len = offsets.last().map(|e| e.tx_end).unwrap_or(0);

    // Step 1: c.* / c.N / c.-N → 0-based transcript position
    let tx_pos: u64 = if is_cds_end {
        if c_pos <= 0 {
            return Err(MappingError::ThreePrimeUtrOverflow(c_pos));
        }
        let pos = cds_tx_end
            .checked_add((c_pos as u64).saturating_sub(1))
            .ok_or(MappingError::ThreePrimeUtrOverflow(c_pos))?;
        if pos >= tx_len {
            return Err(MappingError::ThreePrimeUtrOverflow(c_pos));
        }
        pos
    } else if c_pos > 0 {
        let pos = cds_tx_start + (c_pos as u64) - 1;
        if pos >= cds_tx_end {
            return Err(MappingError::OutsideCds(c_pos));
        }
        pos
    } else if c_pos < 0 {
        // 5' UTR: c.-1 is the base immediately preceding the CDS start.
        let utr = (-c_pos) as u64;
        if utr > cds_tx_start {
            return Err(MappingError::FivePrimeUtrOverflow(c_pos));
        }
        cds_tx_start - utr
    } else {
        // c.0 is not a valid HGVS coordinate.
        return Err(MappingError::OutsideCds(c_pos));
    };

    apply_offset_to_tx_pos(tx, offsets, tx_pos, offset, c_pos)
}

/// Core n. → g. mapper with intronic offset.
fn n_to_g_inner(
    tx: &Transcript,
    offsets: &[ExonOffset],
    n_pos: i64,
    offset: i64,
) -> Result<MappingResult, MappingError> {
    if n_pos <= 0 {
        return Err(MappingError::OutsideTranscript(n_pos));
    }
    let tx_pos = (n_pos as u64) - 1;
    let tx_len = offsets.last().map(|e| e.tx_end).unwrap_or(0);
    if tx_pos >= tx_len {
        return Err(MappingError::OutsideTranscript(n_pos));
    }
    apply_offset_to_tx_pos(tx, offsets, tx_pos, offset, n_pos)
}

/// Apply optional intronic offset to a transcript position and produce genomic.
fn apply_offset_to_tx_pos(
    tx: &Transcript,
    offsets: &[ExonOffset],
    tx_pos: u64,
    offset: i64,
    original_pos: i64,
) -> Result<MappingResult, MappingError> {
    if offset == 0 {
        let g_pos = transcript_to_genomic_fast(tx, tx_pos, offsets)?;
        return Ok(MappingResult {
            position: g_pos,
            chrom_digest: tx.chrom_digest,
        });
    }

    let offset_positive = offset > 0;
    if !is_exon_boundary(tx_pos, offsets, offset_positive) {
        return Err(MappingError::InvalidIntronicOffset {
            pos: original_pos,
            offset,
        });
    }

    // Map tx_pos → genomic at the boundary base, then move into the intron.
    let g_anchor = transcript_to_genomic_fast(tx, tx_pos, offsets)?;

    let g_pos: u64 = match tx.strand {
        Strand::Forward => {
            // Positive tx offset = +genomic; negative tx offset = -genomic.
            if offset_positive {
                g_anchor + offset as u64
            } else {
                let off = (-offset) as u64;
                g_anchor
                    .checked_sub(off)
                    .ok_or(MappingError::InvalidIntronicOffset {
                        pos: original_pos,
                        offset,
                    })?
            }
        }
        Strand::Reverse => {
            // Reverse strand: tx-positive offset = -genomic.
            if offset_positive {
                let off = offset as u64;
                g_anchor
                    .checked_sub(off)
                    .ok_or(MappingError::InvalidIntronicOffset {
                        pos: original_pos,
                        offset,
                    })?
            } else {
                let off = (-offset) as u64;
                g_anchor + off
            }
        }
    };

    Ok(MappingResult {
        position: g_pos,
        chrom_digest: tx.chrom_digest,
    })
}

fn genomic_to_transcript_fast(
    tx: &Transcript,
    g_pos: u64,
    offsets: &[ExonOffset],
) -> Option<u64> {
    for eo in offsets {
        if g_pos >= eo.g_start as u64 && g_pos < eo.g_end as u64 {
            let offset = match tx.strand {
                Strand::Forward => g_pos - eo.g_start as u64,
                Strand::Reverse => eo.g_end as u64 - 1 - g_pos,
            };
            return Some(eo.tx_start + offset);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcripts::models::{Exon, ManeStatus};
    use crate::transcripts::store::{build_reftx_bytes_in_memory, ReadonlyTxStore};

    fn make_tx(acc: &str, strand: Strand, exons: Vec<Exon>) -> Transcript {
        Transcript {
            accession: acc.to_string(),
            gene: "G".to_string(),
            chrom_digest: [0u8; 24],
            strand,
            cds_start: None,
            cds_end: None,
            exons,
            mane: ManeStatus::default(),
        }
    }

    fn store_for(tx: Transcript) -> ReadonlyTxStore {
        let bytes = build_reftx_bytes_in_memory(&[tx]).unwrap();
        ReadonlyTxStore::from_bytes(bytes).unwrap()
    }

    #[test]
    fn g_to_transcript_offset_forward_exonic() {
        // Two forward exons: [10,14) and [20,24). mRNA offsets:
        //   genomic 10..14 -> tx 0..4
        //   genomic 20..24 -> tx 4..8
        let tx = make_tx(
            "NM_F.1",
            Strand::Forward,
            vec![Exon { start: 10, end: 14 }, Exon { start: 20, end: 24 }],
        );
        let store = store_for(tx);
        let mapper = CoordinateMapper::new(&store);

        assert_eq!(mapper.g_to_transcript_offset("NM_F.1", 10).unwrap(), Some(0));
        assert_eq!(mapper.g_to_transcript_offset("NM_F.1", 13).unwrap(), Some(3));
        assert_eq!(mapper.g_to_transcript_offset("NM_F.1", 20).unwrap(), Some(4));
        assert_eq!(mapper.g_to_transcript_offset("NM_F.1", 23).unwrap(), Some(7));
        // Intronic (between the two exons) -> None.
        assert_eq!(mapper.g_to_transcript_offset("NM_F.1", 16).unwrap(), None);
        // Outside all exons -> None.
        assert_eq!(mapper.g_to_transcript_offset("NM_F.1", 0).unwrap(), None);
    }

    #[test]
    fn g_to_transcript_offset_reverse_exonic() {
        // Two exons (genomic order) [10,14) and [20,24) on the reverse strand.
        // Transcript runs 3'->5' in genomic terms, so the last genomic base of
        // the last exon is tx offset 0:
        //   genomic 23 -> tx 0, genomic 20 -> tx 3,
        //   genomic 13 -> tx 4, genomic 10 -> tx 7.
        let tx = make_tx(
            "NM_R.1",
            Strand::Reverse,
            vec![Exon { start: 10, end: 14 }, Exon { start: 20, end: 24 }],
        );
        let store = store_for(tx);
        let mapper = CoordinateMapper::new(&store);

        assert_eq!(mapper.g_to_transcript_offset("NM_R.1", 23).unwrap(), Some(0));
        assert_eq!(mapper.g_to_transcript_offset("NM_R.1", 20).unwrap(), Some(3));
        assert_eq!(mapper.g_to_transcript_offset("NM_R.1", 13).unwrap(), Some(4));
        assert_eq!(mapper.g_to_transcript_offset("NM_R.1", 10).unwrap(), Some(7));
        // Intronic -> None.
        assert_eq!(mapper.g_to_transcript_offset("NM_R.1", 16).unwrap(), None);
    }

    #[test]
    fn g_to_transcript_offset_unknown_accession() {
        let tx = make_tx("NM_F.1", Strand::Forward, vec![Exon { start: 0, end: 4 }]);
        let store = store_for(tx);
        let mapper = CoordinateMapper::new(&store);
        assert!(mapper.g_to_transcript_offset("NM_MISSING.1", 0).is_err());
    }
}
