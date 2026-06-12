//! Mature (spliced) mRNA sequence lookup from a refget sequence store + transcript store.
//!
//! This module provides two layers of API:
//!
//! 1. [`concat_regions`] — a reusable primitive that fetches and stitches an ordered
//!    set of genomic regions from a sequence store, applying reverse-complementation
//!    when the transcript is on the reverse strand.
//!
//! 2. [`mature_mrna_for_transcript`] and [`mature_mrna`] — transcript-level
//!    conveniences that map a [`Transcript`] (or an accession string) to its
//!    mature mRNA reference sequence by calling layer 1 with the transcript's
//!    exons, strand, and chromosome digest.
//!
//! # Notes
//!
//! - The returned sequence is the **unedited reference** mature mRNA.
//!   Applying variants is explicitly out of scope.
//! - There is no new on-disk format; the computation is performed at runtime
//!   using only the existing sequence store and transcript store.

use anyhow::Result;

use crate::store::ReadonlyRefgetStore;
use crate::transcripts::models::{Strand, Transcript};
use crate::transcripts::store::ReadonlyTxStore;

// ============================================================================
// Private helpers
// ============================================================================

/// Reverse-complement a DNA sequence string.
///
/// Handles upper- and lower-case A/C/G/T/N. Any unrecognised byte maps to `N`.
/// Bytes are ASCII per the store's contract; collecting `char` from ASCII bytes
/// is sound.
fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'N' => b'N',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            b'n' => b'n',
            _ => b'N',
        } as char)
        .collect()
}

// ============================================================================
// Layer 1 – region-concatenation primitive
// ============================================================================

/// Fetch and concatenate an ordered set of genomic regions from a sequence
/// store, applying strand correction for the result.
///
/// # Arguments
///
/// * `store` — the refget sequence store to read from.
/// * `chrom_digest` — 24 raw bytes of the SHA-512 hash for the chromosome
///   (i.e. `Transcript::chrom_digest`). These are the same bytes written at
///   ingest; `base64_url::encode` of this array reproduces the exact store key.
/// * `regions` — genomic 0-based half-open intervals `(start, end)` in
///   **genomic 5'->3' order** (as stored in `Transcript::exons`). Empty regions
///   (`start == end`) are silently skipped; a caller should not pass them since
///   real exons are always non-empty.
/// * `strand` — the transcript strand. For `Strand::Reverse`, the full
///   concatenated genomic-order sequence is reverse-complemented in one pass,
///   which simultaneously complements each base AND reverses the segment order —
///   no separate exon re-ordering step is needed.
///
/// # Errors
///
/// * Chromosome digest not in the store → `Err("Sequence not found: <digest>")`.
/// * An exon coordinate exceeds the chromosome length → `Err` with bounds info.
///   This indicates a transcript/assembly mismatch.
pub fn concat_regions(
    store: &ReadonlyRefgetStore,
    chrom_digest: &[u8; 24],
    regions: &[(u32, u32)],
    strand: Strand,
) -> Result<String> {
    if regions.is_empty() {
        return Ok(String::new());
    }

    // Convert the raw 24-byte digest to the base64url string that is the
    // sequence store key. This is the exact inverse of what ingest does:
    //   `base64_url::encode(&sha512_hash[0..24])`
    // so the round-trip is lossless.
    let digest_str = base64_url::encode(chrom_digest);

    // Defensively filter zero-length exons (real exons are non-empty, but be
    // safe). Then convert to `(usize, usize)` for the store API.
    let ranges: Vec<(usize, usize)> = regions
        .iter()
        .filter(|&&(s, e)| s < e)
        .map(|&(s, e)| (s as usize, e as usize))
        .collect();

    if ranges.is_empty() {
        return Ok(String::new());
    }

    // One record resolution, many range reads — O(1) record lookup + one read
    // per range, reusing the open .seq handle.
    let pieces = store.get_substrings(&digest_str, &ranges)?;
    let seq: String = pieces.concat();

    // For reverse-strand transcripts, apply reverse-complement to the
    // genomic-order concatenation. Reverse-complementing the full string
    // simultaneously:
    //   1. Reverses the order of all bases (which also reverses segment order),
    //   2. Complements each base.
    // The result is the correct transcript-order mature sequence — no explicit
    // exon reordering is necessary.
    if strand == Strand::Reverse {
        Ok(reverse_complement(&seq))
    } else {
        Ok(seq)
    }
}

// ============================================================================
// Layer 2 – transcript-level conveniences
// ============================================================================

/// Compose the mature (spliced) mRNA reference sequence for a transcript
/// already in hand.
///
/// This is a thin convenience wrapper over [`concat_regions`] that maps the
/// transcript's exons to `(u32, u32)` region tuples and passes through the
/// chromosome digest and strand.
///
/// # Errors
///
/// Propagates all errors from [`concat_regions`] (sequence not found, out-of-
/// bounds exon, etc.).
pub fn mature_mrna_for_transcript(
    store: &ReadonlyRefgetStore,
    tx: &Transcript,
) -> Result<String> {
    let regions: Vec<(u32, u32)> = tx.exons.iter().map(|e| (e.start, e.end)).collect();
    concat_regions(store, &tx.chrom_digest, &regions, tx.strand)
}

/// Look the accession up in the transcript store, then compose its mature mRNA.
///
/// This is a thin convenience wrapper over [`mature_mrna_for_transcript`] that
/// performs the accession lookup and propagates a descriptive error if the
/// accession is absent.
///
/// # Errors
///
/// * Accession not in the transcript store → `Err("Transcript not found: <acc>")`.
/// * Chromosome digest not in the sequence store → `Err("Sequence not found: <digest>")`.
/// * Out-of-bounds exon coordinate (transcript/assembly mismatch) → descriptive `Err`.
pub fn mature_mrna(
    store: &ReadonlyRefgetStore,
    tx_store: &ReadonlyTxStore,
    accession: &str,
) -> Result<String> {
    let tx = tx_store
        .lookup(accession)
        .ok_or_else(|| anyhow::anyhow!("Transcript not found: {accession}"))?;
    // TranscriptRef derefs to &Transcript.
    mature_mrna_for_transcript(store, &tx)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::digest_sequence;
    use crate::store::ReadonlyRefgetStore;
    use crate::transcripts::models::{Exon, ManeStatus, Strand, Transcript};
    use crate::transcripts::store::{build_reftx_bytes_in_memory, ReadonlyTxStore};

    /// Build a `ReadonlyRefgetStore` that holds a single short chromosome in
    /// Raw (un-encoded) mode so substring retrieval returns the literal ASCII
    /// bytes we put in.
    fn make_seq_store(name: &str, seq: &str) -> (ReadonlyRefgetStore, String) {
        let sr = digest_sequence(name, seq.as_bytes());
        let digest = sr.metadata().sha512t24u.clone();
        let mut store = ReadonlyRefgetStore::new(crate::store::StorageMode::Raw);
        store.add_sequence_record(sr, false).unwrap();
        (store, digest)
    }

    /// Decode a base64url store key back to the raw 24-byte digest.
    fn key_to_chrom_digest(key: &str) -> [u8; 24] {
        let bytes = base64_url::decode(key).expect("valid base64url");
        assert_eq!(bytes.len(), 24, "expected 24 bytes from sha512t24u store key");
        let mut out = [0u8; 24];
        out.copy_from_slice(&bytes);
        out
    }

    fn make_transcript(acc: &str, chrom_digest: [u8; 24], strand: Strand, exons: Vec<Exon>) -> Transcript {
        Transcript {
            accession: acc.to_string(),
            gene: "TEST_GENE".to_string(),
            chrom_digest,
            strand,
            cds_start: None,
            cds_end: None,
            exons,
            mane: ManeStatus::default(),
        }
    }

    #[test]
    fn test_reverse_complement_helper() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("ACGTTTAA"), "TTAAACGT");
        assert_eq!(reverse_complement(""), "");
        // lowercase
        assert_eq!(reverse_complement("acgt"), "acgt");
    }

    #[test]
    fn test_concat_regions_forward() {
        // Chromosome: ATCGATCGATCGATCGATCGATCGATCGATCG (32 bp)
        // Exon 1: [2, 6)  -> "CGAT"
        // Exon 2: [10, 14) -> "CGAT"
        // Expected: "CGATCGAT"
        let chrom = "ATCGATCGATCGATCGATCGATCGATCGATCG";
        let (store, digest) = make_seq_store("chr1", chrom);
        let chrom_digest = key_to_chrom_digest(&digest);

        // Verify round-trip: encode(decode(key)) == key
        let round_tripped = base64_url::encode(&chrom_digest);
        assert_eq!(round_tripped, digest, "chrom_digest round-trip must match store key");

        let regions = vec![(2u32, 6u32), (10u32, 14u32)];
        let result = concat_regions(&store, &chrom_digest, &regions, Strand::Forward).unwrap();
        let expected = format!("{}{}", &chrom[2..6], &chrom[10..14]);
        assert_eq!(result, expected);
        assert_eq!(result, "CGATCGAT");
    }

    #[test]
    fn test_concat_regions_reverse() {
        // Using a known sequence where the reverse-complement result can be
        // verified by hand:
        //   exon 1 substring: "ACGT" (genomic pos 0..4)
        //   exon 2 substring: "TTAA" (genomic pos 4..8)
        //   genomic concat: "ACGTTTAA"
        //   reverse-complement: rev("ACGTTTAA") = "AATTTGCA" then complement...
        //   Actually: complement each base of reversed string:
        //     reversed: "AATTTGCA"
        //     complement: "TTAAACGT"
        //   So: reverse_complement("ACGTTTAA") == "TTAAACGT"
        let chrom = "ACGTTTAAGGCCAACCGGTT"; // 20 bp
        let (store, digest) = make_seq_store("chrRev", chrom);
        let chrom_digest = key_to_chrom_digest(&digest);

        let regions = vec![(0u32, 4u32), (4u32, 8u32)];
        let result = concat_regions(&store, &chrom_digest, &regions, Strand::Reverse).unwrap();
        let genomic_concat = format!("{}{}", &chrom[0..4], &chrom[4..8]);
        assert_eq!(genomic_concat, "ACGTTTAA");
        let expected = reverse_complement(&genomic_concat);
        assert_eq!(result, expected);
        assert_eq!(result, "TTAAACGT");
    }

    #[test]
    fn test_concat_regions_empty() {
        let chrom = "AAAA";
        let (store, digest) = make_seq_store("chrEmpty", chrom);
        let chrom_digest = key_to_chrom_digest(&digest);

        let result = concat_regions(&store, &chrom_digest, &[], Strand::Forward).unwrap();
        assert_eq!(result, "");
    }

    #[test]
    fn test_concat_regions_unknown_digest() {
        let chrom = "AAAA";
        let (store, _) = make_seq_store("chrX", chrom);
        let unknown_digest = [0u8; 24];

        let result = concat_regions(&store, &unknown_digest, &[(0, 2)], Strand::Forward);
        assert!(result.is_err(), "should error for unknown digest");
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("Sequence not found"), "error message should mention 'Sequence not found': {msg}");
    }

    #[test]
    fn test_mature_mrna_for_transcript_forward() {
        let chrom = "ATCGATCGATCGATCGATCGATCGATCGATCG";
        let (store, digest) = make_seq_store("chrFwd", chrom);
        let chrom_digest = key_to_chrom_digest(&digest);

        let tx = make_transcript(
            "NM_FWD.1",
            chrom_digest,
            Strand::Forward,
            vec![Exon { start: 2, end: 6 }, Exon { start: 10, end: 14 }],
        );
        let result = mature_mrna_for_transcript(&store, &tx).unwrap();
        assert_eq!(result, "CGATCGAT");
    }

    #[test]
    fn test_mature_mrna_for_transcript_reverse() {
        let chrom = "ACGTTTAAGGCCAACCGGTT";
        let (store, digest) = make_seq_store("chrRev2", chrom);
        let chrom_digest = key_to_chrom_digest(&digest);

        let tx = make_transcript(
            "NM_REV.1",
            chrom_digest,
            Strand::Reverse,
            vec![Exon { start: 0, end: 4 }, Exon { start: 4, end: 8 }],
        );
        let result = mature_mrna_for_transcript(&store, &tx).unwrap();
        assert_eq!(result, "TTAAACGT");
    }

    #[test]
    fn test_mature_mrna_accession_lookup() {
        let chrom = "ATCGATCGATCGATCGATCGATCGATCGATCG";
        let (store, digest) = make_seq_store("chrLookup", chrom);
        let chrom_digest = key_to_chrom_digest(&digest);

        let tx = make_transcript(
            "NM_LOOKUP.1",
            chrom_digest,
            Strand::Forward,
            vec![Exon { start: 2, end: 6 }, Exon { start: 10, end: 14 }],
        );
        let expected = mature_mrna_for_transcript(&store, &tx).unwrap();

        let bytes = build_reftx_bytes_in_memory(&[tx]).unwrap();
        let tx_store = ReadonlyTxStore::from_bytes(bytes).unwrap();

        let result = mature_mrna(&store, &tx_store, "NM_LOOKUP.1").unwrap();
        assert_eq!(result, expected);
        assert_eq!(result, "CGATCGAT");
    }

    #[test]
    fn test_mature_mrna_accession_not_found() {
        let chrom = "AAAA";
        let (store, digest) = make_seq_store("chrNF", chrom);
        let chrom_digest = key_to_chrom_digest(&digest);

        let tx = make_transcript(
            "NM_REAL.1",
            chrom_digest,
            Strand::Forward,
            vec![Exon { start: 0, end: 2 }],
        );
        let bytes = build_reftx_bytes_in_memory(&[tx]).unwrap();
        let tx_store = ReadonlyTxStore::from_bytes(bytes).unwrap();

        let result = mature_mrna(&store, &tx_store, "NM_MISSING.1");
        assert!(result.is_err(), "missing accession should return Err");
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("Transcript not found"), "error should mention 'Transcript not found': {msg}");
    }

    // Verify that the `_` (Arc) is needed for `make_seq_store` to not drop the
    // backing data. This test uses the `Arc`-backed in-memory store path.
    #[test]
    fn test_concat_regions_single_exon() {
        let chrom = "AAACCCGGGTTT";
        let (store, digest) = make_seq_store("chrSingle", chrom);
        let chrom_digest = key_to_chrom_digest(&digest);

        let result = concat_regions(&store, &chrom_digest, &[(3, 9)], Strand::Forward).unwrap();
        assert_eq!(result, "CCCGGG");
    }
}
