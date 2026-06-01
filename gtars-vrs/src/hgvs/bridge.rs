//! HGVS-to-VRS Allele bridge.
//!
//! Converts a parsed [`HgvsVariant`] into a VRS [`Allele`], resolving
//! transcript coordinates via a [`TranscriptProvider`] and reference
//! sequence via a [`RefgetStore`] / [`ReadonlyRefgetStore`]. The resulting
//! Allele is normalized and digested to the canonical
//! `ga4gh:VA.<digest>` identifier — matching what the equivalent VCF row
//! would produce through `vcf.rs`.
//!
//! # Worked example
//!
//! ```ignore
//! // BRAF V600E -> VRS ID
//! use gtars_vrs::hgvs_str_to_vrs_id;
//! let id = hgvs_str_to_vrs_id(
//!     "NM_004333.6:c.1799T>A",
//!     &provider,
//!     &mut refget,
//!     "ga4gh:SQ.collection_digest",
//! )?;
//! assert!(id.starts_with("ga4gh:VA."));
//! ```
//!
//! # Scope (v1)
//!
//! - **Reference types:** `g.`, `c.`, `n.`
//! - **Edits:** `Sub`, `Del`, `Dup`, `Ins`, `DelIns`, `Identity`
//! - **Rejected:** `m.`/`r.`/`p.` reference types; `Inv` and `Unknown` edits
//!
//! # Note for FFI
//!
//! Python bindings (sibling crate) accept owned `String` and call into the
//! parser's `HgvsVariantOwned` path. This module exposes only borrowed-AST
//! entry points; FFI converts as needed.

use std::collections::HashMap;

use gtars_refget::store::ReadonlyRefgetStore;
// The mutable store is used only by the filesystem-backed bridge entries.
#[cfg(feature = "filesystem")]
use gtars_refget::store::RefgetStore;
use thiserror::Error;

use crate::digest::DigestWriter;
use crate::hgvs::ast::{Datum, Edit, HgvsVariant, LocationRange, Position, ReferenceType};
use crate::hgvs::parser::parse;
use crate::models::{Allele, AlleleState, SequenceLocation, SequenceReference};
use crate::normalize::{
    NormalizeError, RefSeq, RefView, normalize_ref, ref_view_for as ref_view_for_inner,
};
use crate::provider::{ProviderError, TranscriptProvider};

/// Errors raised by the HGVS-to-VRS bridge.
#[derive(Debug, Error)]
pub enum BridgeError {
    #[error("parse error: {0}")]
    Parse(#[from] crate::hgvs::error::HgvsError),
    #[error("provider error: {0}")]
    Provider(#[from] ProviderError),
    #[error("refget error: {0}")]
    Refget(#[source] anyhow::Error),
    #[error("normalize error: {0}")]
    Normalize(#[from] NormalizeError),
    #[error("unsupported reference type: {0:?} (v1 supports g., c., n.)")]
    UnsupportedReferenceType(ReferenceType),
    #[error("unsupported edit: {0}")]
    UnsupportedEdit(&'static str),
    #[error("REF mismatch at {accession}:{pos}: HGVS says {hgvs_ref}, sequence has {actual_ref}")]
    RefMismatch {
        accession: String,
        pos: u64,
        hgvs_ref: String,
        actual_ref: String,
    },
    #[error("chromosome accession {0} not found in collection")]
    UnknownChrom(String),
    #[error("coordinate {pos} out of bounds for {accession} (len={seq_len})")]
    OutOfBounds {
        accession: String,
        pos: u64,
        seq_len: usize,
    },
    #[error("inconsistent edit: {0}")]
    InconsistentEdit(String),
}

// ── Public API ──────────────────────────────────────────────────────────

/// Build a [`RefView`] over a sequence's resident bytes, decoding on the fly
/// (thin wrapper around [`crate::normalize::ref_view_for`] that maps any failure
/// to [`BridgeError::UnknownChrom`] with the chromosome name for context).
fn ref_view_for<'a>(
    store: &'a ReadonlyRefgetStore,
    raw_digest: &str,
    chrom_name: &str,
) -> Result<RefView<'a>, BridgeError> {
    ref_view_for_inner(store, raw_digest)
        .map_err(|_| BridgeError::UnknownChrom(chrom_name.to_string()))
}

/// Convert a parsed HGVS variant into a (non-normalized) VRS [`Allele`].
///
/// Resolves transcript coordinates via `provider` and reads reference bases
/// directly from the encoded store (decoding on the fly). The returned Allele is
/// in genomic orientation, suitable for normalization and digesting.
#[cfg(feature = "filesystem")]
pub fn hgvs_to_allele(
    variant: &HgvsVariant<'_>,
    provider: &dyn TranscriptProvider,
    refget: &mut RefgetStore,
    collection_digest: &str,
) -> Result<Allele, BridgeError> {
    let name_to_digest =
        crate::vcf::build_name_to_digest(refget, collection_digest).map_err(BridgeError::Refget)?;
    let chrom_name = resolve_chrom_for_variant(variant, provider, &name_to_digest)?;
    let raw_digest = name_to_digest
        .get(&chrom_name)
        .ok_or_else(|| BridgeError::UnknownChrom(chrom_name.clone()))?
        .clone();
    crate::vcf::ensure_resident(refget, &raw_digest).map_err(BridgeError::Refget)?;
    let seq = ref_view_for(refget, &raw_digest, &chrom_name)?;
    let parts = build_allele_parts(variant, provider, &chrom_name, &raw_digest, &seq)?;
    Ok(parts.allele)
}

/// Parse + bridge + normalize + digest in one call.
///
/// Returns the canonical `ga4gh:VA.<digest>` identifier matching what the
/// equivalent VCF row would produce through `vcf.rs`.
#[cfg(feature = "filesystem")]
pub fn hgvs_str_to_vrs_id(
    s: &str,
    provider: &dyn TranscriptProvider,
    refget: &mut RefgetStore,
    collection_digest: &str,
) -> Result<String, BridgeError> {
    let variant = parse(s)?;
    let name_to_digest =
        crate::vcf::build_name_to_digest(refget, collection_digest).map_err(BridgeError::Refget)?;
    let chrom_name = resolve_chrom_for_variant(&variant, provider, &name_to_digest)?;
    let raw_digest = name_to_digest
        .get(&chrom_name)
        .ok_or_else(|| BridgeError::UnknownChrom(chrom_name.clone()))?
        .clone();
    crate::vcf::ensure_resident(refget, &raw_digest).map_err(BridgeError::Refget)?;
    let seq = ref_view_for(refget, &raw_digest, &chrom_name)?;
    let parts = build_allele_parts(&variant, provider, &chrom_name, &raw_digest, &seq)?;
    finalize_vrs_id(&parts, &seq)
}

/// Read-only store variant of [`hgvs_to_allele`]. Referenced sequences must be
/// resident (Full/encoded) in the store; reference bases are decoded on the fly.
pub fn hgvs_to_allele_readonly(
    variant: &HgvsVariant<'_>,
    provider: &dyn TranscriptProvider,
    refget: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
) -> Result<Allele, BridgeError> {
    let chrom_name = resolve_chrom_for_variant(variant, provider, name_to_digest)?;
    let raw_digest = name_to_digest
        .get(&chrom_name)
        .ok_or_else(|| BridgeError::UnknownChrom(chrom_name.clone()))?;
    let seq = ref_view_for(refget, raw_digest, &chrom_name)?;
    let parts = build_allele_parts(variant, provider, &chrom_name, raw_digest, &seq)?;
    Ok(parts.allele)
}

/// Read-only store variant of [`hgvs_str_to_vrs_id`].
pub fn hgvs_str_to_vrs_id_readonly(
    s: &str,
    provider: &dyn TranscriptProvider,
    refget: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
) -> Result<String, BridgeError> {
    let variant = parse(s)?;
    let chrom_name = resolve_chrom_for_variant(&variant, provider, name_to_digest)?;
    let raw_digest = name_to_digest
        .get(&chrom_name)
        .ok_or_else(|| BridgeError::UnknownChrom(chrom_name.clone()))?;
    let seq = ref_view_for(refget, raw_digest, &chrom_name)?;
    let parts = build_allele_parts(&variant, provider, &chrom_name, raw_digest, &seq)?;
    finalize_vrs_id(&parts, &seq)
}

// ── Internals ───────────────────────────────────────────────────────────

/// Output of [`build_allele_parts`]: the (non-normalized) Allele plus the
/// owned ref/alt bytes needed for downstream normalization.
struct AlleleParts {
    allele: Allele,
    /// 0-based interbase start used as the normalize input.
    start_ib: u64,
    ref_bytes: Vec<u8>,
    alt_bytes: Vec<u8>,
    /// `SQ.<digest>` for the chromosome.
    refget_accession: String,
}

/// True if the accession looks like a bare gene symbol (not an NCBI / chrom
/// accession), for the gene-symbol → MANE resolution path.
pub(crate) fn looks_like_gene_symbol(accession: &str) -> bool {
    if accession.contains('.') {
        return false;
    }
    let blocked_prefixes = [
        "NC_", "NM_", "NR_", "NG_", "NW_", "NT_", "XM_", "XR_", "ENST", "ENSG", "chr", "GL",
        "KI", "MT",
    ];
    !blocked_prefixes.iter().any(|p| accession.starts_with(p))
}

/// Standard A/T/C/G/N reverse complement. Errors on any other byte.
fn revcomp(bytes: &[u8]) -> Result<Vec<u8>, BridgeError> {
    let mut out = Vec::with_capacity(bytes.len());
    for &b in bytes.iter().rev() {
        let c = match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'N' => b'N',
            _ => {
                return Err(BridgeError::InconsistentEdit(format!(
                    "non-ACGTN base {:?} in HGVS allele",
                    b as char
                )));
            }
        };
        out.push(c);
    }
    Ok(out)
}

fn revcomp_if_neg(bytes: &[u8], strand: i8) -> Result<Vec<u8>, BridgeError> {
    if strand < 0 { revcomp(bytes) } else { Ok(bytes.to_vec()) }
}

/// Resolve the chrom-side refget-name (key into `name_to_digest`) for a
/// variant.
///
/// For `g.` variants whose accession is itself a chromosome accession,
/// the accession is returned unchanged (after MANE fallback for gene
/// symbols). For `c.`/`n.` variants, asks the provider for the
/// chromosome's refget accession (`SQ.<digest>`) and matches it back to
/// a name in `name_to_digest`.
fn resolve_chrom_for_variant(
    variant: &HgvsVariant<'_>,
    provider: &dyn TranscriptProvider,
    name_to_digest: &HashMap<String, String>,
) -> Result<String, BridgeError> {
    let accession = if looks_like_gene_symbol(variant.accession) {
        provider
            .gene_to_mane_accession(variant.accession)
            .ok_or_else(|| ProviderError::NoManeTranscript(variant.accession.to_string()))?
    } else {
        variant.accession.to_string()
    };

    match variant.reference_type {
        ReferenceType::G => {
            if name_to_digest.contains_key(&accession) {
                Ok(accession)
            } else {
                Err(BridgeError::UnknownChrom(accession))
            }
        }
        ReferenceType::C | ReferenceType::N => {
            let sq = provider.get_chrom_accession(&accession)?;
            let raw = sq.strip_prefix("SQ.").unwrap_or(&sq);
            // Find the collection name whose digest matches.
            for (name, digest) in name_to_digest.iter() {
                if digest == raw {
                    return Ok(name.clone());
                }
            }
            // Fall back to direct key lookups (in case caller registered
            // collection by SQ.* form).
            if name_to_digest.contains_key(&sq) {
                Ok(sq)
            } else if name_to_digest.contains_key(raw) {
                Ok(raw.to_string())
            } else {
                Err(BridgeError::UnknownChrom(sq))
            }
        }
        rt => Err(BridgeError::UnsupportedReferenceType(rt)),
    }
}

/// Map a single transcript-relative or genomic [`Position`] to a 0-based
/// interbase coordinate on the chromosome via the provider. Returns
/// `(interbase_position, strand)`.
fn position_to_genomic_interbase(
    pos: &Position,
    accession: &str,
    reference_type: ReferenceType,
    provider: &dyn TranscriptProvider,
) -> Result<(u64, i8), BridgeError> {
    match reference_type {
        ReferenceType::G => {
            if pos.base < 1 {
                return Err(BridgeError::InconsistentEdit(format!(
                    "g. position must be >= 1, got {}",
                    pos.base
                )));
            }
            Ok(((pos.base - 1) as u64, 1))
        }
        ReferenceType::C => {
            let is_cds_end = matches!(pos.datum, Datum::CdsEnd);
            let loc = provider.c_to_genomic_full(accession, pos.base, pos.offset, is_cds_end)?;
            let strand = provider.get_strand(accession)?;
            Ok((loc.start, strand))
        }
        ReferenceType::N => {
            let loc = provider.n_to_genomic_full(accession, pos.base, pos.offset)?;
            let strand = provider.get_strand(accession)?;
            Ok((loc.start, strand))
        }
        rt => Err(BridgeError::UnsupportedReferenceType(rt)),
    }
}

/// Build the (Allele, normalize-input) parts for an HGVS variant. The
/// returned `start_ib`, `ref_bytes`, and `alt_bytes` are all in genomic
/// orientation and ready to feed `normalize()`.
fn build_allele_parts<R: RefSeq + ?Sized>(
    variant: &HgvsVariant<'_>,
    provider: &dyn TranscriptProvider,
    chrom_name: &str,
    raw_digest: &str,
    seq: &R,
) -> Result<AlleleParts, BridgeError> {
    if variant.posedit.uncertain {
        // tracing crate isn't a workspace dep; eprintln keeps us
        // dependency-free without silently swallowing the warning.
        eprintln!(
            "warning: bridging uncertain HGVS expression — VRS digest reflects central coordinate only"
        );
    }

    let accession = if looks_like_gene_symbol(variant.accession) {
        provider
            .gene_to_mane_accession(variant.accession)
            .ok_or_else(|| ProviderError::NoManeTranscript(variant.accession.to_string()))?
    } else {
        variant.accession.to_string()
    };

    let (start_ib, end_ib, strand) = range_and_edit_to_genomic(
        &variant.posedit.pos,
        &variant.posedit.edit,
        &accession,
        variant.reference_type,
        provider,
    )?;

    if end_ib as usize > seq.len() {
        return Err(BridgeError::OutOfBounds {
            accession: chrom_name.to_string(),
            pos: end_ib,
            seq_len: seq.len(),
        });
    }
    if start_ib > end_ib {
        return Err(BridgeError::InconsistentEdit(format!(
            "start_ib {} > end_ib {}",
            start_ib, end_ib
        )));
    }

    let mut actual_ref = Vec::with_capacity((end_ib - start_ib) as usize);
    seq.extend_range(start_ib as usize, end_ib as usize, &mut actual_ref);

    let alt_bytes = compute_alt(&variant.posedit.edit, &actual_ref, strand)?;

    // REF cross-check if the parser supplied a reference.
    if let Some(parsed_ref) = edit_reference(&variant.posedit.edit) {
        let expected = revcomp_if_neg(parsed_ref.as_bytes(), strand)?;
        if expected != actual_ref {
            return Err(BridgeError::RefMismatch {
                accession: chrom_name.to_string(),
                pos: start_ib,
                hgvs_ref: parsed_ref.to_string(),
                actual_ref: String::from_utf8_lossy(&actual_ref).to_string(),
            });
        }
    }

    let refget_accession = format!("SQ.{}", raw_digest);
    let alt_str = std::str::from_utf8(&alt_bytes)
        .map_err(|_| BridgeError::InconsistentEdit("non-UTF8 alt bytes".to_string()))?
        .to_string();
    let allele = Allele {
        location: SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: refget_accession.clone(),
            },
            start: start_ib,
            end: end_ib,
        },
        state: AlleleState::LiteralSequenceExpression { sequence: alt_str },
    };

    Ok(AlleleParts {
        allele,
        start_ib,
        ref_bytes: actual_ref,
        alt_bytes,
        refget_accession,
    })
}

/// Translate (LocationRange, Edit) into genomic-orientation
/// `(start_ib, end_ib, strand)`. Handles single-position edits whose
/// interbase span depends on the edit type (e.g. `Sub` is 1bp, `Ins` is
/// 0bp between two positions, `Del` is the parsed range).
fn range_and_edit_to_genomic(
    pos: &LocationRange,
    edit: &Edit<'_>,
    accession: &str,
    reference_type: ReferenceType,
    provider: &dyn TranscriptProvider,
) -> Result<(u64, u64, i8), BridgeError> {
    match (pos, edit) {
        // Single-position substitution / identity / 1bp del/dup/delins:
        // 1bp interbase span at p.
        (
            LocationRange::Single(p),
            Edit::Sub { .. }
            | Edit::Identity
            | Edit::Del { .. }
            | Edit::Dup { .. }
            | Edit::DelIns { .. },
        ) => {
            let (ib, strand) =
                position_to_genomic_interbase(p, accession, reference_type, provider)?;
            Ok((ib, ib + 1, strand))
        }
        // HGVS `c.123insAT` means insert after position 123, before 124.
        (LocationRange::Single(p), Edit::Ins { .. }) => {
            let (ib, strand) =
                position_to_genomic_interbase(p, accession, reference_type, provider)?;
            // +strand: insertion site is between ib and ib+1 → [ib+1, ib+1).
            // -strand: the genomic insertion site is between ib-1 and ib → [ib, ib).
            if strand >= 0 {
                Ok((ib + 1, ib + 1, strand))
            } else {
                Ok((ib, ib, strand))
            }
        }
        (LocationRange::Range { start, end }, Edit::Ins { .. }) => {
            let (a, strand) =
                position_to_genomic_interbase(start, accession, reference_type, provider)?;
            let (b, _) = position_to_genomic_interbase(end, accession, reference_type, provider)?;
            let lo = a.min(b);
            let hi = a.max(b);
            // Insertion between two positions: HGVS pair must be adjacent.
            if hi - lo != 1 {
                return Err(BridgeError::InconsistentEdit(format!(
                    "ins range positions are not adjacent: {} and {}",
                    a, b
                )));
            }
            Ok((hi, hi, strand))
        }
        // Range over Del/Dup/DelIns/Identity/Sub: span is [min, max+1).
        (
            LocationRange::Range { start, end },
            Edit::Del { .. }
            | Edit::Dup { .. }
            | Edit::DelIns { .. }
            | Edit::Identity
            | Edit::Sub { .. },
        ) => {
            let (a, strand) =
                position_to_genomic_interbase(start, accession, reference_type, provider)?;
            let (b, _) = position_to_genomic_interbase(end, accession, reference_type, provider)?;
            let lo = a.min(b);
            let hi = a.max(b);
            Ok((lo, hi + 1, strand))
        }
        (_, Edit::Inv { .. }) => Err(BridgeError::UnsupportedEdit("inv")),
        (_, Edit::Unknown) => Err(BridgeError::UnsupportedEdit("unknown")),
        (_, Edit::Copy { .. }) => Err(BridgeError::UnsupportedEdit("copy")),
        (_, Edit::Repeat { .. }) => Err(BridgeError::UnsupportedEdit("repeat")),
        (LocationRange::WholeSequence, _) => {
            Err(BridgeError::UnsupportedEdit("whole-sequence edit"))
        }
        (LocationRange::UncertainStart { .. }, _)
        | (LocationRange::UncertainEnd { .. }, _)
        | (LocationRange::UncertainBoth { .. }, _) => {
            Err(BridgeError::UnsupportedEdit("uncertain position range"))
        }
    }
}

/// Compute the genomic-orientation alt bytes for an edit, given the
/// genomic-orientation actual REF bytes already fetched from the store.
fn compute_alt(edit: &Edit<'_>, actual_ref: &[u8], strand: i8) -> Result<Vec<u8>, BridgeError> {
    match edit {
        Edit::Sub { alternate, .. } => revcomp_if_neg(alternate.as_bytes(), strand),
        Edit::Del { .. } => Ok(Vec::new()),
        Edit::Ins { alternate } => revcomp_if_neg(alternate.as_bytes(), strand),
        Edit::Dup { .. } => {
            // Genomic ALT = actual_ref repeated twice.
            let mut v = Vec::with_capacity(actual_ref.len() * 2);
            v.extend_from_slice(actual_ref);
            v.extend_from_slice(actual_ref);
            Ok(v)
        }
        Edit::DelIns { alternate, .. } => revcomp_if_neg(alternate.as_bytes(), strand),
        Edit::Identity => Ok(actual_ref.to_vec()),
        Edit::Inv { .. } => Err(BridgeError::UnsupportedEdit("inv")),
        Edit::Unknown => Err(BridgeError::UnsupportedEdit("unknown")),
        Edit::Copy { .. } => Err(BridgeError::UnsupportedEdit("copy")),
        Edit::Repeat { .. } => Err(BridgeError::UnsupportedEdit("repeat")),
    }
}

/// Extract the parser-supplied reference allele (transcript orientation)
/// for cross-check, if present.
fn edit_reference<'a>(edit: &'a Edit<'_>) -> Option<&'a str> {
    match edit {
        Edit::Sub { reference, .. } => Some(*reference),
        Edit::Del { reference } => *reference,
        Edit::Dup { reference } => *reference,
        Edit::DelIns { reference, .. } => *reference,
        Edit::Identity
        | Edit::Ins { .. }
        | Edit::Inv { .. }
        | Edit::Unknown
        | Edit::Copy { .. }
        | Edit::Repeat { .. } => None,
    }
}

/// Normalize the parts and digest into the final `ga4gh:VA.<digest>` ID,
/// reusing the already-fetched chromosome `seq` slice.
fn finalize_vrs_id<R: RefSeq + ?Sized>(parts: &AlleleParts, seq: &R) -> Result<String, BridgeError> {
    let norm = normalize_ref(seq, parts.start_ib, &parts.ref_bytes, &parts.alt_bytes)?;
    let norm_seq = std::str::from_utf8(&norm.allele).map_err(|_| {
        BridgeError::InconsistentEdit("normalized allele is not valid UTF-8".to_string())
    })?;
    let mut writer = DigestWriter::new();
    Ok(writer.allele_identifier_literal(&parts.refget_accession, norm.start, norm.end, norm_seq))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_looks_like_gene_symbol() {
        assert!(looks_like_gene_symbol("BRAF"));
        assert!(looks_like_gene_symbol("TP53"));
        assert!(!looks_like_gene_symbol("NM_004333.6"));
        assert!(!looks_like_gene_symbol("NC_000007.14"));
        assert!(!looks_like_gene_symbol("chr7"));
        assert!(!looks_like_gene_symbol("MT"));
        assert!(!looks_like_gene_symbol("ENST00000288602"));
    }

    /// Golden VRS id for `chrF:g.6C>T` over the synthetic `chrF` contig. This is
    /// the cross-surface parity anchor: the wasm `hgvs_to_vrs_id` binding
    /// (`gtars-wasm/src/hgvs.rs`) pins the SAME constant, so a divergence between
    /// the native and wasm code paths fails one side or the other.
    pub(crate) const GOLDEN_CHRF_G6CT_VRS_ID: &str = "ga4gh:VA._q-idtHGQxQ4XiEPJ1ExYl_htUeNEkir";

    #[test]
    fn readonly_bridge_matches_golden_vrs_id() {
        use crate::provider::NoTranscriptProvider;
        use gtars_refget::digest::digest_sequence;
        use gtars_refget::store::RefgetStore;
        use std::collections::HashMap;

        // Replicate exactly what the wasm `hgvs_to_vrs_id` wrapper does, so this
        // native result is the reference the wasm binding is asserted against.
        let name = "chrF";
        let bases = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let record = digest_sequence(name, bases.as_bytes());
        let raw_digest = record.metadata().sha512t24u.clone();
        let mut store = RefgetStore::in_memory();
        store.add_sequence_record(record, true).unwrap();
        let mut name_to_digest = HashMap::new();
        name_to_digest.insert(name.to_string(), raw_digest);
        let id = hgvs_str_to_vrs_id_readonly(
            "chrF:g.6C>T",
            &NoTranscriptProvider,
            &store,
            &name_to_digest,
        )
        .unwrap();
        assert_eq!(id, GOLDEN_CHRF_G6CT_VRS_ID);
    }

    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp(b"A").unwrap(), b"T");
        assert_eq!(revcomp(b"ACGT").unwrap(), b"ACGT");
        assert_eq!(revcomp(b"GCG").unwrap(), b"CGC");
        assert_eq!(revcomp(b"N").unwrap(), b"N");
        assert!(revcomp(b"X").is_err());
    }

    #[test]
    fn test_revcomp_if_neg() {
        assert_eq!(revcomp_if_neg(b"AT", 1).unwrap(), b"AT");
        assert_eq!(revcomp_if_neg(b"AT", -1).unwrap(), b"AT");
        assert_eq!(revcomp_if_neg(b"AC", -1).unwrap(), b"GT");
    }
}
