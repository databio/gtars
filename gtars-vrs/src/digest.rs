//! VRS digest computation.
//!
//! Computes GA4GH VRS digests by canonical JSON serialization + SHA-512/24u.
//!
//! The "fast" path (`DigestWriter`) writes canonical JSON directly into a
//! reusable buffer without serde_json::Value allocation, and uses a
//! stack-allocated SHA-512 hasher. For the VRS VCF hot loop.
//!
//! The generic path (`allele_identifier`) uses serde_json for correctness
//! reference and tests.

use sha2::{Digest, Sha512};

use crate::models::{Allele, AlleleState, SequenceLocation};

/// Reusable digest writer that avoids per-call allocations.
///
/// Holds a scratch buffer for canonical JSON bytes and computes SHA-512/24u
/// digests without heap allocation (hasher is stack-allocated, buffer is reused).
pub struct DigestWriter {
    buf: Vec<u8>,
}

impl DigestWriter {
    pub fn new() -> Self {
        Self {
            buf: Vec::with_capacity(512),
        }
    }

    /// Compute `ga4gh:VA.<digest>` for a VRS Allele with a LiteralSequenceExpression state.
    ///
    /// This is the hot-path version that writes canonical JSON directly into a
    /// reusable buffer. It does NOT support `ReferenceLengthExpression`. For that,
    /// use the generic `allele_identifier()` free function.
    ///
    /// # JSON escaping
    ///
    /// This method writes `refget_accession` and `sequence` as raw bytes without
    /// JSON escaping. This is safe by construction: `refget_accession` is always
    /// a GA4GH refget digest (`SQ.` + base64url characters `[A-Za-z0-9_-]`), and
    /// `sequence` is always uppercase nucleotide characters (`[ACGTN]`). Neither
    /// alphabet contains JSON-special characters (`"`, `\`, control chars).
    /// The `test_fast_path_matches_generic` test verifies equivalence with the
    /// serde_json path.
    pub fn allele_identifier_literal(
        &mut self,
        refget_accession: &str,
        start: u64,
        end: u64,
        sequence: &str,
    ) -> String {
        // Step 1: Compute SequenceLocation digest
        // Canonical JSON (keys sorted): {"end":N,"sequenceReference":{"refgetAccession":"...","type":"SequenceReference"},"start":N,"type":"SequenceLocation"}
        self.buf.clear();
        self.buf.extend_from_slice(b"{\"end\":");
        itoa::Buffer::new().format(end).as_bytes().iter().for_each(|&b| self.buf.push(b));
        self.buf.extend_from_slice(b",\"sequenceReference\":{\"refgetAccession\":\"");
        self.buf.extend_from_slice(refget_accession.as_bytes());
        self.buf.extend_from_slice(b"\",\"type\":\"SequenceReference\"},\"start\":");
        itoa::Buffer::new().format(start).as_bytes().iter().for_each(|&b| self.buf.push(b));
        self.buf.extend_from_slice(b",\"type\":\"SequenceLocation\"}");

        let sl_digest = sha512t24u_inline(&self.buf);

        // Step 2: Compute Allele digest
        // Canonical JSON: {"location":"<sl_digest>","state":{"sequence":"...","type":"LiteralSequenceExpression"},"type":"Allele"}
        self.buf.clear();
        self.buf.extend_from_slice(b"{\"location\":\"");
        self.buf.extend_from_slice(sl_digest.as_bytes());
        self.buf.extend_from_slice(b"\",\"state\":{\"sequence\":\"");
        self.buf.extend_from_slice(sequence.as_bytes());
        self.buf.extend_from_slice(b"\",\"type\":\"LiteralSequenceExpression\"},\"type\":\"Allele\"}");

        let allele_digest = sha512t24u_inline(&self.buf);

        format!("ga4gh:VA.{}", allele_digest)
    }
}

/// SHA-512 truncated to 24 bytes, base64url-encoded. Stack-allocated hasher.
#[inline]
fn sha512t24u_inline(data: &[u8]) -> String {
    let mut hasher = Sha512::new();
    hasher.update(data);
    let hash = hasher.finalize();
    base64_url::encode(&hash[..24])
}

// === Generic path (used by tests and non-hot-path callers) ===

use gtars_refget::digest::algorithms::{canonicalize_json, sha512t24u};
use serde_json::json;

/// Compute the GA4GH digest for a SequenceLocation.
pub fn sequence_location_digest(loc: &SequenceLocation) -> String {
    let json_val = json!({
        "end": loc.end,
        "sequenceReference": {
            "refgetAccession": loc.sequence_reference.refget_accession,
            "type": "SequenceReference"
        },
        "start": loc.start,
        "type": "SequenceLocation"
    });
    sha512t24u(canonicalize_json(&json_val).as_bytes())
}

/// Compute the GA4GH digest for an Allele.
pub fn allele_digest(allele: &Allele) -> String {
    let sl_digest = sequence_location_digest(&allele.location);
    let state_json = match &allele.state {
        AlleleState::LiteralSequenceExpression { sequence } => json!({
            "sequence": sequence,
            "type": "LiteralSequenceExpression"
        }),
        AlleleState::ReferenceLengthExpression {
            length,
            repeat_subunit_length,
            sequence,
        } => {
            let mut obj = json!({
                "length": length,
                "repeatSubunitLength": repeat_subunit_length,
                "type": "ReferenceLengthExpression"
            });
            if let Some(seq) = sequence {
                obj["sequence"] = serde_json::Value::String(seq.clone());
            }
            obj
        }
    };
    let json_val = json!({
        "location": sl_digest,
        "state": state_json,
        "type": "Allele"
    });
    sha512t24u(canonicalize_json(&json_val).as_bytes())
}

/// Compute the full GA4GH VRS identifier for an Allele.
pub fn allele_identifier(allele: &Allele) -> String {
    format!("ga4gh:VA.{}", allele_digest(allele))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::{SequenceLocation, SequenceReference};

    #[test]
    fn test_sequence_location_digest_deterministic() {
        let loc = SequenceLocation {
            sequence_reference: SequenceReference {
                refget_accession: "SQ.F-LrLnMKIjgbR1HECnsl_VGjXfs3QHDE".to_string(),
            },
            start: 55181319,
            end: 55181320,
        };
        let d1 = sequence_location_digest(&loc);
        let d2 = sequence_location_digest(&loc);
        assert_eq!(d1, d2);
        assert_eq!(d1.len(), 32);
    }

    #[test]
    fn test_allele_identifier_format() {
        let allele = Allele {
            location: SequenceLocation {
                sequence_reference: SequenceReference {
                    refget_accession: "SQ.F-LrLnMKIjgbR1HECnsl_VGjXfs3QHDE".to_string(),
                },
                start: 55181319,
                end: 55181320,
            },
            state: AlleleState::LiteralSequenceExpression {
                sequence: "T".to_string(),
            },
        };
        let id = allele_identifier(&allele);
        assert!(id.starts_with("ga4gh:VA."));
        assert_eq!(id.len(), 9 + 32);
    }

    /// Verify the fast path produces identical results to the generic path.
    #[test]
    fn test_fast_path_matches_generic() {
        let allele = Allele {
            location: SequenceLocation {
                sequence_reference: SequenceReference {
                    refget_accession: "SQ.F-LrLnMKIjgbR1HECnsl_VGjXfs3QHDE".to_string(),
                },
                start: 55181319,
                end: 55181320,
            },
            state: AlleleState::LiteralSequenceExpression {
                sequence: "T".to_string(),
            },
        };

        let generic = allele_identifier(&allele);

        let mut writer = DigestWriter::new();
        let fast = writer.allele_identifier_literal(
            "SQ.F-LrLnMKIjgbR1HECnsl_VGjXfs3QHDE",
            55181319,
            55181320,
            "T",
        );

        assert_eq!(generic, fast);
    }
}
