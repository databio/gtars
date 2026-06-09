//! WASM bindings for HGVS → VRS allele identifier computation.
//!
//! Exposes browser-callable entry points that convert an HGVS string into a
//! canonical GA4GH VRS allele identifier (`ga4gh:VA.<digest>`).
//!
//! Reference bases come from JavaScript (there is no filesystem in the browser):
//! JS fetches the bases for the referenced sequence and hands them in. The Rust
//! entry point is synchronous over already-resident bytes — JS owns the async
//! HTTP (fetch / `gloo-net`).
//!
//! # Reference types
//!
//! - Genomic (`g.`) variants resolve with **no transcript store**: use the bare
//!   [`hgvs_to_vrs_id`], which passes a [`NoTranscriptProvider`].
//! - Transcript (`c.` / `n.`) variants — and gene-symbol → MANE lookups —
//!   resolve via the transcript-aware entries
//!   ([`hgvs_to_vrs_id_with_transcripts`] and
//!   [`RefgetStore::hgvs_to_vrs_id_with_transcripts`](crate::refget::RefgetStore)),
//!   given a [`TranscriptStore`](crate::transcripts::TranscriptStore) built from
//!   `.reftx` bytes. The transcript store is the WASM-safe in-memory backend
//!   wrapped in the real `TxProvider`; nothing here touches `memmap2` or the
//!   filesystem.
//!
//! # Example
//!
//! ```javascript
//! import init, {
//!     hgvs_to_vrs_id, hgvs_to_vrs_id_with_transcripts, parse_hgvs, TranscriptStore,
//! } from "gtars-js";
//! await init();
//!
//! // Fetch the full bases for the referenced sequence (small contig or
//! // a correctly-positioned window covering the whole sequence).
//! const bases = await (await fetch(seqUrl)).text();
//! const id = hgvs_to_vrs_id("NC_000001.11:g.100A>T", "NC_000001.11", bases);
//! // -> "ga4gh:VA.<digest>"
//!
//! // c./n. variants: load a transcript store once, then resolve against it.
//! const tx = new TranscriptStore(reftxBytes);
//! const cid = hgvs_to_vrs_id_with_transcripts(
//!     "NM_004333.6:c.1799T>A", "chr7", chr7Bases, tx);
//! ```

use std::collections::HashMap;

use gtars_refget::digest::digest_sequence;
use gtars_refget::store::RefgetStore;
use gtars_vrs::hgvs::ast::ReferenceType;
use gtars_vrs::hgvs::bridge::hgvs_str_to_vrs_id_readonly;
use gtars_vrs::hgvs::parser::parse;
use gtars_vrs::provider::NoTranscriptProvider;

use serde::Serialize;
use wasm_bindgen::prelude::*;

use crate::transcripts::TranscriptStore;

/// Build a one-off single-sequence in-memory `RefgetStore` from raw bases plus
/// the matching `name → raw-digest` map, exactly as the one-shot entries need.
/// Shared by [`hgvs_to_vrs_id`] and [`hgvs_to_vrs_id_with_transcripts`].
fn build_single_seq_store(
    sequence_name: &str,
    sequence_bases: &str,
) -> Result<(RefgetStore, HashMap<String, String>), JsValue> {
    // `digest_sequence` computes the sha512t24u/md5 digests and stores the raw
    // bases, which the readonly bridge decodes on the fly.
    let record = digest_sequence(sequence_name, sequence_bases.as_bytes());
    let raw_digest = record.metadata().sha512t24u.clone();

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_record(record, true)
        .map_err(|e| JsValue::from_str(&format!("failed to add sequence: {e}")))?;

    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    name_to_digest.insert(sequence_name.to_string(), raw_digest);
    Ok((store, name_to_digest))
}

/// Convert a genomic (`g.`) HGVS string into a GA4GH VRS allele identifier.
///
/// This is the zero-transcript convenience path: it uses a
/// [`NoTranscriptProvider`], so only `g.` resolves. For transcript-relative
/// (`c.`/`n.`) variants use [`hgvs_to_vrs_id_with_transcripts`] with a
/// [`TranscriptStore`](crate::transcripts::TranscriptStore) built from `.reftx`
/// bytes.
///
/// # Arguments
/// * `hgvs` - The HGVS expression, e.g. `"NC_000001.11:g.100A>T"`. Must be a
///   genomic (`g.`) variant; `c.`/`n.` are rejected here (no transcript store).
/// * `sequence_name` - The name/accession the HGVS expression references (e.g.
///   `"NC_000001.11"`). This MUST match the accession in `hgvs`.
/// * `sequence_bases` - The reference bases for `sequence_name`. To produce a
///   canonical VRS id these must be the full bases of the referenced sequence,
///   so the computed refget digest matches the true reference.
///
/// # Returns
/// The canonical `ga4gh:VA.<digest>` identifier, or a `JsValue` error string.
#[wasm_bindgen]
pub fn hgvs_to_vrs_id(
    hgvs: &str,
    sequence_name: &str,
    sequence_bases: &str,
) -> Result<String, JsValue> {
    let (store, name_to_digest) = build_single_seq_store(sequence_name, sequence_bases)?;
    // The readonly bridge resolves a `g.` variant's accession against this map
    // and looks the resulting raw digest up in the store.
    hgvs_str_to_vrs_id_readonly(
        hgvs,
        &NoTranscriptProvider,
        &store,
        &name_to_digest,
        &mut Vec::new(),
    )
    .map_err(|e| JsValue::from_str(&format!("{e}")))
}

/// Convert an HGVS string into a GA4GH VRS allele identifier, resolving
/// transcript-relative (`c.`/`n.`) coordinates and gene-symbol → MANE lookups
/// against a [`TranscriptStore`](crate::transcripts::TranscriptStore).
///
/// Builds the in-memory `RefgetStore` from `sequence_bases` exactly like
/// [`hgvs_to_vrs_id`], but passes the real [`TxProvider`] from `transcripts`
/// instead of [`NoTranscriptProvider`]. The transcript provider maps the
/// `c.`/`n.` coordinate to a genomic interbase position on the chromosome the
/// transcript lives on; the bridge then reads the reference bases for that
/// chromosome from the `RefgetStore`.
///
/// # Resident-chromosome requirement
/// The chromosome the transcript maps onto (the genome contig the `c.`/`n.`
/// variant lands on) MUST be resident in the `RefgetStore` for the variant to
/// resolve. Here that means `sequence_name`/`sequence_bases` must be that
/// chromosome. For multi-chromosome / genome-scale use, prefer
/// [`RefgetStore::hgvs_to_vrs_id_with_transcripts`](crate::refget::RefgetStore),
/// which resolves against a resident genome holding every contig.
///
/// # Arguments
/// * `hgvs` - The HGVS expression, e.g. `"NM_004333.6:c.1799T>A"` or a `g.`/`n.`.
/// * `sequence_name` - The chromosome name the transcript maps onto (e.g.
///   `"chr7"`).
/// * `sequence_bases` - The reference bases for that chromosome.
/// * `transcripts` - A loaded [`TranscriptStore`](crate::transcripts::TranscriptStore).
#[wasm_bindgen]
pub fn hgvs_to_vrs_id_with_transcripts(
    hgvs: &str,
    sequence_name: &str,
    sequence_bases: &str,
    transcripts: &TranscriptStore,
) -> Result<String, JsValue> {
    let (store, name_to_digest) = build_single_seq_store(sequence_name, sequence_bases)?;
    hgvs_str_to_vrs_id_readonly(
        hgvs,
        transcripts.provider(),
        &store,
        &name_to_digest,
        &mut Vec::new(),
    )
    .map_err(|e| JsValue::from_str(&format!("{e}")))
}

/// Parse-only summary of an HGVS expression, mirroring Python's `parse_hgvs`.
#[derive(Serialize)]
struct ParsedHgvs {
    /// The reference accession (e.g. `"NC_000001.11"`).
    accession: String,
    /// Gene symbol, if present in the expression.
    gene: Option<String>,
    /// Reference type: `"g"`, `"c"`, `"n"`, `"m"`, `"r"`, or `"p"`.
    reference_type: String,
}

/// Parse an HGVS string, returning a small structured summary (accession, gene,
/// reference type) or throwing on a parse error. Pure; no network or store.
#[wasm_bindgen]
pub fn parse_hgvs(hgvs: &str) -> Result<JsValue, JsValue> {
    let variant = parse(hgvs).map_err(|e| JsValue::from_str(&format!("{e}")))?;
    let parsed = ParsedHgvs {
        accession: variant.accession.to_string(),
        gene: variant.gene.map(|g| g.to_string()),
        reference_type: reference_type_str(variant.reference_type).to_string(),
    };
    serde_wasm_bindgen::to_value(&parsed).map_err(|e| JsValue::from_str(&format!("{e}")))
}

fn reference_type_str(rt: ReferenceType) -> &'static str {
    match rt {
        ReferenceType::G => "g",
        ReferenceType::C => "c",
        ReferenceType::N => "n",
        ReferenceType::M => "m",
        ReferenceType::R => "r",
        ReferenceType::P => "p",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;
    use wasm_bindgen_test::*;

    // Forward-strand synthetic contig; "chr"-prefixed so the bridge treats the
    // name as an accession rather than a gene symbol.
    const CHR_F_NAME: &str = "chrF";
    const CHR_F_BASES: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    #[derive(Deserialize)]
    struct ParsedOut {
        accession: String,
        gene: Option<String>,
        reference_type: String,
    }

    // Cross-surface parity anchor. This MUST equal
    // `gtars_vrs::hgvs::bridge::tests::GOLDEN_CHRF_G6CT_VRS_ID`, the id the native
    // (non-wasm) readonly bridge produces for the same inputs. If the wasm and
    // native code paths ever diverge, exactly one of the two tests fails.
    const GOLDEN_CHRF_G6CT_VRS_ID: &str = "ga4gh:VA._q-idtHGQxQ4XiEPJ1ExYl_htUeNEkir";

    #[wasm_bindgen_test]
    fn g_substitution_matches_native_vrs_id() {
        // chrF[5] (1-based pos 6) = 'C'. Sub C>T.
        let id = hgvs_to_vrs_id("chrF:g.6C>T", CHR_F_NAME, CHR_F_BASES)
            .expect("conversion should succeed");
        // Not just well-formed (`ga4gh:VA.`) — byte-identical to the native result.
        assert_eq!(id, GOLDEN_CHRF_G6CT_VRS_ID);
    }

    #[wasm_bindgen_test]
    fn bare_entry_rejects_transcript_variant() {
        // The bare `hgvs_to_vrs_id` uses NoTranscriptProvider, so c. is rejected
        // (no transcript store). The transcript-aware path is covered by
        // `transcripts::tests::c_variant_resolves_via_transcript_store`, which
        // RESOLVES a c. input through `hgvs_to_vrs_id_with_transcripts`.
        let res = hgvs_to_vrs_id("NM_000546.6:c.215C>G", CHR_F_NAME, CHR_F_BASES);
        assert!(
            res.is_err(),
            "c. variant should be rejected without a transcript store"
        );
    }

    #[wasm_bindgen_test]
    fn parse_hgvs_extracts_fields() {
        let val = parse_hgvs("chrF:g.6C>T").expect("parse should succeed");
        let parsed: ParsedOut = serde_wasm_bindgen::from_value(val).expect("deserialize");
        assert_eq!(parsed.accession, "chrF");
        assert_eq!(parsed.reference_type, "g");
        assert!(parsed.gene.is_none());
    }

    #[wasm_bindgen_test]
    fn parse_hgvs_rejects_garbage() {
        assert!(parse_hgvs("not an hgvs string").is_err());
    }
}
