//! WASM batch entry point: VCF → VRS allele identifiers.
//!
//! Given a persistent [`RefgetStore`](crate::refget::RefgetStore) holding the
//! reference genome and the text of a dropped VCF, streams a VRS allele id for
//! every variant. The genome is reused across the whole file (no per-variant
//! store rebuild), and the reference is decoded on the fly from its encoded form.
//!
//! This is the browser analog of the native `compute_vrs_ids_streaming_readonly`:
//! instead of opening a path, it runs the shared reader-generic core
//! ([`gtars_vrs::vcf_core::compute_vrs_ids_streaming_readonly_from_reader`]) over
//! a `Cursor` of the in-memory VCF bytes.
//!
//! # Example (JavaScript)
//! ```javascript
//! import init, { RefgetStore, vcf_to_vrs_ids } from "gtars-js";
//! await init();
//!
//! const store = new RefgetStore();
//! store.add_sequence("chr1", chr1Bases);
//!
//! const vcfText = await file.text();
//! const count = vcf_to_vrs_ids(store, vcfText, (chrom, pos, ref, alt, vrsId) => {
//!     // pos is 0-based interbase; one call per ALT allele
//!     console.log(`${chrom}:${pos} ${ref}>${alt} -> ${vrsId}`);
//! });
//! console.log(`${count} VRS ids`);
//! ```

use std::io::Cursor;

use gtars_vrs::vcf_core::{VrsResult, compute_vrs_ids_streaming_readonly_from_reader};
use wasm_bindgen::prelude::*;

use crate::refget::RefgetStore;

/// Compute a VRS allele id for every variant in `vcf_text`, invoking `on_result`
/// once per ALT allele with `(chrom, pos, ref, alt, vrs_id)`.
///
/// `pos` is the 0-based interbase coordinate (1-based VCF POS minus one), passed
/// as a JS number. Returns the number of results produced.
///
/// # Arguments
/// * `store` - The resident genome handle (see [`RefgetStore`]).
/// * `vcf_text` - The VCF contents as text (already decompressed; the caller
///   inflates `.vcf.gz` in JS).
/// * `on_result` - A JS callback `(chrom, pos, ref, alt, vrs_id)`.
///
/// # Errors
/// Returns a `JsValue` error string if reference normalization fails for a
/// variant (e.g. a chromosome present in the VCF but missing from the store),
/// or if the callback itself throws.
#[wasm_bindgen]
pub fn vcf_to_vrs_ids(
    store: &RefgetStore,
    vcf_text: &str,
    on_result: &js_sys::Function,
) -> Result<usize, JsValue> {
    let this = JsValue::null();
    // First JS callback error is captured and surfaced after the loop unwinds.
    let mut callback_err: Option<JsValue> = None;

    let reader = Cursor::new(vcf_text.as_bytes());
    let count = compute_vrs_ids_streaming_readonly_from_reader(
        store.inner_readonly(),
        store.name_to_digest(),
        reader,
        |r: VrsResult| {
            if callback_err.is_some() {
                return; // already failed; stop calling out to JS
            }
            let args = js_sys::Array::new();
            args.push(&JsValue::from_str(&r.chrom));
            args.push(&JsValue::from_f64(r.pos as f64));
            args.push(&JsValue::from_str(&r.ref_allele));
            args.push(&JsValue::from_str(&r.alt_allele));
            args.push(&JsValue::from_str(&r.vrs_id));
            if let Err(e) = on_result.apply(&this, &args) {
                callback_err = Some(e);
            }
        },
    )
    .map_err(|e| JsValue::from_str(&format!("{e}")))?;

    if let Some(e) = callback_err {
        return Err(e);
    }
    Ok(count)
}

// ── VCF + transcript seam (design) ───────────────────────────────────────────
//
// Standard VCF rows are GENOMIC (CHROM/POS), so the streaming core
// (`compute_vrs_ids_streaming_readonly_from_reader`) needs no transcript
// provider and `vcf_to_vrs_ids` above passes none. The transcript provider only
// matters for transcript-coordinate inputs (a future VCF dialect or an HGVS
// column carrying `c.`/`n.`).
//
// The seam is designed to mirror the HGVS path exactly: the SAME
// `TranscriptStore` holder that serves `hgvs_to_vrs_id_with_transcripts` would
// serve VCF too — one transcript store, one genome store, both reused across the
// whole file (just as one `RefgetStore` backs the genomic VCF path today).
//
// FOLLOW-UP: a transcript-aware VCF batch entry
//
//     pub fn vcf_to_vrs_ids_with_transcripts(
//         store: &RefgetStore,
//         transcripts: &TranscriptStore,
//         vcf_text: &str,
//         on_result: &js_sys::Function,
//     ) -> Result<usize, JsValue>
//
// would thread `transcripts.provider()` into the streaming core. That requires a
// provider-accepting variant of
// `compute_vrs_ids_streaming_readonly_from_reader` in gtars-vrs
// (`gtars-vrs/src/vcf_core.rs`): add a `provider: &dyn TranscriptProvider`
// parameter, used only where a row's reference type is `c.`/`n.`, and have the
// existing genomic entry pass `&NoTranscriptProvider`. The streaming core does
// not yet parse transcript-coordinate references, so that is a non-trivial core
// change tracked here rather than landed in this gtars-wasm-scoped change; the
// `TranscriptStore` holder and `provider()` accessor designed in
// `transcripts.rs` are exactly what that follow-up plugs into.

#[cfg(test)]
#[cfg(target_arch = "wasm32")]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use std::rc::Rc;
    use wasm_bindgen::closure::Closure;
    use wasm_bindgen::JsCast;
    use wasm_bindgen_test::*;

    use crate::refget::RefgetStore;

    // Shares hgvs.rs's synthetic contig + golden id so the batch path is proven
    // to land on the same VRS id as the single-variant entry.
    const CHR_F_BASES: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    const GOLDEN_CHRF_G6CT_VRS_ID: &str = "ga4gh:VA._q-idtHGQxQ4XiEPJ1ExYl_htUeNEkir";

    #[wasm_bindgen_test]
    fn batch_vcf_matches_golden_id() {
        let mut store = RefgetStore::new();
        store
            .add_sequence("chrF", CHR_F_BASES.as_bytes())
            .expect("add_sequence");

        // One real SNV at chrF pos 6 (C>T), a header line, and a no-ALT row that
        // must be skipped.
        let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chrF\t6\t.\tC\tT\t.\t.\t.
chrF\t10\t.\tA\t.\t.\t.\t.
";

        let ids: Rc<RefCell<Vec<String>>> = Rc::new(RefCell::new(Vec::new()));
        let ids_cb = Rc::clone(&ids);
        let cb = Closure::wrap(Box::new(
            move |_chrom: String, _pos: f64, _ref: String, _alt: String, vrs_id: String| {
                ids_cb.borrow_mut().push(vrs_id);
            },
        )
            as Box<dyn FnMut(String, f64, String, String, String)>);
        let func: &js_sys::Function = cb.as_ref().unchecked_ref();

        let count = vcf_to_vrs_ids(&store, vcf, func).expect("vcf_to_vrs_ids");

        assert_eq!(count, 1, "only the C>T row yields a result");
        assert_eq!(ids.borrow().len(), 1);
        assert_eq!(ids.borrow()[0], GOLDEN_CHRF_G6CT_VRS_ID);
    }

    #[wasm_bindgen_test]
    fn batch_vcf_skips_unknown_chrom() {
        let mut store = RefgetStore::new();
        store
            .add_sequence("chrF", CHR_F_BASES.as_bytes())
            .expect("add_sequence");

        // chrZ is not in the store -> row is skipped, no error, no results.
        let vcf = "chrZ\t6\t.\tC\tT\t.\t.\t.\n";
        let noop = Closure::wrap(Box::new(
            move |_c: String, _p: f64, _r: String, _a: String, _v: String| {},
        )
            as Box<dyn FnMut(String, f64, String, String, String)>);
        let func: &js_sys::Function = noop.as_ref().unchecked_ref();

        let count = vcf_to_vrs_ids(&store, vcf, func).expect("vcf_to_vrs_ids");
        assert_eq!(count, 0);
    }
}
