//! AUDIT TEST — finding #4 (MEDIUM): VCF→VRS does not validate the VCF REF
//! allele against the actual reference sequence.
//!
//! The readonly VCF streaming core
//! (`gtars_vrs::compute_vrs_ids_streaming_readonly_from_reader`, in
//! `vcf_core.rs`) hands the VCF's REF allele straight into
//! `normalize::normalize_ref(view, pos, ref_allele, alt_allele)`. `normalize_ref`
//! uses the supplied REF only for its LENGTH (to compute the end coordinate and
//! to trim shared prefix/suffix between REF and ALT); it then ROLLS against the
//! *real* reference bases pulled from the store. At no point does it compare the
//! supplied REF bytes to the actual reference bases at that position.
//!
//! Contrast: the HGVS bridge DOES cross-check the parsed REF against the real
//! reference and raises `BridgeError::RefMismatch` on disagreement
//! (`hgvs/bridge.rs` ~line 372 / 730). The VCF path has no equivalent guard.
//!
//! Consequence: a VCF record whose REF disagrees with the loaded reference (a
//! common real-world symptom of wrong-assembly / liftover / data-entry errors)
//! is silently processed. Because trimming is driven by the (wrong) REF while
//! rolling uses the (right) reference, the emitted VRS id can be computed from a
//! REF the reference never had — a wrong id, emitted with no error and no flag.
//!
//! This test loads a reference where position (1-based) 6 is 'C'. It computes a
//! VRS id from a record claiming `REF=G` (WRONG) and compares it to the id from
//! the correct `REF=C`. The asserts encode the EXPECTED-CORRECT behavior — a
//! REF/reference mismatch should be rejected (error) or at minimum should not
//! silently yield the same id as a correct REF — so the test FAILS while the
//! bug is present, CONFIRMING it. We also print the concrete ids so the
//! silent-acceptance is visible in `--nocapture` output.

use std::collections::HashMap;
use std::io::Cursor;

use gtars_refget::store::{FastaImportOptions, RefgetStore};
use gtars_vrs::{VrsResult, compute_vrs_ids_streaming_readonly_from_reader};

const CHR_NAME: &str = "chrF";
// 1-based pos:  1234567...
//               ACGTACGT...
// pos 6 (1-based) == interbase 5 == 'C'.
const CHR_SEQ: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // len 40

/// Build an in-memory readonly RefgetStore holding `chrF`, plus the
/// name->raw-digest map the streaming core needs.
fn build_readonly_store() -> (
    gtars_refget::store::ReadonlyRefgetStore,
    HashMap<String, String>,
    tempfile::TempDir,
) {
    let temp = tempfile::tempdir().unwrap();
    let fasta_path = temp.path().join("fix.fa");
    std::fs::write(
        &fasta_path,
        format!(">{}\n{}\n", CHR_NAME, std::str::from_utf8(CHR_SEQ).unwrap()),
    )
    .unwrap();

    let mut store = RefgetStore::in_memory();
    store.disable_encoding();
    store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
        .unwrap();
    store.load_all_sequences().unwrap();

    let collection_digest = store
        .iter_collections()
        .next()
        .map(|c| c.metadata.digest.clone())
        .expect("collection");
    let coll = store.get_collection(&collection_digest).unwrap();
    let mut name_to_digest = HashMap::new();
    for r in &coll.sequences {
        let m = r.metadata();
        name_to_digest.insert(m.name.clone(), m.sha512t24u.clone());
    }

    (store.into_readonly(), name_to_digest, temp)
}

/// Run a single one-line VCF through the readonly streaming core and collect
/// any `VrsResult`s. Returns the full `Result` so we can observe whether the
/// pipeline ERRORS on a bad REF (it should; it does not).
fn run_one_record(
    store: &gtars_refget::store::ReadonlyRefgetStore,
    n2d: &HashMap<String, String>,
    vcf_line: &str,
) -> anyhow::Result<Vec<VrsResult>> {
    let vcf = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{vcf_line}\n"
    );
    let mut out = Vec::new();
    compute_vrs_ids_streaming_readonly_from_reader(
        store,
        n2d,
        Cursor::new(vcf.into_bytes()),
        |r| out.push(r),
    )?;
    Ok(out)
}

/// BUG B. A VCF record at chrF:6 with the CORRECT ref ('C') and one with a
/// WRONG ref ('G', the reference is actually 'C') are both fed to the VCF→VRS
/// path.
///
/// EXPECTED-CORRECT behavior (encoded by the asserts): the wrong-REF record is
/// rejected (the call errors, or yields no result, or at the very least yields
/// a DIFFERENT id than the correct-REF record — never silently the same id as
/// if the REF had been right).
///
/// ACTUAL behavior: no REF validation exists. The wrong REF 'G' has the same
/// length as the correct 'C', so trimming behaves identically and the variant
/// is rolled against the real reference. The path emits an id with NO error.
/// The assert below therefore FAILS, confirming the bug.
#[ignore = "audit reproduction (PR #256) — run with --ignored"]
#[test]
fn vcf_wrong_ref_must_be_rejected() {
    let (store, n2d, _temp) = build_readonly_store();

    // Correct REF: chrF[ib 5] (1-based 6) == 'C'. C>T SNV.
    let correct = run_one_record(&store, &n2d, "chrF\t6\t.\tC\tT\t.\tPASS\t.")
        .expect("correct-REF record should process");
    assert_eq!(correct.len(), 1, "correct REF should yield exactly one result");
    let correct_id = correct[0].vrs_id.clone();
    eprintln!("correct REF=C -> {correct_id}");

    // WRONG REF: claims 'G' at a position whose reference base is 'C'.
    let wrong = run_one_record(&store, &n2d, "chrF\t6\t.\tG\tT\t.\tPASS\t.");

    match &wrong {
        Err(e) => {
            // Desired outcome: pipeline errored on the REF/reference mismatch.
            eprintln!("wrong REF=G -> Err: {e:#}  (validation present — good)");
        }
        Ok(results) => {
            let wrong_id = results
                .first()
                .map(|r| r.vrs_id.clone())
                .unwrap_or_else(|| "<no result>".to_string());
            eprintln!(
                "wrong REF=G -> Ok, vrs_id = {wrong_id}  \
                 (NO error raised on a REF that disagrees with the reference)"
            );
        }
    }

    // EXPECTED-CORRECT contract: a REF that contradicts the loaded reference
    // must NOT be silently accepted and turned into a VRS id. Either the call
    // errors, or it must at least not produce a result that pretends the REF
    // was valid.
    assert!(
        wrong.is_err() || wrong.as_ref().unwrap().is_empty(),
        "BUG #4: VCF→VRS accepted REF=G at a position whose reference base is \
         'C' and emitted a VRS id with no error/flag. The VCF REF allele is \
         never validated against the actual reference sequence. \
         correct(REF=C) id = {correct_id}; wrong(REF=G) outcome = {wrong:?}"
    );
}

/// Companion that makes the silent-corruption concrete: with a WRONG ref whose
/// trimming differs from the correct ref, the emitted id can differ from the
/// correct one while STILL being emitted without error. Here REF="GT" ALT="G"
/// (a claimed 1bp deletion of the base after pos 6) is processed against a
/// reference whose bytes at 6_7 are actually "CG", so the claimed REF is wrong.
/// A correct pipeline would reject the mismatch; this one emits an id.
#[ignore = "audit reproduction (PR #256) — run with --ignored"]
#[test]
fn vcf_wrong_ref_indel_emits_without_validation() {
    let (store, n2d, _temp) = build_readonly_store();

    // chrF[ib5..ib7] (1-based 6..7) == "CG". A truthful deletion record would be
    // REF="CG" ALT="C". We instead claim REF="GG" ALT="G" — first base wrong.
    let wrong = run_one_record(&store, &n2d, "chrF\t6\t.\tGG\tG\t.\tPASS\t.");
    eprintln!("wrong indel REF=GG -> {wrong:?}");

    assert!(
        wrong.is_err() || wrong.as_ref().unwrap().is_empty(),
        "BUG #4: VCF→VRS accepted an indel REF='GG' that disagrees with the \
         reference ('CG') and emitted a VRS id with no error. got {wrong:?}"
    );
}
