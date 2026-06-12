//! REGRESSION GUARDS: VCF→VRS must validate the VCF REF allele against the
//! actual reference sequence.
//!
//! The readonly VCF streaming core
//! (`gtars_vrs::compute_vrs_ids_streaming_readonly_from_reader`, in
//! `vcf_core.rs`) hands the VCF's REF allele into
//! `normalize::normalize_ref(view, pos, ref_allele, alt_allele)`, which now
//! compares the supplied REF bytes to the actual reference bases at that
//! position and raises `NormalizeError::RefMismatch` on disagreement. This
//! mirrors the HGVS bridge's `BridgeError::RefMismatch` guard.
//!
//! A VCF record whose REF disagrees with the loaded reference (a common
//! real-world symptom of wrong-assembly / liftover / data-entry errors) must
//! fail loudly rather than emit a VRS id computed from a REF the reference
//! never had. These tests pin that behavior for both a same-length SNV
//! mismatch and an indel mismatch.

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

/// A VCF record at chrF:6 with the CORRECT ref ('C') processes, while one with
/// a WRONG ref ('G', the reference is actually 'C') must be rejected with a
/// `RefMismatch` error. The wrong REF 'G' has the same length as the correct
/// 'C', so trimming behaves identically — only an explicit REF-vs-reference
/// comparison can catch it.
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

    let err = wrong.expect_err(
        "VCF→VRS must reject REF=G at a position whose reference base is 'C'",
    );
    let msg = format!("{err:#}");
    eprintln!("wrong REF=G -> Err: {msg}");
    assert!(
        msg.contains("ref allele mismatch"),
        "expected a RefMismatch error, got: {msg}"
    );
}

/// Indel companion: REF="GG" ALT="G" at chrF:6, where the reference bytes at
/// interbase 5..7 are actually "CG". The first REF base disagrees, so the
/// record must be rejected with a `RefMismatch` error rather than emitting an id.
#[test]
fn vcf_wrong_ref_indel_must_be_rejected() {
    let (store, n2d, _temp) = build_readonly_store();

    // chrF[ib5..ib7] (1-based 6..7) == "CG". A truthful deletion record would be
    // REF="CG" ALT="C". We instead claim REF="GG" ALT="G" — first base wrong.
    let wrong = run_one_record(&store, &n2d, "chrF\t6\t.\tGG\tG\t.\tPASS\t.");

    let err = wrong.expect_err(
        "VCF→VRS must reject indel REF='GG' that disagrees with reference 'CG'",
    );
    let msg = format!("{err:#}");
    eprintln!("wrong indel REF=GG -> Err: {msg}");
    assert!(
        msg.contains("ref allele mismatch"),
        "expected a RefMismatch error, got: {msg}"
    );
}
