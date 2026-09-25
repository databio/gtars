//! Regression tests for residues that did not round-trip through Encoded mode.
//!
//! Each test compares decoded output against the ORIGINAL input, not against
//! another decoder. Written test-first: these fail until the fixes land.
//!
//! 1. DnaIupac decoding table: D came back as H, H came back as V.
//! 2. DnaIupac encoding table: U shared T's code, so U came back as T.
//! 3. Protein alphabet: U (selenocysteine) and B (Asx) had no code, were
//!    encoded as 0 (Alanine), and came back as A.

use gtars_refget::digest::{
    decode_string_from_bytes, encode_sequence, guess_alphabet, lookup_alphabet, sha512t24u,
    AlphabetType,
};
use gtars_refget::store::{FastaImportOptions, RefgetStore};
use std::io::Write;
use tempfile::NamedTempFile;

/// Guess the alphabet, encode, decode, and return (alphabet, decoded).
fn roundtrip(seq: &[u8]) -> (AlphabetType, Vec<u8>) {
    let alphabet_type = guess_alphabet(seq);
    let alphabet = lookup_alphabet(&alphabet_type);
    let encoded = encode_sequence(seq, alphabet);
    let decoded = decode_string_from_bytes(&encoded, seq.len(), alphabet);
    (alphabet_type, decoded)
}

fn assert_roundtrips(seq: &str) {
    let (alphabet_type, decoded) = roundtrip(seq.as_bytes());
    assert_eq!(
        String::from_utf8_lossy(&decoded),
        seq,
        "{seq} did not round-trip through {alphabet_type:?}"
    );
}

// --- 1. DnaIupac D and H (decoding table) ---

#[test]
fn iupac_d_roundtrips() {
    // Was: ACGTACGTHACGT
    assert_roundtrips("ACGTACGTDACGT");
}

#[test]
fn iupac_h_roundtrips() {
    // Was: ACGTACGTVACGT
    assert_roundtrips("ACGTACGTHACGT");
}

#[test]
fn iupac_all_codes_roundtrip() {
    assert_roundtrips("ACGTRYSWKMBDHVN");
}

#[test]
fn short_protein_assigned_iupac_roundtrips() {
    // ENSP00000499040.1 (NOTCH2). Only IUPAC letters, so it lands in DnaIupac.
    // Was: MCVTYVNGTGYC
    assert_roundtrips("MCVTYHNGTGYC");
}

// --- 2. DnaIupac U (encoding table) ---

#[test]
fn iupac_u_roundtrips() {
    // Was: ACGTACGTACGT
    assert_roundtrips("ACGUACGUACGU");
}

#[test]
fn iupac_u_and_t_stay_distinct() {
    assert_roundtrips("ACGTUACGTU");
}

// --- 3. Protein U and B (missing codes + alphabet selection) ---

#[test]
fn protein_selenocysteine_roundtrips() {
    // Was: MKLEEFWQPAIIAA
    assert_roundtrips("MKLEEFWQPUIIAA");
}

#[test]
fn protein_asx_roundtrips() {
    // Was: MKLEEFWQPAIIAA
    assert_roundtrips("MKLEEFWQPBIIAA");
}

// --- End to end: returned bytes must hash to the stored digest ---

#[test]
fn store_returns_bytes_matching_digest() {
    let cases = [
        ("iupac_d", "ACGTACGTDACGT"),
        ("iupac_h", "ACGTACGTHACGT"),
        ("notch2", "MCVTYHNGTGYC"),
        ("rna_u", "ACGUACGUACGU"),
        ("selenoprotein", "MKLEEFWQPUIIAA"),
        ("asx", "MKLEEFWQPBIIAA"),
    ];

    let mut fasta = NamedTempFile::new().unwrap();
    for (name, seq) in &cases {
        writeln!(fasta, ">{name}\n{seq}").unwrap();
    }

    let mut store = RefgetStore::in_memory();
    store
        .add_sequence_collection_from_fasta(fasta.path(), FastaImportOptions::new())
        .unwrap();
    let collection = store.list_collections(0, usize::MAX, &[]).unwrap().results[0]
        .digest
        .clone();

    let mut failures = Vec::new();
    for (name, seq) in &cases {
        let record = store.get_sequence_by_name(&collection, name).unwrap();
        let decoded = record.decode().unwrap();
        if decoded != *seq {
            failures.push(format!(
                "{name}: stored digest {} but returned {decoded} (digest {})",
                sha512t24u(seq.as_bytes()),
                sha512t24u(decoded.as_bytes()),
            ));
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}
