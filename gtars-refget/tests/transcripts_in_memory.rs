//! WASM-config proof: the in-memory transcript backend stands alone without
//! the file backends.
//!
//! This test compiles and runs with `transcripts` but WITHOUT `filesystem`
//! (the wasm configuration), proving `ReadonlyTxStore::from_bytes` +
//! `TxBytes::InMemory` is fully usable on its own. It builds a valid `.reftx`
//! byte image entirely in memory (no `std::fs`, no `memmap2`) and runs lookups.
//!
//! Run with: `cargo test -p gtars-refget --no-default-features \
//!            --features transcripts --test transcripts_in_memory`
#![cfg(feature = "transcripts")]

use gtars_refget::{
    build_reftx_bytes_in_memory, Exon, ManeStatus, ReadonlyTxStore, Strand, Transcript,
};

#[test]
fn test_in_memory_backend_standalone() {
    let braf = Transcript {
        accession: "NM_004333.6".to_string(),
        gene: "BRAF".to_string(),
        chrom_digest: [1u8; 24],
        strand: Strand::Forward,
        cds_start: Some(100),
        cds_end: Some(400),
        exons: vec![Exon { start: 50, end: 500 }],
        mane: ManeStatus { mane_select: true, mane_clinical: false },
    };
    let tp53 = Transcript {
        accession: "NM_000546.6".to_string(),
        gene: "TP53".to_string(),
        chrom_digest: [2u8; 24],
        strand: Strand::Reverse,
        cds_start: None,
        cds_end: None,
        exons: vec![Exon { start: 10, end: 20 }],
        mane: Default::default(),
    };

    let bytes = build_reftx_bytes_in_memory(&[braf, tp53]).unwrap();
    let store = ReadonlyTxStore::from_bytes(bytes).unwrap();
    assert_eq!(store.len(), 2);

    let tx = store.lookup("NM_004333.6").expect("must find BRAF");
    assert_eq!(tx.gene, "BRAF");
    assert!(store.lookup("NM_MISSING.1").is_none());

    assert!(store.has_mane_index());
    assert_eq!(store.lookup_mane("braf").unwrap().accession, "NM_004333.6");
    assert!(store.lookup_mane("TP53").is_none()); // not MANE-select
}
