//! Integration tests for gtars-reftx.

use std::sync::Arc;
use std::thread;

use gtars_reftx::{Exon, Strand, Transcript, TxStore, TxStoreBuilder};

#[test]
fn test_concurrent_access() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("concurrent.reftx");

    let mut builder = TxStoreBuilder::new();
    for i in 0..1000 {
        builder.transcripts.push(Transcript {
            accession: format!("NM_{:06}.1", i),
            gene: format!("GENE{}", i),
            chrom_digest: [0u8; 24],
            strand: Strand::Forward,
            cds_start: Some(100),
            cds_end: Some(200),
            exons: vec![Exon { start: 50, end: 250 }],
            mane: Default::default(),
        });
    }
    builder.build(&path).unwrap();

    let store = Arc::new(TxStore::open(&path).unwrap().into_readonly());

    let handles: Vec<_> = (0..8)
        .map(|i| {
            let store = Arc::clone(&store);
            thread::spawn(move || {
                for j in 0..100 {
                    let acc = format!("NM_{:06}.1", (i * 100 + j) % 1000);
                    let tx = store.lookup(&acc);
                    assert!(tx.is_some(), "lookup failed for {}", acc);
                }
            })
        })
        .collect();

    for h in handles {
        h.join().unwrap();
    }
}
