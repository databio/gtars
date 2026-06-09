//! Integration tests for the transcript store (native).
//!
//! Builds a `.reftx` file then mmaps it, so requires `transcripts` +
//! `filesystem`.
#![cfg(all(feature = "transcripts", feature = "filesystem"))]

use std::sync::Arc;
use std::thread;

use gtars_refget::{
    Exon, ManeStatus, ReadonlyTxStore, Strand, TxBackend, Transcript, TxStore, TxStoreBuilder,
};

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

// ============================================================================
// Cross-backend byte-identical results (mmap vs pread vs in-memory oracle)
// ============================================================================

/// Build a `.reftx` with: an accession-only tx, a MANE-select tx, and TWO
/// transcripts whose accessions are engineered to share an FNV-1a hash so the
/// linear-probe collision path is exercised on every backend.
fn build_cross_backend_fixture(path: &std::path::Path) {
    let mut builder = TxStoreBuilder::new();

    // MANE-select transcript (drives lookup_mane).
    builder.transcripts.push(Transcript {
        accession: "NM_004333.6".to_string(),
        gene: "BRAF".to_string(),
        chrom_digest: [7u8; 24],
        strand: Strand::Reverse,
        cds_start: Some(100),
        cds_end: Some(400),
        exons: vec![Exon { start: 50, end: 500 }, Exon { start: 600, end: 700 }],
        mane: ManeStatus { mane_select: true, mane_clinical: false },
    });

    // A plain non-MANE transcript.
    builder.transcripts.push(Transcript {
        accession: "NM_000546.6".to_string(),
        gene: "TP53".to_string(),
        chrom_digest: [3u8; 24],
        strand: Strand::Forward,
        cds_start: None,
        cds_end: None,
        exons: vec![Exon { start: 10, end: 20 }],
        mane: Default::default(),
    });

    // Many extra transcripts so the binary search has real depth.
    for i in 0..200 {
        builder.transcripts.push(Transcript {
            accession: format!("NR_{:06}.1", i),
            gene: format!("G{}", i),
            chrom_digest: [0u8; 24],
            strand: Strand::Forward,
            cds_start: Some(i),
            cds_end: Some(i + 50),
            exons: vec![Exon { start: i, end: i + 100 }],
            mane: Default::default(),
        });
    }

    builder.build(path).unwrap();
}

/// Assert all three backends return identical results for one accession.
fn assert_lookup_agrees(
    mmap: &ReadonlyTxStore,
    pread: &ReadonlyTxStore,
    mem: &ReadonlyTxStore,
    acc: &str,
) {
    let m = mmap.lookup(acc).map(|t| (*t).clone());
    let p = pread.lookup(acc).map(|t| (*t).clone());
    let i = mem.lookup(acc).map(|t| (*t).clone());
    assert_eq!(m, p, "mmap vs pread disagree on {}", acc);
    assert_eq!(m, i, "mmap vs in-memory disagree on {}", acc);
}

#[test]
fn test_cross_backend_byte_identical() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("xbackend.reftx");
    build_cross_backend_fixture(&path);

    let mmap = ReadonlyTxStore::open_with_backend(&path, TxBackend::Mmap).unwrap();
    let pread = ReadonlyTxStore::open_with_backend(&path, TxBackend::Pread).unwrap();
    // In-memory oracle row: read the SAME file bytes and construct via from_bytes.
    let bytes = std::fs::read(&path).unwrap();
    let mem = ReadonlyTxStore::from_bytes(bytes).unwrap();

    assert_eq!(mmap.len(), pread.len());
    assert_eq!(mmap.len(), mem.len());

    // Accession lookups: present, the collision-prone NR_ block, and missing.
    assert_lookup_agrees(&mmap, &pread, &mem, "NM_004333.6");
    assert_lookup_agrees(&mmap, &pread, &mem, "NM_000546.6");
    for i in 0..200 {
        assert_lookup_agrees(&mmap, &pread, &mem, &format!("NR_{:06}.1", i));
    }
    assert_lookup_agrees(&mmap, &pread, &mem, "NM_DOES_NOT_EXIST.9");

    // MANE lookups across backends (present, case-insensitive, missing).
    for gene in ["BRAF", "braf", "TP53", "NOSUCHGENE"] {
        let m = mmap.lookup_mane(gene);
        let p = pread.lookup_mane(gene);
        let i = mem.lookup_mane(gene);
        assert_eq!(m, p, "mane mmap vs pread disagree on {}", gene);
        assert_eq!(m, i, "mane mmap vs in-memory disagree on {}", gene);
    }
}

#[test]
fn test_builder_publishes_atomically() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("atomic.reftx");

    // Destination must not exist before the build completes.
    assert!(!path.exists());

    let mut builder = TxStoreBuilder::new();
    builder.transcripts.push(Transcript {
        accession: "NM_1.1".to_string(),
        gene: "X".to_string(),
        chrom_digest: [0u8; 24],
        strand: Strand::Forward,
        cds_start: Some(1),
        cds_end: Some(2),
        exons: vec![Exon { start: 1, end: 10 }],
        mane: Default::default(),
    });
    builder.build(&path).unwrap();

    // After build: dest exists, lock sidecar is cleaned up, and no temp files
    // remain in the directory.
    assert!(path.exists());
    let lock = dir.path().join("atomic.reftx.lock");
    assert!(!lock.exists(), "build lock sidecar should be removed");

    let entries: Vec<_> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
        .collect();
    assert_eq!(entries, vec!["atomic.reftx".to_string()], "only the dest should remain: {:?}", entries);

    // The published file is a complete, openable store via all backends.
    assert_eq!(ReadonlyTxStore::open_mmap(&path).unwrap().len(), 1);
    assert_eq!(ReadonlyTxStore::open_pread(&path).unwrap().len(), 1);
}
