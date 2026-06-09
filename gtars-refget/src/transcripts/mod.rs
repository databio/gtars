//! Binary transcript store and HGVS coordinate mapper.
//!
//! Mirrors gtars-refget's own core/filesystem split:
//! - CORE (this `transcripts` feature, all targets, no mmap/fs): `models`,
//!   `mapper`, the `store` lookup logic, and an in-memory byte-backed transcript
//!   store built from owned `.reftx` bytes (`ReadonlyTxStore::from_bytes`).
//! - NATIVE-ONLY (`filesystem`): the file-backed byte backends in `mmap`
//!   (`TxBytes::Mmap` / `TxBytes::Pread`, opened via `open_mmap` / `open_pread`
//!   / `open_with_backend`, plus the mmap-backed mutable `TxStore`) and the
//!   atomic file-writing `builder`. The byte source (`store::TxBytes`) is the
//!   single place the three backends differ; all lookup logic is
//!   backend-agnostic and compiles on wasm with only `TxBytes::InMemory`.

pub mod mapper;
pub mod models;
pub mod store;

#[cfg(feature = "filesystem")]
pub mod builder;
#[cfg(feature = "filesystem")]
pub mod mmap;

pub use mapper::{CoordinateMapper, CoordinateMapperWriter, MappingError, MappingResult};
pub use models::{Exon, ManeStatus, Strand, Transcript};
pub use store::{build_reftx_bytes_in_memory, ReadonlyTxStore, TranscriptRef};

// Native-only: the mmap-backed mutable store, the file-backed backend selector,
// and the atomic file-writing builder.
#[cfg(feature = "filesystem")]
pub use builder::TxStoreBuilder;
#[cfg(feature = "filesystem")]
pub use mmap::{TxBackend, TxStore};

// ============================================================================
// Native integration-style unit tests (build a file, mmap it).
// ============================================================================

#[cfg(all(test, feature = "filesystem"))]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn sample_transcript() -> Transcript {
        Transcript {
            accession: "NM_004333.6".to_string(),
            gene: "BRAF".to_string(),
            chrom_digest: [0u8; 24],
            strand: Strand::Reverse,
            cds_start: Some(140719327),
            cds_end: Some(140924929),
            exons: vec![
                Exon { start: 140719327, end: 140719536 },
                Exon { start: 140726493, end: 140726516 },
                Exon { start: 140753274, end: 140753393 },
            ],
            mane: ManeStatus { mane_select: true, mane_clinical: false },
        }
    }

    #[test]
    fn test_exon_len() {
        let exon = Exon { start: 100, end: 200 };
        assert_eq!(exon.len(), 100);
    }

    #[test]
    fn test_transcript_length() {
        let tx = sample_transcript();
        let expected = (140719536 - 140719327) + (140726516 - 140726493) + (140753393 - 140753274);
        assert_eq!(tx.transcript_length(), expected);
    }

    #[test]
    fn test_transcript_is_coding() {
        let mut tx = sample_transcript();
        assert!(tx.is_coding());

        tx.cds_start = None;
        assert!(!tx.is_coding());
    }

    #[test]
    fn test_build_and_lookup() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.reftx");

        let mut builder = TxStoreBuilder::new();
        builder.transcripts.push(sample_transcript());
        builder.build(&path).unwrap();

        let store = TxStore::open(&path).unwrap();
        assert_eq!(store.len(), 1);

        let tx = store.lookup("NM_004333.6").unwrap();
        assert_eq!(tx.gene, "BRAF");
        assert_eq!(tx.exons.len(), 3);
        assert!(matches!(tx.strand, Strand::Reverse));
    }

    #[test]
    fn test_lookup_not_found() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.reftx");

        let mut builder = TxStoreBuilder::new();
        builder.transcripts.push(sample_transcript());
        builder.build(&path).unwrap();

        let store = TxStore::open(&path).unwrap();
        assert!(store.lookup("NM_NONEXISTENT.1").is_none());
    }

    #[test]
    fn test_readonly_store() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.reftx");

        let mut builder = TxStoreBuilder::new();
        builder.transcripts.push(sample_transcript());
        builder.build(&path).unwrap();

        let store = TxStore::open(&path).unwrap();
        let readonly = store.into_readonly();

        let tx = readonly.lookup("NM_004333.6").unwrap();
        assert_eq!(tx.gene, "BRAF");
    }

    #[test]
    fn test_coordinate_mapping_forward_strand() {
        let tx = Transcript {
            accession: "NM_TEST.1".to_string(),
            gene: "TEST".to_string(),
            chrom_digest: [0u8; 24],
            strand: Strand::Forward,
            cds_start: Some(100),
            cds_end: Some(200),
            exons: vec![
                Exon { start: 50, end: 150 },
                Exon { start: 170, end: 220 },
            ],
            mane: Default::default(),
        };

        let dir = tempdir().unwrap();
        let path = dir.path().join("test.reftx");

        let mut builder = TxStoreBuilder::new();
        builder.transcripts.push(tx);
        builder.build(&path).unwrap();

        let store = TxStore::open(&path).unwrap().into_readonly();
        let mapper = CoordinateMapper::new(&store);

        // c.1 should map to genomic 100 (CDS start)
        let result = mapper.c_to_g("NM_TEST.1", 1).unwrap();
        assert_eq!(result.position, 100);
    }

    fn test_mane_select_transcript() -> Transcript {
        Transcript {
            accession: "NM_004333.6".to_string(),
            gene: "BRAF".to_string(),
            chrom_digest: [1u8; 24],
            strand: Strand::Forward,
            cds_start: Some(100),
            cds_end: Some(400),
            exons: vec![Exon { start: 50, end: 500 }],
            mane: ManeStatus { mane_select: true, mane_clinical: false },
        }
    }

    #[test]
    fn test_mane_status_byte_roundtrip() {
        for byte in 0u8..4 {
            let m = ManeStatus::from_flags_byte(byte);
            assert_eq!(m.to_flags_byte(), byte & 0x03);
        }
        let m = ManeStatus { mane_select: true, mane_clinical: true };
        assert_eq!(m.to_flags_byte(), 0x03);
    }

    #[test]
    fn test_build_with_mane_and_lookup_mane() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("mane.reftx");

        let mut builder = TxStoreBuilder::new();
        builder.transcripts.push(test_mane_select_transcript());
        // Also a non-MANE transcript to confirm filtering.
        let mut other = test_mane_select_transcript();
        other.accession = "NM_OTHER.1".to_string();
        other.gene = "OTHER".to_string();
        other.mane = Default::default();
        builder.transcripts.push(other);

        builder.build(&path).unwrap();

        let store = TxStore::open(&path).unwrap().into_readonly();
        assert!(store.has_mane_index());

        let tx = store.lookup_mane("BRAF").expect("MANE BRAF must resolve");
        assert_eq!(tx.accession, "NM_004333.6");
        assert!(tx.mane.mane_select);

        // Case insensitive.
        let tx_lc = store.lookup_mane("braf").unwrap();
        assert_eq!(tx_lc.accession, "NM_004333.6");

        // Non-MANE gene returns None.
        assert!(store.lookup_mane("OTHER").is_none());
        assert!(store.lookup_mane("MISSING").is_none());
    }

    #[test]
    fn test_build_with_no_mane_has_zero_offset() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("nomane.reftx");
        let mut builder = TxStoreBuilder::new();
        builder.transcripts.push(sample_transcript()); // mane_select=true from helper
        // Override mane to false:
        builder.transcripts[0].mane = Default::default();
        builder.build(&path).unwrap();

        let store = TxStore::open(&path).unwrap().into_readonly();
        assert!(!store.has_mane_index());
        assert!(store.lookup_mane("BRAF").is_none());
    }

    #[test]
    fn test_c_to_g_full_5prime_utr_forward() {
        // Forward strand: CDS [100, 200), single exon [50, 250).
        let tx = Transcript {
            accession: "NM_TEST.1".to_string(),
            gene: "TEST".to_string(),
            chrom_digest: [0u8; 24],
            strand: Strand::Forward,
            cds_start: Some(100),
            cds_end: Some(200),
            exons: vec![Exon { start: 50, end: 250 }],
            mane: Default::default(),
        };
        let dir = tempdir().unwrap();
        let path = dir.path().join("t.reftx");
        let mut builder = TxStoreBuilder::new();
        builder.transcripts.push(tx);
        builder.build(&path).unwrap();
        let store = TxStore::open(&path).unwrap().into_readonly();
        let mapper = CoordinateMapper::new(&store);

        // c.1 → genomic 100, c.-1 → 99, c.-10 → 90.
        assert_eq!(mapper.c_to_g_full("NM_TEST.1", 1, 0, false).unwrap().position, 100);
        assert_eq!(mapper.c_to_g_full("NM_TEST.1", -1, 0, false).unwrap().position, 99);
        assert_eq!(mapper.c_to_g_full("NM_TEST.1", -10, 0, false).unwrap().position, 90);

        // c.*1 → CDS end position (200, since cds_end is exclusive).
        assert_eq!(mapper.c_to_g_full("NM_TEST.1", 1, 0, true).unwrap().position, 200);

        // 5' UTR underflow.
        assert!(mapper.c_to_g_full("NM_TEST.1", -100, 0, false).is_err());
    }

    #[test]
    fn test_c_to_g_full_intronic_forward() {
        // Two exons [50, 100) and [150, 250). CDS = [60, 220).
        let tx = Transcript {
            accession: "NM_T2.1".to_string(),
            gene: "T2".to_string(),
            chrom_digest: [0u8; 24],
            strand: Strand::Forward,
            cds_start: Some(60),
            cds_end: Some(220),
            exons: vec![
                Exon { start: 50, end: 100 },
                Exon { start: 150, end: 250 },
            ],
            mane: Default::default(),
        };
        let dir = tempdir().unwrap();
        let path = dir.path().join("t2.reftx");
        let mut builder = TxStoreBuilder::new();
        builder.transcripts.push(tx);
        builder.build(&path).unwrap();
        let store = TxStore::open(&path).unwrap().into_readonly();
        let mapper = CoordinateMapper::new(&store);

        // c.40 → genomic 99 (last base of exon 1).
        assert_eq!(mapper.c_to_g_full("NM_T2.1", 40, 0, false).unwrap().position, 99);
        // c.40+5 → 5 bases downstream in intron → 99 + 5 = 104.
        assert_eq!(mapper.c_to_g_full("NM_T2.1", 40, 5, false).unwrap().position, 104);
        // c.41-3 → 3 bases upstream of first base of exon 2 (genomic 150) → 147.
        assert_eq!(mapper.c_to_g_full("NM_T2.1", 41, -3, false).unwrap().position, 147);
        // c.30+1 (not at exon boundary in tx) → error.
        assert!(mapper.c_to_g_full("NM_T2.1", 30, 1, false).is_err());
    }

    #[test]
    fn test_c_to_g_by_gene_uses_mane() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("mane2.reftx");
        let mut builder = TxStoreBuilder::new();
        builder.transcripts.push(test_mane_select_transcript());
        builder.build(&path).unwrap();
        let store = TxStore::open(&path).unwrap().into_readonly();
        let mapper = CoordinateMapper::new(&store);

        let (acc, result) = mapper.c_to_g_by_gene("BRAF", 1, 0, false).unwrap();
        assert_eq!(acc, "NM_004333.6");
        assert_eq!(result.position, 100);

        assert!(mapper.c_to_g_by_gene("UNKNOWN_GENE", 1, 0, false).is_err());
    }

    #[test]
    fn test_mapper_writer_reuses_buffer() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.reftx");

        let tx = Transcript {
            accession: "NM_TEST.1".to_string(),
            gene: "TEST".to_string(),
            chrom_digest: [0u8; 24],
            strand: Strand::Forward,
            cds_start: Some(100),
            cds_end: Some(200),
            exons: vec![Exon { start: 50, end: 250 }],
            mane: Default::default(),
        };

        let mut builder = TxStoreBuilder::new();
        builder.transcripts.push(tx);
        builder.build(&path).unwrap();

        let store = TxStore::open(&path).unwrap().into_readonly();
        let mut mapper = CoordinateMapperWriter::new(&store);

        for i in 1..10 {
            let result = mapper.c_to_g("NM_TEST.1", i).unwrap();
            assert_eq!(result.position, 99 + i as u64);
        }
    }
}
