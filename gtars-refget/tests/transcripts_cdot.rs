//! `TxStoreBuilder::ingest_cdot` against real cdot files.
//!
//! The fixtures in `tests/data/cdot/` are copied unchanged from cdot's own
//! test data, so these tests fail if the parser drifts from the real layout.
#![cfg(all(feature = "transcripts", feature = "filesystem"))]

use std::path::PathBuf;

use gtars_refget::{
    build_reftx_bytes_in_memory, CoordinateMapper, ReadonlyTxStore, Strand, TxStoreBuilder,
};

fn cdot_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/cdot")
}

fn builder_with(contigs: &[&str]) -> TxStoreBuilder {
    let mut b = TxStoreBuilder::new();
    for (i, c) in contigs.iter().enumerate() {
        b.add_chrom_mapping(c, [i as u8 + 1; 24]);
    }
    b
}

#[test]
fn ingests_real_refseq_file() {
    let mut b = builder_with(&["NC_000007.13", "NC_000002.11"]);
    let n = b
        .ingest_cdot(cdot_dir().join("cdot.refseq.grch37.json"), None)
        .unwrap();
    assert_eq!(n, 2);

    let tx = b
        .transcripts
        .iter()
        .find(|t| t.accession == "NM_001637.3")
        .unwrap();
    assert_eq!(tx.gene, "AOAH");
    assert_eq!(tx.strand, Strand::Reverse);
    assert_eq!(tx.cds_start, Some(36552857));
    assert_eq!(tx.cds_end, Some(36763753));
    assert_eq!(tx.exons.len(), 21);
    assert_eq!((tx.exons[0].start, tx.exons[0].end), (36552548, 36552986));
    assert!(tx.exons.windows(2).all(|w| w[0].start < w[1].start));

    let nc = b
        .transcripts
        .iter()
        .find(|t| t.accession == "NR_023343.1")
        .unwrap();
    assert_eq!(nc.strand, Strand::Forward);
    assert_eq!(nc.cds_start, None);
    assert_eq!(nc.exons.len(), 1);
}

#[test]
fn skips_transcripts_on_unmapped_contigs() {
    let mut b = builder_with(&["NC_000007.13"]);
    let n = b
        .ingest_cdot(cdot_dir().join("cdot.refseq.grch37.json"), None)
        .unwrap();
    assert_eq!(n, 1);
    assert_eq!(b.transcripts[0].accession, "NM_001637.3");
}

#[test]
fn rejects_build_not_in_file() {
    let mut b = builder_with(&["NC_000007.14"]);
    let err = b
        .ingest_cdot(cdot_dir().join("cdot.ensembl.grch38.json"), Some("GRCh37"))
        .unwrap_err();
    assert!(err.to_string().contains("GRCh37"), "{err}");
}

/// c.1 must land on the first base of the start codon. On a minus-strand
/// transcript that is the last base of the CDS, i.e. `cds_end - 1`. This
/// checks that cdot's coordinates mean what the mapper assumes (0-based,
/// half-open genomic).
#[test]
fn maps_c1_on_real_minus_strand_transcript() {
    let mut b = builder_with(&["NC_000007.14"]);
    b.ingest_cdot(cdot_dir().join("cdot.ensembl.grch38.json"), Some("GRCh38"))
        .unwrap();
    let store =
        ReadonlyTxStore::from_bytes(build_reftx_bytes_in_memory(&b.transcripts).unwrap()).unwrap();
    let mapper = CoordinateMapper::new(&store);

    let c1 = mapper.c_to_g("ENST00000617537.5", 1).unwrap();
    assert_eq!(c1.position, 36724148 - 1);
    // c.-1 is one base upstream on the transcript, so one base higher on the genome.
    let c_minus1 = mapper.c_to_g("ENST00000617537.5", -1).unwrap();
    assert_eq!(c_minus1.position, 36724148);
}

#[test]
fn multi_build_file_needs_a_build() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("multi.json");
    std::fs::write(
        &path,
        r#"{
          "genome_builds": ["GRCh37", "GRCh38"],
          "transcripts": {
            "NM_X.1": {
              "id": "NM_X.1",
              "gene_name": "X",
              "genome_builds": {
                "GRCh37": {"contig": "NC_A", "strand": "+", "exons": [[100, 200, 0, 1, 100, null]]},
                "GRCh38": {"contig": "NC_B", "strand": "+", "exons": [[500, 600, 0, 1, 100, null]]}
              }
            }
          }
        }"#,
    )
    .unwrap();

    let mut b = builder_with(&["NC_A", "NC_B"]);
    assert!(b.ingest_cdot(&path, None).is_err());

    assert_eq!(b.ingest_cdot(&path, Some("GRCh38")).unwrap(), 1);
    assert_eq!(b.transcripts[0].exons[0].start, 500);
}
