//! End-to-end integration test: FASTA + VCF → VRS IDs

use std::io::Write;

use gtars_refget::store::{FastaImportOptions, RefgetStore};
use gtars_vrs::vcf::{compute_vrs_ids_from_vcf, unique_chroms_from_vcf};
use tempfile::tempdir;

/// Create a test FASTA with a known sequence and a matching VCF with variants.
/// The sequence is designed to include repeat regions for normalization testing.
#[test]
fn test_vcf_to_vrs_ids_end_to_end() {
    let dir = tempdir().unwrap();

    // chr1: sequence with an A-repeat region (positions 50-59)
    let chr1_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACAAAAAAAAAAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCTTTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
    let chr2_seq = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG";

    let fasta_path = dir.path().join("test.fa");
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        write!(f, ">chr1\n{}\n>chr2\n{}\n", chr1_seq, chr2_seq).unwrap();
    }

    // Create a VCF with various variant types
    let vcf_path = dir.path().join("test.vcf");
    {
        let mut f = std::fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        // SNV at chr1:5 (1-based), 0-based pos 4
        writeln!(f, "chr1\t5\t.\tA\tT\t.\tPASS\t.").unwrap();
        // SNV at chr1:10 (1-based)
        writeln!(f, "chr1\t10\t.\tG\tC\t.\tPASS\t.").unwrap();
        // Insertion in A-repeat: chr1:51 (1-based = 0-based 50), ref=A, alt=AA
        writeln!(f, "chr1\t51\t.\tA\tAA\t.\tPASS\t.").unwrap();
        // Deletion in A-repeat: chr1:51, ref=AA, alt=A
        writeln!(f, "chr1\t51\t.\tAA\tA\t.\tPASS\t.").unwrap();
        // Multi-allelic: chr2:5, ref=G, alt=A,T
        writeln!(f, "chr2\t5\t.\tG\tA,T\t.\tPASS\t.").unwrap();
        // Symbolic allele (should be skipped)
        writeln!(f, "chr1\t20\t.\tA\t<DEL>\t.\tPASS\t.").unwrap();
    }

    // Load FASTA into RefgetStore
    let mut store = RefgetStore::in_memory();
    store.disable_encoding(); // Use raw mode for simplicity
    store
        .add_sequence_collection_from_fasta(&fasta_path, FastaImportOptions::new())
        .unwrap();

    // Get collection digest
    let collection_digest = {
        for meta in store.sequence_metadata() {
            eprintln!(
                "Sequence: {} (len={}, digest={})",
                meta.name, meta.length, meta.sha512t24u
            );
        }
        let paged = store.list_collections(0, 10, &[]).unwrap();
        assert_eq!(paged.results.len(), 1);
        paged.results[0].digest.clone()
    };

    // Run VRS ID computation
    let results = compute_vrs_ids_from_vcf(
        &mut store,
        &collection_digest,
        vcf_path.to_str().unwrap(),
    )
    .unwrap();

    // Print results for debugging
    eprintln!("\nVRS Results:");
    for r in &results {
        eprintln!(
            "  {}:{} {}>{} → {}",
            r.chrom, r.pos, r.ref_allele, r.alt_allele, r.vrs_id
        );
    }

    // Should have 6 results: 2 SNVs + 1 insertion + 1 deletion + 2 from multi-allelic
    // (symbolic <DEL> should be skipped)
    assert_eq!(results.len(), 6, "Expected 6 results, got {}", results.len());

    // All VRS IDs should have the right format
    for r in &results {
        assert!(
            r.vrs_id.starts_with("ga4gh:VA."),
            "Bad VRS ID format: {}",
            r.vrs_id
        );
        assert_eq!(r.vrs_id.len(), 9 + 32, "Bad VRS ID length: {}", r.vrs_id);
    }

    // Same variant should always produce same ID (determinism)
    let results2 = compute_vrs_ids_from_vcf(
        &mut store,
        &collection_digest,
        vcf_path.to_str().unwrap(),
    )
    .unwrap();
    for (a, b) in results.iter().zip(results2.iter()) {
        assert_eq!(a.vrs_id, b.vrs_id, "Non-deterministic VRS ID");
    }

    // The two multi-allelic results (A and T) should have different IDs
    let multi_a = &results[4];
    let multi_t = &results[5];
    assert_eq!(multi_a.alt_allele, "A");
    assert_eq!(multi_t.alt_allele, "T");
    assert_ne!(multi_a.vrs_id, multi_t.vrs_id);
}

/// Regression test for unique_chroms_from_vcf: validates that it correctly
/// returns all unique chromosomes from a multi-chromosome VCF.
///
/// This test protects against the line buffer accumulation bug described in:
/// https://github.com/databio/gtars/issues/..._
/// If the string buffer passed to read_line is not cleared between iterations,
/// read_line appends to the existing buffer, causing split_once('\t') to always
/// return the first chromosome encountered instead of extracting from each line.
#[test]
fn test_unique_chroms_from_vcf_multi_chrom_correctness() {
    let dir = tempdir().unwrap();

    // Create a VCF with at least 3 different chromosomes
    let vcf_path = dir.path().join("multi_chrom.vcf");
    {
        let mut f = std::fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        // Variants on chr1
        writeln!(f, "chr1\t5\t.\tA\tT\t.\tPASS\t.").unwrap();
        writeln!(f, "chr1\t10\t.\tG\tC\t.\tPASS\t.").unwrap();
        // Variants on chr2
        writeln!(f, "chr2\t15\t.\tT\tA\t.\tPASS\t.").unwrap();
        // Variants on chr3
        writeln!(f, "chr3\t20\t.\tC\tG\t.\tPASS\t.").unwrap();
        writeln!(f, "chr3\t25\t.\tG\tA\t.\tPASS\t.").unwrap();
    }

    // Call unique_chroms_from_vcf
    let chroms = unique_chroms_from_vcf(vcf_path.to_str().unwrap()).unwrap();

    // Sort for deterministic comparison
    let mut chroms = chroms;
    chroms.sort();

    // Should contain all 3 chromosomes
    assert_eq!(
        chroms.len(),
        3,
        "Expected exactly 3 unique chromosomes, got {}. Chroms: {:?}",
        chroms.len(),
        chroms
    );

    // Verify the expected chromosomes are present
    assert_eq!(chroms[0], "chr1");
    assert_eq!(chroms[1], "chr2");
    assert_eq!(chroms[2], "chr3");
}
