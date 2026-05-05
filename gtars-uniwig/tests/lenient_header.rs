use std::path::PathBuf;

use gtars_uniwig::uniwig_main;

#[test]
fn uniwig_processes_bam_with_vn_less_pg_record() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest_dir.join("tests/data/vn_less_pg.bam");
    let chrom_sizes = manifest_dir.join("tests/data/vn_less_pg.chrom.sizes");

    assert!(bam.exists(), "fixture missing: {}", bam.display());
    assert!(chrom_sizes.exists(), "fixture missing: {}", chrom_sizes.display());

    let tempdir = tempfile::tempdir().unwrap();
    let prefix = tempdir.path().join("out").to_string_lossy().into_owned();

    uniwig_main(
        vec!["start", "end", "core"],
        1,
        bam.to_str().unwrap(),
        chrom_sizes.to_str().unwrap(),
        &prefix,
        "wig",
        "bam",
        2,
        false,
        1,
        0,
        false,
        true,
        1.0,
        "fixed",
    )
    .expect("uniwig_main panicked or returned an error on a BAM with a VN-less @PG record");
}
