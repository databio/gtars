//! Roundtrip test: decode the tiny fixture produced by the Python POC,
//! compare to expected records.

use gtars_rgvcf::{Record, RgvcfReader};

#[test]
#[ignore] // requires /tmp/plan9-work/tiny fixture
fn tiny_fixture_iteration() {
    let path = "/tmp/plan9-work/tiny/t.rgvcf";
    if !std::path::Path::new(path).exists() {
        eprintln!("skipping: {} not present", path);
        return;
    }
    let reader = RgvcfReader::open(path).unwrap();
    assert_eq!(reader.chromosomes().len(), 2);
    // chr1 — expect 3 inline SNVs + 2 escape (indel + multi-alt)
    let chr1: Vec<_> = reader
        .iter_chrom("chr1")
        .unwrap()
        .map(|r| r.unwrap())
        .collect();
    assert_eq!(chr1.len(), 5);
    let n_complex = chr1
        .iter()
        .filter(|r| matches!(r, Record::Complex { .. }))
        .count();
    assert_eq!(n_complex, 2);
}
