//! Transcript store lookup benchmark.
//!
//! Native: builds and mmaps a `.reftx` store, so requires `transcripts` +
//! `filesystem`. Run with `--features "transcripts filesystem"`.
#![cfg(all(feature = "transcripts", feature = "filesystem"))]

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gtars_refget::{
    CoordinateMapper, CoordinateMapperWriter, Exon, Strand, Transcript, TxStore, TxStoreBuilder,
};
use tempfile::tempdir;

fn create_test_store(count: usize) -> tempfile::TempDir {
    let dir = tempdir().unwrap();
    let path = dir.path().join("bench.reftx");

    let mut builder = TxStoreBuilder::new();
    for i in 0..count {
        builder.transcripts.push(Transcript {
            accession: format!("NM_{:06}.1", i),
            gene: format!("GENE{}", i),
            chrom_digest: [0u8; 24],
            strand: if i % 2 == 0 {
                Strand::Forward
            } else {
                Strand::Reverse
            },
            cds_start: Some(1000),
            cds_end: Some(5000),
            exons: vec![
                Exon { start: 500, end: 1500 },
                Exon { start: 2000, end: 3000 },
                Exon { start: 4000, end: 5500 },
            ],
            mane: Default::default(),
        });
    }
    builder.build(&path).unwrap();

    dir
}

fn bench_lookup(c: &mut Criterion) {
    let dir = create_test_store(100_000);
    let path = dir.path().join("bench.reftx");
    let store = TxStore::open(&path).unwrap().into_readonly_lazy();

    c.bench_function("lookup_single_cold", |b| {
        b.iter(|| black_box(store.lookup("NM_050000.1")))
    });

    let _ = store.lookup("NM_050000.1");

    c.bench_function("lookup_single_warm", |b| {
        b.iter(|| black_box(store.lookup("NM_050000.1")))
    });
}

fn bench_lookup_batch(c: &mut Criterion) {
    let dir = create_test_store(100_000);
    let path = dir.path().join("bench.reftx");
    let store = TxStore::open(&path).unwrap().into_readonly();

    let accessions: Vec<String> = (0..1000).map(|i| format!("NM_{:06}.1", i * 100)).collect();

    c.bench_function("lookup_batch_1000", |b| {
        b.iter(|| {
            for acc in &accessions {
                black_box(store.lookup(acc));
            }
        })
    });
}

fn bench_coordinate_mapping(c: &mut Criterion) {
    let dir = create_test_store(10_000);
    let path = dir.path().join("bench.reftx");
    let store = TxStore::open(&path).unwrap().into_readonly();

    let mapper = CoordinateMapper::new(&store);

    c.bench_function("c_to_g_single", |b| {
        b.iter(|| black_box(mapper.c_to_g("NM_005000.1", 100)))
    });
}

fn bench_coordinate_mapping_batch(c: &mut Criterion) {
    let dir = create_test_store(10_000);
    let path = dir.path().join("bench.reftx");
    let store = TxStore::open(&path).unwrap().into_readonly();

    let mut mapper = CoordinateMapperWriter::new(&store);

    c.bench_function("c_to_g_batch_1000_zero_alloc", |b| {
        b.iter(|| {
            for i in 1..=1000 {
                black_box(mapper.c_to_g("NM_005000.1", i));
            }
        })
    });
}

criterion_group!(
    benches,
    bench_lookup,
    bench_lookup_batch,
    bench_coordinate_mapping,
    bench_coordinate_mapping_batch
);
criterion_main!(benches);
