//! Audit tests for `gtars_uniwig::bamqc`.
//!
//! These tests synthesize tiny coordinate-sorted, BGZF-compressed BAM files
//! WITH a real BAI index (built via noodles' `csi::binning_index::Indexer`),
//! because `compute_bam_qc` queries the BAM by region per-chromosome and
//! therefore REQUIRES a valid index plus `@SQ` header entries.
//!
//! They are intended to CONFIRM two suspected bugs:
//!
//!   BUG 1 (bamqc.rs:212 / :286): the REPORTED `m2` field is forced to
//!     `m2.max(1)`. `m2` is the count of positions/read-pairs observed
//!     EXACTLY twice. The `.max(1)` is only needed for the pbc2 *denominator*
//!     (already handled separately via `m2_f`). Forcing the reported value
//!     to >= 1 means the "Two_read_pairs" column shows 1 even when the true
//!     count is 0.
//!
//!   BUG 2 (bamqc.rs ~133-139, ~200-211): the NRF denominator
//!     (`total_pairs = paired_read_count / 2`) counts ALL segmented reads,
//!     including reads whose mate lies on another chromosome, while the
//!     M1 / M_distinct numerators are built only from within-chromosome
//!     read1+read2 joins. An inter-chromosomal pair therefore inflates the
//!     denominator without contributing to M1/M_distinct, understating NRF.

use std::num::NonZeroUsize;
use std::path::PathBuf;

use noodles::bam;
use noodles::bam::bai;
use noodles::core::Position;
use noodles::csi::binning_index::Indexer;
use noodles::sam::alignment::record::cigar::{op::Kind, Op};
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::io::Write as AlignmentWrite;
// `Record` trait brings `alignment_end()`/`alignment_start()`/etc. into scope
// for the concrete `bam::Record` used while building the BAI index.
use noodles::sam::alignment::Record as _;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::map::ReferenceSequence;
use noodles::sam::header::record::value::Map;
use noodles::sam::Header;

use gtars_uniwig::bamqc::{compute_bam_qc, compute_bam_qc_parallel};

// --- SAM flag bit constants (from the SAM spec) ---
const FLAG_PAIRED: u16 = 0x1;
const FLAG_READ1: u16 = 0x40; // first segment
const FLAG_READ2: u16 = 0x80; // last segment

/// A minimal record description used by the test harness to build records.
struct ReadSpec {
    name: &'static str,
    ref_id: usize,
    pos: usize,
    flags: u16,
    /// (mate_ref_id, mate_pos) when this is a paired read.
    mate: Option<(usize, usize)>,
    template_length: i32,
}

/// Build a `RecordBuf` from a `ReadSpec`. All records get a 4bp Match CIGAR
/// and a 4bp sequence so `alignment_end` is well-defined for indexing.
fn build_record(spec: &ReadSpec) -> RecordBuf {
    let mut builder = RecordBuf::builder()
        .set_name(spec.name)
        .set_flags(Flags::from(spec.flags))
        .set_reference_sequence_id(spec.ref_id)
        .set_alignment_start(Position::try_from(spec.pos).unwrap())
        .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
        .set_sequence(b"ACGT".as_slice().into())
        .set_template_length(spec.template_length);

    if let Some((mref, mpos)) = spec.mate {
        builder = builder
            .set_mate_reference_sequence_id(mref)
            .set_mate_alignment_start(Position::try_from(mpos).unwrap());
    }

    builder.build()
}

/// Serialize a header + records to an in-memory BGZF BAM byte buffer.
fn write_bam_bytes(header: &Header, records: &[RecordBuf]) -> Vec<u8> {
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_header(header).expect("write header");
    for record in records {
        writer
            .write_alignment_record(header, record)
            .expect("write record");
    }
    // `into_inner()` returns the bgzf::Writer; `finish()` flushes the BGZF EOF.
    writer.into_inner().finish().expect("finish bgzf")
}

/// Build a BAI index for an in-memory BAM byte buffer by re-reading it and
/// recording each record's BGZF virtual-position chunk. This mirrors noodles'
/// own internal `index()` test helper.
fn build_bai_index(bam_bytes: &[u8]) -> bai::Index {
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;

    let mut reader = bam::io::Reader::new(bam_bytes);
    let header = reader.read_header().expect("read header for indexing");

    let mut indexer = Indexer::default();
    let mut chunk_start = reader.get_ref().virtual_position();

    let mut record = bam::Record::default();
    while reader.read_record(&mut record).expect("read record") != 0 {
        let chunk_end = reader.get_ref().virtual_position();

        let alignment_context = match (
            record.reference_sequence_id().transpose().unwrap(),
            record.alignment_start().transpose().unwrap(),
            record.alignment_end().transpose().unwrap(),
        ) {
            (Some(id), Some(start), Some(end)) => {
                let is_mapped = !record.flags().is_unmapped();
                Some((id, start, end, is_mapped))
            }
            _ => None,
        };

        let chunk = Chunk::new(chunk_start, chunk_end);
        indexer.add_record(alignment_context, chunk).expect("add record");

        chunk_start = chunk_end;
    }

    indexer.build(header.reference_sequences().len())
}

/// Write a coordinate-sorted BGZF BAM + its `.bai` index into a temp dir and
/// return the path to the BAM (the temp dir is returned too so it lives long
/// enough). Records MUST already be sorted by (ref_id, pos).
fn write_indexed_bam(header: &Header, records: &[RecordBuf]) -> (tempfile::TempDir, PathBuf) {
    let dir = tempfile::tempdir().expect("tempdir");
    let bam_path = dir.path().join("audit.bam");
    let bai_path = dir.path().join("audit.bam.bai");

    let bytes = write_bam_bytes(header, records);
    std::fs::write(&bam_path, &bytes).expect("write bam file");

    let index = build_bai_index(&bytes);
    bai::write(&bai_path, &index).expect("write bai");

    (dir, bam_path)
}

/// Two-chromosome header (sq0, sq1), each 1000bp long.
fn two_chrom_header() -> Header {
    Header::builder()
        .add_reference_sequence(
            "sq0",
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1000).unwrap()),
        )
        .add_reference_sequence(
            "sq1",
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(1000).unwrap()),
        )
        .build()
}

/// Sanity check that we really did write a queryable, indexed BAM and the
/// public API can read it. If THIS fails, the index-writing approach is
/// blocked and the bug verdicts below are unreliable.
#[ignore = "audit reproduction (PR #256) — run with --ignored"]
#[test]
fn test_harness_indexed_bam_is_queryable() {
    let header = two_chrom_header();
    // Two distinct single-end reads on sq0 -> m1 == 2, m2 == 0 expected
    // (but we only check the harness reads *something* here).
    let records = vec![
        build_record(&ReadSpec {
            name: "r1",
            ref_id: 0,
            pos: 10,
            flags: 0,
            mate: None,
            template_length: 0,
        }),
        build_record(&ReadSpec {
            name: "r2",
            ref_id: 0,
            pos: 200,
            flags: 0,
            mate: None,
            template_length: 0,
        }),
    ];

    let (_dir, bam_path) = write_indexed_bam(&header, &records);
    let result = compute_bam_qc(&bam_path).expect("compute_bam_qc on indexed BAM");
    assert!(
        result.total_reads > 0,
        "harness produced an unreadable/empty BAM: {result:?}"
    );
    eprintln!("[harness] result = {result:?}");
}

/// Helper: build a BAM where NO position is observed exactly twice, using
/// single-end reads only.
///
/// Position key for an unpaired read is `(pos, qlen, 0, 0)` where qlen is the
/// sequence length (always 4 here). So multiplicity is controlled purely by
/// `pos`. We create:
///   - one distinct position  (multiplicity 1)
///   - one position seen 3 times (multiplicity 3)
/// No position has multiplicity 2 -> true m2 == 0.
fn no_multiplicity_two_records() -> (Header, Vec<RecordBuf>) {
    let header = two_chrom_header();
    let mut specs: Vec<ReadSpec> = Vec::new();

    // multiplicity 1 at pos 10
    specs.push(ReadSpec { name: "a", ref_id: 0, pos: 10, flags: 0, mate: None, template_length: 0 });

    // multiplicity 3 at pos 100 (three reads share identical (pos,qlen))
    specs.push(ReadSpec { name: "b1", ref_id: 0, pos: 100, flags: 0, mate: None, template_length: 0 });
    specs.push(ReadSpec { name: "b2", ref_id: 0, pos: 100, flags: 0, mate: None, template_length: 0 });
    specs.push(ReadSpec { name: "b3", ref_id: 0, pos: 100, flags: 0, mate: None, template_length: 0 });

    // records must be coordinate-sorted by (ref_id, pos): 10 then 100,100,100
    let records: Vec<RecordBuf> = specs.iter().map(build_record).collect();
    (header, records)
}

/// BUG 1: reported `m2` should be the TRUE count of multiplicity-2 positions
/// (here 0), but the code returns `m2.max(1)`.
#[ignore = "audit reproduction (PR #256) — run with --ignored"]
#[test]
fn test_bug1_m2_forced_to_one_serial() {
    let (header, records) = no_multiplicity_two_records();
    let (_dir, bam_path) = write_indexed_bam(&header, &records);

    let result = compute_bam_qc(&bam_path).expect("compute_bam_qc");
    eprintln!("[bug1 serial] result = {result:?}");

    // Sanity: we engineered distinct positions {pos=10 (x1), pos=100 (x3)}.
    // distinct == 2, m1 == 1 (only pos=10 has count 1), and NO position has
    // count == 2.
    assert_eq!(result.distinct, 2, "expected 2 distinct positions");
    assert_eq!(result.m1, 1, "expected exactly one multiplicity-1 position");

    // The TRUE number of positions seen exactly twice is 0.
    assert_eq!(
        result.m2, 0,
        "BUG 1 CONFIRMED if this fails: reported m2 should be 0 (no \
         multiplicity-2 positions) but code returns m2.max(1) == 1"
    );
}

/// BUG 1, parallel path (bamqc.rs:212).
#[ignore = "audit reproduction (PR #256) — run with --ignored"]
#[test]
fn test_bug1_m2_forced_to_one_parallel() {
    let (header, records) = no_multiplicity_two_records();
    let (_dir, bam_path) = write_indexed_bam(&header, &records);

    let result = compute_bam_qc_parallel(&bam_path, 2).expect("compute_bam_qc_parallel");
    eprintln!("[bug1 parallel] result = {result:?}");

    assert_eq!(result.distinct, 2, "expected 2 distinct positions");
    assert_eq!(result.m1, 1, "expected exactly one multiplicity-1 position");
    assert_eq!(
        result.m2, 0,
        "BUG 1 CONFIRMED if this fails (parallel path, bamqc.rs:212): \
         reported m2 should be 0 but code returns m2.max(1) == 1"
    );
}

/// Build a control BAM with N within-chromosome proper pairs (all on sq0),
/// each pair at a distinct position so every pair is multiplicity-1.
/// Returns (header, sorted records).
///
/// Each pair contributes:
///   - read1 (FLAG_PAIRED|FLAG_READ1) at pos p
///   - read2 (FLAG_PAIRED|FLAG_READ2) at pos p+50
/// joined by qname into key (p, tlen1, p+50, tlen2) with multiplicity 1.
fn paired_records(within_chrom_pairs: usize, split_pairs: usize) -> (Header, Vec<RecordBuf>) {
    let header = two_chrom_header();
    let mut specs: Vec<ReadSpec> = Vec::new();

    for i in 0..within_chrom_pairs {
        let p1 = 10 + i * 100;
        let p2 = p1 + 50;
        let tlen = 54i32; // arbitrary but consistent
        specs.push(ReadSpec {
            name: leak_name(&format!("pair{i}")),
            ref_id: 0,
            pos: p1,
            flags: FLAG_PAIRED | FLAG_READ1,
            mate: Some((0, p2)),
            template_length: tlen,
        });
        specs.push(ReadSpec {
            name: leak_name(&format!("pair{i}")),
            ref_id: 0,
            pos: p2,
            flags: FLAG_PAIRED | FLAG_READ2,
            mate: Some((0, p1)),
            template_length: -tlen,
        });
    }

    // Inter-chromosomal pairs: read1 on sq0, read2 (mate) on sq1. Both reads
    // are mapped and segmented, so each increments `paired_read_count` on ITS
    // OWN chromosome -> total_pairs grows. But the read1/read2 qname join only
    // happens WITHIN a chromosome, so these never match and contribute NOTHING
    // to m1/m_distinct.
    //
    // NOTE: total_pairs is summed per-chromosome as `paired_read_count / 2`
    // (integer division). A *single* split pair adds one orphan segmented read
    // per chromosome, which truncation hides. We therefore add split pairs in
    // even multiples per chromosome so the inflation actually surfaces.
    for i in 0..split_pairs {
        let p = 800 + i * 20;
        specs.push(ReadSpec {
            name: leak_name(&format!("split{i}")),
            ref_id: 0,
            pos: p,
            flags: FLAG_PAIRED | FLAG_READ1,
            mate: Some((1, p)),
            template_length: 0,
        });
        specs.push(ReadSpec {
            name: leak_name(&format!("split{i}")),
            ref_id: 1,
            pos: p,
            flags: FLAG_PAIRED | FLAG_READ2,
            mate: Some((0, p)),
            template_length: 0,
        });
    }

    // Sort by (ref_id, pos) for coordinate-sorted output.
    specs.sort_by_key(|s| (s.ref_id, s.pos));
    let records: Vec<RecordBuf> = specs.iter().map(build_record).collect();
    (header, records)
}

/// `set_name` takes anything Into<Vec<u8>> but our ReadSpec holds a &'static
/// str. We synthesize unique names at runtime, so leak them to get 'static.
/// (Only used in tests; the small leak is acceptable.)
fn leak_name(s: &str) -> &'static str {
    Box::leak(s.to_string().into_boxed_str())
}

/// BUG 2: inter-chromosomal pair inflates the NRF denominator (total_pairs)
/// without contributing to m1/distinct, so NRF drops even though the read
/// population that builds m1 is unchanged.
#[ignore = "audit reproduction (PR #256) — run with --ignored"]
#[test]
fn test_bug2_interchromosomal_pair_inflates_denominator() {
    // Control: 2 within-chromosome pairs, no split pairs.
    let (h_ctrl, r_ctrl) = paired_records(2, 0);
    let (_d1, p_ctrl) = write_indexed_bam(&h_ctrl, &r_ctrl);
    let ctrl = compute_bam_qc(&p_ctrl).expect("ctrl qc");
    eprintln!("[bug2 control]   = {ctrl:?}");

    // Treatment: same 2 within-chromosome pairs PLUS two inter-chromosomal
    // pairs (even per-chromosome so the `/2` truncation doesn't mask them).
    let (h_tx, r_tx) = paired_records(2, 2);
    let (_d2, p_tx) = write_indexed_bam(&h_tx, &r_tx);
    let tx = compute_bam_qc(&p_tx).expect("tx qc");
    eprintln!("[bug2 treatment] = {tx:?}");

    // The within-chromosome joins are identical between control and treatment,
    // so m1 and distinct must be unchanged.
    assert_eq!(ctrl.m1, tx.m1, "m1 should be unchanged by the split pair");
    assert_eq!(
        ctrl.distinct, tx.distinct,
        "distinct should be unchanged by the split pair"
    );
    assert_eq!(ctrl.m1, 2, "expected 2 multiplicity-1 within-chrom pairs");

    // total_reads (== effective_total == total_pairs for paired data) is the
    // NRF denominator. The two split pairs add 2 orphan segmented reads to sq0
    // (paired_read_count: 4 -> 6 -> num_pairs 2 -> 3) and 2 to sq1
    // (paired_read_count 0 -> 2 -> num_pairs 0 -> 1), so total_pairs 2 -> 4,
    // even though none of these reads join into an m1/distinct entry.
    assert!(
        tx.total_reads > ctrl.total_reads,
        "BUG 2: inter-chromosomal pairs should inflate total_pairs (denominator). \
         ctrl.total_reads={}, tx.total_reads={}",
        ctrl.total_reads,
        tx.total_reads
    );

    // Because the numerator (m1) is unchanged but the denominator grew, NRF
    // must DECREASE solely due to the inter-chromosomal pair: numerator and
    // denominator populations differ.
    assert!(
        tx.nrf < ctrl.nrf,
        "BUG 2 CONFIRMED if NRF drops with no change to m1/distinct: \
         ctrl.nrf={}, tx.nrf={} (denominator inflated by cross-chrom pair)",
        ctrl.nrf,
        tx.nrf
    );
}
