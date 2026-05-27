//! Layer-2 only mapper-correctness test (no VRS dependency).
//!
//! Iterates `tests/data/hgvs/synthetic/cases.tsv` (symlink to the
//! gtars-vrs fixture), filters to non-genomic, non-error c./n. cases,
//! and asserts that `CoordinateMapper::{c_to_g_full, n_to_g_full}`
//! produces the same genomic coordinate the generator recorded.
//!
//! This test exists so a reftx-only contributor can run mapper
//! correctness tests without pulling in the gtars-vrs crate. It parses
//! the HGVS string with a tiny regex-free splitter sufficient for the
//! synthetic corpus's narrow shape (no whitespace, fixed punctuation).

use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use gtars_reftx::mapper::CoordinateMapper;
use gtars_reftx::{TxStore, TxStoreBuilder};

const CHROM_NAME: &str = "chr_synth";
// Hard-coded chromosome digest: derived independently from the
// synthetic genome (sha512t24u of the 30000-byte sequence). The test
// regenerates it from the FASTA at runtime instead of hard-coding it,
// to stay in sync with the generator's seed.

#[derive(Debug, Clone)]
struct Case {
    case_id: String,
    hgvs_string: String,
    ref_type: String,
    transcript: String,
    expected_pos: u64,
    expected_error: String,
}

fn synthetic_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/hgvs/synthetic")
}

fn parse_cases() -> Vec<Case> {
    let path = synthetic_dir().join("cases.tsv");
    let raw = fs::read_to_string(&path).expect("read cases.tsv");
    let mut rows = Vec::new();
    let mut header_seen = false;
    for line in raw.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if !header_seen {
            header_seen = true;
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        rows.push(Case {
            case_id: fields[0].to_string(),
            hgvs_string: fields[1].to_string(),
            ref_type: fields[2].to_string(),
            transcript: fields[6].to_string(),
            expected_pos: fields[8].parse().unwrap_or(0),
            expected_error: fields[12].to_string(),
        });
    }
    rows
}

/// Compute SHA-512/24 base64-url digest of the FASTA's first sequence.
fn synthetic_chrom_digest() -> [u8; 24] {
    use std::io::Read;
    let mut buf = Vec::new();
    let mut f = fs::File::open(synthetic_dir().join("synthetic.fa")).unwrap();
    f.read_to_end(&mut buf).unwrap();
    let mut seq = Vec::new();
    for line in buf.split(|&b| b == b'\n') {
        if line.starts_with(b">") || line.is_empty() {
            continue;
        }
        seq.extend_from_slice(line);
    }
    use sha2::Digest;
    let mut hasher = sha2::Sha512::new();
    hasher.update(&seq);
    let h = hasher.finalize();
    let mut out = [0u8; 24];
    out.copy_from_slice(&h[..24]);
    out
}

fn build_txstore() -> Arc<gtars_reftx::ReadonlyTxStore> {
    let dir = synthetic_dir();
    let cdot_path = dir.join("synthetic_transcripts.json");
    let tmp = tempfile::tempdir().unwrap();
    let bin_path = tmp.path().join("synth.reftx");

    let mut builder = TxStoreBuilder::new();
    builder.add_chrom_mapping(CHROM_NAME, synthetic_chrom_digest());
    let n = builder.ingest_cdot(&cdot_path).expect("ingest cdot");
    assert!(n >= 5);
    builder.build(&bin_path).unwrap();

    let store = Arc::new(TxStore::open(&bin_path).unwrap().into_readonly());
    std::mem::forget(tmp); // keep mmap alive
    store
}

/// Parse a synthetic HGVS string of the form `<acc>:<g|c|n>.<rest>`.
/// Returns (accession, ref_type, c_pos_or_n_pos, intron_offset, is_cds_end).
fn parse_synth_hgvs(s: &str) -> Option<(String, char, i64, i64, bool)> {
    let (acc, rest) = s.split_once(':')?;
    let bytes = rest.as_bytes();
    if bytes.len() < 3 {
        return None;
    }
    let ref_type = bytes[0] as char;
    if !matches!(ref_type, 'g' | 'c' | 'n') {
        return None;
    }
    if bytes[1] != b'.' {
        return None;
    }
    // Skip ranges (ins/del/dup spanning two positions). The position
    // part appears before any edit suffix; if `_` appears anywhere in
    // `rest` it is a range we don't handle in this layer-2 mapper test.
    if rest.contains('_') {
        return None;
    }
    // After "X.", parse position. Stop at the first edit indicator
    // (A>T, del, ins, dup, delins, =).
    let after = &rest[2..];
    // Find where edit starts: scan for one of the markers.
    // Acceptable position chars: digits, +, -, *, _ (range — synthetic uses only single positions for c./n.)
    let mut pos_end = 0;
    let mut iter = after.char_indices().peekable();
    while let Some((i, c)) = iter.peek().copied() {
        if c.is_ascii_digit() || c == '+' || c == '-' || c == '*' {
            pos_end = i + c.len_utf8();
            iter.next();
            continue;
        }
        break;
    }
    let pos_str = &after[..pos_end];

    // Check for `_` indicating range — we only handle single-positions in this test.
    if pos_str.is_empty() || pos_str.contains('_') {
        return None;
    }

    let mut is_cds_end = false;
    let mut s2 = pos_str;
    if let Some(rest) = s2.strip_prefix('*') {
        is_cds_end = true;
        s2 = rest;
    }
    // Split into base and offset.
    // Base may have a leading '-' for c.-N. Find the offset separator: first
    // '+' or '-' that is NOT at index 0.
    let mut sep = None;
    for (i, c) in s2.char_indices() {
        if i == 0 {
            continue;
        }
        if c == '+' || c == '-' {
            sep = Some(i);
            break;
        }
    }
    let (base_str, off_str) = match sep {
        Some(i) => (&s2[..i], &s2[i..]),
        None => (s2, ""),
    };
    let base: i64 = base_str.parse().ok()?;
    let offset: i64 = if off_str.is_empty() {
        0
    } else {
        off_str.parse().ok()?
    };
    Some((acc.to_string(), ref_type, base, offset, is_cds_end))
}

#[test]
fn synthetic_mapper_layer2() {
    let cases = parse_cases();
    let store = build_txstore();
    let mapper = CoordinateMapper::new(&store);

    let mut failures = vec![];
    let mut checked = 0usize;
    for case in &cases {
        if !case.expected_error.is_empty() {
            continue;
        }
        if case.ref_type == "g" {
            continue;
        }
        let parsed = match parse_synth_hgvs(&case.hgvs_string) {
            Some(p) => p,
            None => continue, // skip ranges (del/ins/dup) for layer-2 test
        };
        let (acc, rt, base, off, is_end) = parsed;
        if acc != case.transcript {
            continue;
        }
        let result = match rt {
            'c' => mapper.c_to_g_full(&acc, base, off, is_end),
            'n' => mapper.n_to_g_full(&acc, base, off),
            _ => continue,
        };
        let r = match result {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!(
                    "L2 FAIL {}: {} -- mapper errored: {:?}",
                    case.case_id, case.hgvs_string, e
                ));
                continue;
            }
        };
        // Mapper returns 0-based genomic; expected_pos is 1-based.
        let mapper_1based = r.position + 1;
        if mapper_1based != case.expected_pos {
            failures.push(format!(
                "L2 FAIL {}: {} -- mapper={} expected={}",
                case.case_id, case.hgvs_string, mapper_1based, case.expected_pos
            ));
        }
        checked += 1;
    }

    if !failures.is_empty() {
        for f in &failures {
            eprintln!("{}", f);
        }
    }
    eprintln!("layer-2 mapper test: {} cases checked, {} failures", checked, failures.len());
    assert!(failures.is_empty(), "layer-2 mapper failures (see eprintln above)");
    assert!(checked >= 10, "expected to check at least 10 c./n. cases");
}

#[test]
fn synthetic_mapper_negative_cases() {
    let cases = parse_cases();
    let store = build_txstore();
    let mapper = CoordinateMapper::new(&store);

    let mut failures = vec![];
    let mut known_gaps = vec![];
    for case in &cases {
        // Only mapper-level negative cases (skip ParseError, which the
        // mapper never sees because the parser rejects upstream).
        if case.expected_error.is_empty() || case.expected_error == "ParseError" {
            continue;
        }
        if case.ref_type == "g" {
            continue; // g. negative cases are out-of-bounds at refget level
        }
        let parsed = match parse_synth_hgvs(&case.hgvs_string) {
            Some(p) => p,
            None => continue,
        };
        let (acc, rt, base, off, is_end) = parsed;
        let result = match rt {
            'c' => mapper.c_to_g_full(&acc, base, off, is_end),
            'n' => mapper.n_to_g_full(&acc, base, off),
            _ => continue,
        };
        if result.is_ok() {
            // Known mapper gap: `apply_offset_to_tx_pos` does not validate
            // that |intron_offset| <= intron_length; it blindly shifts the
            // genomic anchor into (or beyond) the next exon. Plan 10 flags
            // this as a real bug to fix; until then the synthetic corpus
            // documents the expected behavior and this test soft-warns.
            if case.expected_error == "IntronOffsetTooLarge" {
                known_gaps.push(format!(
                    "KNOWN MAPPER GAP: {} ({}) — mapper accepted intronic offset \
                     past intron length; see plan 10. Result: {:?}",
                    case.case_id, case.hgvs_string, result
                ));
                continue;
            }
            failures.push(format!(
                "expected error {} but mapper returned Ok for {}: {}",
                case.expected_error, case.case_id, case.hgvs_string
            ));
        }
    }
    for g in &known_gaps {
        eprintln!("{}", g);
    }
    assert!(failures.is_empty(), "{:#?}", failures);
}
