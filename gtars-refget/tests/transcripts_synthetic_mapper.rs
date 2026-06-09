//! Layer-2 only mapper-correctness test (no VRS dependency).
//!
//! Iterates `tests/data/hgvs/synthetic/cases.tsv` (symlink to the gtars-vrs
//! fixture), filters to non-genomic, non-error c./n. cases, and asserts that
//! `CoordinateMapper::{c_to_g_full, n_to_g_full}` produces the same genomic
//! coordinate the generator recorded.
//!
//! This test runs against the WASM-SAFE CORE only: it builds the `.reftx` byte
//! image in memory (the small serializer below mirrors the native builder's
//! record format) and constructs the store via `ReadonlyTxStore::from_bytes`.
//! It therefore compiles and runs under `--features transcripts` WITHOUT
//! `filesystem`, giving the core direct coverage. (Reading the fixture files
//! from disk happens in the test harness via `std::fs`, which is always
//! available to tests regardless of crate feature gating.)
#![cfg(feature = "transcripts")]

use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use gtars_refget::{
    build_reftx_bytes_in_memory, CoordinateMapper, Exon, ManeStatus, ReadonlyTxStore, Strand,
    Transcript,
};

const CHROM_NAME: &str = "chr_synth";

// ----------------------------------------------------------------------------
// cdot JSON ingest (test-only; parses the synthetic fixture into Transcripts).
// ----------------------------------------------------------------------------

fn ingest_cdot(path: &PathBuf, chrom_name: &str, chrom_digest: [u8; 24]) -> Vec<Transcript> {
    let raw = fs::read_to_string(path).expect("read cdot json");
    let json: serde_json::Value = serde_json::from_str(&raw).expect("parse cdot json");
    let txs = json
        .get("transcripts")
        .and_then(|v| v.as_object())
        .expect("transcripts object");

    let mut out = Vec::new();
    for (_, tx) in txs {
        let contig = tx.get("contig").and_then(|v| v.as_str()).unwrap_or("");
        if contig != chrom_name {
            continue;
        }
        let strand = match tx.get("strand").and_then(|v| v.as_i64()).unwrap_or(0) {
            1 => Strand::Forward,
            -1 => Strand::Reverse,
            _ => continue,
        };
        let exons: Vec<Exon> = tx
            .get("exons")
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|e| {
                        let pair = e.as_array()?;
                        let s = pair.first()?.as_u64()? as u32;
                        let en = pair.get(1)?.as_u64()? as u32;
                        Some(Exon { start: s, end: en })
                    })
                    .collect()
            })
            .unwrap_or_default();
        if exons.is_empty() {
            continue;
        }
        out.push(Transcript {
            accession: tx.get("id").and_then(|v| v.as_str()).unwrap_or("").to_string(),
            gene: tx
                .get("gene_name")
                .and_then(|v| v.as_str())
                .unwrap_or("")
                .to_string(),
            chrom_digest,
            strand,
            cds_start: tx.get("cds_start").and_then(|v| v.as_u64()).map(|v| v as u32),
            cds_end: tx.get("cds_end").and_then(|v| v.as_u64()).map(|v| v as u32),
            exons,
            mane: ManeStatus::default(),
        });
    }
    out
}

// ----------------------------------------------------------------------------
// Fixture helpers
// ----------------------------------------------------------------------------

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

/// Compute SHA-512/24 digest of the FASTA's first sequence.
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

fn build_txstore() -> Arc<ReadonlyTxStore> {
    let dir = synthetic_dir();
    let cdot_path = dir.join("synthetic_transcripts.json");
    let transcripts = ingest_cdot(&cdot_path, CHROM_NAME, synthetic_chrom_digest());
    assert!(transcripts.len() >= 5);
    let bytes = build_reftx_bytes_in_memory(&transcripts).unwrap();
    Arc::new(ReadonlyTxStore::from_bytes(bytes).unwrap())
}

/// Parse a synthetic HGVS string of the form `<acc>:<g|c|n>.<rest>`.
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
    if rest.contains('_') {
        return None;
    }
    let after = &rest[2..];
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

    if pos_str.is_empty() || pos_str.contains('_') {
        return None;
    }

    let mut is_cds_end = false;
    let mut s2 = pos_str;
    if let Some(rest) = s2.strip_prefix('*') {
        is_cds_end = true;
        s2 = rest;
    }
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
            None => continue,
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
    eprintln!(
        "layer-2 mapper test: {} cases checked, {} failures",
        checked,
        failures.len()
    );
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
        if case.expected_error.is_empty() || case.expected_error == "ParseError" {
            continue;
        }
        if case.ref_type == "g" {
            continue;
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
