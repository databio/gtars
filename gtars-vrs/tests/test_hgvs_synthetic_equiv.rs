//! HGVS equivalence-group harness — layer 3 (plan 12).
//!
//! Each group in `fixtures/equivalence_groups.json` lists multiple representations
//! of the same biological variant (g., c., n., or VCF tuple). All members of a
//! group MUST collapse to the same VRS ID — and equal the group's
//! `expected_vrs_id` (computed once per group at fixture-generation time
//! via the canonical VCF→VRS path).
//!
//! Failure-attribution semantics: on mismatch within a group, the harness
//! does NOT panic on the first diverging member. It buckets all members of
//! the group by computed VRS ID (or error message) and prints a structured
//! report so the developer immediately sees which forms agree and which
//! diverge. All groups are run; one panic at the end summarizes everything.
//!
//! Cross-link to plan 10 (`cases.tsv`): each `expected_vrs_id` in
//! `equivalence_groups.json` must also appear somewhere in `cases.tsv`,
//! since groups are built from biological tuples that already exist as
//! per-row cases. Drift between the files trips a sanity assertion.

use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use gtars_refget::store::RefgetStore;
use gtars_reftx::provider::ReftxProvider;
use gtars_reftx::{TxStore, TxStoreBuilder};
use gtars_vrs::digest::DigestWriter;
use gtars_vrs::hgvs::bridge::hgvs_str_to_vrs_id;
use gtars_vrs::normalize::normalize;
use serde::Deserialize;

const CHROM_NAME: &str = "chr_synth";

#[derive(Debug, Clone, Deserialize)]
#[allow(dead_code)] // `transcript` and `notes` are loaded for diagnostics, not asserted on
struct EquivRow {
    group_id: String,
    member_kind: String,
    expression: String,
    #[serde(default)]
    transcript: String,
    expected_vrs_id: String,
    #[serde(default)]
    notes: String,
}

#[derive(Debug, Clone)]
struct EquivGroup {
    group_id: String,
    expected_vrs_id: String,
    members: Vec<EquivRow>,
}

fn synthetic_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/hgvs/synthetic")
}

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn parse_equivalence_groups() -> Vec<EquivGroup> {
    let path = fixtures_dir().join("equivalence_groups.json");
    let raw = fs::read_to_string(&path).expect("read equivalence_groups.json");
    let rows: Vec<EquivRow> = serde_json::from_str(&raw).expect("parse equivalence_groups.json");

    // Group rows by group_id while preserving first-seen order.
    let mut order: Vec<String> = Vec::new();
    let mut by_id: BTreeMap<String, Vec<EquivRow>> = BTreeMap::new();
    for r in rows {
        if !by_id.contains_key(&r.group_id) {
            order.push(r.group_id.clone());
        }
        by_id.entry(r.group_id.clone()).or_default().push(r);
    }

    let mut out = Vec::new();
    for gid in order {
        let members = by_id.remove(&gid).unwrap();
        // Validate: ≥2 members, single shared expected_vrs_id, valid kinds.
        assert!(
            members.len() >= 2,
            "group {} has fewer than 2 members",
            gid
        );
        let expected = members[0].expected_vrs_id.clone();
        assert!(
            !expected.is_empty(),
            "group {}: expected_vrs_id is empty (regenerate fixtures with --no-vrs OFF)",
            gid
        );
        for m in &members {
            assert_eq!(
                m.expected_vrs_id, expected,
                "group {} has inconsistent expected_vrs_id across rows",
                gid
            );
            match m.member_kind.as_str() {
                "hgvs_g" | "hgvs_c" | "hgvs_n" | "vcf" => {}
                other => panic!("group {}: unknown member_kind {:?}", gid, other),
            }
        }
        out.push(EquivGroup {
            group_id: gid,
            expected_vrs_id: expected,
            members,
        });
    }
    out
}

/// Cross-link sanity check: every `expected_vrs_id` in
/// `equivalence_groups.tsv` should appear in `cases.tsv` too. (Groups are
/// built from biological tuples already covered by per-row cases.)
fn assert_cases_cross_link(groups: &[EquivGroup]) {
    let cases_path = synthetic_dir().join("cases.tsv");
    let raw = fs::read_to_string(&cases_path).expect("read cases.tsv");
    let mut case_ids: HashSet<String> = HashSet::new();
    let mut header_seen = false;
    for line in raw.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if !header_seen {
            header_seen = true;
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 12 {
            let id = f[11].trim();
            if !id.is_empty() {
                case_ids.insert(id.to_string());
            }
        }
    }
    let mut missing: Vec<(String, String)> = Vec::new();
    for g in groups {
        if !case_ids.contains(&g.expected_vrs_id) {
            missing.push((g.group_id.clone(), g.expected_vrs_id.clone()));
        }
    }
    assert!(
        missing.is_empty(),
        "equivalence groups reference VRS IDs not in cases.tsv: {:?}",
        missing
    );
}

fn build_fixture() -> (RefgetStore, Arc<gtars_reftx::ReadonlyTxStore>, ReftxProvider, String) {
    let dir = synthetic_dir();
    let fasta = dir.join("synthetic.fa");

    let store_tmpdir = tempfile::tempdir().unwrap();
    let store_path = store_tmpdir.path().join("store");
    let mut store = RefgetStore::on_disk(&store_path).expect("create store");
    store.disable_encoding();
    store.set_quiet(true);
    store
        .add_sequence_collection_from_fasta(&fasta, gtars_refget::store::FastaImportOptions::new())
        .expect("import synthetic.fa");
    store.load_all_sequences().expect("load sequences");
    std::mem::forget(store_tmpdir);

    let collection_digest = store
        .iter_collections()
        .next()
        .map(|c| c.metadata.digest.clone())
        .expect("collection");

    let coll = store.get_collection(&collection_digest).unwrap();
    let mut chrom_digest_b64 = String::new();
    for r in &coll.sequences {
        let m = r.metadata();
        if m.name == CHROM_NAME {
            chrom_digest_b64 = m.sha512t24u.clone();
        }
    }
    assert!(!chrom_digest_b64.is_empty());

    let digest_bytes = base64_url::decode(&chrom_digest_b64).expect("decode digest");
    let mut digest_arr = [0u8; 24];
    digest_arr.copy_from_slice(&digest_bytes);

    let cdot_path = dir.join("synthetic_transcripts.json");
    let tmpdir = tempfile::tempdir().unwrap();
    let bin_path = tmpdir.path().join("synth.reftx");

    let mut builder = TxStoreBuilder::new();
    builder.add_chrom_mapping(CHROM_NAME, digest_arr);
    builder.ingest_cdot(&cdot_path).expect("ingest cdot");
    builder.build(&bin_path).expect("build reftx");

    let txstore = Arc::new(TxStore::open(&bin_path).unwrap().into_readonly());
    std::mem::forget(tmpdir);
    let provider = ReftxProvider::new(Arc::clone(&txstore));
    (store, txstore, provider, collection_digest)
}

fn resolve_seq_bytes(store: &mut RefgetStore, coll: &str, name: &str) -> Option<Vec<u8>> {
    let c = store.get_collection(coll).ok()?;
    let sha = c
        .sequences
        .iter()
        .find(|r| r.metadata().name == name)
        .map(|r| r.metadata().sha512t24u.clone())?;
    store.ensure_decoded(&sha).ok()?;
    store.sequence_bytes(&sha).map(|b| b.to_vec())
}

fn resolve_raw_digest(store: &mut RefgetStore, coll: &str, name: &str) -> Option<String> {
    let c = store.get_collection(coll).ok()?;
    c.sequences
        .iter()
        .find(|r| r.metadata().name == name)
        .map(|r| r.metadata().sha512t24u.clone())
}

/// Compute VRS ID for a single member. For HGVS forms, this routes through
/// `hgvs_str_to_vrs_id`. For `vcf`, it parses `chrom:pos:ref:alt` and
/// applies the canonical normalize+digest path directly (mirroring the
/// generator).
fn compute_member_vrs_id(
    member: &EquivRow,
    store: &mut RefgetStore,
    provider: &ReftxProvider,
    coll: &str,
) -> Result<String, String> {
    match member.member_kind.as_str() {
        "hgvs_g" | "hgvs_c" | "hgvs_n" => {
            hgvs_str_to_vrs_id(&member.expression, provider, store, coll)
                .map_err(|e| format!("{:?}", e))
        }
        "vcf" => {
            // Parse "chrom:pos:ref:alt" — note that ref or alt may be empty
            // (ins / del). Use splitn(4, ':') so we don't mis-split if a
            // ref/alt contains a colon (it shouldn't, but defensive).
            let parts: Vec<&str> = member.expression.splitn(4, ':').collect();
            if parts.len() != 4 {
                return Err(format!(
                    "malformed vcf expression {:?}",
                    member.expression
                ));
            }
            let chrom = parts[0];
            let pos: u64 = parts[1].parse().map_err(|e| format!("bad pos: {:?}", e))?;
            let ref_str = parts[2];
            let alt_str = parts[3];

            let raw = resolve_seq_bytes(store, coll, chrom)
                .ok_or_else(|| format!("missing seq bytes for {}", chrom))?;
            let raw_digest = resolve_raw_digest(store, coll, chrom)
                .ok_or_else(|| format!("missing raw digest for {}", chrom))?;
            let sq = format!("SQ.{raw_digest}");
            let pos_ib = pos.saturating_sub(1);
            let norm = normalize(&raw, pos_ib, ref_str.as_bytes(), alt_str.as_bytes())
                .map_err(|e| format!("normalize errored: {:?}", e))?;
            let mut writer = DigestWriter::new();
            Ok(writer.allele_identifier_literal(
                &sq,
                norm.start,
                norm.end,
                std::str::from_utf8(&norm.allele).unwrap(),
            ))
        }
        other => Err(format!("unknown member_kind {:?}", other)),
    }
}

#[derive(Debug)]
struct GroupResult {
    group_id: String,
    expected_vrs_id: String,
    /// Per-member: (member_kind, expression, computed_id_or_err)
    member_results: Vec<(String, String, Result<String, String>)>,
}

impl GroupResult {
    fn ok(&self) -> bool {
        for (_, _, r) in &self.member_results {
            match r {
                Err(_) => return false,
                Ok(id) if id != &self.expected_vrs_id => return false,
                _ => {}
            }
        }
        true
    }

    /// Bucket members by computed VRS ID (or error message). Returns
    /// (bucket_key, member_descriptions) pairs sorted for stable output.
    fn buckets(&self) -> Vec<(String, Vec<String>)> {
        let mut by_key: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for (kind, expr, r) in &self.member_results {
            let key = match r {
                Ok(id) => format!("id={id}"),
                Err(e) => format!("ERROR: {e}"),
            };
            by_key
                .entry(key)
                .or_default()
                .push(format!("{kind}[{expr}]"));
        }
        by_key.into_iter().collect()
    }

    fn render_failure(&self) -> String {
        let mut s = String::new();
        s.push_str(&format!(
            "EQUIVALENCE FAILURE: group_id = {}\n",
            self.group_id
        ));
        s.push_str(&format!(
            "  expected_vrs_id (from TSV): {}\n",
            self.expected_vrs_id
        ));
        for (i, (key, members)) in self.buckets().into_iter().enumerate() {
            s.push_str(&format!("  Bucket {} ({}):\n", i + 1, key));
            for m in members {
                s.push_str(&format!("    - {}\n", m));
            }
        }
        s
    }
}

#[test]
fn synthetic_equivalence_groups() {
    let groups = parse_equivalence_groups();
    assert!(
        groups.len() >= 10,
        "expected >= 10 equivalence groups, got {}",
        groups.len()
    );

    assert_cases_cross_link(&groups);

    let (mut store, _txstore, provider, coll) = build_fixture();

    let mut results: Vec<GroupResult> = Vec::with_capacity(groups.len());
    for g in &groups {
        let mut member_results = Vec::with_capacity(g.members.len());
        for m in &g.members {
            let r = compute_member_vrs_id(m, &mut store, &provider, &coll);
            member_results.push((m.member_kind.clone(), m.expression.clone(), r));
        }
        results.push(GroupResult {
            group_id: g.group_id.clone(),
            expected_vrs_id: g.expected_vrs_id.clone(),
            member_results,
        });
    }

    let failed: Vec<&GroupResult> = results.iter().filter(|r| !r.ok()).collect();
    let pass = results.len() - failed.len();
    eprintln!(
        "\nsynthetic equivalence: {} groups, {} pass, {} fail",
        results.len(),
        pass,
        failed.len(),
    );

    if !failed.is_empty() {
        eprintln!("\n=== Equivalence-group failures ({}) ===", failed.len());
        for r in &failed {
            eprintln!("{}", r.render_failure());
        }
        panic!(
            "{} equivalence group(s) failed; see structured reports above",
            failed.len()
        );
    }
}

#[test]
fn equivalence_groups_tsv_well_formed() {
    let groups = parse_equivalence_groups();
    assert!(groups.len() >= 10);
    let mut ids = HashSet::new();
    for g in &groups {
        assert!(
            ids.insert(g.group_id.clone()),
            "duplicate group_id: {}",
            g.group_id
        );
        // At least one non-vcf and one vcf member, OR at least two HGVS-form
        // members, or at least g + vcf — i.e. the group must be non-trivially
        // multi-form. Easiest invariant: not all members of identical kind.
        let kinds: HashSet<_> = g.members.iter().map(|m| m.member_kind.clone()).collect();
        assert!(
            kinds.len() >= 2,
            "group {} has only one member_kind: {:?}",
            g.group_id,
            kinds
        );
    }
}
