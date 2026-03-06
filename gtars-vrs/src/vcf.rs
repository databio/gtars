//! VCF parsing and VRS ID computation.
//!
//! Reads a VCF file (plain text or gzipped/bgzf), normalizes each variant,
//! and computes GA4GH VRS Allele identifiers using a fast zero-allocation
//! digest path (no serde_json per variant).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore};

use crate::digest::DigestWriter;
use crate::normalize::normalize;

/// Result of computing a VRS identifier for a single VCF variant.
#[derive(Debug, Clone)]
pub struct VrsResult {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub vrs_id: String,
}

/// Open a VCF file, auto-detecting gzip/bgzf compression.
fn open_vcf(path: &str) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).context(format!("Failed to open VCF: {}", path))?;
    let capacity = 256 * 1024; // 256KB buffer for large VCF files
    if path.ends_with(".gz") || path.ends_with(".bgz") {
        Ok(Box::new(BufReader::with_capacity(
            capacity,
            MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(capacity, file)))
    }
}

/// Build the name→digest mapping for a collection.
fn build_name_to_digest(
    store: &mut RefgetStore,
    collection_digest: &str,
) -> Result<HashMap<String, String>> {
    let collection = store
        .get_collection(collection_digest)
        .context("Failed to get sequence collection")?;

    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    for seq_record in &collection.sequences {
        let meta = seq_record.metadata();
        name_to_digest.insert(meta.name.clone(), meta.sha512t24u.clone());
    }
    Ok(name_to_digest)
}

/// Read the next non-empty line from a VCF reader. Returns false at EOF.
fn read_vcf_line(reader: &mut dyn BufRead, buf: &mut String) -> Result<bool> {
    buf.clear();
    match reader.read_line(buf) {
        Ok(0) => Ok(false),
        Ok(_) => Ok(true),
        Err(e) => {
            // BGZF files end with an empty gzip block that MultiGzDecoder
            // may interpret as an invalid header. Treat as EOF.
            if e.kind() == std::io::ErrorKind::InvalidInput {
                Ok(false)
            } else {
                Err(e.into())
            }
        }
    }
}

// ── Streaming APIs (primary) ────────────────────────────────────────────

/// Stream VRS results via callback. Lazily decodes chromosomes, clears cache when done.
/// Returns the number of results processed.
pub fn compute_vrs_ids_streaming(
    store: &mut RefgetStore,
    collection_digest: &str,
    vcf_path: &str,
    mut on_result: impl FnMut(VrsResult),
) -> Result<usize> {
    let name_to_digest = build_name_to_digest(store, collection_digest)?;
    let mut chrom_accessions: HashMap<String, String> = HashMap::new();
    let mut digest_writer = DigestWriter::new();
    let mut line_buf = String::new();
    let mut count = 0;

    let mut reader = open_vcf(vcf_path)?;

    while read_vcf_line(&mut *reader, &mut line_buf)? {
        let line = line_buf.trim_end_matches('\n').trim_end_matches('\r');
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.splitn(10, '\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let chrom = fields[0];
        let pos: u64 = fields[1]
            .parse::<u64>()
            .context("Invalid POS field")?
            .saturating_sub(1);
        let ref_allele = fields[3];
        let alt_field = fields[4];

        if !chrom_accessions.contains_key(chrom) {
            if let Some(digest) = name_to_digest.get(chrom) {
                store.ensure_decoded(digest.as_str())?;
                chrom_accessions.insert(chrom.to_string(), format!("SQ.{}", digest));
            }
        }

        let seq_accession = match chrom_accessions.get(chrom) {
            Some(acc) => acc,
            None => continue,
        };

        let raw_digest = &name_to_digest[chrom];
        let sequence = store
            .sequence_bytes(raw_digest.as_str())
            .context(format!("Chromosome {} not decoded in store", chrom))?;

        for alt in alt_field.split(',') {
            if alt.starts_with('<') || alt == "*" || alt == "." {
                continue;
            }

            let norm = normalize(sequence, pos, ref_allele.as_bytes(), alt.as_bytes())
                .context(format!("Failed to normalize variant at {}:{}", chrom, pos + 1))?;
            let norm_seq = std::str::from_utf8(&norm.allele)
                .context(format!("Normalized allele is not valid UTF-8 at {}:{}", chrom, pos + 1))?;

            let vrs_id =
                digest_writer.allele_identifier_literal(seq_accession, norm.start, norm.end, norm_seq);

            on_result(VrsResult {
                chrom: chrom.to_string(),
                pos,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt.to_string(),
                vrs_id,
            });
            count += 1;
        }
    }

    store.clear_decoded_cache();
    Ok(count)
}

/// Stream VRS results via callback using a pre-loaded read-only store.
/// All referenced sequences must already be decoded.
/// Returns the number of results processed.
pub fn compute_vrs_ids_streaming_readonly(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
    mut on_result: impl FnMut(VrsResult),
) -> Result<usize> {
    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(name, digest)| (name.clone(), format!("SQ.{}", digest)))
        .collect();

    let mut digest_writer = DigestWriter::new();
    let mut line_buf = String::new();
    let mut count = 0;

    let mut reader = open_vcf(vcf_path)?;

    while read_vcf_line(&mut *reader, &mut line_buf)? {
        let line = line_buf.trim_end_matches('\n').trim_end_matches('\r');
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.splitn(10, '\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let chrom = fields[0];
        let pos: u64 = fields[1]
            .parse::<u64>()
            .context("Invalid POS field")?
            .saturating_sub(1);
        let ref_allele = fields[3];
        let alt_field = fields[4];

        let seq_accession = match chrom_accessions.get(chrom) {
            Some(acc) => acc,
            None => continue,
        };

        let raw_digest = match name_to_digest.get(chrom) {
            Some(d) => d,
            None => continue,
        };
        let sequence = store
            .sequence_bytes(raw_digest.as_str())
            .context(format!("Chromosome {} not decoded in store", chrom))?;

        for alt in alt_field.split(',') {
            if alt.starts_with('<') || alt == "*" || alt == "." {
                continue;
            }

            let norm = normalize(sequence, pos, ref_allele.as_bytes(), alt.as_bytes())
                .context(format!("Failed to normalize variant at {}:{}", chrom, pos + 1))?;
            let norm_seq = std::str::from_utf8(&norm.allele)
                .context(format!("Normalized allele is not valid UTF-8 at {}:{}", chrom, pos + 1))?;

            let vrs_id =
                digest_writer.allele_identifier_literal(seq_accession, norm.start, norm.end, norm_seq);

            on_result(VrsResult {
                chrom: chrom.to_string(),
                pos,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt.to_string(),
                vrs_id,
            });
            count += 1;
        }
    }

    Ok(count)
}

// ── Vec-collecting APIs (convenience wrappers) ──────────────────────────

/// Convenience wrapper: collects all VRS results into a Vec.
pub fn compute_vrs_ids_from_vcf(
    store: &mut RefgetStore,
    collection_digest: &str,
    vcf_path: &str,
) -> Result<Vec<VrsResult>> {
    let mut results = Vec::new();
    compute_vrs_ids_streaming(store, collection_digest, vcf_path, |r| results.push(r))?;
    Ok(results)
}

/// Collect all VRS results from a pre-loaded read-only store.
pub fn compute_vrs_ids_from_vcf_readonly(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    vcf_path: &str,
) -> Result<Vec<VrsResult>> {
    let mut results = Vec::new();
    compute_vrs_ids_streaming_readonly(store, name_to_digest, vcf_path, |r| results.push(r))?;
    Ok(results)
}
