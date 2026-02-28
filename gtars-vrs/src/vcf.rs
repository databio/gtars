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
use gtars_refget::store::RefgetStore;

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

/// Compute VRS identifiers for all variants in a VCF file.
pub fn compute_vrs_ids_from_vcf(
    store: &mut RefgetStore,
    collection_digest: &str,
    vcf_path: &str,
) -> Result<Vec<VrsResult>> {
    // Get the collection to resolve chromosome names → sequence digests
    let collection = store
        .get_collection(collection_digest)
        .context("Failed to get sequence collection")?;

    // Build name → digest mapping from collection metadata (no decoding yet)
    let mut name_to_digest: HashMap<String, String> = HashMap::new();
    for seq_record in &collection.sequences {
        let meta = seq_record.metadata();
        name_to_digest.insert(meta.name.clone(), meta.sha512t24u.clone());
    }

    // Track which chromosomes have been decoded, mapping name → "SQ.xxx" accession
    let mut chrom_digests: HashMap<String, String> = HashMap::new();
    let mut results = Vec::new();
    let mut digest_writer = DigestWriter::new();
    let mut line_buf = String::new();

    let mut reader = open_vcf(vcf_path)?;

    loop {
        line_buf.clear();
        if reader.read_line(&mut line_buf)? == 0 {
            break;
        }

        let line = line_buf.trim_end_matches('\n').trim_end_matches('\r');
        if line.starts_with('#') {
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
            .saturating_sub(1); // VCF is 1-based → 0-based
        let ref_allele = fields[3];
        let alt_field = fields[4];

        // Lazily decode chromosome in the store on first encounter
        if !chrom_digests.contains_key(chrom) {
            if let Some(digest) = name_to_digest.get(chrom) {
                store.ensure_decoded(digest.as_str())?;
                chrom_digests.insert(chrom.to_string(), format!("SQ.{}", digest));
            }
        }

        let seq_accession = match chrom_digests.get(chrom) {
            Some(acc) => acc,
            None => continue, // Skip chromosomes not in the collection
        };

        // Look up decoded sequence bytes from store
        let raw_digest = &name_to_digest[chrom];
        let sequence = store
            .sequence_bytes(raw_digest.as_str())
            .context(format!("Chromosome {} not decoded in store", chrom))?;

        // Each ALT allele gets its own VRS ID
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

            results.push(VrsResult {
                chrom: chrom.to_string(),
                pos,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt.to_string(),
                vrs_id,
            });
        }
    }

    Ok(results)
}
