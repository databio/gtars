//! WASM-safe VCF→VRS primitives.
//!
//! This module holds the parts of the VCF pipeline that have **no** filesystem,
//! threading, or BGZF dependencies, so they compile for both the native
//! (`filesystem`) build and the `wasm` build. The heavy machinery — file
//! opening, BGZF detection, and the thread/`crossbeam` parallel pipeline —
//! lives in [`crate::vcf`] behind `#[cfg(feature = "filesystem")]`.
//!
//! The key entry point is [`compute_vrs_ids_streaming_readonly_from_reader`]:
//! the per-record streaming loop over a [`ReadonlyRefgetStore`], generic over
//! any [`std::io::BufRead`] source. Native code feeds it `open_vcf(path)?`;
//! wasm feeds it a `std::io::Cursor` over the dropped VCF bytes.

use std::collections::HashMap;
use std::io::BufRead;

use anyhow::{Context, Result};
use gtars_refget::store::ReadonlyRefgetStore;

use crate::digest::DigestWriter;
use crate::normalize::{RefView, normalize_ref, ref_view_for as ref_view_for_inner};

/// Build a [`RefView`] over a resident sequence's bytes (thin `anyhow` wrapper
/// around [`crate::normalize::ref_view_for`]).
pub(crate) fn ref_view_for<'a>(
    store: &'a ReadonlyRefgetStore,
    raw_digest: &str,
) -> Result<RefView<'a>> {
    ref_view_for_inner(store, raw_digest).map_err(|e| anyhow::anyhow!(e))
}

/// True for an ALT allele we actually process: not symbolic (`<DEL>` etc.),
/// not the `*` overlapping-deletion marker, not the `.` missing marker, and
/// not empty.
pub fn is_real_alt(alt: &str) -> bool {
    !(alt.is_empty() || alt.starts_with('<') || alt == "*" || alt == ".")
}

/// One parsed VCF data record. Borrows slices out of the caller's line buffer.
/// `pos` is already 0-based interbase (1-based POS minus one; POS `0` is
/// rejected at parse time, never saturated).
/// `alts` is the raw ALT field; callers split it on `,` and filter with
/// [`is_real_alt`] (or use [`ParsedRecord::real_alts`]).
pub struct ParsedRecord<'a> {
    pub chrom: &'a str,
    pub pos: u64,
    pub ref_allele: &'a str,
    pub alts: &'a str,
}

impl<'a> ParsedRecord<'a> {
    /// Iterator over the real (non-symbolic, non-missing) ALT alleles.
    pub fn real_alts(&self) -> impl Iterator<Item = &'a str> {
        self.alts.split(',').filter(|a| is_real_alt(a))
    }
}

/// Parse one VCF line into CHROM/POS/REF/ALT, returning `None` for header/blank
/// lines, lines with fewer than 5 fields, or an unparseable / out-of-spec POS.
/// POS is converted to 0-based interbase. Borrows from `line`.
///
/// Per the VCF spec POS is 1-based and MUST be `>= 1`; a POS of `0` is rejected
/// (returns `None`) rather than silently saturating to interbase `0`, which
/// would corrupt the variant's coordinate.
pub fn parse_vcf_record(line: &str) -> Option<ParsedRecord<'_>> {
    let line = line.trim_end_matches(['\n', '\r']);
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    let mut it = line.splitn(6, '\t');
    let chrom = it.next()?;
    let pos_s = it.next()?;
    let _id = it.next()?;
    let ref_allele = it.next()?;
    let alts = it.next()?;
    let pos_1based = pos_s.parse::<u64>().ok()?;
    if pos_1based < 1 {
        return None;
    }
    let pos = pos_1based - 1;
    Some(ParsedRecord {
        chrom,
        pos,
        ref_allele,
        alts,
    })
}

/// Result of computing a VRS identifier for a single VCF variant.
#[derive(Debug, Clone)]
pub struct VrsResult {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    pub vrs_id: String,
}

/// Read the next non-empty line from a VCF reader. Returns false at EOF.
pub(crate) fn read_vcf_line(reader: &mut dyn BufRead, buf: &mut String) -> Result<bool> {
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

/// Stream VRS results via callback from any line source, using a pre-loaded
/// read-only store. All referenced sequences must already be decodable from the
/// store (resident or on-the-fly encoded). Returns the number of results.
///
/// This is the reader-generic core of the readonly VCF→VRS path: it carries no
/// filesystem dependency, so the native path-based wrapper
/// (`compute_vrs_ids_streaming_readonly` in [`crate::vcf`]) feeds it
/// `open_vcf(path)?`, and the wasm path feeds it a `std::io::Cursor` over the
/// dropped VCF bytes.
pub fn compute_vrs_ids_streaming_readonly_from_reader<R: BufRead>(
    store: &ReadonlyRefgetStore,
    name_to_digest: &HashMap<String, String>,
    mut reader: R,
    mut on_result: impl FnMut(VrsResult),
) -> Result<usize> {
    let chrom_accessions: HashMap<String, String> = name_to_digest
        .iter()
        .map(|(name, digest)| (name.clone(), format!("SQ.{}", digest)))
        .collect();

    let mut digest_writer = DigestWriter::new();
    let mut line_buf = String::new();
    let mut count = 0;

    while read_vcf_line(&mut reader, &mut line_buf)? {
        let Some(rec) = parse_vcf_record(&line_buf) else {
            continue;
        };
        let chrom = rec.chrom;
        let pos = rec.pos;
        let ref_allele = rec.ref_allele;

        let seq_accession = match chrom_accessions.get(chrom) {
            Some(acc) => acc,
            None => continue,
        };

        let raw_digest = match name_to_digest.get(chrom) {
            Some(d) => d,
            None => continue,
        };
        let view = ref_view_for(store, raw_digest.as_str())
            .context(format!("Chromosome {} not available in store", chrom))?;

        for alt in rec.real_alts() {
            let norm = normalize_ref(&view, pos, ref_allele.as_bytes(), alt.as_bytes())
                .context(format!("Failed to normalize variant at {}:{}", chrom, pos + 1))?;
            let norm_seq = std::str::from_utf8(&norm.allele).context(format!(
                "Normalized allele is not valid UTF-8 at {}:{}",
                chrom,
                pos + 1
            ))?;

            let vrs_id = digest_writer.allele_identifier_literal(
                seq_accession,
                norm.start,
                norm.end,
                norm_seq,
            );

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
