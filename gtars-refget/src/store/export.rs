//! FASTA export and BED region extraction for RefgetStore.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::ffi::OsStr;

use indexmap::IndexMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use flate2::read::MultiGzDecoder;
use gtars_core::utils::get_file_info;

use crate::digest::{
    SequenceMetadata, SequenceRecord,
    decode_substring_from_bytes, lookup_alphabet,
};
use crate::hashkeyable::{DigestKey, HashKeyable};


// ============================================================================
// Free function
// ============================================================================

/// Helper function to decode a sequence record and write it as a FASTA entry.
///
/// Handles decoding (encoded or raw storage modes), header formatting (with optional
/// description), and line-wrapped sequence output.
pub(crate) fn write_fasta_record(
    writer: &mut dyn Write,
    metadata: &SequenceMetadata,
    sequence_data: &[u8],
    mode: StorageMode,
    line_width: usize,
) -> Result<()> {
    let decoded_sequence = match mode {
        StorageMode::Encoded => {
            let alphabet = lookup_alphabet(&metadata.alphabet);
            let decoded =
                decode_substring_from_bytes(sequence_data, 0, metadata.length, alphabet);
            String::from_utf8(decoded).context("Failed to decode sequence as UTF-8")?
        }
        StorageMode::Raw => String::from_utf8(sequence_data.to_vec())
            .context("Failed to decode raw sequence as UTF-8")?,
    };

    let header = match &metadata.description {
        Some(desc) => format!(">{} {}", metadata.name, desc),
        None => format!(">{}", metadata.name),
    };
    writeln!(writer, "{}", header)?;

    for chunk in decoded_sequence.as_bytes().chunks(line_width) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }

    Ok(())
}

// ============================================================================
// Shared per-region substring retrieval
// ============================================================================

/// Retrieve the substring for a single region, shared by the file-based and
/// vectors-based iterators.
///
/// Mutates the caller-owned `previous_parsed_chr` / `current_seq_digest` cache
/// so the consecutive-same-chrom optimization works across either path.
/// `region_label` (e.g. `"Line 5"` or `"Region 5"`) is interpolated into error
/// messages.
#[allow(clippy::too_many_arguments)]
fn retrieve_substring_for_region<K: AsRef<[u8]>>(
    store: &ReadonlyRefgetStore,
    collection_digest: &K,
    parsed_chr: String,
    parsed_start: i64,
    parsed_end: i64,
    previous_parsed_chr: &mut String,
    current_seq_digest: &mut String,
    region_label: &str,
) -> Result<RetrievedSequence> {
    if parsed_start < 0 || parsed_end < 0 {
        return Err(anyhow!(
            "{} has invalid start or end coordinates: start={}, end={}",
            region_label,
            parsed_start,
            parsed_end
        ));
    }

    if *previous_parsed_chr != parsed_chr {
        *previous_parsed_chr = parsed_chr.clone();

        let result = store
            .get_sequence_by_name(collection_digest, &parsed_chr)
            .map_err(|e| {
                anyhow!(
                    "{}: sequence '{}' not found in collection '{}': {}",
                    region_label,
                    parsed_chr,
                    String::from_utf8_lossy(collection_digest.as_ref()),
                    e
                )
            })?;

        *current_seq_digest = result.metadata().sha512t24u.clone();
    }

    let retrieved_substring = store
        .get_substring(&*current_seq_digest, parsed_start as usize, parsed_end as usize)
        .map_err(|e| {
            anyhow!(
                "{}: failed to get substring for digest '{}' from {} to {}: {}",
                region_label,
                current_seq_digest,
                parsed_start,
                parsed_end,
                e
            )
        })?;

    Ok(RetrievedSequence {
        sequence: retrieved_substring,
        chrom_name: parsed_chr,
        start: parsed_start as u32,
        end: parsed_end as u32,
    })
}

// ============================================================================
// SubstringsFromRegions Iterator impl
// ============================================================================

impl<K> Iterator for SubstringsFromRegions<'_, K>
where
    K: AsRef<[u8]>,
{
    type Item = Result<RetrievedSequence, anyhow::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        use gtars_core::utils::parse_bedlike_file;

        let mut line_string = String::new();

        let num_bytes = self.reader.read_line(&mut line_string);
        match num_bytes {
            Ok(bytes) => {
                if bytes == 0 {
                    return None;
                }
            }
            Err(err) => return Some(Err(err.into())),
        };

        self.line_num += 1;

        let (parsed_chr, parsed_start, parsed_end) = match parse_bedlike_file(line_string.trim()) {
            Some(coords) => coords,
            None => {
                let err_str = format!(
                    "Error reading line {} because it could not be parsed as a BED-like entry: '{}'",
                    self.line_num + 1,
                    line_string
                );
                return Some(Err(anyhow!(err_str)));
            }
        };

        Some(retrieve_substring_for_region(
            self.store,
            &self.collection_digest,
            parsed_chr,
            parsed_start as i64,
            parsed_end as i64,
            &mut self.previous_parsed_chr,
            &mut self.current_seq_digest,
            &format!("Line {}", self.line_num + 1),
        ))
    }
}

// ============================================================================
// SubstringsFromRegionVectors Iterator impl
// ============================================================================

impl<K> Iterator for SubstringsFromRegionVectors<'_, K>
where
    K: AsRef<[u8]>,
{
    type Item = Result<RetrievedSequence, anyhow::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.chroms.len() {
            return None;
        }
        let i = self.index;
        self.index += 1;

        let parsed_chr = self.chroms[i].clone();
        let parsed_start = self.starts[i] as i64;
        let parsed_end = self.ends[i] as i64;

        Some(retrieve_substring_for_region(
            self.store,
            &self.collection_digest,
            parsed_chr,
            parsed_start,
            parsed_end,
            &mut self.previous_parsed_chr,
            &mut self.current_seq_digest,
            &format!("Region {}", i + 1),
        ))
    }
}

// ============================================================================
// ReadonlyRefgetStore export methods
// ============================================================================

impl ReadonlyRefgetStore {
    /// Get an iterator over substrings defined by BED file regions.
    pub fn substrings_from_regions<'a, K: AsRef<[u8]>>(
        &'a self,
        collection_digest: K,
        bed_file_path: &str,
    ) -> Result<SubstringsFromRegions<'a, K>> {
        let path = Path::new(bed_file_path);
        let file_info = get_file_info(path);
        let is_gzipped = file_info.is_gzipped;

        let opened_bed_file = File::open(path)?;

        let reader: Box<dyn std::io::Read> = match is_gzipped {
            true => Box::new(MultiGzDecoder::new(std::io::BufReader::new(opened_bed_file))),
            false => Box::new(opened_bed_file),
        };
        let reader = std::io::BufReader::new(reader);

        Ok(SubstringsFromRegions {
            store: self,
            reader,
            collection_digest,
            previous_parsed_chr: String::new(),
            current_seq_digest: String::new(),
            line_num: 0,
        })
    }

    /// Get an iterator over substrings defined by in-memory region vectors.
    ///
    /// `chroms`, `starts`, and `ends` are parallel slices; element `i` defines
    /// region `chroms[i]:starts[i]-ends[i]`. Yields the same `RetrievedSequence`
    /// results as `substrings_from_regions` without any temp-file round-trip.
    pub fn substrings_from_region_vectors<'a, K, S>(
        &'a self,
        collection_digest: K,
        chroms: &[S],
        starts: &[u32],
        ends: &[u32],
    ) -> Result<SubstringsFromRegionVectors<'a, K>>
    where
        K: AsRef<[u8]>,
        S: AsRef<str>,
    {
        if chroms.len() != starts.len() || chroms.len() != ends.len() {
            return Err(anyhow!(
                "Mismatched region vector lengths: chroms={}, starts={}, ends={}",
                chroms.len(),
                starts.len(),
                ends.len()
            ));
        }

        Ok(SubstringsFromRegionVectors {
            store: self,
            collection_digest,
            chroms: chroms.iter().map(|c| c.as_ref().to_string()).collect(),
            starts: starts.to_vec(),
            ends: ends.to_vec(),
            index: 0,
            previous_parsed_chr: String::new(),
            current_seq_digest: String::new(),
        })
    }

    /// Export sequences from BED file regions to a FASTA file.
    pub fn export_fasta_from_regions<K: AsRef<[u8]>>(
        &self,
        collection_digest: K,
        bed_file_path: &str,
        output_file_path: &str,
    ) -> Result<()> {
        let output_path_obj = Path::new(output_file_path);
        if let Some(parent) = output_path_obj.parent() {
            std::fs::create_dir_all(parent)?;
        }

        let file = File::create(output_file_path)?;

        let mut writer: Box<dyn Write> = if output_path_obj.extension() == Some(OsStr::new("gz")) {
            Box::new(GzEncoder::new(file, Compression::default()))
        } else {
            Box::new(file)
        };

        let seq_iter = self.substrings_from_regions(&collection_digest, bed_file_path)?;

        for rs in seq_iter.into_iter() {
            let rs = rs?;

            // Write one header per region with coordinates
            let header = format!(">{}:{}-{}\n", rs.chrom_name, rs.start, rs.end);
            writer.write_all(header.as_bytes())?;
            writer.write_all(rs.sequence.as_bytes())?;
            writer.write_all(b"\n")?;
        }

        writer.flush()?;

        Ok(())
    }

    /// Export sequences from a collection to a FASTA file.
    pub fn export_fasta<K: AsRef<[u8]>, P: AsRef<Path>>(
        &self,
        collection_digest: K,
        output_path: P,
        sequence_names: Option<Vec<&str>>,
        line_width: Option<usize>,
    ) -> Result<()> {
        let line_width = line_width.unwrap_or(80);
        let output_path = output_path.as_ref();
        let collection_key = collection_digest.as_ref().to_key();

        let name_to_digest: IndexMap<String, DigestKey> = self
            .name_lookup
            .get(&collection_key)
            .ok_or_else(|| {
                anyhow!(
                    "Collection not found: {:?}",
                    String::from_utf8_lossy(collection_digest.as_ref())
                )
            })?
            .clone();

        let names_to_export: Vec<String> = if let Some(names) = sequence_names {
            names.iter().map(|s| s.to_string()).collect()
        } else {
            name_to_digest.keys().cloned().collect()
        };

        let file = File::create(output_path).context(format!(
            "Failed to create output file: {}",
            output_path.display()
        ))?;

        let mut writer: Box<dyn Write> = if output_path.extension() == Some(OsStr::new("gz")) {
            Box::new(GzEncoder::new(BufWriter::new(file), Compression::default()))
        } else {
            Box::new(BufWriter::new(file))
        };

        for seq_name in names_to_export {
            let seq_digest = name_to_digest
                .get(&seq_name)
                .ok_or_else(|| anyhow!("Sequence '{}' not found in collection", seq_name))?;

            let record = self
                .sequence_store
                .get(seq_digest)
                .ok_or_else(|| anyhow!("Sequence record not found for digest: {:?}", seq_digest))?;

            let (metadata, sequence_data): (&SequenceMetadata, &[u8]) = match record {
                SequenceRecord::Stub(_) => {
                    return Err(anyhow!("Sequence data not loaded for '{}'. Call load_sequence() or load_all_sequences() first.", seq_name));
                }
                SequenceRecord::Full { metadata, sequence } => (metadata, sequence.as_slice()),
            };

            write_fasta_record(&mut *writer, metadata, sequence_data, self.mode, line_width)?;
        }

        writer.flush()?;

        Ok(())
    }

    /// Export sequences by their sequence digests to a FASTA file.
    pub fn export_fasta_by_digests<P: AsRef<Path>>(
        &self,
        seq_digests: Vec<&str>,
        output_path: P,
        line_width: Option<usize>,
    ) -> Result<()> {
        let line_width = line_width.unwrap_or(80);
        let output_path = output_path.as_ref();

        let file = File::create(output_path).context(format!(
            "Failed to create output file: {}",
            output_path.display()
        ))?;

        let mut writer: Box<dyn Write> = if output_path.extension() == Some(OsStr::new("gz")) {
            Box::new(GzEncoder::new(BufWriter::new(file), Compression::default()))
        } else {
            Box::new(BufWriter::new(file))
        };

        for digest_str in seq_digests {
            let digest_key = digest_str.as_bytes().to_key();

            let record = self
                .sequence_store
                .get(&digest_key)
                .ok_or_else(|| anyhow!("Sequence record not found for digest: {}", digest_str))?;

            let (metadata, sequence_data): (&SequenceMetadata, &[u8]) = match record {
                SequenceRecord::Stub(_) => {
                    return Err(anyhow!(
                        "Sequence data not loaded for digest: {}",
                        digest_str
                    ));
                }
                SequenceRecord::Full { metadata, sequence } => (metadata, sequence.as_slice()),
            };

            write_fasta_record(&mut *writer, metadata, sequence_data, self.mode, line_width)?;
        }

        writer.flush()?;

        Ok(())
    }
}
