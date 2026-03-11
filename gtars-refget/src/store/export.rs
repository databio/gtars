//! FASTA export and BED region extraction for RefgetStore.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::collections::HashMap;
use std::ffi::OsStr;
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
use crate::hashkeyable::HashKeyable;


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

        if parsed_start < 0 || parsed_end < 0 {
            let err_str = format!(
                "Error reading line {} due to invalid start or end coordinates: '{}'",
                self.line_num + 1,
                line_string
            );
            return Some(Err(anyhow!(err_str)));
        }

        if self.previous_parsed_chr != parsed_chr {
            self.previous_parsed_chr = parsed_chr.clone();

            let result = match self
                .store
                .get_sequence_by_name(&self.collection_digest, &parsed_chr)
            {
                Ok(seq_record) => seq_record,
                Err(e) => {
                    let err_str = format!(
                        "Line {}: sequence '{}' not found in collection '{}': {}",
                        self.line_num + 1,
                        parsed_chr,
                        String::from_utf8_lossy(self.collection_digest.as_ref()),
                        e
                    );
                    return Some(Err(anyhow!(err_str)));
                }
            };

            self.current_seq_digest = result.metadata().sha512t24u.clone();
        }

        let retrieved_substring = match self.store.get_substring(
            &self.current_seq_digest,
            parsed_start as usize,
            parsed_end as usize,
        ) {
            Ok(substring) => substring,
            Err(e) => {
                let err_str = format!(
                    "Line {}: failed to get substring for digest '{}' from {} to {}: {}",
                    self.line_num + 1,
                    self.current_seq_digest,
                    parsed_start,
                    parsed_end,
                    e
                );
                return Some(Err(anyhow!(err_str)));
            }
        };

        Some(Ok(RetrievedSequence {
            sequence: retrieved_substring,
            chrom_name: parsed_chr,
            start: parsed_start as u32,
            end: parsed_end as u32,
        }))
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

        let collection_key = collection_digest.as_ref().to_key();

        let name_to_metadata: HashMap<String, SequenceMetadata> = self
            .name_lookup
            .get(&collection_key)
            .map(|name_map| {
                name_map
                    .iter()
                    .filter_map(|(name, seq_digest)| {
                        self.sequence_store.get(seq_digest).map(|record| {
                            (name.clone(), record.metadata().clone())
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();

        let seq_iter = self.substrings_from_regions(&collection_digest, bed_file_path)?;

        let mut previous_parsed_chr = String::new();
        let mut current_header: String = String::new();
        let mut previous_header: String = String::new();

        for rs in seq_iter.into_iter() {
            let rs = rs?;

            if previous_parsed_chr != rs.chrom_name {
                previous_parsed_chr = rs.chrom_name.clone();

                if let Some(meta) = name_to_metadata.get(&rs.chrom_name) {
                    current_header =
                        format!(">{} {} {} {} {}", meta.name, meta.length, meta.alphabet, meta.sha512t24u, meta.md5);
                }
            }

            let retrieved_substring = rs.sequence;

            if previous_header != current_header {
                let prefix = if previous_header.is_empty() { "" } else { "\n" };

                previous_header = current_header.clone();

                let header_to_be_written = format!("{}{}\n", prefix, current_header);
                writer.write_all(header_to_be_written.as_bytes())?;
            }

            writer.write_all(retrieved_substring.as_ref())?;
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

        let name_to_digest: HashMap<String, [u8; 32]> = self
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

            let (metadata, sequence_data) = match record {
                SequenceRecord::Stub(_) => {
                    return Err(anyhow!("Sequence data not loaded for '{}'. Call load_sequence() or load_all_sequences() first.", seq_name));
                }
                SequenceRecord::Full { metadata, sequence } => (metadata, sequence),
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

            let (metadata, sequence_data) = match record {
                SequenceRecord::Stub(_) => {
                    return Err(anyhow!(
                        "Sequence data not loaded for digest: {}",
                        digest_str
                    ));
                }
                SequenceRecord::Full { metadata, sequence } => (metadata, sequence),
            };

            write_fasta_record(&mut *writer, metadata, sequence_data, self.mode, line_width)?;
        }

        writer.flush()?;

        Ok(())
    }
}
