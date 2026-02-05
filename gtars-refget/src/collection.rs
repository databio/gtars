//! Extended SequenceCollection with filesystem operations.
//!
//! This module adds file I/O and caching methods to the core types
//! defined in `crate::digest::types`.

use anyhow::Result;
use std::io::Write;
use std::path::Path;

use crate::fasta::digest_fasta;
use crate::utils::PathExtension;

// Re-export core types from digest::types for backward compatibility
pub use crate::digest::types::{
    FaiMetadata, SeqColDigestLvl1, SequenceCollection, SequenceCollectionMetadata,
    SequenceCollectionRecord, SequenceMetadata, SequenceRecord, digest_sequence,
    digest_sequence_with_description, parse_rgsi_line,
};

/// Shared implementation for writing RGSI (Refget Sequence Index) files.
fn write_rgsi_impl<P: AsRef<Path>>(
    file_path: P,
    metadata: &SequenceCollectionMetadata,
    sequences: Option<&[SequenceRecord]>,
) -> Result<()> {
    let file_path = file_path.as_ref();
    let mut file = std::fs::File::create(file_path)?;

    // Write collection digest headers
    writeln!(file, "##seqcol_digest={}", metadata.digest)?;
    writeln!(file, "##names_digest={}", metadata.names_digest)?;
    writeln!(file, "##sequences_digest={}", metadata.sequences_digest)?;
    writeln!(file, "##lengths_digest={}", metadata.lengths_digest)?;
    writeln!(
        file,
        "#name\tlength\talphabet\tsha512t24u\tmd5\tdescription"
    )?;

    // Write sequence metadata if available
    if let Some(seqs) = sequences {
        for seq_record in seqs {
            let seq_meta = seq_record.metadata();
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}",
                seq_meta.name,
                seq_meta.length,
                seq_meta.alphabet,
                seq_meta.sha512t24u,
                seq_meta.md5,
                seq_meta.description.as_deref().unwrap_or("")
            )?;
        }
    }
    Ok(())
}

// ============================================================================
// Extensions to SequenceMetadata for filesystem operations
// ============================================================================

/// Extension trait for SequenceMetadata with filesystem-dependent methods.
pub trait SequenceMetadataExt {
    fn disk_size(&self, mode: &crate::store::StorageMode) -> usize;
}

impl SequenceMetadataExt for SequenceMetadata {
    /// Calculate the disk size in bytes for this sequence.
    fn disk_size(&self, mode: &crate::store::StorageMode) -> usize {
        match mode {
            crate::store::StorageMode::Raw => self.length,
            crate::store::StorageMode::Encoded => {
                let bits_per_symbol = self.alphabet.bits_per_symbol();
                let total_bits = self.length * bits_per_symbol;
                total_bits.div_ceil(8)
            }
        }
    }
}

// ============================================================================
// Extensions to SequenceRecord for filesystem operations
// ============================================================================

/// Extension trait for SequenceRecord with filesystem-dependent methods.
pub trait SequenceRecordExt {
    fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<()>;
}

impl SequenceRecordExt for SequenceRecord {
    /// Write a single sequence to a file.
    fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        use std::fs::{self, File};
        use std::io::Write;

        let data = match self {
            SequenceRecord::Stub(_) => {
                return Err(anyhow::anyhow!(
                    "Cannot write file: sequence data not loaded"
                ));
            }
            SequenceRecord::Full { sequence, .. } => sequence,
        };

        // Create parent directories if they don't exist
        if let Some(parent) = path.as_ref().parent() {
            fs::create_dir_all(parent)?;
        }

        let mut file = File::create(path)?;
        file.write_all(data)?;
        Ok(())
    }
}

// ============================================================================
// Extensions to SequenceCollectionRecord for filesystem operations
// ============================================================================

/// Extension trait for SequenceCollectionRecord with filesystem-dependent methods.
pub trait SequenceCollectionRecordExt {
    fn write_collection_rgsi<P: AsRef<Path>>(&self, file_path: P) -> Result<()>;
}

impl SequenceCollectionRecordExt for SequenceCollectionRecord {
    /// Write the collection to an RGSI file.
    fn write_collection_rgsi<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        write_rgsi_impl(file_path, self.metadata(), self.sequences())
    }
}

// ============================================================================
// Extensions to SequenceCollection for filesystem operations
// ============================================================================

/// Extension trait for SequenceCollection with filesystem-dependent methods.
pub trait SequenceCollectionExt {
    fn from_fasta<P: AsRef<Path>>(file_path: P) -> Result<SequenceCollection>;
    fn from_rgsi<P: AsRef<Path>>(file_path: P) -> Result<SequenceCollection>;
    fn from_path_no_cache<P: AsRef<Path>>(file_path: P) -> Result<SequenceCollection>;
    fn from_path_with_cache<P: AsRef<Path>>(
        file_path: P,
        read_cache: bool,
        write_cache: bool,
    ) -> Result<SequenceCollection>;
    fn write_collection_rgsi<P: AsRef<Path>>(&self, file_path: P) -> Result<()>;
    fn write_rgsi(&self) -> Result<()>;
    fn to_record(&self) -> SequenceCollectionRecord;
    fn write_fasta<P: AsRef<Path>>(&self, file_path: P, line_width: Option<usize>) -> Result<()>;
}

impl SequenceCollectionExt for SequenceCollection {
    /// Default behavior: read and write cache
    fn from_fasta<P: AsRef<Path>>(file_path: P) -> Result<SequenceCollection> {
        Self::from_path_with_cache(file_path, true, true)
    }

    fn from_rgsi<P: AsRef<Path>>(file_path: P) -> Result<SequenceCollection> {
        let rgsi_file_path = file_path.as_ref().replace_exts_with("rgsi");

        if rgsi_file_path.exists() {
            read_rgsi_file(&rgsi_file_path)
        } else {
            Err(anyhow::anyhow!(
                "RGSI file does not exist at {:?}",
                rgsi_file_path
            ))
        }
    }

    /// No caching at all
    fn from_path_no_cache<P: AsRef<Path>>(file_path: P) -> Result<SequenceCollection> {
        Self::from_path_with_cache(file_path, false, false)
    }

    fn from_path_with_cache<P: AsRef<Path>>(
        file_path: P,
        read_cache: bool,
        write_cache: bool,
    ) -> Result<SequenceCollection> {
        let fa_file_path = file_path.as_ref();
        let rgsi_file_path = fa_file_path.replace_exts_with("rgsi");

        // Check if the cache file exists and is valid
        if read_cache && rgsi_file_path.exists() {
            let seqcol = read_rgsi_file(&rgsi_file_path)?;
            if !seqcol.sequences.is_empty() {
                return Ok(seqcol);
            }
            // Cache is empty/stale - fall through to re-digest
            let _ = std::fs::remove_file(&rgsi_file_path);
        }

        // Digest the fasta file
        let seqcol: SequenceCollection = digest_fasta(file_path.as_ref())?;

        // Write the SequenceCollection to the RGSI file
        if write_cache && !rgsi_file_path.exists() {
            seqcol.write_rgsi()?;
        }
        Ok(seqcol)
    }

    /// Write the SequenceCollection to an RGSI file.
    fn write_collection_rgsi<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        write_rgsi_impl(file_path, &self.metadata, Some(&self.sequences))
    }

    /// Write the SequenceCollection to an RGSI file using the stored file path.
    fn write_rgsi(&self) -> Result<()> {
        if let Some(file_path) = &self.metadata.file_path {
            let rgsi_file_path = file_path.replace_exts_with("rgsi");
            self.write_collection_rgsi(rgsi_file_path)
        } else {
            Err(anyhow::anyhow!(
                "No file path specified for RGSI output. Use `write_collection_rgsi` to specify a file path."
            ))
        }
    }

    /// Convert to a SequenceCollectionRecord for storage.
    fn to_record(&self) -> SequenceCollectionRecord {
        SequenceCollectionRecord::from(self.clone())
    }

    /// Write the SequenceCollection to a FASTA file.
    fn write_fasta<P: AsRef<Path>>(&self, file_path: P, line_width: Option<usize>) -> Result<()> {
        use std::fs::File;

        let line_width = line_width.unwrap_or(70);
        let mut output_file = File::create(file_path)?;

        for record in &self.sequences {
            if !record.is_loaded() {
                return Err(anyhow::anyhow!(
                    "Cannot write FASTA: sequence '{}' does not have data loaded",
                    record.metadata().name
                ));
            }

            let decoded_sequence = record.decode().ok_or_else(|| {
                anyhow::anyhow!("Failed to decode sequence '{}'", record.metadata().name)
            })?;

            // Write FASTA header
            let metadata = record.metadata();
            let header = match &metadata.description {
                Some(desc) => format!(">{} {}", metadata.name, desc),
                None => format!(">{}", metadata.name),
            };
            writeln!(output_file, "{}", header)?;

            // Write sequence with line wrapping
            for chunk in decoded_sequence.as_bytes().chunks(line_width) {
                output_file.write_all(chunk)?;
                output_file.write_all(b"\n")?;
            }
        }

        Ok(())
    }
}

// ============================================================================
// RGSI File I/O
// ============================================================================

/// Read an RGSI file and return a SequenceCollection.
pub fn read_rgsi_file<T: AsRef<Path>>(file_path: T) -> Result<SequenceCollection> {
    use std::io::BufRead;

    let file = std::fs::File::open(&file_path)?;
    let reader = std::io::BufReader::new(file);
    let mut results = Vec::new();

    // Variables to store header metadata
    let mut seqcol_digest = String::new();
    let mut names_digest = String::new();
    let mut sequences_digest = String::new();
    let mut lengths_digest = String::new();

    for line in reader.lines() {
        let line = line?;

        // Parse header metadata lines
        if line.starts_with("##") {
            if let Some(stripped) = line.strip_prefix("##") {
                if let Some((key, value)) = stripped.split_once('=') {
                    match key {
                        "seqcol_digest" => seqcol_digest = value.to_string(),
                        "names_digest" => names_digest = value.to_string(),
                        "sequences_digest" => sequences_digest = value.to_string(),
                        "lengths_digest" => lengths_digest = value.to_string(),
                        _ => {} // Ignore unknown metadata keys
                    }
                }
            }
            continue;
        }

        // Skip other comment lines
        if line.starts_with('#') {
            continue;
        }

        // Parse sequence data line
        if let Some(metadata) = parse_rgsi_line(&line) {
            results.push(SequenceRecord::Stub(metadata));
        }
    }

    // Build collection metadata
    let collection_metadata = if sequences_digest.is_empty()
        || names_digest.is_empty()
        || lengths_digest.is_empty()
    {
        // Compute from sequence records
        SequenceCollectionMetadata::from_sequences(&results, Some(file_path.as_ref().to_path_buf()))
    } else {
        // Use digests from file header
        let lvl1 = SeqColDigestLvl1 {
            sequences_digest,
            names_digest,
            lengths_digest,
        };

        let digest = if seqcol_digest.is_empty() {
            lvl1.to_digest()
        } else {
            seqcol_digest
        };

        SequenceCollectionMetadata {
            digest,
            n_sequences: results.len(),
            names_digest: lvl1.names_digest,
            sequences_digest: lvl1.sequences_digest,
            lengths_digest: lvl1.lengths_digest,
            file_path: Some(file_path.as_ref().to_path_buf()),
        }
    };

    Ok(SequenceCollection {
        metadata: collection_metadata,
        sequences: results,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::encode_sequence;
    use crate::digest::{AlphabetType, lookup_alphabet};
    use crate::fasta::{digest_fasta, load_fasta};

    #[test]
    fn test_decode_returns_none_when_no_data() {
        let seqcol =
            digest_fasta("../tests/data/fasta/base.fa").expect("Failed to load test FASTA file");

        for seq_record in &seqcol.sequences {
            assert!(!seq_record.is_loaded());
            assert_eq!(seq_record.decode(), None);
        }
    }

    #[test]
    fn test_decode_with_loaded_data() {
        let seqcol =
            load_fasta("../tests/data/fasta/base.fa").expect("Failed to load test FASTA file");

        let expected_sequences = vec![("chrX", "TTGGGGAA"), ("chr1", "GGAA"), ("chr2", "GCGC")];

        for (seq_record, (expected_name, expected_seq)) in
            seqcol.sequences.iter().zip(expected_sequences.iter())
        {
            assert_eq!(seq_record.metadata().name, *expected_name);
            assert!(seq_record.is_loaded());
            let decoded = seq_record.decode().expect("decode() should return Some");
            assert_eq!(decoded, *expected_seq);
        }
    }

    #[test]
    fn test_decode_handles_encoded_data() {
        let sequence = b"ACGT";
        let alphabet = lookup_alphabet(&AlphabetType::Dna2bit);
        let encoded_data = encode_sequence(sequence, alphabet);

        let record = SequenceRecord::Full {
            metadata: SequenceMetadata {
                name: "test_seq".to_string(),
                description: None,
                length: 4,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Dna2bit,
                fai: None,
            },
            sequence: encoded_data,
        };

        let decoded = record.decode().expect("Should decode encoded data");
        assert_eq!(decoded, "ACGT");
    }

    #[test]
    fn test_sequence_collection_from_fasta() {
        let seqcol = SequenceCollection::from_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load FASTA");

        assert_eq!(seqcol.sequences.len(), 3);
        assert!(!seqcol.metadata.digest.is_empty());
    }

    #[test]
    fn test_sequence_collection_write_and_read_rgsi() {
        use tempfile::tempdir;

        let seqcol =
            digest_fasta("../tests/data/fasta/base.fa").expect("Failed to load test FASTA file");

        let dir = tempdir().expect("Failed to create temp dir");
        let rgsi_path = dir.path().join("test.rgsi");

        seqcol
            .write_collection_rgsi(&rgsi_path)
            .expect("Failed to write RGSI");

        let loaded = read_rgsi_file(&rgsi_path).expect("Failed to read RGSI");

        assert_eq!(loaded.metadata.digest, seqcol.metadata.digest);
        assert_eq!(loaded.sequences.len(), seqcol.sequences.len());
    }

    #[test]
    fn test_sequence_collection_write_fasta() {
        use std::fs;
        use tempfile::tempdir;

        let seqcol =
            load_fasta("../tests/data/fasta/base.fa").expect("Failed to load test FASTA file");

        let dir = tempdir().expect("Failed to create temp dir");
        let fasta_path = dir.path().join("output.fa");

        seqcol
            .write_fasta(&fasta_path, None)
            .expect("Failed to write FASTA");

        let content = fs::read_to_string(&fasta_path).expect("Failed to read file");
        assert!(content.contains(">chrX"));
        assert!(content.contains("TTGGGGAA"));
    }

    #[test]
    fn test_digest_sequence_basic() {
        let seq = digest_sequence("test_seq", b"ACGTACGT");
        assert_eq!(seq.metadata().name, "test_seq");
        assert_eq!(seq.metadata().length, 8);
        assert!(!seq.metadata().sha512t24u.is_empty());
        assert!(seq.is_loaded());
        assert_eq!(seq.sequence().unwrap(), b"ACGTACGT");
    }

    #[test]
    fn test_digest_sequence_matches_fasta_digest() {
        let seqcol =
            load_fasta("../tests/data/fasta/base.fa").expect("Failed to load test FASTA file");
        let fasta_seq = &seqcol.sequences[0]; // chrX: TTGGGGAA

        let prog_seq = digest_sequence("chrX", b"TTGGGGAA");

        assert_eq!(
            prog_seq.metadata().sha512t24u,
            fasta_seq.metadata().sha512t24u
        );
        assert_eq!(prog_seq.metadata().md5, fasta_seq.metadata().md5);
    }

    #[test]
    fn test_sequence_metadata_disk_size() {
        use crate::store::StorageMode;

        let metadata = SequenceMetadata {
            name: "test".to_string(),
            description: None,
            length: 1000,
            sha512t24u: "test".to_string(),
            md5: "test".to_string(),
            alphabet: AlphabetType::Dna2bit,
            fai: None,
        };

        assert_eq!(metadata.disk_size(&StorageMode::Raw), 1000);
        assert_eq!(metadata.disk_size(&StorageMode::Encoded), 250);
    }

    #[test]
    fn test_sequence_record_to_file() {
        use std::fs;
        use tempfile::tempdir;

        let seqcol =
            load_fasta("../tests/data/fasta/base.fa").expect("Failed to load test FASTA file");

        let dir = tempdir().expect("Failed to create temp dir");
        let file_path = dir.path().join("test_seq.txt");

        seqcol.sequences[0]
            .to_file(&file_path)
            .expect("Failed to write file");

        let content = fs::read(&file_path).expect("Failed to read file");
        assert!(!content.is_empty());
    }
}
