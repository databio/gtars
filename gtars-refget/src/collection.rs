use crate::alphabet::AlphabetType;
use crate::digest::{canonicalize_json, sha512t24u};
use crate::fasta::{digest_fasta, read_fasta_refget_file};
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::fmt::Display;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::utils::PathExtension;

/// A single Sequence Collection, which may or may not hold data.
#[derive(Clone, Debug)]
pub struct SequenceCollection {
    /// Vector of SequenceRecords, which contain metadata (name, length, digests, alphabet)
    /// and optionally the actual sequence data.
    pub sequences: Vec<SequenceRecord>,

    /// The SHA512t24u collection digest identifier as a string.
    pub digest: String,

    /// Level 1 digest components. Contains separate digests for
    /// sequences, names, and lengths arrays.
    pub lvl1: SeqColDigestLvl1,

    /// Optional path to the source file from which this collection was loaded.
    /// Used for caching and reference purposes.
    pub file_path: Option<PathBuf>,
}

/// A struct representing the first level of digests for a refget sequence collection.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SeqColDigestLvl1 {
    pub sequences_digest: String,
    pub names_digest: String,
    pub lengths_digest: String,
}

impl SeqColDigestLvl1 {
    /// Compute collection digest from lvl1 digests
    pub fn to_digest(&self) -> String {
        // Create JSON object with the lvl1 digest strings
        let mut lvl1_object = serde_json::Map::new();
        lvl1_object.insert(
            "names".to_string(),
            serde_json::Value::String(self.names_digest.clone()),
        );
        lvl1_object.insert(
            "sequences".to_string(),
            serde_json::Value::String(self.sequences_digest.clone()),
        );

        let lvl1_json = serde_json::Value::Object(lvl1_object);

        // Canonicalize the JSON object and compute collection digest
        let lvl1_canonical = canonicalize_json(&lvl1_json);
        let digest = sha512t24u(lvl1_canonical.as_bytes());
        println!("lvl1 digest: {}", digest);
        digest
    }

    /// Compute lvl1 digests from a collection of SequenceMetadata
    pub fn from_metadata(metadata_vec: &[&SequenceMetadata]) -> Self {
        use serde_json::Value;

        // Extract arrays for each field
        let sequences: Vec<String> = metadata_vec
            .iter()
            .map(|md| format!("SQ.{}", md.sha512t24u))
            .collect();
        let names: Vec<&str> = metadata_vec.iter().map(|md| md.name.as_str()).collect();
        let lengths: Vec<usize> = metadata_vec.iter().map(|md| md.length).collect();

        // Convert to JSON Values and canonicalize
        let sequences_json = Value::Array(
            sequences
                .iter()
                .map(|s| Value::String(s.to_string()))
                .collect(),
        );
        let names_json = Value::Array(names.iter().map(|s| Value::String(s.to_string())).collect());
        let lengths_json = Value::Array(
            lengths
                .iter()
                .map(|l| Value::Number(serde_json::Number::from(*l)))
                .collect(),
        );

        // Canonicalize to JCS format
        let sequences_canonical = canonicalize_json(&sequences_json);
        let names_canonical = canonicalize_json(&names_json);
        let lengths_canonical = canonicalize_json(&lengths_json);

        // Hash the canonicalized arrays
        SeqColDigestLvl1 {
            sequences_digest: sha512t24u(sequences_canonical.as_bytes()),
            names_digest: sha512t24u(names_canonical.as_bytes()),
            lengths_digest: sha512t24u(lengths_canonical.as_bytes()),
        }
    }
}

/// A representation of a single sequence that includes metadata and optionally data.
/// Combines sequence metadata with optional raw/encoded data
#[derive(Clone, Debug)]
pub struct SequenceRecord {
    pub metadata: SequenceMetadata,
    pub fai: Option<FaiMetadata>,
    pub data: Option<Vec<u8>>,
}

use std::fs::{self, File};
/// Metadata for a single sequence, including its name, length, digests, and alphabet type.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SequenceMetadata {
    pub name: String,
    pub length: usize,
    pub sha512t24u: String,
    pub md5: String,
    pub alphabet: AlphabetType,
}

/// FASTA index (FAI) metadata for a sequence.
/// This data is only present when a sequence was loaded from a FASTA file.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FaiMetadata {
    pub offset: u64,       // byte offset to first base of sequence data
    pub line_bases: u32,   // number of bases per line
    pub line_bytes: u32,   // number of bytes per line (including newline chars)
}

impl SequenceRecord {
    /// Utility function to write a single sequence to a file
    pub fn to_file<P: AsRef<Path>>(&self, path: P) -> anyhow::Result<()> {
        // Create parent directories if they don't exist
        if let Some(parent) = path.as_ref().parent() {
            fs::create_dir_all(parent)?;
        }

        let mut file = File::create(path)?;
        if let Some(data) = &self.data {
            file.write_all(data)?;
        }
        Ok(())
    }

    /// Decodes the sequence data to a string.
    ///
    /// This method attempts to decode the sequence data stored in this record.
    /// It handles both raw (uncompressed UTF-8) and encoded (bit-packed) data.
    /// The decoding strategy depends on the alphabet type:
    /// - For ASCII alphabet: data is already in raw form, just convert to string
    /// - For other alphabets: attempt encoded decoding first, fall back to raw
    ///
    /// # Returns
    ///
    /// * `Some(String)` - The decoded sequence if data is loaded
    /// * `None` - If no data is loaded in this record
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use gtars_refget::collection::SequenceRecord;
    /// # let record: SequenceRecord = todo!();
    /// let sequence = record.decode();
    /// if let Some(seq) = sequence {
    ///     println!("Sequence: {}", seq);
    /// }
    /// ```
    pub fn decode(&self) -> Option<String> {
        use crate::alphabet::lookup_alphabet;
        use crate::encoder::decode_substring_from_bytes;

        let data = self.data.as_ref()?;

        // For ASCII alphabet (8 bits per symbol), the data is always stored raw
        if self.metadata.alphabet == crate::alphabet::AlphabetType::Ascii {
            return String::from_utf8(data.clone()).ok();
        }

        // Try to detect if data is raw or encoded
        // Heuristic: for encoded data, the size should be approximately length * bits_per_symbol / 8
        // For raw data, the size should be approximately equal to length
        let alphabet = lookup_alphabet(&self.metadata.alphabet);

        // If data size matches the expected length (not the encoded size), it's probably raw
        if data.len() == self.metadata.length {
            // Try to decode as UTF-8
            if let Ok(raw_string) = String::from_utf8(data.clone()) {
                // Data appears to be raw UTF-8
                return Some(raw_string);
            }
        }

        // Data is probably encoded (size matches expected encoded size), try to decode it
        let decoded_bytes = decode_substring_from_bytes(
            data,
            0,
            self.metadata.length,
            alphabet
        );

        // Convert to string
        String::from_utf8(decoded_bytes).ok()
    }
}

impl SequenceCollection {
    /// Default behavior: read and write cache
    pub fn from_fasta<P: AsRef<Path>>(file_path: P) -> Result<Self> {
        Self::from_path_with_cache(file_path, true, true)
    }

    pub fn from_farg<P: AsRef<Path>>(file_path: P) -> Result<Self> {
        let farg_file_path = file_path.as_ref().replace_exts_with("farg");
        println!("From_farg - Reading from file: {:?}", file_path.as_ref());
        println!("Farg file path: {:?}", farg_file_path);

        if farg_file_path.exists() {
            println!("Reading from existing farg file: {:?}", farg_file_path);
            read_fasta_refget_file(&farg_file_path)
        } else {
            Err(anyhow::anyhow!(
                "FARG file does not exist at {:?}",
                farg_file_path
            ))
        }
    }

    /// Create a SequenceCollection from a vector of SequenceRecords.
    pub fn from_records(records: Vec<SequenceRecord>) -> Self {
        // Compute lvl1 digests from the metadata
        let metadata_refs: Vec<&SequenceMetadata> = records.iter().map(|r| &r.metadata).collect();
        let lvl1 = SeqColDigestLvl1::from_metadata(&metadata_refs);

        // Compute collection digest from lvl1 digests
        let collection_digest = lvl1.to_digest();

        SequenceCollection {
            sequences: records,
            digest: collection_digest,
            lvl1,
            file_path: None,
        }
    }

    /// No caching at all
    pub fn from_path_no_cache<P: AsRef<Path>>(file_path: P) -> Result<Self> {
        Self::from_path_with_cache(file_path, false, false)
    }

    pub fn from_path_with_cache<P: AsRef<Path>>(
        file_path: P,
        read_cache: bool,
        write_cache: bool,
    ) -> Result<Self> {
        // If the farg file exists, just use that.
        let fa_file_path = file_path.as_ref();
        let farg_file_path = fa_file_path.replace_exts_with("farg");
        println!(
            "from path with cache: reading from file: {:?}",
            file_path.as_ref()
        );
        println!("Farg file path: {:?}", farg_file_path);
        // Check if the file already exists
        if read_cache && farg_file_path.exists() {
            println!("Reading from existing farg file: {:?}", farg_file_path);
            // Read the existing farg file
            let seqcol = read_fasta_refget_file(&farg_file_path)?;

            // seqcol is already a SequenceCollection, just return it
            return Ok(seqcol);
        }
        println!("Computing digests...: {:?}", farg_file_path);

        // If the farg file does not exist, compute the digests
        // Digest the fasta file (your function)
        let seqcol: SequenceCollection = digest_fasta(file_path.as_ref())?;

        // Write the SequenceCollection to the FARG file
        if write_cache && !farg_file_path.exists() {
            seqcol.write_farg()?;
            println!("Farg file written to {:?}", farg_file_path);
        } else {
            println!(
                "Farg file already exists, not writing: {:?}",
                farg_file_path
            );
        }
        Ok(seqcol)
    }

    /// Write the SequenceCollection to a collection FARG file.
    ///
    /// Creates a FARG file with collection-level digest headers followed by
    /// sequence metadata for all sequences in this collection.
    ///
    /// # Arguments
    /// * `file_path` - The path to the FARG file to be written
    ///
    /// # Returns
    /// Result indicating success or error
    ///
    /// # Format
    /// The file includes:
    /// - Collection digest headers (##seqcol_digest, ##names_digest, etc.)
    /// - Column header (#name, length, alphabet, sha512t24u, md5)
    /// - One line per sequence with metadata
    pub fn write_collection_farg<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        let file_path = file_path.as_ref();
        println!("Writing collection farg file: {:?}", file_path);
        let mut file = std::fs::File::create(file_path)?;

        // Write collection digest headers
        writeln!(file, "##seqcol_digest={}", self.digest)?;
        writeln!(file, "##names_digest={}", self.lvl1.names_digest)?;
        writeln!(file, "##sequences_digest={}", self.lvl1.sequences_digest)?;
        writeln!(file, "##lengths_digest={}", self.lvl1.lengths_digest)?;
        writeln!(file, "#name\tlength\talphabet\tsha512t24u\tmd5")?;

        // Write sequence metadata
        for result_sr in &self.sequences {
            let result = result_sr.metadata.clone();
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}",
                result.name, result.length, result.alphabet, result.sha512t24u, result.md5
            )?;
        }
        Ok(())
    }

    /// Write the SequenceCollection to a FARG file, using the file path stored in the struct.
    pub fn write_farg(&self) -> Result<()> {
        if let Some(ref file_path) = self.file_path {
            let farg_file_path = file_path.replace_exts_with("farg");
            self.write_collection_farg(farg_file_path)
        } else {
            Err(anyhow::anyhow!(
                "No file path specified for FARG output. Use `write_collection_farg` to specify a file path."
            ))
        }
    }
}

impl Display for SequenceCollection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SequenceCollection with {} sequences, digest: {}",
            self.sequences.len(),
            self.digest
        )?;
        write!(f, "\nFirst 3 sequences:")?;
        for seqrec in self.sequences.iter().take(3) {
            write!(f, "\n- {}", seqrec)?;
        }
        Ok(())
    }
}

impl Display for SequenceRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SequenceRecord: {} (length: {}, alphabet: {}, ga4gh: {:02x?}, md5: {:02x?})",
            &self.metadata.name,
            &self.metadata.length,
            &self.metadata.alphabet,
            &self.metadata.sha512t24u,
            &self.metadata.md5
        )?;
        Ok(())
    }
}

// Iterator implementations for SequenceCollection
// Allows: for seq in &collection { ... }
impl<'a> IntoIterator for &'a SequenceCollection {
    type Item = &'a SequenceRecord;
    type IntoIter = std::slice::Iter<'a, SequenceRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.sequences.iter()
    }
}

// Consuming iterator
// Allows: for seq in collection { ... } (consumes the collection)
impl IntoIterator for SequenceCollection {
    type Item = SequenceRecord;
    type IntoIter = std::vec::IntoIter<SequenceRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.sequences.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::{digest_fasta, load_fasta};
    use crate::encoder::encode_sequence;

    #[test]
    fn test_decode_returns_none_when_no_data() {
        // Test that decode() returns None when data is None
        let seqcol = digest_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        for seq_record in &seqcol.sequences {
            assert!(seq_record.data.is_none(), "digest_fasta should not load sequence data");
            assert_eq!(seq_record.decode(), None, "decode() should return None when data is None");
        }
    }

    #[test]
    fn test_decode_with_loaded_data() {
        // Test that decode() returns the correct sequence when data is loaded
        let seqcol = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        // Expected sequences from base.fa
        let expected_sequences = vec![
            ("chrX", "TTGGGGAA"),
            ("chr1", "GGAA"),
            ("chr2", "GCGC"),
        ];

        for (seq_record, (expected_name, expected_seq)) in seqcol.sequences.iter().zip(expected_sequences.iter()) {
            assert_eq!(seq_record.metadata.name, *expected_name);
            assert!(seq_record.data.is_some(), "load_fasta should load sequence data");

            let decoded = seq_record.decode().expect("decode() should return Some when data is present");
            assert_eq!(decoded, *expected_seq,
                "Decoded sequence for {} should match expected sequence", expected_name);
        }
    }

    #[test]
    fn test_decode_handles_encoded_data() {
        // Test that decode() correctly handles bit-packed encoded data
        use crate::alphabet::{lookup_alphabet, AlphabetType};

        let sequence = b"ACGT";
        let alphabet = lookup_alphabet(&AlphabetType::Dna2bit);
        let encoded_data = encode_sequence(sequence, alphabet);

        let record = SequenceRecord {
            metadata: SequenceMetadata {
                name: "test_seq".to_string(),
                length: 4,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Dna2bit,
            },
            fai: None,
            data: Some(encoded_data),
        };

        let decoded = record.decode().expect("Should decode encoded data");
        assert_eq!(decoded, "ACGT", "Should correctly decode bit-packed DNA sequence");
    }

    #[test]
    fn test_decode_handles_raw_utf8_data() {
        // Test that decode() handles raw UTF-8 data (not encoded)
        let raw_sequence = b"ACGTACGT".to_vec();

        let record = SequenceRecord {
            metadata: SequenceMetadata {
                name: "test_seq".to_string(),
                length: 8,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Ascii,
            },
            fai: None,
            data: Some(raw_sequence),
        };

        let decoded = record.decode().expect("Should decode raw UTF-8 data");
        assert_eq!(decoded, "ACGTACGT", "Should correctly decode raw sequence data");
    }

    #[test]
    fn test_decode_with_iupac_alphabet() {
        // Test decode with IUPAC DNA alphabet (4-bit encoding)
        use crate::alphabet::{lookup_alphabet, AlphabetType};

        let sequence = b"ACGTRYMK";
        let alphabet = lookup_alphabet(&AlphabetType::DnaIupac);
        let encoded_data = encode_sequence(sequence, alphabet);

        let record = SequenceRecord {
            metadata: SequenceMetadata {
                name: "iupac_test".to_string(),
                length: 8,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::DnaIupac,
            },
            fai: None,
            data: Some(encoded_data),
        };

        let decoded = record.decode().expect("Should decode IUPAC encoded data");
        assert_eq!(decoded, "ACGTRYMK", "Should correctly decode IUPAC DNA sequence");
    }

    #[test]
    fn test_decode_with_protein_alphabet() {
        // Test decode with protein alphabet (5-bit encoding)
        use crate::alphabet::{lookup_alphabet, AlphabetType};

        let sequence = b"ACDEFGHIKLMNPQRSTVWY";
        let alphabet = lookup_alphabet(&AlphabetType::Protein);
        let encoded_data = encode_sequence(sequence, alphabet);

        let record = SequenceRecord {
            metadata: SequenceMetadata {
                name: "protein_test".to_string(),
                length: 20,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Protein,
            },
            fai: None,
            data: Some(encoded_data),
        };

        let decoded = record.decode().expect("Should decode protein encoded data");
        assert_eq!(decoded, "ACDEFGHIKLMNPQRSTVWY", "Should correctly decode protein sequence");
    }

    #[test]
    fn test_decode_empty_sequence() {
        // Test decode with empty data
        let record = SequenceRecord {
            metadata: SequenceMetadata {
                name: "empty_seq".to_string(),
                length: 0,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Dna2bit,
            },
            fai: None,
            data: Some(Vec::new()),
        };

        let decoded = record.decode().expect("Should handle empty sequence");
        assert_eq!(decoded, "", "Empty sequence should decode to empty string");
    }

    #[test]
    fn test_decode_detects_raw_vs_encoded() {
        // Test that decode() correctly distinguishes between raw and encoded data
        use crate::alphabet::{lookup_alphabet, AlphabetType};

        // Create raw UTF-8 data that matches the sequence length (heuristic for raw data)
        let raw_sequence = b"GGGGGGGG";  // 8 bytes for 8-symbol sequence
        let record_raw = SequenceRecord {
            metadata: SequenceMetadata {
                name: "raw_test".to_string(),
                length: 8,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Dna2bit,
            },
            fai: None,
            data: Some(raw_sequence.to_vec()),
        };

        let decoded_raw = record_raw.decode().expect("Should decode raw data");
        assert_eq!(decoded_raw, "GGGGGGGG");

        // Create encoded data (2 bits per symbol, so 8 symbols = 2 bytes)
        let alphabet = lookup_alphabet(&AlphabetType::Dna2bit);
        let encoded_data = encode_sequence(b"GGGGGGGG", alphabet);
        assert_eq!(encoded_data.len(), 2, "Encoded 2-bit DNA should be 2 bytes for 8 symbols");

        let record_encoded = SequenceRecord {
            metadata: SequenceMetadata {
                name: "encoded_test".to_string(),
                length: 8,
                sha512t24u: "test_digest".to_string(),
                md5: "test_md5".to_string(),
                alphabet: AlphabetType::Dna2bit,
            },
            fai: None,
            data: Some(encoded_data),
        };

        let decoded_encoded = record_encoded.decode().expect("Should decode encoded data");
        assert_eq!(decoded_encoded, "GGGGGGGG");
    }

    #[test]
    fn test_sequence_collection_iterator() {
        // Test that SequenceCollection can be iterated over
        let collection = load_fasta("../tests/data/fasta/base.fa")
            .expect("Failed to load test FASTA file");

        // Test borrowing iterator (&collection)
        let mut count = 0;
        for seq in &collection {
            assert!(seq.metadata.length > 0, "Sequence should have length");
            count += 1;
        }
        assert_eq!(count, 3, "base.fa should have 3 sequences");

        // Collection should still be usable after borrowing iteration
        assert_eq!(collection.sequences.len(), 3);

        // Test consuming iterator (collection)
        let names: Vec<String> = collection
            .into_iter()
            .map(|seq| seq.metadata.name)
            .collect();

        assert_eq!(names.len(), 3);
        assert!(names.contains(&"chrX".to_string()));
        assert!(names.contains(&"chr1".to_string()));
        assert!(names.contains(&"chr2".to_string()));
    }
}
