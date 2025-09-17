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

    /// Flag indicating whether sequence data is actually loaded in memory.
    /// When false, only metadata is available.
    pub has_data: bool,
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
            has_data: true,
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
            seqcol.to_farg()?;
            println!("Farg file written to {:?}", farg_file_path);
        } else {
            println!(
                "Farg file already exists, not writing: {:?}",
                farg_file_path
            );
        }
        Ok(seqcol)
    }

    /// Write the SequenceCollection to a FASTA refget (FARG) file
    /// * `file_path` - The path to the FARG file to be written.
    pub fn to_farg_path<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        // Write the FARGI file
        let file_path = file_path.as_ref();
        println!("Writing farg file: {:?}", file_path);
        let mut file = std::fs::File::create(file_path)?;

        // Write header with digest metadata
        writeln!(file, "##seqcol_digest={}", self.digest)?;
        writeln!(file, "##names_digest={}", self.lvl1.names_digest)?;
        writeln!(file, "##sequences_digest={}", self.lvl1.sequences_digest)?;
        writeln!(file, "##lengths_digest={}", self.lvl1.lengths_digest)?;
        writeln!(file, "#name\tlength\talphabet\tsha512t24u\tmd5")?;

        // Write sequence data
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

    /// Write the SeqColDigest to a FARG file, using the file path stored in the struct.
    pub fn to_farg(&self) -> Result<()> {
        if let Some(ref file_path) = self.file_path {
            let farg_file_path = file_path.replace_exts_with("farg");
            self.to_farg_path(farg_file_path)
        } else {
            Err(anyhow::anyhow!(
                "No file path specified for FARG output. Use `to_farg_path` to specify a file path."
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
