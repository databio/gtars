use super::alphabet::{AlphabetType, lookup_alphabet};
use seq_io::fasta::{Reader, Record};
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::path::{Path, PathBuf};

use super::encoder::SequenceEncoder;
use super::encoder::decode_substring_from_bytes;
use crate::fasta::read_fasta_refget_file;
use crate::hashkeyable::HashKeyable;
use anyhow::anyhow;
use anyhow::{Context, Result};
use chrono::Utc;
use flate2::read::GzDecoder;
use gtars_core::utils::{get_dynamic_reader, get_file_info, parse_bedlike_file};
use serde::{Deserialize, Serialize};
use std::fs::{self, File, OpenOptions, create_dir_all};
use std::io::{BufRead, BufReader, Read, Write};
use std::{io, str};
// Import the HashKeyable trait for converting types to a 32-byte key

// Import collection types
use super::collection::{SequenceCollection, SequenceMetadata, SequenceRecord};

// const DEFAULT_COLLECTION_ID: [u8; 32] = [0u8; 32]; // Default collection ID for the name lookup table

const DEFAULT_COLLECTION_ID: &str = "DEFAULT_REFGET_SEQUENCE_COLLECTION"; // Default collection ID for the name lookup table

/// Enum storing whether sequences will be stored in Raw or Encoded form
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum StorageMode {
    Raw,
    Encoded,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RetrievedSequence {
    pub sequence: String,
    pub chrom_name: String,
    pub start: u32,
    pub end: u32,
}

/// Global store handling cross-collection sequence management
/// Holds a global sequence_store, which holds all sequences (across collections) so that
/// sequences are deduplicated.
/// This allows lookup by sequence digest directly (bypassing collection information).
/// The GlobalRefgetStore also holds a collections hashmap, to provide lookup by collection+name
#[derive(Debug)]
pub struct GlobalRefgetStore {
    /// SHA512t24u digest -> SequenceRecord (metadata + optional data)
    sequence_store: HashMap<[u8; 32], SequenceRecord>,
    /// MD5 digest -> SHA512t24u digest lookup
    md5_lookup: HashMap<[u8; 32], [u8; 32]>,

    /// Collection digest -> {name -> SHA512t24u digest}
    name_lookup: HashMap<[u8; 32], HashMap<String, [u8; 32]>>,
    /// Active sequence collections
    collections: HashMap<[u8; 32], SequenceCollection>,
    /// Storage strategy for sequences
    mode: StorageMode,
    /// Where the store lives on disk (local store or cache directory)
    local_path: Option<PathBuf>,
    /// Where to pull sequences from (if remote-backed)
    remote_source: Option<String>,
    /// Template for sequence file paths (e.g., "sequences/%s2/%s.seq")
    seqdata_path_template: Option<String>,
}

/// Metadata for the entire store.
/// This is used to serialize metadata to a `index.json`, which can be loaded by the application.
#[derive(Serialize, Deserialize, Debug)]
struct StoreMetadata {
    /// Version of the metadata format
    version: u32,
    /// Template for sequence file paths
    seqdata_path_template: String,
    /// Template for collection file paths
    collections_path_template: String,
    /// Path to the sequence metadata index file
    sequence_index: String,
    /// Storage mode (Raw or Encoded)
    mode: StorageMode,
    /// Creation timestamp
    created_at: String,
}

pub struct GetSeqsBedFileIter<'a, K>
where
    K: AsRef<[u8]>,
{
    store: &'a mut GlobalRefgetStore,
    reader: BufReader<Box<dyn Read>>,
    collection_digest: K,
    previous_parsed_chr: String,
    current_seq_digest: String,
    line_num: usize,
}

impl<K> Iterator for GetSeqsBedFileIter<'_, K>
where
    K: AsRef<[u8]>,
{
    type Item = Result<RetrievedSequence, Box<dyn std::error::Error>>;

    fn next(&mut self) -> Option<Self::Item> {
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
                return Some(Err(err_str.into()));
            }
        };

        if parsed_start == -1 || parsed_end == -1 {
            let err_str = format!(
                "Error reading line {} due to invalid start or end coordinates: '{}'",
                self.line_num + 1,
                line_string
            );
            return Some(Err(err_str.into()));
        }

        if self.previous_parsed_chr != parsed_chr {
            self.previous_parsed_chr = parsed_chr.clone();

            let result = match self
                .store
                .get_sequence_by_collection_and_name(&self.collection_digest, &parsed_chr)
            {
                Some(seq_record) => seq_record,
                None => {
                    let err_str = format!(
                        "Warning: Skipping line {} because sequence '{}' not found in collection '{}'.",
                        self.line_num + 1,
                        parsed_chr,
                        String::from_utf8_lossy(self.collection_digest.as_ref())
                    );
                    return Some(Err(err_str.into()));
                }
            };

            self.current_seq_digest = result.metadata.sha512t24u.clone();
        }

        let retrieved_substring = match self.store.get_substring(
            &self.current_seq_digest,
            parsed_start as usize,
            parsed_end as usize,
        ) {
            Some(substring) => substring,
            None => {
                let err_str = format!(
                    "Warning: Skipping line {} because substring for digest '{}' from {} to {} not found or invalid.",
                    self.line_num + 1,
                    self.current_seq_digest,
                    parsed_start,
                    parsed_end
                );
                return Some(Err(err_str.into()));
            }
        };

        Some(Ok(RetrievedSequence {
            sequence: retrieved_substring,
            chrom_name: parsed_chr,
            start: parsed_start as u32, // Convert i32 to u32
            end: parsed_end as u32,     // Convert i32 to u32
        }))
    }
}

impl GlobalRefgetStore {
    /// Generic constructor. Creates a new, empty `GlobalRefgetStore`.
    pub fn new(mode: StorageMode) -> Self {
        // Initialize the name lookup with a default collection
        let mut name_lookup = HashMap::new();
        name_lookup.insert(DEFAULT_COLLECTION_ID.to_key(), HashMap::new());

        GlobalRefgetStore {
            sequence_store: HashMap::new(),
            md5_lookup: HashMap::new(),
            name_lookup,
            collections: HashMap::new(),
            mode,
            local_path: None,
            remote_source: None,
            seqdata_path_template: None,
        }
    }

    /// Adds a sequence to the Store
    /// Ensure that it is added to the appropriate collection.
    /// If no collection is specified, it will be added to the default collection.
    // Using Into here  instead of the Option direction allows us to accept
    // either None or [u8; 32], without having to wrap it in Some().
    pub fn add_sequence<T: Into<Option<[u8; 32]>>>(
        &mut self,
        sequence_record: SequenceRecord,
        collection_digest: T,
    ) -> Result<()> {
        // Ensure collection exists; otherwise use the default collection
        let collection_digest = collection_digest
            .into()
            .unwrap_or(DEFAULT_COLLECTION_ID.to_key());
        self.collections.get(&collection_digest).ok_or_else(|| {
            anyhow::anyhow!("Collection not found for digest: {:?}", collection_digest)
        })?;

        // Add to name lookup for the collection
        self.name_lookup
            .entry(collection_digest)
            .or_default()
            .insert(
                sequence_record.metadata.name.clone(),
                sequence_record.metadata.sha512t24u.to_key(),
            );

        // Finally, add SequenceRecord to store (consuming the object)
        self.add_sequence_record(sequence_record)?;

        Ok(())
    }

    /// Adds a collection, and all sequences in it, to the store.
    pub fn add_sequence_collection(&mut self, collection: SequenceCollection) -> Result<()> {
        let coll_digest = collection.digest.to_key();

        // Register the collection
        self.collections.insert(coll_digest, collection.clone());

        // Add all sequences in the collection to the store
        for sequence_record in collection.sequences {
            self.add_sequence(sequence_record, coll_digest)?;
        }

        Ok(())
    }

    // Adds SequenceRecord to the store.
    // Should only be used internally, via `add_sequence`, which ensures sequences are added to collections.
    fn add_sequence_record(&mut self, sr: SequenceRecord) -> Result<()> {
        self.md5_lookup
            .insert(sr.metadata.md5.to_key(), sr.metadata.sha512t24u.to_key());
        self.sequence_store
            .insert(sr.metadata.sha512t24u.to_key(), sr);
        Ok(())
    }

    // Loading the sequence data requires 2 passes through the FASTA file, because we
    // have to guess the alphabet before we can encode. So, the first pass is the digesting
    // function, which digests and guesses the alphabet, and calculates length, to produce the
    // SequenceDigest object. Then, using this object, we can go through the sequence a second time and encode it.
    pub fn import_fasta<P: AsRef<Path>>(&mut self, file_path: P) -> Result<()> {
        println!("Loading farg index...");
        let seqcol = SequenceCollection::from_fasta(&file_path)?;

        // Register the collection
        self.add_sequence_collection(seqcol.clone())?;

        // Local hashmap to store SequenceMetadata (digests)
        let mut seqmeta_hashmap: HashMap<String, SequenceMetadata> = HashMap::new();
        let seqcol_sequences = seqcol.sequences.clone(); // Clone to avoid partial move
        for record in seqcol_sequences {
            let seqmeta = record.metadata;
            seqmeta_hashmap.insert(seqmeta.name.clone(), seqmeta);
        }

        let file_reader = get_dynamic_reader(file_path.as_ref())?;
        let mut fasta_reader = Reader::new(file_reader);

        println!("Preparing to load sequences into refget SeqColStore...");

        while let Some(record) = fasta_reader.next() {
            let record = record?;
            let id = record.id()?;
            let dr = seqmeta_hashmap[id].clone();
            println!("Digest result: {:?}", dr);

            match self.mode {
                StorageMode::Raw => {
                    let mut raw_sequence = Vec::with_capacity(dr.length);
                    // For raw, just extend with the line content.
                    for seq_line in record.seq_lines() {
                        raw_sequence.extend(seq_line);
                    }
                    println!(
                        "Storing raw sequence. Name: {}; Alphabet: {}; Digest: {}",
                        id, dr.alphabet, dr.sha512t24u
                    );

                    self.add_sequence(
                        SequenceRecord {
                            metadata: dr,
                            fai: None,  // Store doesn't preserve FAI data
                            data: Some(raw_sequence),
                        },
                        seqcol.digest.to_key(),
                    )?;
                }
                StorageMode::Encoded => {
                    // Create a SequenceEncoder to handle the encoding of the sequence.
                    let mut encoder = SequenceEncoder::new(dr.alphabet, dr.length);
                    for seq_line in record.seq_lines() {
                        encoder.update(seq_line);
                    }
                    // let encoded_sequence = BitVec::<u8, Msb0>::from_vec(encoder.finalize());
                    let encoded_sequence = encoder.finalize();

                    println!(
                        "Storing encoded sequence. Name: {}; Alphabet: {}; Digest: {}",
                        id, dr.alphabet, dr.sha512t24u
                    );
                    self.add_sequence(
                        SequenceRecord {
                            metadata: dr,
                            fai: None,  // Store doesn't preserve FAI data
                            data: Some(encoded_sequence),
                        },
                        seqcol.digest.to_key(),
                    )?;
                }
            }
        }

        println!("Finished loading sequences into refget SeqColStore.");
        Ok(())
    }

    /// Returns a list of all sequence digests in the store
    pub fn list_sequence_digests(&self) -> Vec<[u8; 32]> {
        self.sequence_store.keys().cloned().collect()
    }

    /// Returns a list of all sequences with their metadata
    pub fn list_sequences(&self) -> Vec<&SequenceMetadata> {
        self.sequence_store.values().map(|rec| &rec.metadata).collect()
    }

    /// Returns a list of all collection digests
    pub fn list_collection_digests(&self) -> Vec<[u8; 32]> {
        self.collections.keys().cloned().collect()
    }

    /// Returns a list of all collections
    pub fn list_collections(&self) -> Vec<&SequenceCollection> {
        self.collections.values().collect()
    }

    /// Returns the local path where the store is located (if any)
    pub fn local_path(&self) -> Option<&PathBuf> {
        self.local_path.as_ref()
    }

    /// Returns the remote source URL (if any)
    pub fn remote_source(&self) -> Option<&str> {
        self.remote_source.as_deref()
    }

    /// Retrieve a SequenceRecord from the store by its SHA512t24u digest
    pub fn get_sequence_by_id<K: AsRef<[u8]>>(&mut self, seq_digest: K) -> Option<&SequenceRecord> {
        let digest_key = seq_digest.to_key();

        // Ensure the sequence data is loaded
        if let Err(e) = self.ensure_sequence_loaded(&digest_key) {
            eprintln!("Failed to load sequence: {}", e);
            return None;
        }

        self.sequence_store.get(&digest_key)
    }

    /// Retrieve a SequenceRecord from the store by its collection digest and name
    pub fn get_sequence_by_collection_and_name<K: AsRef<[u8]>>(
        &mut self,
        collection_digest: K,
        sequence_name: &str,
    ) -> Option<&SequenceRecord> {
        let collection_key = collection_digest.to_key();

        // Try to ensure collection is loaded (lazy-load from remote if needed)
        if let Err(e) = self.ensure_collection_loaded(&collection_key) {
            eprintln!("Failed to load collection: {}", e);
            return None;
        }

        // Look up the collection by digest
        let digest_key = if let Some(name_map) = self.name_lookup.get(&collection_key) {
            // Look up the sequence name in the collection's name map
            name_map.get(sequence_name).cloned()?
        } else {
            return None;
        };

        // Ensure the sequence data is loaded
        if let Err(e) = self.ensure_sequence_loaded(&digest_key) {
            eprintln!("Failed to load sequence: {}", e);
            return None;
        }

        // Retrieve the sequence record from the store
        self.sequence_store.get(&digest_key)
    }

    pub fn get_seqs_in_bed_file_iter<'a, K: AsRef<[u8]>>(
        &'a mut self,
        collection_digest: K,
        bed_file_path: &str,
    ) -> Result<GetSeqsBedFileIter<'a, K>, Box<dyn std::error::Error>> {
        let path = Path::new(bed_file_path);
        let file_info = get_file_info(path);
        let is_gzipped = file_info.is_gzipped;

        let opened_bed_file = File::open(path)?;

        let reader: Box<dyn Read> = match is_gzipped {
            true => Box::new(GzDecoder::new(BufReader::new(opened_bed_file))),
            false => Box::new(opened_bed_file),
        };
        let reader = BufReader::new(reader);

        Ok(GetSeqsBedFileIter {
            store: self,
            reader,
            collection_digest,
            previous_parsed_chr: String::new(),
            current_seq_digest: String::new(),
            line_num: 0,
        })
    }

    /// Given a Sequence Collection digest and a bed file path
    /// retrieve sequences and write to a file
    pub fn get_seqs_bed_file<K: AsRef<[u8]>>(
        &mut self,
        collection_digest: K,
        bed_file_path: &str,
        output_file_path: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Set up the output path and create directories if they don't exist
        let output_path = Path::new(output_file_path).parent().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "Invalid output file path: parent directory not found",
            )
        })?;

        create_dir_all(output_path)?;

        // Open file for writing to
        let mut output_file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(output_file_path)?;

        // Pre-fetch all sequence metadata from the collection to avoid borrowing issues
        let collection_key = collection_digest.as_ref().to_key();
        let name_to_metadata: HashMap<String, (String, usize, AlphabetType, String, String)> = self
            .name_lookup
            .get(&collection_key)
            .map(|name_map| {
                name_map
                    .iter()
                    .filter_map(|(name, seq_digest)| {
                        self.sequence_store.get(seq_digest).map(|record| {
                            (
                                name.clone(),
                                (
                                    record.metadata.name.clone(),
                                    record.metadata.length,
                                    record.metadata.alphabet,
                                    record.metadata.sha512t24u.clone(),
                                    record.metadata.md5.clone(),
                                ),
                            )
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();

        let seq_iter = self.get_seqs_in_bed_file_iter(&collection_digest, bed_file_path)?;

        let mut previous_parsed_chr = String::new();
        let mut current_header: String = String::new();
        let mut previous_header: String = String::new();

        for rs in seq_iter.into_iter() {
            if let Err(err) = rs {
                eprintln!("{err}");
                continue;
            }
            let rs = rs.unwrap();

            if previous_parsed_chr != rs.chrom_name {
                previous_parsed_chr = rs.chrom_name.clone();

                // Look up metadata from our pre-fetched map
                if let Some((name, length, alphabet, sha512, md5)) =
                    name_to_metadata.get(&rs.chrom_name)
                {
                    current_header = format!(
                        ">{} {} {} {} {}",
                        name, length, alphabet, sha512, md5
                    );
                }
            }

            let retrieved_substring = rs.sequence;

            if previous_header != current_header {
                let prefix = if previous_header.is_empty() { "" } else { "\n" };

                previous_header = current_header.clone();

                // Combine the prefix, current_header, and a trailing newline
                let header_to_be_written = format!("{}{}\n", prefix, current_header);
                output_file.write_all(header_to_be_written.as_bytes())?;
            }

            output_file.write_all(retrieved_substring.as_ref())?;
        }

        Ok(())
    }

    /// Given a Sequence Collection digest and a bed file path,
    /// retrieve sequences and return them as a Vec<RetrievedSequence>.
    /// Lines that cannot be parsed or sequences that cannot be retrieved will be skipped
    pub fn get_seqs_bed_file_to_vec<K: AsRef<[u8]>>(
        &mut self,
        collection_digest: K,
        bed_file_path: &str,
    ) -> Result<Vec<RetrievedSequence>, Box<dyn std::error::Error>> {
        let seq_iter = self.get_seqs_in_bed_file_iter(collection_digest, bed_file_path)?;

        let retrieved_sequences = seq_iter
            .into_iter()
            .filter_map(|rs| match rs {
                Ok(rs) => Some(rs),
                Err(err) => {
                    eprintln!("{err}");
                    None
                }
            })
            .collect::<Vec<RetrievedSequence>>();

        Ok(retrieved_sequences)
    }

    /// Retrieve a SequenceRecord from the store by its MD5 digest
    pub fn get_sequence_by_md5<K: AsRef<[u8]>>(&mut self, seq_md5: K) -> Option<&SequenceRecord> {
        // Look up the SHA512t24u digest using the MD5 lookup
        let sha512_digest = self.md5_lookup.get(&seq_md5.to_key()).cloned()?;

        // Ensure the sequence data is loaded
        if let Err(e) = self.ensure_sequence_loaded(&sha512_digest) {
            eprintln!("Failed to load sequence: {}", e);
            return None;
        }

        // Retrieve the sequence record from the store
        self.sequence_store.get(&sha512_digest)
    }

    /// Retrieves a substring from an encoded sequence by its SHA512t24u digest.
    ///
    /// # Arguments
    ///
    /// * `sha512_digest` - The SHA512t24u digest of the sequence
    /// * `start` - The start index of the substring (inclusive)
    /// * `end` - The end index of the substring (exclusive)
    ///
    /// # Returns
    ///
    /// The substring if the sequence is found, or None if not found
    pub fn get_substring<K: AsRef<[u8]>>(
        &mut self,
        sha512_digest: K,
        start: usize,
        end: usize,
    ) -> Option<String> {
        let digest_key = sha512_digest.to_key();

        // Ensure the sequence data is loaded
        if let Err(e) = self.ensure_sequence_loaded(&digest_key) {
            eprintln!("Failed to load sequence: {}", e);
            return None;
        }

        let record = self.sequence_store.get(&digest_key)?;
        let sequence = record.data.as_ref()?;

        if start >= record.metadata.length || end > record.metadata.length || start >= end {
            println!(
                "Invalid substring range: start={}, end={}, sequence length={}",
                start, end, record.metadata.length
            );
            return None;
        }

        match self.mode {
            StorageMode::Encoded => {
                let alphabet = lookup_alphabet(&record.metadata.alphabet);
                let decoded_sequence = decode_substring_from_bytes(sequence, start, end, alphabet);
                Some(String::from_utf8(decoded_sequence).expect("Invalid UTF-8"))
            }
            StorageMode::Raw => {
                let raw_slice: &[u8] = &sequence[start..end];
                println!("Raw sequence slice: {:?}", raw_slice);
                match String::from_utf8(raw_slice.to_vec()) {
                    Ok(raw_string) => Some(raw_string),
                    Err(e) => {
                        eprintln!("Failed to decode UTF-8 sequence: {}", e);
                        None
                    }
                }
            }
        }
    }

    /// Export sequences from a collection to a FASTA file
    ///
    /// # Arguments
    /// * `collection_digest` - The digest of the collection to export from
    /// * `output_path` - Path to write the FASTA file
    /// * `sequence_names` - Optional list of sequence names to export.
    ///                      If None, exports all sequences in the collection.
    /// * `line_width` - Optional line width for wrapping sequences (default: 80)
    ///
    /// # Returns
    /// Result indicating success or error
    pub fn export_fasta<K: AsRef<[u8]>, P: AsRef<Path>>(
        &mut self,
        collection_digest: K,
        output_path: P,
        sequence_names: Option<Vec<&str>>,
        line_width: Option<usize>,
    ) -> Result<()> {
        let line_width = line_width.unwrap_or(80);
        let output_path = output_path.as_ref();
        let collection_key = collection_digest.as_ref().to_key();

        // Get the name map for this collection and build a map of name -> digest
        let name_to_digest: HashMap<String, [u8; 32]> = self
            .name_lookup
            .get(&collection_key)
            .ok_or_else(|| {
                anyhow!("Collection not found: {:?}", String::from_utf8_lossy(collection_digest.as_ref()))
            })?
            .clone();

        // Determine which sequences to export
        let names_to_export: Vec<String> = if let Some(names) = sequence_names {
            // Filter to only requested names
            names.iter().map(|s| s.to_string()).collect()
        } else {
            // Export all sequences in the collection
            name_to_digest.keys().cloned().collect()
        };

        // Create output file
        let mut output_file = File::create(output_path)
            .context(format!("Failed to create output file: {}", output_path.display()))?;

        // Export each sequence
        for seq_name in names_to_export {
            // Get the sequence digest from the name map
            let seq_digest = name_to_digest.get(&seq_name).ok_or_else(|| {
                anyhow!("Sequence '{}' not found in collection", seq_name)
            })?;

            // Ensure sequence is loaded
            self.ensure_sequence_loaded(seq_digest)?;

            // Get the sequence record
            let record = self.sequence_store.get(seq_digest).ok_or_else(|| {
                anyhow!("Sequence record not found for digest: {:?}", seq_digest)
            })?;

            // Get the sequence data
            let sequence_data = record.data.as_ref().ok_or_else(|| {
                anyhow!("Sequence data not loaded for '{}'", seq_name)
            })?;

            // Decode the sequence based on storage mode
            let decoded_sequence = match self.mode {
                StorageMode::Encoded => {
                    let alphabet = lookup_alphabet(&record.metadata.alphabet);
                    let decoded = decode_substring_from_bytes(
                        sequence_data,
                        0,
                        record.metadata.length,
                        alphabet,
                    );
                    String::from_utf8(decoded).context("Failed to decode sequence as UTF-8")?
                }
                StorageMode::Raw => {
                    String::from_utf8(sequence_data.clone())
                        .context("Failed to decode raw sequence as UTF-8")?
                }
            };

            // Write FASTA header
            writeln!(output_file, ">{}", seq_name)?;

            // Write sequence with line wrapping
            for chunk in decoded_sequence.as_bytes().chunks(line_width) {
                output_file.write_all(chunk)?;
                output_file.write_all(b"\n")?;
            }
        }

        Ok(())
    }

    /// Export sequences by their digests to a FASTA file
    ///
    /// # Arguments
    /// * `digests` - List of SHA512t24u digests to export
    /// * `output_path` - Path to write the FASTA file
    /// * `line_width` - Optional line width for wrapping sequences (default: 80)
    ///
    /// # Returns
    /// Result indicating success or error
    pub fn export_fasta_by_digests<P: AsRef<Path>>(
        &mut self,
        digests: Vec<&str>,
        output_path: P,
        line_width: Option<usize>,
    ) -> Result<()> {
        let line_width = line_width.unwrap_or(80);
        let output_path = output_path.as_ref();

        // Create output file
        let mut output_file = File::create(output_path)
            .context(format!("Failed to create output file: {}", output_path.display()))?;

        // Export each sequence
        for digest_str in digests {
            let digest_key = digest_str.as_bytes().to_key();

            // Ensure sequence is loaded
            self.ensure_sequence_loaded(&digest_key)?;

            // Get the sequence record
            let record = self.sequence_store.get(&digest_key).ok_or_else(|| {
                anyhow!("Sequence record not found for digest: {}", digest_str)
            })?;

            // Get the sequence data
            let sequence_data = record.data.as_ref().ok_or_else(|| {
                anyhow!("Sequence data not loaded for digest: {}", digest_str)
            })?;

            // Decode the sequence based on storage mode
            let decoded_sequence = match self.mode {
                StorageMode::Encoded => {
                    let alphabet = lookup_alphabet(&record.metadata.alphabet);
                    let decoded = decode_substring_from_bytes(
                        sequence_data,
                        0,
                        record.metadata.length,
                        alphabet,
                    );
                    String::from_utf8(decoded).context("Failed to decode sequence as UTF-8")?
                }
                StorageMode::Raw => {
                    String::from_utf8(sequence_data.clone())
                        .context("Failed to decode raw sequence as UTF-8")?
                }
            };

            // Write FASTA header with sequence name
            writeln!(output_file, ">{}", record.metadata.name)?;

            // Write sequence with line wrapping
            for chunk in decoded_sequence.as_bytes().chunks(line_width) {
                output_file.write_all(chunk)?;
                output_file.write_all(b"\n")?;
            }
        }

        Ok(())
    }

    /// Helper function to get the relative path for a sequence based on its SHA512t24u digest string
    fn get_sequence_path(digest_str: &str, template: &str) -> PathBuf {
        let path_str = template
            .replace("%s2", &digest_str[0..2])
            .replace("%s", digest_str);

        PathBuf::from(path_str)
    }

    /// Helper function to fetch a file from local path or remote source
    /// Returns the file contents as Vec<u8>
    fn fetch_file(
        local_path: &Path,
        remote_source: &Option<String>,
        relative_path: &str,
    ) -> Result<Vec<u8>> {
        let full_local_path = local_path.join(relative_path);

        // Check if file exists locally
        if full_local_path.exists() {
            return fs::read(&full_local_path)
                .context(format!("Failed to read local file: {}", full_local_path.display()));
        }

        // If not local and we have a remote source, fetch from remote
        if let Some(remote_url) = remote_source {
            let full_remote_url = if remote_url.ends_with('/') {
                format!("{}{}", remote_url, relative_path)
            } else {
                format!("{}/{}", remote_url, relative_path)
            };

            let response = ureq::get(&full_remote_url)
                .call()
                .map_err(|e| anyhow!("Failed to fetch from remote: {}", e))?;

            let mut data = Vec::new();
            response
                .into_reader()
                .read_to_end(&mut data)
                .context("Failed to read response body")?;

            // Create parent directory if needed
            if let Some(parent) = full_local_path.parent() {
                create_dir_all(parent)?;
            }

            // Save to local cache
            fs::write(&full_local_path, &data)
                .context(format!("Failed to cache file to: {}", full_local_path.display()))?;

            Ok(data)
        } else {
            Err(anyhow!(
                "File not found locally and no remote source configured: {}",
                relative_path
            ))
        }
    }

    /// Load a local RefgetStore from a directory
    /// This loads metadata only; sequence data is loaded on-demand
    pub fn load_local<P: AsRef<Path>>(cache_path: P) -> Result<Self> {
        let root_path = cache_path.as_ref();

        // Read and parse index.json
        let index_path = root_path.join("index.json");
        let json = fs::read_to_string(&index_path).context(format!(
            "Failed to read index.json from {}",
            index_path.display()
        ))?;

        let metadata: StoreMetadata =
            serde_json::from_str(&json).context("Failed to parse index.json")?;

        // Create a new empty store with the correct mode
        let mut store = GlobalRefgetStore::new(metadata.mode);
        store.local_path = Some(root_path.to_path_buf());
        store.seqdata_path_template = Some(metadata.seqdata_path_template.clone());

        // Load sequence metadata from the sequence index file (metadata only, no data)
        let sequence_index_path = root_path.join(&metadata.sequence_index);
        if sequence_index_path.exists() {
            let file = std::fs::File::open(&sequence_index_path)?;
            let reader = std::io::BufReader::new(file);

            for line in reader.lines() {
                let line = line?;

                // Skip comment lines
                if line.starts_with('#') {
                    continue;
                }

                // Parse sequence metadata lines
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() != 5 {
                    continue; // Skip lines that don't have exactly 5 columns
                }

                let seq_metadata = SequenceMetadata {
                    name: parts[0].to_string(),
                    length: parts[1].parse().unwrap_or(0),
                    alphabet: parts[2].parse().unwrap_or(AlphabetType::Unknown),
                    sha512t24u: parts[3].to_string(),
                    md5: parts[4].to_string(),
                };

                // Create a SequenceRecord with no data (lazy loading)
                let record = SequenceRecord {
                    metadata: seq_metadata.clone(),
                    fai: None,
                    data: None,  // Data will be loaded on-demand
                };

                // Add to store
                let sha512_key = seq_metadata.sha512t24u.to_key();
                store.sequence_store.insert(sha512_key, record);

                // Add to MD5 lookup
                let md5_key = seq_metadata.md5.to_key();
                store.md5_lookup.insert(md5_key, sha512_key);
            }
        }

        // Load collections from the collections directory
        let collections_dir = root_path.join("collections");
        if collections_dir.exists() {
            for entry in fs::read_dir(&collections_dir)? {
                let entry = entry?;
                let path = entry.path();

                if path.is_file() && path.extension() == Some(std::ffi::OsStr::new("farg")) {
                    // Load the collection from the .farg file
                    let collection = read_fasta_refget_file(&path)?;
                    let collection_digest = collection.digest.to_key();

                    // Add collection to store
                    store
                        .collections
                        .insert(collection_digest, collection.clone());

                    // Build name lookup for this collection
                    let mut name_map = HashMap::new();
                    for sequence_record in &collection.sequences {
                        let sha512_key = sequence_record.metadata.sha512t24u.to_key();
                        name_map.insert(sequence_record.metadata.name.clone(), sha512_key);
                    }
                    store.name_lookup.insert(collection_digest, name_map);
                }
            }
        }

        Ok(store)
    }

    /// Load a remote-backed RefgetStore
    /// This loads metadata from remote and caches sequence data on-demand
    pub fn load_remote<P: AsRef<Path>, S: AsRef<str>>(
        cache_path: P,
        remote_url: S,
    ) -> Result<Self> {
        let cache_path = cache_path.as_ref();
        let remote_url = remote_url.as_ref().to_string();

        // Create cache directory if it doesn't exist
        create_dir_all(cache_path)?;

        // Fetch index.json from remote
        let index_data = Self::fetch_file(cache_path, &Some(remote_url.clone()), "index.json")?;
        let json = String::from_utf8(index_data)
            .context("index.json contains invalid UTF-8")?;

        let metadata: StoreMetadata =
            serde_json::from_str(&json).context("Failed to parse index.json")?;

        // Create a new empty store with the correct mode
        let mut store = GlobalRefgetStore::new(metadata.mode);
        store.local_path = Some(cache_path.to_path_buf());
        store.remote_source = Some(remote_url.clone());
        store.seqdata_path_template = Some(metadata.seqdata_path_template.clone());

        // Fetch sequence index from remote
        let sequence_index_data = Self::fetch_file(
            cache_path,
            &Some(remote_url.clone()),
            &metadata.sequence_index,
        )?;
        let sequence_index_str = String::from_utf8(sequence_index_data)
            .context("sequence index contains invalid UTF-8")?;

        // Parse sequence metadata (metadata only, no data)
        for line in sequence_index_str.lines() {
            // Skip comment lines
            if line.starts_with('#') {
                continue;
            }

            // Parse sequence metadata lines
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() != 5 {
                continue;
            }

            let seq_metadata = SequenceMetadata {
                name: parts[0].to_string(),
                length: parts[1].parse().unwrap_or(0),
                alphabet: parts[2].parse().unwrap_or(AlphabetType::Unknown),
                sha512t24u: parts[3].to_string(),
                md5: parts[4].to_string(),
            };

            // Create a SequenceRecord with no data (lazy loading)
            let record = SequenceRecord {
                metadata: seq_metadata.clone(),
                fai: None,
                data: None,
            };

            // Add to store
            let sha512_key = seq_metadata.sha512t24u.to_key();
            store.sequence_store.insert(sha512_key, record);

            // Add to MD5 lookup
            let md5_key = seq_metadata.md5.to_key();
            store.md5_lookup.insert(md5_key, sha512_key);
        }

        // Load collections from the collections directory
        let collections_dir_path = "collections";
        let local_collections_dir = cache_path.join(collections_dir_path);

        // Try to read from local cache first, or create if doesn't exist
        if !local_collections_dir.exists() {
            create_dir_all(&local_collections_dir)?;
        }

        // Note: We'll try to fetch collection files on-demand when needed
        // For now, just check if any exist locally
        if local_collections_dir.exists() {
            for entry in fs::read_dir(&local_collections_dir)? {
                let entry = entry?;
                let path = entry.path();

                if path.is_file() && path.extension() == Some(std::ffi::OsStr::new("farg")) {
                    // Load the collection from the .farg file
                    let collection = read_fasta_refget_file(&path)?;
                    let collection_digest = collection.digest.to_key();

                    // Add collection to store
                    store
                        .collections
                        .insert(collection_digest, collection.clone());

                    // Build name lookup for this collection
                    let mut name_map = HashMap::new();
                    for sequence_record in &collection.sequences {
                        let sha512_key = sequence_record.metadata.sha512t24u.to_key();
                        name_map.insert(sequence_record.metadata.name.clone(), sha512_key);
                    }
                    store.name_lookup.insert(collection_digest, name_map);
                }
            }
        }

        Ok(store)
    }

    /// Ensure a collection is loaded into the store
    /// If the collection is not already loaded, try to fetch it from local or remote
    fn ensure_collection_loaded(&mut self, collection_digest: &[u8; 32]) -> Result<()> {
        // Check if collection is already loaded
        if self.name_lookup.contains_key(collection_digest) {
            return Ok(());
        }

        // Collection not found - try to fetch it if we have a remote source
        // Convert the collection digest to a string for the path
        let digest_str = String::from_utf8_lossy(collection_digest);

        // Build the collection file path using the template
        // Default template is "collections/%s.farg"
        let relative_path = format!("collections/{}.farg", digest_str);

        // Try to fetch the collection file
        let local_path = self
            .local_path
            .as_ref()
            .ok_or_else(|| anyhow!("No local path configured"))?;

        // Try to fetch (this will use cache if available)
        // fetch_file will save the file to local_path/relative_path
        let _collection_data = Self::fetch_file(local_path, &self.remote_source, &relative_path)?;

        // Read the collection from the cached file
        let collection_file_path = local_path.join(&relative_path);
        let collection = read_fasta_refget_file(&collection_file_path)?;

        // Verify the collection digest matches what we requested
        let loaded_digest = collection.digest.to_key();
        if loaded_digest != *collection_digest {
            eprintln!(
                "Warning: Collection digest mismatch. Expected {:?}, got {:?}",
                String::from_utf8_lossy(collection_digest),
                String::from_utf8_lossy(&loaded_digest)
            );
        }

        // Add collection to store
        self.collections.insert(*collection_digest, collection.clone());

        // Build name lookup for this collection
        let mut name_map = HashMap::new();
        for sequence_record in &collection.sequences {
            let sha512_key = sequence_record.metadata.sha512t24u.to_key();
            name_map.insert(sequence_record.metadata.name.clone(), sha512_key);
        }
        self.name_lookup.insert(*collection_digest, name_map);

        Ok(())
    }

    /// Ensure a sequence is loaded into memory
    /// If the sequence data is not already loaded, fetch it from local or remote
    fn ensure_sequence_loaded(&mut self, digest: &[u8; 32]) -> Result<()> {
        // Check if sequence exists
        let record = self
            .sequence_store
            .get(digest)
            .ok_or_else(|| anyhow!("Sequence not found in store"))?;

        // If data is already loaded, return early
        if record.data.is_some() {
            return Ok(());
        }

        // Get the necessary information before borrowing mutably
        let digest_str = &record.metadata.sha512t24u;
        let template = self
            .seqdata_path_template
            .as_ref()
            .ok_or_else(|| anyhow!("No sequence data path template configured"))?;

        // Build the relative path using the template
        let relative_path = template
            .replace("%s2", &digest_str[0..2])
            .replace("%s4", &digest_str[0..4])
            .replace("%s", digest_str);

        // Fetch the sequence data
        let local_path = self
            .local_path
            .as_ref()
            .ok_or_else(|| anyhow!("No local path configured"))?;

        let data = Self::fetch_file(local_path, &self.remote_source, &relative_path)?;

        // Update the record with the loaded data
        if let Some(record) = self.sequence_store.get_mut(digest) {
            record.data = Some(data);
        }

        Ok(())
    }

    /// Write the sequence_store component to a FARG file (without collection headers).
    /// * `file_path` - The path to the FARG file to be written.
    pub fn to_farg<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
        // Write the FARGI file
        let file_path = file_path.as_ref();
        println!("Writing farg file: {:?}", file_path);
        let mut file = std::fs::File::create(file_path)?;

        // Write header with digest metadata
        writeln!(file, "#name\tlength\talphabet\tsha512t24u\tmd5")?;

        // Write sequence data
        for result_sr in self.sequence_store.values() {
            let result = result_sr.metadata.clone();
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}",
                result.name, result.length, result.alphabet, result.sha512t24u, result.md5
            )?;
        }
        Ok(())
    }

    /// Write a GlobalRefgetStore object to a directory
    pub fn write_store_to_directory<P: AsRef<Path>>(
        &self,
        root_path: P,
        seqdata_path_template: &str,
    ) -> Result<()> {
        let root_path = root_path.as_ref();
        println!(
            "Writing store to directory: {}; Using seqdata path template: {}",
            root_path.display(),
            seqdata_path_template
        );

        // Create the root directory if it doesn't exist
        fs::create_dir_all(root_path)?;

        // Create sequences directory
        let sequences_dir = root_path.join("sequences");
        fs::create_dir_all(&sequences_dir)?;

        // Create collections directory
        let collections_dir = root_path.join("collections");
        fs::create_dir_all(&collections_dir)?;

        // Write each sequence to its own file
        for record in self.sequence_store.values() {
            if record.data.is_some() {
                // Get the path for this sequence using the template and base64url-encoded digest
                let rel_path =
                    Self::get_sequence_path(&record.metadata.sha512t24u, seqdata_path_template);
                let full_path = root_path.join(&rel_path);

                // Write the sequence data to file
                record.to_file(full_path)?;
            } else {
                return Err(anyhow!(
                    "Sequence data is missing for digest: {}",
                    &record.metadata.sha512t24u
                ));
            }
        }

        // Write each collection to its own .farg file
        for collection in self.collections.values() {
            let collection_file_path =
                root_path.join(format!("collections/{}.farg", collection.digest));
            collection.to_farg_path(&collection_file_path)?;
        }

        // Write the sequence metadata index file
        let sequence_index_path = root_path.join("sequences.farg");
        self.to_farg(&sequence_index_path)?;

        // Create the metadata structure
        let metadata = StoreMetadata {
            version: 1,
            seqdata_path_template: seqdata_path_template.to_string(),
            collections_path_template: "collections/%s.farg".to_string(),
            sequence_index: "sequences.farg".to_string(),
            mode: self.mode,
            created_at: Utc::now().to_rfc3339(),
        };

        // Write metadata to index.json
        let json = serde_json::to_string_pretty(&metadata)
            .context("Failed to serialize metadata to JSON")?;
        fs::write(root_path.join("index.json"), json).context("Failed to write index.json")?;

        Ok(())
    }
}

impl Display for GlobalRefgetStore {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "SeqColStore object:")?;
        writeln!(f, ">Sequences (n={}):", self.sequence_store.len())?;
        // Print out the sequences in the store
        for (i, (sha512_digest, sequence_record)) in self.sequence_store.iter().take(10).enumerate()
        {
            let seq = sequence_record.data.as_ref().unwrap();
            // Extract the first 3 characters of the sequence (or fewer if the sequence is shorter)
            let first_8_chars = match self.mode {
                StorageMode::Encoded => {
                    let alphabet = lookup_alphabet(&sequence_record.metadata.alphabet);
                    let decoded = decode_substring_from_bytes(
                        seq,
                        0,
                        8.min(sequence_record.metadata.length),
                        alphabet,
                    );
                    String::from_utf8(decoded).unwrap_or_else(|_| "???".to_string())
                }
                StorageMode::Raw => String::from_utf8(seq[0..8.min(seq.len())].to_vec())
                    .unwrap_or_else(|_| "???".to_string()),
            };

            writeln!(
                f,
                "   - {}. {:02x?}, MD5: {:02x?}, Length: {}, Alphabet: {:?}, Start: {}",
                i + 1,
                std::str::from_utf8(sha512_digest).unwrap(),
                &sequence_record.metadata.md5,
                &sequence_record.metadata.length,
                &sequence_record.metadata.alphabet,
                first_8_chars
            )?;
        }
        writeln!(f, ">Collections (n={:?}):", self.name_lookup.len())?;
        // Print out the collections in the store
        for (i, (digest, name_map)) in self.name_lookup.iter().enumerate() {
            // Convert the digest to a hex string
            let seqcol_digest_str = String::from_utf8_lossy(digest);
            writeln!(
                f,
                "  {}. Collection Digest: {:02x?}",
                i + 1,
                seqcol_digest_str
            )?;
            for (name, sha512_digest) in name_map {
                // Convert the sha512_digest to a hex string
                let sha512_str = String::from_utf8_lossy(sha512_digest);
                writeln!(f, "   - Name: {}, SHA512: {:02x?}", name, sha512_str)?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use std::time::Instant;
    use crate::collection::{
        SeqColDigestLvl1, SequenceCollection, SequenceMetadata, SequenceRecord,
    };
    use crate::digest::{md5, sha512t24u};
    use tempfile::tempdir;

    #[test]
    fn store_fa_to_farg() {
        // Create temporary directory
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Copy test FASTA file to temp directory
        let test_fa = "../tests/data/fasta/base.fa";
        let temp_fa = temp_path.join("base.fa");
        let temp_farg = temp_path.join("base.farg");

        std::fs::copy(test_fa, &temp_fa).expect("Failed to copy test FASTA file");

        // Create sequence collection from temporary file
        let seqcol = SequenceCollection::from_fasta(&temp_fa)
            .expect("Failed to create SeqColDigest from FASTA file");

        // Write FARG to temporary directory
        seqcol.to_farg().expect("Failed to write farg file");

        // Load and verify
        let loaded_seqcol = read_fasta_refget_file(&temp_farg).expect("Failed to read refget file");

        // Test round-trip integrity
        for (original, loaded) in seqcol.sequences.iter().zip(loaded_seqcol.sequences.iter()) {
            assert_eq!(original.metadata.name, loaded.metadata.name);
            assert_eq!(original.metadata.length, loaded.metadata.length);
            assert_eq!(original.metadata.sha512t24u, loaded.metadata.sha512t24u);
            assert_eq!(original.metadata.md5, loaded.metadata.md5);
            assert_eq!(original.metadata.alphabet, loaded.metadata.alphabet);
        }
    }

    // Helper function to calculate actual digests for testing
    fn calculate_test_digests(sequence: &[u8]) -> (String, String) {
        (sha512t24u(sequence), md5(sequence))
    }

    #[test]
    fn test_refget_store_retrieve_seq_and_vec() {
        // Create temporary directory for all test files
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // --- 1. Prepare Test FASTA Data ---
        let fasta_content = "\
>chr1
ATGCATGCATGC
>chr2
GGGGAAAA
";
        let temp_fasta_path = temp_path.join("test.fa");

        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        // --- 2. Initialize GlobalRefgetStore and import FASTA ---
        let mut store = GlobalRefgetStore::new(StorageMode::Encoded);
        store.import_fasta(&temp_fasta_path).unwrap();

        let sequence_keys: Vec<[u8; 32]> = store.sequence_store.keys().cloned().collect();

        let _ = sequence_keys[0]; //ww1QMyfFm1f4qa3fRLqqJGafIeEuZR1V
        let _ = sequence_keys[1]; //OyXJErGtjgcIVSdobGkHE3sBdQ5faDTf
        let collection_digest_ref: &str = "uC_UorBNf3YUu1YIDainBhI94CedlNeH";

        // Calculate expected SHA512t24u and MD5 for test sequences
        let (chr1_sha, chr1_md5) = calculate_test_digests(b"ATGCATGCATGC");
        let (chr2_sha, chr2_md5) = calculate_test_digests(b"GGGGAAAA");
        println!("chr1 values: {}  {}", chr1_sha, chr1_md5);
        println!("chr2 values: {}  {}", chr2_sha, chr2_md5);

        // --- 3. Prepare Test BED Data ---
        let bed_content = "\
chr1\t0\t5
chr1\t8\t12
chr2\t0\t4
chr_nonexistent\t10\t20
chr1\t-5\t100
";
        let temp_bed_path = temp_path.join("test.bed");

        fs::write(&temp_bed_path, bed_content).expect("Failed to write test BED file");

        let temp_output_fa_path = temp_path.join("output.fa");

        store
            .get_seqs_bed_file(
                collection_digest_ref,
                temp_bed_path.to_str().unwrap(),
                temp_output_fa_path.to_str().unwrap(),
            )
            .expect("get_seqs_bed_file failed");

        // Read the output FASTA file and verify its content
        let output_fa_content =
            fs::read_to_string(&temp_output_fa_path).expect("Failed to read output FASTA file");

        // Expected output content (headers and sequences should match the logic of the function)
        let expected_fa_content = format!(
            ">chr1 12 dna2bit {} {}\nATGCAATGC\n>chr2 8 dna2bit {} {}\nGGGG\n",
            chr1_sha, chr1_md5, chr2_sha, chr2_md5
        );
        assert_eq!(
            output_fa_content.trim(),
            expected_fa_content.trim(),
            "Output FASTA file content mismatch"
        );
        println!("✓ get_seqs_bed_file test passed.");

        // --- Test get_seqs_bed_file_to_vec (returns Vec<RetrievedSequence>) ---
        let vec_result = store
            .get_seqs_bed_file_to_vec(collection_digest_ref, temp_bed_path.to_str().unwrap())
            .expect("get_seqs_bed_file_to_vec failed");

        // Define the expected vector of RetrievedSequence structs
        let expected_vec = vec![
            RetrievedSequence {
                sequence: "ATGCA".to_string(),
                chrom_name: "chr1".to_string(),
                start: 0,
                end: 5,
            },
            RetrievedSequence {
                sequence: "ATGC".to_string(),
                chrom_name: "chr1".to_string(),
                start: 8,
                end: 12,
            },
            RetrievedSequence {
                sequence: "GGGG".to_string(),
                chrom_name: "chr2".to_string(),
                start: 0,
                end: 4,
            },
        ];

        // Assert that the returned vector matches the expected vector
        assert_eq!(
            vec_result, expected_vec,
            "Retrieved sequence vector mismatch"
        );
        println!("✓ get_seqs_bed_file_to_vec test passed.");
    }

    #[test]
    fn test_global_refget_store() {
        let sequence = b"ACGT";
        let name = "test_seq";
        println!("Testing GlobalRefgetStore with sequence: {}", name);

        // Create a sequence collection
        let mut collection = SequenceCollection {
            sequences: Vec::new(),
            digest: "test_collection".to_string(),
            lvl1: SeqColDigestLvl1 {
                names_digest: "test".to_string(),
                sequences_digest: "test".to_string(),
                lengths_digest: "test".to_string(),
            },
            file_path: None,
        };

        // Create a sequence record
        let seq_metadata = SequenceMetadata {
            name: name.to_string(),
            length: sequence.len(),
            sha512t24u: sha512t24u(sequence),
            md5: md5(sequence),
            alphabet: AlphabetType::Dna2bit,
        };

        let record = SequenceRecord {
            metadata: seq_metadata.clone(),
            fai: None,  // Test data doesn't need FAI
            data: Some(sequence.to_vec()),
        };

        collection.sequences.push(record);

        // Add the sequence to the store
        let mut store = GlobalRefgetStore::new(StorageMode::Encoded);
        store.add_sequence_collection(collection.clone()).unwrap();

        // Verify the store has the sequence
        assert!(!store.sequence_store.is_empty());

        // Test sequence lookup by collection+name (using string digest)
        let retrieved_by_name_str =
            store.get_sequence_by_collection_and_name(&collection.digest, name);
        assert!(retrieved_by_name_str.is_some());
        let retrieved_record = retrieved_by_name_str.unwrap();
        assert_eq!(retrieved_record.metadata.name, name);
        assert_eq!(retrieved_record.data.as_ref().unwrap(), sequence);

        // Test sequence lookup by collection+name (using [u8; 32] digest)
        let retrieved_by_name_key =
            store.get_sequence_by_collection_and_name(collection.digest.to_key(), name);
        assert!(retrieved_by_name_key.is_some());
        let retrieved_record = retrieved_by_name_key.unwrap();
        assert_eq!(retrieved_record.metadata.name, name);
        assert_eq!(retrieved_record.data.as_ref().unwrap(), sequence);

        // Test sequence lookup by SHA512 digest (using string)
        let retrieved_by_sha512_str = store.get_sequence_by_id(&seq_metadata.sha512t24u);
        assert!(retrieved_by_sha512_str.is_some());
        let retrieved_record = retrieved_by_sha512_str.unwrap();
        assert_eq!(retrieved_record.metadata.name, name);
        assert_eq!(retrieved_record.data.as_ref().unwrap(), sequence);

        // Test sequence lookup by SHA512 digest (using [u8; 32])
        let retrieved_by_sha512_key = store.get_sequence_by_id(seq_metadata.sha512t24u.to_key());
        assert!(retrieved_by_sha512_key.is_some());
        let retrieved_record = retrieved_by_sha512_key.unwrap();
        assert_eq!(retrieved_record.metadata.name, name);
        assert_eq!(retrieved_record.data.as_ref().unwrap(), sequence);

        // Test sequence lookup by MD5 digest (using string)
        let retrieved_by_md5_str = store.get_sequence_by_md5(&seq_metadata.md5);
        assert!(retrieved_by_md5_str.is_some());
        let retrieved_record = retrieved_by_md5_str.unwrap();
        assert_eq!(retrieved_record.metadata.name, name);
        assert_eq!(retrieved_record.data.as_ref().unwrap(), sequence);

        // Test sequence lookup by MD5 digest (using [u8; 32])
        let retrieved_by_md5_key = store.get_sequence_by_md5(seq_metadata.md5.to_key());
        assert!(retrieved_by_md5_key.is_some());
        let retrieved_record = retrieved_by_md5_key.unwrap();
        assert_eq!(retrieved_record.metadata.name, name);
        assert_eq!(retrieved_record.data.as_ref().unwrap(), sequence);
    }

    #[test]
    fn test_import_fasta() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Copy test FASTA file to temp directory
        let test_fa = "../tests/data/fasta/base.fa";
        let temp_fa = temp_path.join("base.fa");

        std::fs::copy(test_fa, &temp_fa).expect("Failed to copy test FASTA file");

        let mut store = GlobalRefgetStore::new(StorageMode::Encoded);

        // Import the FASTA file
        store.import_fasta(temp_fa).unwrap();

        // Check that the store has sequences
        assert!(!store.sequence_store.is_empty());

        // Try writing to a file
        let seq_template = "sequences/%s2/%s.seq";
        // let col_template = "collections/%s.farg";
        store
            .write_store_to_directory(temp_path.to_str().unwrap(), seq_template)
            .unwrap();
    }

    #[test]
    fn test_disk_persistence() {
        // Create a temporary directory for the test
        let temp_dir = tempdir().unwrap();
        let temp_path = temp_dir.path();
        let temp_fasta = temp_path.join("base.fa.gz");
        std::fs::copy("../tests/data/fasta/base.fa.gz", &temp_fasta)
            .expect("Failed to copy base.fa.gz to tempdir");

        // Create a new sequence store
        let mut store = GlobalRefgetStore::new(StorageMode::Encoded);

        // Import a FASTA file into the store
        // store.import_fasta("../tests/data/subset.fa.gz").unwrap();
        store.import_fasta(&temp_fasta).unwrap();

        // Get the sequence keys for verification (assuming we know the test file contains 3 sequences)
        let sequence_keys: Vec<[u8; 32]> = store.sequence_store.keys().cloned().collect();
        assert_eq!(
            sequence_keys.len(),
            3,
            "Test file should contain exactly 3 sequences"
        );

        let sha512_key1 = sequence_keys[0];
        let sha512_key2 = sequence_keys[1];

        // Store original sequences for comparison
        let original_seq1 = store.sequence_store.get(&sha512_key1).unwrap().clone();
        let original_seq2 = store.sequence_store.get(&sha512_key2).unwrap().clone();

        // Write the store to the temporary directory
        let seq_template = "sequences/%s2/%s.seq";
        store
            .write_store_to_directory(temp_path, seq_template)
            .unwrap();

        // Verify that the files were created
        assert!(temp_path.join("sequences").exists());
        assert!(temp_path.join("sequences").read_dir().unwrap().count() > 0);
        assert!(temp_path.join("index.json").exists());
        assert!(temp_path.join("sequences.farg").exists());
        assert!(temp_path.join("collections").exists());

        // Load the store from disk
        let mut loaded_store = GlobalRefgetStore::load_local(temp_path).unwrap();

        // Verify that the loaded store has the same sequences
        assert_eq!(loaded_store.sequence_store.len(), 3);

        // Verify that we can retrieve sequences by their keys
        assert!(loaded_store.sequence_store.contains_key(&sha512_key1));
        assert!(loaded_store.sequence_store.contains_key(&sha512_key2));

        // Verify the content of the sequences
        let loaded_seq1 = loaded_store.sequence_store.get(&sha512_key1).unwrap();
        let loaded_seq2 = loaded_store.sequence_store.get(&sha512_key2).unwrap();

        // Check metadata equality
        assert_eq!(original_seq1.metadata.name, loaded_seq1.metadata.name);
        assert_eq!(original_seq1.metadata.length, loaded_seq1.metadata.length);
        assert_eq!(
            original_seq1.metadata.sha512t24u,
            loaded_seq1.metadata.sha512t24u
        );
        assert_eq!(original_seq1.metadata.md5, loaded_seq1.metadata.md5);

        assert_eq!(original_seq2.metadata.name, loaded_seq2.metadata.name);
        assert_eq!(original_seq2.metadata.length, loaded_seq2.metadata.length);
        assert_eq!(
            original_seq2.metadata.sha512t24u,
            loaded_seq2.metadata.sha512t24u
        );
        assert_eq!(original_seq2.metadata.md5, loaded_seq2.metadata.md5);

        // Check data is not loaded initially (lazy loading)
        assert_eq!(loaded_seq1.data, None, "Data should not be loaded initially with lazy loading");
        assert_eq!(loaded_seq2.data, None, "Data should not be loaded initially with lazy loading");

        // Verify MD5 lookup is preserved
        assert_eq!(loaded_store.md5_lookup.len(), 3);

        // Verify collections are preserved
        assert_eq!(loaded_store.collections.len(), store.collections.len());

        // Test sequence retrieval functionality
        for (digest, original_record) in &store.sequence_store {
            let loaded_record = loaded_store.get_sequence_by_id(*digest).unwrap();
            assert_eq!(original_record.metadata.name, loaded_record.metadata.name);
            assert_eq!(
                original_record.metadata.length,
                loaded_record.metadata.length
            );

            // Test substring retrieval works on loaded store
            if original_record.metadata.length > 0 {
                let substring_len = std::cmp::min(5, original_record.metadata.length);
                let substring = loaded_store.get_substring(digest, 0, substring_len);
                assert!(
                    substring.is_some(),
                    "Should be able to retrieve substring from loaded sequence"
                );
                println!("Do we ever get here?");
            }
        }

        println!("✓ Disk persistence test passed - all data preserved correctly");
    }

    #[test]
    fn test_export_fasta_all_sequences() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Create test FASTA
        let fasta_content = "\
>chr1
ATGCATGCATGC
>chr2
GGGGAAAA
>chr3
TTTTCCCC
";
        let temp_fasta_path = temp_path.join("test.fa");
        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        // Import into store
        let mut store = GlobalRefgetStore::new(StorageMode::Encoded);
        store.import_fasta(&temp_fasta_path).unwrap();

        // Get the collection digest
        let collections: Vec<_> = store.collections.keys().cloned().collect();
        assert_eq!(collections.len(), 1, "Should have exactly one collection");
        let collection_digest = collections[0];

        // Export all sequences
        let output_path = temp_path.join("exported_all.fa");
        store
            .export_fasta(&collection_digest, &output_path, None, Some(80))
            .expect("Failed to export FASTA");

        // Read and verify the exported file
        let exported_content = fs::read_to_string(&output_path).expect("Failed to read exported file");

        // Check that all sequences are present
        assert!(exported_content.contains(">chr1"), "Should contain chr1");
        assert!(exported_content.contains(">chr2"), "Should contain chr2");
        assert!(exported_content.contains(">chr3"), "Should contain chr3");
        assert!(exported_content.contains("ATGCATGCATGC"), "Should contain chr1 sequence");
        assert!(exported_content.contains("GGGGAAAA"), "Should contain chr2 sequence");
        assert!(exported_content.contains("TTTTCCCC"), "Should contain chr3 sequence");

        println!("✓ Export all sequences test passed");
    }

    #[test]
    fn test_export_fasta_subset_sequences() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Create test FASTA
        let fasta_content = "\
>chr1
ATGCATGCATGC
>chr2
GGGGAAAA
>chr3
TTTTCCCC
";
        let temp_fasta_path = temp_path.join("test.fa");
        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        // Import into store
        let mut store = GlobalRefgetStore::new(StorageMode::Encoded);
        store.import_fasta(&temp_fasta_path).unwrap();

        // Get the collection digest
        let collections: Vec<_> = store.collections.keys().cloned().collect();
        let collection_digest = collections[0];

        // Export only chr1 and chr3
        let output_path = temp_path.join("exported_subset.fa");
        let subset_names = vec!["chr1", "chr3"];
        store
            .export_fasta(&collection_digest, &output_path, Some(subset_names), Some(80))
            .expect("Failed to export subset FASTA");

        // Read and verify the exported file
        let exported_content = fs::read_to_string(&output_path).expect("Failed to read exported file");

        // Check that only selected sequences are present
        assert!(exported_content.contains(">chr1"), "Should contain chr1");
        assert!(!exported_content.contains(">chr2"), "Should NOT contain chr2");
        assert!(exported_content.contains(">chr3"), "Should contain chr3");
        assert!(exported_content.contains("ATGCATGCATGC"), "Should contain chr1 sequence");
        assert!(!exported_content.contains("GGGGAAAA"), "Should NOT contain chr2 sequence");
        assert!(exported_content.contains("TTTTCCCC"), "Should contain chr3 sequence");

        println!("✓ Export subset sequences test passed");
    }

    #[test]
    fn test_export_fasta_roundtrip() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Create test FASTA with longer sequences
        let fasta_content = "\
>seq1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>seq2
GGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCCTTTTGGGGAAAACCCC
";
        let temp_fasta_path = temp_path.join("original.fa");
        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        // Import into store
        let mut store1 = GlobalRefgetStore::new(StorageMode::Encoded);
        store1.import_fasta(&temp_fasta_path).unwrap();

        // Get original digests
        let original_digests: Vec<String> = store1
            .sequence_store
            .values()
            .map(|r| r.metadata.sha512t24u.clone())
            .collect();

        // Export to new FASTA
        let collections: Vec<_> = store1.collections.keys().cloned().collect();
        let collection_digest = collections[0];
        let exported_path = temp_path.join("exported.fa");
        store1
            .export_fasta(&collection_digest, &exported_path, None, Some(60))
            .expect("Failed to export FASTA");

        // Re-import the exported FASTA
        let mut store2 = GlobalRefgetStore::new(StorageMode::Encoded);
        store2.import_fasta(&exported_path).unwrap();

        // Verify digests match (same sequences)
        let new_digests: Vec<String> = store2
            .sequence_store
            .values()
            .map(|r| r.metadata.sha512t24u.clone())
            .collect();

        assert_eq!(original_digests.len(), new_digests.len(), "Should have same number of sequences");
        for digest in original_digests {
            assert!(
                new_digests.contains(&digest),
                "Digest {} should be present after roundtrip",
                digest
            );
        }

        println!("✓ Export/import roundtrip test passed");
    }

    #[test]
    fn test_export_fasta_by_digests() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Create test FASTA
        let fasta_content = "\
>chr1
ATGCATGCATGC
>chr2
GGGGAAAA
";
        let temp_fasta_path = temp_path.join("test.fa");
        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        // Import into store
        let mut store = GlobalRefgetStore::new(StorageMode::Encoded);
        store.import_fasta(&temp_fasta_path).unwrap();

        // Get digests
        let digests: Vec<String> = store
            .sequence_store
            .values()
            .map(|r| r.metadata.sha512t24u.clone())
            .collect();

        // Export by digests
        let output_path = temp_path.join("exported_by_digests.fa");
        let digest_refs: Vec<&str> = digests.iter().map(|s| s.as_str()).collect();
        store
            .export_fasta_by_digests(digest_refs, &output_path, Some(80))
            .expect("Failed to export FASTA by digests");

        // Read and verify the exported file
        let exported_content = fs::read_to_string(&output_path).expect("Failed to read exported file");

        // Check that all sequences are present
        assert!(exported_content.contains(">chr1"), "Should contain chr1");
        assert!(exported_content.contains(">chr2"), "Should contain chr2");
        assert!(exported_content.contains("ATGCATGCATGC"), "Should contain chr1 sequence");
        assert!(exported_content.contains("GGGGAAAA"), "Should contain chr2 sequence");

        println!("✓ Export by digests test passed");
    }

    #[test]
    fn test_export_fasta_error_handling() {
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Create test FASTA
        let fasta_content = "\
>chr1
ATGCATGCATGC
";
        let temp_fasta_path = temp_path.join("test.fa");
        fs::write(&temp_fasta_path, fasta_content).expect("Failed to write test FASTA file");

        // Import into store
        let mut store = GlobalRefgetStore::new(StorageMode::Encoded);
        store.import_fasta(&temp_fasta_path).unwrap();

        // Test with non-existent collection
        let output_path = temp_path.join("should_fail.fa");
        let fake_collection = b"fake_collection_digest_12345678";
        let result = store.export_fasta(fake_collection, &output_path, None, Some(80));
        assert!(result.is_err(), "Should fail with non-existent collection");

        // Test with non-existent sequence name
        let collections: Vec<_> = store.collections.keys().cloned().collect();
        let collection_digest = collections[0];
        let result = store.export_fasta(
            &collection_digest,
            &output_path,
            Some(vec!["nonexistent_chr"]),
            Some(80),
        );
        assert!(result.is_err(), "Should fail with non-existent sequence name");

        println!("✓ Error handling test passed");
    }
}
