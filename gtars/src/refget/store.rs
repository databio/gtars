use super::alphabet::{lookup_alphabet, AlphabetType};
use seq_io::fasta::{Reader, Record};
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::path::{Path, PathBuf};

use anyhow::anyhow;
use anyhow::{Context, Result};
use chrono::Utc;
use memmap2::Mmap;
use serde::{Deserialize, Serialize};
use std::fs::{self, create_dir_all, File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::str;
use flate2::read::GzDecoder;
use super::encoder::decode_substring_from_bytes;
use super::encoder::SequenceEncoder;
use crate::common::utils::get_dynamic_reader;
use crate::refget::fasta::read_fasta_refget_file;
use crate::refget::hashkeyable::HashKeyable;
use crate::uniwig::reading::parse_bedlike_file;
use crate::uniwig::utils::get_file_info;
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
        self.collections
            .get(&collection_digest)
            .ok_or_else(|| anyhow::anyhow!("Collection not found"))?;

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
                            data: Some(raw_sequence),
                        },
                        seqcol.digest.to_key(),
                    )?;
                }
                StorageMode::Encoded => {
                    // Create a SequenceEncoder to handle the encoding of the sequence.
                    let mut encoder = SequenceEncoder::new(dr.alphabet, dr.length);
                    for seq_line in record.seq_lines() {
                        encoder.update(&seq_line);
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

    /// Returns a list of all sequence names in the store
    pub fn list_sequence_digests(&self) -> Vec<[u8; 32]> {
        self.sequence_store.keys().cloned().collect()
    }

    /// Retrieve a SequenceRecord from the store by its SHA512t24u digest
    pub fn get_sequence_by_id<K: AsRef<[u8]>>(&self, seq_digest: K) -> Option<&SequenceRecord> {
        self.sequence_store.get(&seq_digest.to_key())
    }

    /// Retrieve a SequenceRecord from the store by its collection digest and name
    pub fn get_sequence_by_collection_and_name<K: AsRef<[u8]>>(
        &self,
        collection_digest: K,
        sequence_name: &str,
    ) -> Option<&SequenceRecord> {
        // Look up the collection by digest
        if let Some(name_map) = self.name_lookup.get(&collection_digest.to_key()) {
            // Look up the sequence name in the collection's name map
            if let Some(sha512t24u) = name_map.get(sequence_name) {
                // Retrieve the sequence record from the store
                return self.sequence_store.get(sha512t24u);
            }
        }
        None
    }


    /// Given a Sequence Collection digest and a bed file path
    /// retrieve sequences and write to a file
    pub fn get_seqs_bed_file<K: AsRef<[u8]>>(        &self,
                                                     collection_digest: K,
                                                     bed_file_path: &str,
                                                     output_file_path: &str) -> Result<()>{

        // Read each line of bed file
        // for each line in bed file retrieve sequences
        // write them to fasta file


        // Set up the output path and the directories along the way
        let output_path = std::path::Path::new(&output_file_path).parent().unwrap();
        let _ = create_dir_all(output_path);

        // Open file for writing to
        let mut output_file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(output_file_path)
            .unwrap();

        // Read in Bed file and get the sequences, writing to output file along the way
        let path = Path::new(bed_file_path);
        let pathbuf = PathBuf::from(bed_file_path);
        let file_info = get_file_info(&pathbuf);
        let is_gzipped = file_info.is_gzipped;
        let opened_bed_file = File::open(path).unwrap();

        let reader: Box<dyn Read> = match is_gzipped {
            true => Box::new(GzDecoder::new(opened_bed_file)),
            false => Box::new(opened_bed_file),
        };
        let reader = BufReader::new(reader);

        let mut parsed_chr: String = String::new();
        let mut previous_parsed_chr: String = String::new();

        // let mut seq_record: &SequenceRecord = SequenceRecord::new();

        let mut current_seq_digest: String = String::new();
        let mut current_header: String = String::new();
        let mut previous_header: String = String::new();

        for line in reader.lines() {
            //println!("Here is line{:?}", line);

            // Must use a 2nd let statement to appease the borrow-checker
            let line_string = line.unwrap();
            let s = line_string.as_str();

            let (parsed_chr, parsed_start, parsed_end) = parse_bedlike_file(s).unwrap();

            if previous_parsed_chr != parsed_chr{


                previous_parsed_chr = parsed_chr.clone();

                let result = self.get_sequence_by_collection_and_name(&collection_digest, &*parsed_chr).unwrap();
                current_seq_digest = result.metadata.sha512t24u.clone();
                current_header = format!(">{} {} {} {} {}", result.metadata.name, result.metadata.length, result.metadata.alphabet, result.metadata.sha512t24u, result.metadata.md5);

            }


            if previous_header != current_header{

                previous_header = current_header.clone();

                // WRITE NEW HEADER
                let mut header_to_be_written = current_header.clone();
                header_to_be_written.push('\n');
                output_file.write_all(header_to_be_written.as_ref()).unwrap();

            }

            // NOW WRITE SUBSTRING IF IT IS FOUND
            let retrieved_substring = self.get_substring(&current_seq_digest, parsed_start as usize, parsed_end as usize).unwrap();

            output_file.write_all(retrieved_substring.as_ref()).unwrap();

            //println!("{}", retrieved_substring);

        }


        Ok(())



    }

    /// Retrieve a SequenceRecord from the store by its MD5 digest
    pub fn get_sequence_by_md5<K: AsRef<[u8]>>(&self, seq_md5: K) -> Option<&SequenceRecord> {
        // Look up the SHA512t24u digest using the MD5 lookup
        if let Some(sha512_digest) = self.md5_lookup.get(&seq_md5.to_key()) {
            // Retrieve the sequence record from the store
            return self.sequence_store.get(sha512_digest);
        }
        None
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
        &self,
        sha512_digest: K,
        start: usize,
        end: usize,
    ) -> Option<String> {
        let record = self.sequence_store.get(&sha512_digest.to_key())?;
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
                let decoded_sequence = decode_substring_from_bytes(&sequence, start, end, alphabet);
                Some(String::from_utf8(decoded_sequence).expect("Invalid UTF-8"))
            }
            StorageMode::Raw => {
                let raw_slice: &[u8] = &sequence[start..end];
                println!("Raw sequence slice: {:?}", raw_slice);
                Some(String::from_utf8(raw_slice.to_vec()).expect("Invalid UTF-8"))
            }
        }
    }



    /// Helper function to get the relative path for a sequence based on its SHA512t24u digest string
    fn get_sequence_path(digest_str: &str, template: &str) -> PathBuf {
        let path_str = template
            .replace("%s2", &digest_str[0..2])
            .replace("%s", digest_str);

        PathBuf::from(path_str)
    }

    /// Load a GlobalRefgetStore from a directory
    pub fn load_from_directory<P: AsRef<Path>>(root_path: P) -> Result<Self> {
        let root_path = root_path.as_ref();

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

        // Load sequence metadata from the sequence index file
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

                // Generate the path using the template
                let path_str = metadata
                    .seqdata_path_template
                    .replace("%s2", &seq_metadata.sha512t24u[0..2])
                    .replace("%s4", &seq_metadata.sha512t24u[0..4])
                    .replace("%s", &seq_metadata.sha512t24u);
                let seq_path = root_path.join(path_str);

                // Open and memory-map the sequence file
                let file = File::open(&seq_path).context(format!(
                    "Failed to open sequence file: {}",
                    seq_path.display()
                ))?;

                let mmap =
                    unsafe { Mmap::map(&file) }.context("Failed to memory-map sequence file")?;

                // Convert to Vec<u8> for storage
                let data = mmap.to_vec();

                // Create a SequenceRecord
                let record = SequenceRecord {
                    metadata: seq_metadata.clone(),
                    data: Some(data),
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
        for (_digest, result_sr) in &self.sequence_store {
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
        for (_digest, record) in &self.sequence_store {
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
        for (_digest, collection) in &self.collections {
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
    use crate::refget::collection::{
        SeqColDigestLvl1, SequenceCollection, SequenceMetadata, SequenceRecord,
    };
    use crate::refget::digest::{md5, sha512t24u};
    use tempfile::tempdir;

    #[test]
    fn store_fa_to_farg() {
        // Create temporary directory
        let temp_dir = tempdir().expect("Failed to create temporary directory");
        let temp_path = temp_dir.path();

        // Copy test FASTA file to temp directory
        let test_fa = "tests/data/fasta/base.fa";
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


    #[test]
    fn test_refget_store_retrieve_seq(){

        println!("All tests passing.");

        let query_bed_str ="/home/drc/Downloads/test_adding_fastas/INPUT_FILE/testbed.bed";
        let mut store = GlobalRefgetStore::new(StorageMode::Encoded);
        store.import_fasta("/home/drc/Downloads/test_adding_fastas/INPUT_FILE/100linesGRCH38.fa").unwrap();

        let known_seq_col_digest = "VRAhrMWmnghggvNFd4qNoRg-J3AqugDh";

        let output_file = "/home/drc/Downloads/test_adding_fastas/OUTPUT_FILES/mytest.fa";

        let result  = store.get_seqs_bed_file(known_seq_col_digest,query_bed_str,output_file).unwrap();
        // let know_name = "1";
        //
        // let result = store.get_sequence_by_collection_and_name(known_seq_col_digest, know_name).unwrap();
        //
        // let seq_digest = result.metadata.sha512t24u.clone();
        // println!("{}", seq_digest);
        //
        // let retrieved_substring = store.get_substring(seq_digest, 100, 105).unwrap();
        //
        // println!("{}", retrieved_substring);
        println!("ok");

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
            has_data: false,
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
        let test_fa = "tests/data/fasta/base.fa";
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
            .write_store_to_directory("tests/store_test", seq_template)
            .unwrap();
    }

    #[test]
    fn test_disk_persistence() {
        // Create a temporary directory for the test
        let temp_dir = tempdir().unwrap();
        let temp_path = temp_dir.path();
        let temp_fasta = temp_path.join("base.fa.gz");
        std::fs::copy("tests/data/fasta/base.fa.gz", &temp_fasta)
            .expect("Failed to copy base.fa.gz to tempdir");

        // Create a new sequence store
        let mut store = GlobalRefgetStore::new(StorageMode::Encoded);

        // Import a FASTA file into the store
        // store.import_fasta("tests/data/subset.fa.gz").unwrap();
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
        let loaded_store = GlobalRefgetStore::load_from_directory(temp_path).unwrap();

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

        // Check data equality
        assert_eq!(original_seq1.data, loaded_seq1.data);
        assert_eq!(original_seq2.data, loaded_seq2.data);

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
}
