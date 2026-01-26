//! WASM bindings for gtars-refget functionality.
//!
//! This module exposes refget digest computation functions to JavaScript/TypeScript.
//! It provides functionality to compute sequence collection (seqcol) digests from FASTA data
//! directly in the browser without requiring filesystem access.
//!
//! # Streaming API
//!
//! For large files, use the streaming API to process data in chunks:
//! ```javascript
//! const handle = fasta_hasher_new();
//! await fetch(url).then(response => {
//!     const reader = response.body.getReader();
//!     while (true) {
//!         const { done, value } = await reader.read();
//!         if (done) break;
//!         fasta_hasher_update(handle, value);
//!     }
//! });
//! const result = fasta_hasher_finish(handle);
//! ```

use wasm_bindgen::prelude::*;
use gtars_refget::digest::{
    sha512t24u, md5, canonicalize_json,
    digest_fasta_bytes, FastaStreamHasher, SequenceCollection,
};
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::sync::Mutex;

// ============================================================================
// Global storage for streaming hasher instances
// ============================================================================

/// Global storage for hasher instances, protected by a mutex.
/// Uses a simple counter for handle IDs.
static HASHER_STORAGE: Mutex<Option<HasherStorage>> = Mutex::new(None);

struct HasherStorage {
    hashers: HashMap<u32, FastaStreamHasher>,
    next_id: u32,
}

impl HasherStorage {
    fn new() -> Self {
        Self {
            hashers: HashMap::new(),
            next_id: 1,
        }
    }

    fn insert(&mut self, hasher: FastaStreamHasher) -> u32 {
        let id = self.next_id;
        self.next_id = self.next_id.wrapping_add(1);
        if self.next_id == 0 {
            self.next_id = 1; // Skip 0 as it's reserved for errors
        }
        self.hashers.insert(id, hasher);
        id
    }

    fn get_mut(&mut self, id: u32) -> Option<&mut FastaStreamHasher> {
        self.hashers.get_mut(&id)
    }

    fn remove(&mut self, id: u32) -> Option<FastaStreamHasher> {
        self.hashers.remove(&id)
    }
}

fn with_storage<F, R>(f: F) -> R
where
    F: FnOnce(&mut HasherStorage) -> R,
{
    let mut guard = HASHER_STORAGE.lock().unwrap();
    if guard.is_none() {
        *guard = Some(HasherStorage::new());
    }
    f(guard.as_mut().unwrap())
}

// ============================================================================
// Batch API (existing functions)
// ============================================================================

/// Compute a seqcol (sequence collection) digest from FASTA content.
///
/// This function takes FASTA file content as bytes and returns a JSON object
/// containing the sequence collection digest and metadata. Automatically handles
/// gzip-compressed FASTA files.
///
/// # Arguments
/// * `fasta_content` - The FASTA file content as a Uint8Array (supports gzip)
///
/// # Returns
/// A JavaScript object containing:
/// - `digest`: The top-level seqcol digest
/// - `names_digest`: Level 1 digest of sequence names
/// - `sequences_digest`: Level 1 digest of sequence digests
/// - `lengths_digest`: Level 1 digest of sequence lengths
/// - `n_sequences`: Number of sequences
/// - `sequences`: Array of sequence metadata objects
///
/// # Example (JavaScript)
/// ```javascript
/// const fastaContent = new TextEncoder().encode(">chr1\nACGT\n");
/// const result = digest_seqcol(fastaContent);
/// console.log(result.digest);  // Top-level seqcol digest
/// ```
#[wasm_bindgen]
pub fn digest_seqcol(fasta_content: &[u8]) -> Result<JsValue, JsError> {
    let collection = digest_fasta_bytes(fasta_content)
        .map_err(|e| JsError::new(&format!("Failed to digest FASTA: {}", e)))?;

    // Convert to a JS-friendly format
    let result = SeqColResult::from_collection(&collection);
    serde_wasm_bindgen::to_value(&result)
        .map_err(|e| JsError::new(&format!("Serialization error: {}", e)))
}

/// Compute a sequence digest (sha512t24u) from raw sequence data.
///
/// This is the GA4GH standard digest for sequences. The sequence is uppercased
/// before computing the digest.
///
/// # Arguments
/// * `sequence` - The sequence string (A, C, G, T, etc.)
///
/// # Returns
/// The sha512t24u digest string (base64url encoded)
///
/// # Example (JavaScript)
/// ```javascript
/// const digest = sequence_digest("ACGT");
/// console.log(digest);  // "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
/// ```
#[wasm_bindgen]
pub fn sequence_digest(sequence: &str) -> String {
    sha512t24u(sequence.to_ascii_uppercase())
}

/// Compute an MD5 digest from raw sequence data.
///
/// The sequence is uppercased before computing the digest.
///
/// # Arguments
/// * `sequence` - The sequence string
///
/// # Returns
/// The MD5 digest as a hex string
#[wasm_bindgen]
pub fn sequence_md5(sequence: &str) -> String {
    md5(sequence.to_ascii_uppercase())
}

/// Compute the sha512t24u digest of a string.
///
/// This is a low-level function that computes the GA4GH standard digest
/// of arbitrary input. Unlike `sequence_digest`, this does NOT uppercase
/// the input.
///
/// # Arguments
/// * `input` - The input string
///
/// # Returns
/// The sha512t24u digest string (base64url encoded)
#[wasm_bindgen]
pub fn compute_sha512t24u(input: &str) -> String {
    sha512t24u(input)
}

/// Compute the MD5 digest of a string.
///
/// Low-level function that does NOT uppercase the input.
///
/// # Arguments
/// * `input` - The input string
///
/// # Returns
/// The MD5 digest as a hex string
#[wasm_bindgen]
pub fn compute_md5(input: &str) -> String {
    md5(input)
}

/// Canonicalize a JSON string according to RFC-8785.
///
/// This is useful for computing digests of JSON objects in a deterministic way.
///
/// # Arguments
/// * `json_str` - A JSON string
///
/// # Returns
/// The canonicalized JSON string
#[wasm_bindgen]
pub fn canonicalize_json_string(json_str: &str) -> Result<String, JsError> {
    let value: serde_json::Value = serde_json::from_str(json_str)
        .map_err(|e| JsError::new(&format!("Invalid JSON: {}", e)))?;
    Ok(canonicalize_json(&value))
}

// ============================================================================
// Streaming API
// ============================================================================

/// Create a new streaming FASTA hasher.
///
/// Returns a handle (u32) that must be used in subsequent calls.
/// The handle must be freed with `fasta_hasher_finish()` or `fasta_hasher_free()`.
///
/// # Returns
/// A handle ID (> 0) on success, or 0 on error.
///
/// # Example (JavaScript)
/// ```javascript
/// const handle = fasta_hasher_new();
/// // ... use handle with fasta_hasher_update ...
/// const result = fasta_hasher_finish(handle);
/// ```
#[wasm_bindgen]
pub fn fasta_hasher_new() -> u32 {
    with_storage(|storage| storage.insert(FastaStreamHasher::new()))
}

/// Process a chunk of FASTA data.
///
/// This can be called multiple times with successive chunks.
/// Handles both plain text and gzip-compressed FASTA.
///
/// # Arguments
/// * `handle` - The hasher handle from `fasta_hasher_new()`
/// * `chunk` - A chunk of FASTA data as Uint8Array
///
/// # Returns
/// `true` on success, `false` on error.
///
/// # Example (JavaScript)
/// ```javascript
/// const handle = fasta_hasher_new();
/// fasta_hasher_update(handle, new TextEncoder().encode(">chr1\nACGT"));
/// fasta_hasher_update(handle, new TextEncoder().encode("\n>chr2\nTGCA\n"));
/// ```
#[wasm_bindgen]
pub fn fasta_hasher_update(handle: u32, chunk: &[u8]) -> bool {
    with_storage(|storage| {
        if let Some(hasher) = storage.get_mut(handle) {
            hasher.update(chunk).is_ok()
        } else {
            false
        }
    })
}

/// Finalize the hasher and return the results.
///
/// This consumes the hasher and frees its resources.
/// After calling this, the handle is no longer valid.
///
/// # Arguments
/// * `handle` - The hasher handle from `fasta_hasher_new()`
///
/// # Returns
/// A JavaScript object with the sequence collection result,
/// or an error if the handle is invalid.
///
/// # Example (JavaScript)
/// ```javascript
/// const handle = fasta_hasher_new();
/// fasta_hasher_update(handle, fastaData);
/// const result = fasta_hasher_finish(handle);
/// console.log(result.digest);
/// ```
#[wasm_bindgen]
pub fn fasta_hasher_finish(handle: u32) -> Result<JsValue, JsError> {
    let hasher = with_storage(|storage| storage.remove(handle))
        .ok_or_else(|| JsError::new("Invalid hasher handle"))?;

    let collection = hasher
        .finish()
        .map_err(|e| JsError::new(&format!("Failed to finalize hasher: {}", e)))?;

    let result = SeqColResult::from_collection(&collection);
    serde_wasm_bindgen::to_value(&result)
        .map_err(|e| JsError::new(&format!("Serialization error: {}", e)))
}

/// Free a hasher without getting results.
///
/// Use this if you need to cancel processing early.
/// After calling this, the handle is no longer valid.
///
/// # Arguments
/// * `handle` - The hasher handle from `fasta_hasher_new()`
///
/// # Returns
/// `true` if the handle was valid and freed, `false` otherwise.
#[wasm_bindgen]
pub fn fasta_hasher_free(handle: u32) -> bool {
    with_storage(|storage| storage.remove(handle).is_some())
}

/// Get the current progress of a streaming hasher.
///
/// # Arguments
/// * `handle` - The hasher handle from `fasta_hasher_new()`
///
/// # Returns
/// A JavaScript object with progress information, or null if handle is invalid.
/// The object contains:
/// - `completed_sequences`: Number of fully processed sequences
/// - `current_sequence_name`: Name of sequence being processed (if any)
/// - `current_sequence_length`: Length of sequence being processed
#[wasm_bindgen]
pub fn fasta_hasher_progress(handle: u32) -> Result<JsValue, JsError> {
    let progress = with_storage(|storage| {
        storage.get_mut(handle).map(|hasher| HasherProgress {
            completed_sequences: hasher.sequence_count(),
            current_sequence_name: hasher.current_sequence_name().map(String::from),
            current_sequence_length: hasher.current_sequence_length(),
        })
    });

    match progress {
        Some(p) => serde_wasm_bindgen::to_value(&p)
            .map_err(|e| JsError::new(&format!("Serialization error: {}", e))),
        None => Err(JsError::new("Invalid hasher handle")),
    }
}

// ============================================================================
// Helper types for JavaScript interop
// ============================================================================

/// JavaScript-friendly representation of a sequence collection result
#[derive(Serialize, Deserialize)]
pub struct SeqColResult {
    /// Top-level seqcol digest
    pub digest: String,
    /// Number of sequences
    pub n_sequences: usize,
    /// Level 1 digest of names array
    pub names_digest: String,
    /// Level 1 digest of sequences array
    pub sequences_digest: String,
    /// Level 1 digest of lengths array
    pub lengths_digest: String,
    /// Sequence metadata
    pub sequences: Vec<SeqMetadata>,
}

/// JavaScript-friendly representation of sequence metadata
#[derive(Serialize, Deserialize)]
pub struct SeqMetadata {
    /// Sequence name (e.g., "chr1")
    pub name: String,
    /// Sequence length in bases
    pub length: usize,
    /// sha512t24u digest of the sequence
    pub sha512t24u: String,
    /// MD5 digest of the sequence
    pub md5: String,
    /// Guessed alphabet type (e.g., "dna2bit")
    pub alphabet: String,
    /// Optional description from FASTA header
    pub description: Option<String>,
}

/// Progress information for streaming hasher
#[derive(Serialize, Deserialize)]
pub struct HasherProgress {
    /// Number of fully processed sequences
    pub completed_sequences: usize,
    /// Name of sequence currently being processed (if any)
    pub current_sequence_name: Option<String>,
    /// Current length of sequence being processed
    pub current_sequence_length: usize,
}

impl SeqColResult {
    fn from_collection(collection: &SequenceCollection) -> Self {
        let sequences: Vec<SeqMetadata> = collection
            .sequences
            .iter()
            .map(|r| {
                let meta = r.metadata();
                SeqMetadata {
                    name: meta.name.clone(),
                    length: meta.length,
                    sha512t24u: meta.sha512t24u.clone(),
                    md5: meta.md5.clone(),
                    alphabet: format!("{}", meta.alphabet),
                    description: meta.description.clone(),
                }
            })
            .collect();

        SeqColResult {
            digest: collection.metadata.digest.clone(),
            n_sequences: collection.metadata.n_sequences,
            names_digest: collection.metadata.names_digest.clone(),
            sequences_digest: collection.metadata.sequences_digest.clone(),
            lengths_digest: collection.metadata.lengths_digest.clone(),
            sequences,
        }
    }
}
