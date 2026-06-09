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
//! const handle = fastaHasherNew();
//! try {
//!     const response = await fetch(url);
//!     const reader = response.body.getReader();
//!     while (true) {
//!         const { done, value } = await reader.read();
//!         if (done) break;
//!         fastaHasherUpdate(handle, value);  // throws on error
//!     }
//!     const result = fastaHasherFinish(handle);
//! } catch (err) {
//!     fastaHasherFree(handle);  // cleanup on error
//!     throw err;
//! }
//! ```

use gtars_refget::digest::{
    canonicalize_json, digest_fasta_bytes, digest_sequence, lookup_alphabet, md5, sha512t24u,
    AlphabetType, FastaStreamHasher, SequenceCollection, SequenceMetadata, SequenceRecord,
};
use gtars_refget::store::{ReadonlyRefgetStore, RefgetStore as CoreRefgetStore};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use wasm_bindgen::prelude::*;

// ============================================================================
// Global storage for streaming hasher instances
// ============================================================================

/// Global storage for hasher instances, protected by a mutex.
/// The Mutex ensures only one caller can access the storage at a time,
/// preventing data races (no two mutable references can exist simultaneously).
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
        // Find an unused ID (handles wrap-around after ~4 billion allocations)
        let mut id = self.next_id;
        while self.hashers.contains_key(&id) || id == 0 {
            id = id.wrapping_add(1);
            if id == 0 {
                id = 1; // Skip 0
            }
        }
        self.next_id = id.wrapping_add(1);
        if self.next_id == 0 {
            self.next_id = 1;
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
    let mut guard = HASHER_STORAGE
        .lock()
        .expect("HASHER_STORAGE mutex poisoned");
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
/// const result = digestSeqcol(fastaContent);
/// console.log(result.digest);  // Top-level seqcol digest
/// ```
#[wasm_bindgen(js_name = "digestSeqcol")]
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
/// const digest = sequenceDigest("ACGT");
/// console.log(digest);  // "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
/// ```
#[wasm_bindgen(js_name = "sequenceDigest")]
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
#[wasm_bindgen(js_name = "sequenceMd5")]
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
#[wasm_bindgen(js_name = "computeSha512t24u")]
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
#[wasm_bindgen(js_name = "computeMd5")]
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
#[wasm_bindgen(js_name = "canonicalizeJsonString")]
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
/// The handle must be freed with `fastaHasherFinish()` or `fastaHasherFree()`.
///
/// # Returns
/// A handle ID (always > 0).
///
/// # Example (JavaScript)
/// ```javascript
/// const handle = fastaHasherNew();
/// // ... use handle with fastaHasherUpdate ...
/// const result = fastaHasherFinish(handle);
/// ```
#[wasm_bindgen(js_name = "fastaHasherNew")]
pub fn fasta_hasher_new() -> u32 {
    with_storage(|storage| storage.insert(FastaStreamHasher::new()))
}

/// Process a chunk of FASTA data.
///
/// This can be called multiple times with successive chunks.
/// Handles both plain text and gzip-compressed FASTA.
///
/// # Arguments
/// * `handle` - The hasher handle from `fastaHasherNew()`
/// * `chunk` - A chunk of FASTA data as Uint8Array
///
/// # Returns
/// Ok on success, or an error with details about what went wrong.
///
/// # Example (JavaScript)
/// ```javascript
/// const handle = fastaHasherNew();
/// fastaHasherUpdate(handle, new TextEncoder().encode(">chr1\nACGT"));
/// fastaHasherUpdate(handle, new TextEncoder().encode("\n>chr2\nTGCA\n"));
/// ```
#[wasm_bindgen(js_name = "fastaHasherUpdate")]
pub fn fasta_hasher_update(handle: u32, chunk: &[u8]) -> Result<(), JsError> {
    with_storage(|storage| {
        if let Some(hasher) = storage.get_mut(handle) {
            hasher
                .update(chunk)
                .map_err(|e| JsError::new(&format!("Failed to process chunk: {}", e)))
        } else {
            Err(JsError::new("Invalid hasher handle"))
        }
    })
}

/// Finalize the hasher and return the results.
///
/// This consumes the hasher and frees its resources.
/// After calling this, the handle is no longer valid.
///
/// # Arguments
/// * `handle` - The hasher handle from `fastaHasherNew()`
///
/// # Returns
/// A JavaScript object with the sequence collection result,
/// or an error if the handle is invalid.
///
/// # Example (JavaScript)
/// ```javascript
/// const handle = fastaHasherNew();
/// fastaHasherUpdate(handle, fastaData);
/// const result = fastaHasherFinish(handle);
/// console.log(result.digest);
/// ```
#[wasm_bindgen(js_name = "fastaHasherFinish")]
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
/// * `handle` - The hasher handle from `fastaHasherNew()`
///
/// # Returns
/// `true` if the handle was valid and freed, `false` otherwise.
#[wasm_bindgen(js_name = "fastaHasherFree")]
pub fn fasta_hasher_free(handle: u32) -> bool {
    with_storage(|storage| storage.remove(handle).is_some())
}

/// Get the current progress of a streaming hasher.
///
/// # Arguments
/// * `handle` - The hasher handle from `fastaHasherNew()`
///
/// # Returns
/// A JavaScript object with progress information, or null if handle is invalid.
/// The object contains:
/// - `completed_sequences`: Number of fully processed sequences
/// - `current_sequence_name`: Name of sequence being processed (if any)
/// - `current_sequence_length`: Length of sequence being processed
#[wasm_bindgen(js_name = "fastaHasherProgress")]
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

// ============================================================================
// Persistent, reusable refget store (genome-scale, encoded)
// ============================================================================

/// A persistent, reusable refget store for the browser.
///
/// Unlike the one-shot `hgvs_to_vrs_id` (which rebuilds a fresh single-sequence
/// store on every call), this holds the reference genome **once** so a whole VCF
/// of variants can reuse it. The genome is held in an encoded in-memory store and
/// decoded on the fly per variant via the VRS `RefView` path; the decoded genome
/// (~3 GB for a human reference) is never materialized.
///
/// Two ingest paths:
/// - [`add_sequence`](RefgetStore::add_sequence): hand in raw bases; the store
///   digests them and stores them as a resident record (good for small contigs
///   or tests).
/// - [`add_encoded_sequence`](RefgetStore::add_encoded_sequence): hand in bytes
///   that are **already** 2-/3-bit encoded (the genome blob downloaded once into
///   OPFS), skipping re-encoding of the ~750 MB reference.
///
/// Chromosome-name aliasing: VCFs reference sequences as `chr1`, `1`, or an
/// accession like `NC_000001.11`. Each ingest registers the given name plus its
/// `chr`-prefix toggle against the same digest; use [`add_alias`](RefgetStore::add_alias)
/// to map an accession (or any other name) the VCF might use to a known digest.
///
/// # Example (JavaScript)
/// ```javascript
/// import init, { RefgetStore } from "gtars-js";
/// await init();
/// const store = new RefgetStore();
/// const digest = store.add_sequence("chr1", new TextEncoder().encode("ACGT..."));
/// store.add_alias("NC_000001.11", digest);
/// // ... hand `store` to vcf_to_vrs_ids(store, vcfText, onResult) ...
/// ```
///
/// The same resident store serves both `g.` (via
/// [`hgvs_to_vrs_id`](RefgetStore::hgvs_to_vrs_id)) and `c.`/`n.` (via
/// [`hgvs_to_vrs_id_with_transcripts`](RefgetStore::hgvs_to_vrs_id_with_transcripts),
/// given a [`TranscriptStore`](crate::transcripts::TranscriptStore)) — exactly
/// one genome holder backs every reference type.
#[wasm_bindgen]
pub struct RefgetStore {
    inner: CoreRefgetStore,
    /// Maps every name a VCF might use (canonical, `chr`-toggled, explicit
    /// aliases) to the raw sha512t24u digest of the underlying sequence.
    name_to_digest: HashMap<String, String>,
}

#[wasm_bindgen]
impl RefgetStore {
    /// Create a new, empty encoded in-memory store.
    #[wasm_bindgen(constructor)]
    pub fn new() -> RefgetStore {
        RefgetStore {
            inner: CoreRefgetStore::in_memory(),
            name_to_digest: HashMap::new(),
        }
    }

    /// Ingest a sequence from **raw bases**. Computes the GA4GH digests and
    /// stores the record resident; the VRS path decodes on the fly. Returns the
    /// raw sha512t24u digest, which is also the VRS sequence accession.
    ///
    /// # Arguments
    /// * `name` - The sequence name (e.g. `"chr1"`).
    /// * `bases` - The reference bases as a Uint8Array.
    #[wasm_bindgen]
    pub fn add_sequence(&mut self, name: &str, bases: &[u8]) -> Result<String, JsError> {
        let record = digest_sequence(name, bases);
        let digest = record.metadata().sha512t24u.clone();
        self.inner
            .add_sequence_record(record, true)
            .map_err(|e| JsError::new(&format!("failed to add sequence '{name}': {e}")))?;
        self.register_name(name, &digest);
        Ok(digest)
    }

    /// Ingest a sequence whose bytes are **already encoded** (2-/3-bit packed),
    /// as shipped in a prebuilt genome blob. Skips re-encoding the reference.
    ///
    /// The GA4GH digest is computed over the *raw* bases at genome-build time and
    /// cannot be recovered from the packed bytes, so it (and the alphabet) must be
    /// supplied from the blob's seqcol metadata.
    ///
    /// # Arguments
    /// * `name` - The sequence name (e.g. `"chr1"`).
    /// * `sha512t24u` - The raw-sequence GA4GH digest (the VRS accession).
    /// * `md5` - The raw-sequence MD5 (may be empty; used only for md5 lookups).
    /// * `length` - The sequence length in **bases** (not packed bytes).
    /// * `encoded` - The packed bytes; must be exactly `ceil(length * bits/8)`.
    /// * `alphabet` - One of `"dna2bit"`, `"dna3bit"`, `"dnaio"` (IUPAC).
    #[wasm_bindgen]
    pub fn add_encoded_sequence(
        &mut self,
        name: &str,
        sha512t24u: &str,
        md5: &str,
        length: usize,
        encoded: &[u8],
        alphabet: &str,
    ) -> Result<(), JsError> {
        let alphabet_type = parse_alphabet(alphabet)
            .ok_or_else(|| JsError::new(&format!("unsupported alphabet: {alphabet}")))?;
        let bps = lookup_alphabet(&alphabet_type).bits_per_symbol;
        let expected = length.saturating_mul(bps).div_ceil(8);
        if encoded.len() != expected {
            return Err(JsError::new(&format!(
                "encoded length {} does not match expected {} for length={} alphabet={} ({} bits/symbol)",
                encoded.len(),
                expected,
                length,
                alphabet,
                bps
            )));
        }
        let metadata = SequenceMetadata {
            name: name.to_string(),
            description: None,
            length,
            sha512t24u: sha512t24u.to_string(),
            md5: md5.to_string(),
            alphabet: alphabet_type,
            fai: None,
        };
        let record = SequenceRecord::Full {
            metadata,
            sequence: Arc::new(encoded.to_vec()),
        };
        self.inner
            .add_sequence_record(record, true)
            .map_err(|e| JsError::new(&format!("failed to add encoded sequence '{name}': {e}")))?;
        self.register_name(name, sha512t24u);
        Ok(())
    }

    /// Map an additional `name` (e.g. an accession like `"NC_000001.11"`) the VCF
    /// might use to an already-ingested sequence's `digest`.
    #[wasm_bindgen]
    pub fn add_alias(&mut self, name: &str, digest: &str) {
        self.name_to_digest.insert(name.to_string(), digest.to_string());
    }

    /// Number of distinct name→digest entries registered (including aliases).
    #[wasm_bindgen(js_name = "nameCount")]
    pub fn name_count(&self) -> usize {
        self.name_to_digest.len()
    }

    /// Convert a single genomic (`g.`) HGVS string into a GA4GH VRS allele id,
    /// resolving the reference against THIS resident store (reusing the loaded
    /// genome — unlike the standalone `hgvs_to_vrs_id`, which builds a one-off
    /// store per call). The referenced sequence must already be ingested
    /// (`add_sequence` / `add_encoded_sequence`). Uses a `NoTranscriptProvider`,
    /// so only `g.` resolves; for `c.`/`n.` use
    /// [`hgvs_to_vrs_id_with_transcripts`](RefgetStore::hgvs_to_vrs_id_with_transcripts).
    #[wasm_bindgen]
    pub fn hgvs_to_vrs_id(&self, hgvs: &str) -> Result<String, JsValue> {
        gtars_vrs::hgvs::bridge::hgvs_str_to_vrs_id_readonly(
            hgvs,
            &gtars_vrs::provider::NoTranscriptProvider,
            self.inner_readonly(),
            self.name_to_digest(),
        )
        .map(|b| b.value)
        .map_err(|e| JsValue::from_str(&format!("{e}")))
    }

    /// Convert an HGVS string into a GA4GH VRS allele id, resolving
    /// transcript-relative (`c.`/`n.`) coordinates and gene-symbol → MANE
    /// lookups against `transcripts`, while reading reference bases from THIS
    /// resident genome.
    ///
    /// This is the genome-scale path: a caller holding a resident genome (every
    /// contig ingested) AND a [`TranscriptStore`](crate::transcripts::TranscriptStore)
    /// can resolve `c.`/`n.` without rebuilding either. The chromosome the
    /// transcript maps onto must already be resident here (`add_sequence` /
    /// `add_encoded_sequence`); otherwise the variant cannot resolve.
    #[wasm_bindgen]
    pub fn hgvs_to_vrs_id_with_transcripts(
        &self,
        hgvs: &str,
        transcripts: &crate::transcripts::TranscriptStore,
    ) -> Result<String, JsValue> {
        gtars_vrs::hgvs::bridge::hgvs_str_to_vrs_id_readonly(
            hgvs,
            transcripts.provider(),
            self.inner_readonly(),
            self.name_to_digest(),
        )
        .map(|b| b.value)
        .map_err(|e| JsValue::from_str(&format!("{e}")))
    }
}

impl RefgetStore {
    /// Register a name plus its `chr`-prefix toggle against `digest`, so VCFs
    /// that drop or add the `chr` prefix still resolve. Does not overwrite an
    /// explicit alias already pointing elsewhere for the toggled form.
    fn register_name(&mut self, name: &str, digest: &str) {
        self.name_to_digest.insert(name.to_string(), digest.to_string());
        let toggled = match name.strip_prefix("chr") {
            Some(rest) => rest.to_string(),
            None => format!("chr{name}"),
        };
        self.name_to_digest
            .entry(toggled)
            .or_insert_with(|| digest.to_string());
    }

    /// The underlying read-only store, for the VCF→VRS streaming core.
    pub(crate) fn inner_readonly(&self) -> &ReadonlyRefgetStore {
        &self.inner
    }

    /// The name→digest map built from ingested sequences and aliases.
    pub(crate) fn name_to_digest(&self) -> &HashMap<String, String> {
        &self.name_to_digest
    }
}

impl Default for RefgetStore {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse a JS-supplied alphabet name into an [`AlphabetType`]. Accepts the
/// `Display` spellings used elsewhere in gtars (`"dna2bit"`, `"dna3bit"`,
/// `"dnaio"`), plus the friendlier `"dnaiupac"`.
fn parse_alphabet(alphabet: &str) -> Option<AlphabetType> {
    match alphabet.to_ascii_lowercase().as_str() {
        "dna2bit" => Some(AlphabetType::Dna2bit),
        "dna3bit" => Some(AlphabetType::Dna3bit),
        "dnaio" | "dnaiupac" => Some(AlphabetType::DnaIupac),
        _ => None,
    }
}

// Tests use wasm-bindgen-test since JsError/JsValue require WASM runtime
#[cfg(test)]
#[cfg(target_arch = "wasm32")]
mod tests {
    use super::*;
    use wasm_bindgen_test::*;

    // Persistent-store tests share the hgvs.rs synthetic contig and golden id, so
    // the store class is proven to drive the SAME readonly VRS path the one-shot
    // `hgvs_to_vrs_id` does. If the store ingest ever corrupts bytes/metadata,
    // these diverge from the golden id.
    const CHR_F_NAME: &str = "chrF";
    const CHR_F_BASES: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    const GOLDEN_CHRF_G6CT_VRS_ID: &str = "ga4gh:VA._q-idtHGQxQ4XiEPJ1ExYl_htUeNEkir";

    fn vrs_id_via_store(store: &RefgetStore) -> String {
        gtars_vrs::hgvs::bridge::hgvs_str_to_vrs_id_readonly(
            "chrF:g.6C>T",
            &gtars_vrs::provider::NoTranscriptProvider,
            store.inner_readonly(),
            store.name_to_digest(),
        )
        .expect("readonly bridge should succeed")
        .value
    }

    #[wasm_bindgen_test]
    fn store_add_sequence_drives_readonly_vrs() {
        let mut store = RefgetStore::new();
        let digest = store
            .add_sequence(CHR_F_NAME, CHR_F_BASES.as_bytes())
            .expect("add_sequence should succeed");
        // Raw-bases ingest resolves through the same path as the one-shot entry.
        assert_eq!(vrs_id_via_store(&store), GOLDEN_CHRF_G6CT_VRS_ID);
        // chr-toggle alias was registered ("chrF" -> also "F").
        assert_eq!(store.name_to_digest().get("F"), Some(&digest));
    }

    #[wasm_bindgen_test]
    fn store_add_encoded_sequence_drives_readonly_vrs() {
        // Encode the contig exactly as the store would, then ingest the packed
        // bytes via the encoded path; the on-the-fly decode must reproduce the id.
        let meta = digest_sequence(CHR_F_NAME, CHR_F_BASES.as_bytes());
        let digest = meta.metadata().sha512t24u.clone();
        let md5 = meta.metadata().md5.clone();
        let encoded = gtars_refget::digest::encode_sequence(
            CHR_F_BASES.as_bytes(),
            &gtars_refget::digest::DNA_2BIT_ALPHABET,
        );

        let mut store = RefgetStore::new();
        store
            .add_encoded_sequence(
                CHR_F_NAME,
                &digest,
                &md5,
                CHR_F_BASES.len(),
                &encoded,
                "dna2bit",
            )
            .expect("add_encoded_sequence should succeed");
        assert_eq!(vrs_id_via_store(&store), GOLDEN_CHRF_G6CT_VRS_ID);
    }

    #[wasm_bindgen_test]
    fn store_hgvs_to_vrs_id_uses_resident_genome() {
        let mut store = RefgetStore::new();
        store
            .add_sequence(CHR_F_NAME, CHR_F_BASES.as_bytes())
            .expect("add_sequence");
        // The store method resolves against the resident genome (no per-call
        // store rebuild) and must match the same golden id.
        let id = store.hgvs_to_vrs_id("chrF:g.6C>T").expect("hgvs_to_vrs_id");
        assert_eq!(id, GOLDEN_CHRF_G6CT_VRS_ID);
    }

    #[wasm_bindgen_test]
    fn store_add_encoded_sequence_rejects_bad_length() {
        let mut store = RefgetStore::new();
        // 40 bases at 2 bits/symbol needs 10 packed bytes; give it 3.
        let res = store.add_encoded_sequence("chrF", "deadbeef", "", 40, &[0u8; 3], "dna2bit");
        assert!(res.is_err(), "mismatched encoded length should be rejected");
    }

    #[wasm_bindgen_test]
    fn test_streaming_lifecycle() {
        // Basic new -> update -> finish cycle
        let handle = fasta_hasher_new();
        assert!(handle > 0);

        // Update with FASTA data
        let fasta = b">chr1\nACGT\n>chr2\nTGCA\n";
        fasta_hasher_update(handle, fasta).expect("update should succeed");

        // Finish and verify result
        let result = fasta_hasher_finish(handle).expect("finish should succeed");
        let seqcol: SeqColResult = serde_wasm_bindgen::from_value(result).expect("deserialize");
        assert_eq!(seqcol.n_sequences, 2);
        assert_eq!(seqcol.sequences[0].name, "chr1");
        assert_eq!(seqcol.sequences[1].name, "chr2");
    }

    #[wasm_bindgen_test]
    fn test_invalid_handle_error() {
        // Using handle 9999 which doesn't exist
        let result = fasta_hasher_update(9999, b"data");
        assert!(result.is_err());

        let result = fasta_hasher_finish(9999);
        assert!(result.is_err());
    }
}
