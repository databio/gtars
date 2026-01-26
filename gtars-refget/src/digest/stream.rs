//! Streaming FASTA hasher - WASM-safe, constant memory usage.
//!
//! This module provides a streaming FASTA parser that can process data in chunks,
//! maintaining constant memory usage regardless of file size. Ideal for WASM
//! environments where large files are fetched in chunks.

use md5::Md5;
use sha2::{Digest, Sha512};
use std::io::{self, Write};

use super::alphabet::AlphabetGuesser;
use super::fasta::parse_fasta_header;
use super::types::{SequenceCollection, SequenceCollectionMetadata, SequenceMetadata, SequenceRecord};

/// Streaming state for FASTA processing.
#[derive(Clone, Copy, Debug, PartialEq)]
enum ParserState {
    /// Waiting for a header line starting with '>'
    AwaitingHeader,
    /// Currently in sequence data lines
    InSequence,
}

/// Inner FASTA processor that implements Write.
/// Decompressed data flows into this via the Write trait.
struct FastaProcessor {
    state: ParserState,
    line_buffer: Vec<u8>,
    current_name: Option<String>,
    current_description: Option<String>,
    current_length: usize,
    sha512_hasher: Sha512,
    md5_hasher: Md5,
    alphabet_guesser: AlphabetGuesser,
    sequences: Vec<SequenceRecord>,
    /// Store any processing error (Write trait can't return custom errors)
    processing_error: Option<String>,
}

impl FastaProcessor {
    fn new() -> Self {
        Self {
            state: ParserState::AwaitingHeader,
            line_buffer: Vec::with_capacity(8192),
            current_name: None,
            current_description: None,
            current_length: 0,
            sha512_hasher: Sha512::new(),
            md5_hasher: Md5::new(),
            alphabet_guesser: AlphabetGuesser::new(),
            sequences: Vec::new(),
            processing_error: None,
        }
    }

    fn process_byte(&mut self, byte: u8) {
        // Stop processing if we already have an error
        if self.processing_error.is_some() {
            return;
        }

        if byte == b'\n' || byte == b'\r' {
            if let Err(e) = self.process_line() {
                self.processing_error = Some(e.to_string());
            }
            self.line_buffer.clear();
        } else {
            self.line_buffer.push(byte);
        }
    }

    fn process_line(&mut self) -> anyhow::Result<()> {
        if self.line_buffer.is_empty() {
            return Ok(());
        }

        if self.line_buffer[0] == b'>' {
            // Header line - finalize previous sequence if any
            if self.current_name.is_some() {
                self.finalize_current_sequence();
            }

            // Parse header
            let header = std::str::from_utf8(&self.line_buffer[1..])
                .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in header: {}", e))?;
            let (name, description) = parse_fasta_header(header);

            // Start new sequence
            self.current_name = Some(name);
            self.current_description = description;
            self.current_length = 0;
            self.sha512_hasher = Sha512::new();
            self.md5_hasher = Md5::new();
            self.alphabet_guesser = AlphabetGuesser::new();
            self.state = ParserState::InSequence;
        } else if self.state == ParserState::InSequence && self.current_name.is_some() {
            // Sequence line - uppercase and hash
            let uppercased: Vec<u8> = self.line_buffer.iter()
                .filter(|&&b| !b.is_ascii_whitespace())
                .map(|b| b.to_ascii_uppercase())
                .collect();

            if !uppercased.is_empty() {
                self.sha512_hasher.update(&uppercased);
                self.md5_hasher.update(&uppercased);
                self.alphabet_guesser.update(&uppercased);
                self.current_length += uppercased.len();
            }
        }

        Ok(())
    }

    fn finalize_current_sequence(&mut self) {
        if let Some(name) = self.current_name.take() {
            let sha512 = base64_url::encode(&self.sha512_hasher.clone().finalize()[0..24]);
            let md5 = format!("{:x}", self.md5_hasher.clone().finalize());
            let alphabet = self.alphabet_guesser.guess();

            let metadata = SequenceMetadata {
                name,
                description: self.current_description.take(),
                length: self.current_length,
                sha512t24u: sha512,
                md5,
                alphabet,
                fai: None,
            };

            self.sequences.push(SequenceRecord::Stub(metadata));
        }
    }

    fn finish(mut self) -> anyhow::Result<SequenceCollection> {
        // Check for any stored processing error
        if let Some(err) = self.processing_error {
            return Err(anyhow::anyhow!("Processing error: {}", err));
        }

        // Process any remaining data in the line buffer
        if !self.line_buffer.is_empty() {
            self.process_line()?;
        }

        // Finalize the last sequence
        if self.current_name.is_some() {
            self.finalize_current_sequence();
        }

        // Build the collection
        let metadata = SequenceCollectionMetadata::from_sequences(&self.sequences, None);

        Ok(SequenceCollection {
            metadata,
            sequences: self.sequences,
        })
    }
}

impl Write for FastaProcessor {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        for &byte in buf {
            self.process_byte(byte);
            // Check for processing errors after each byte
            if let Some(ref err) = self.processing_error {
                return Err(io::Error::new(io::ErrorKind::InvalidData, err.clone()));
            }
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// State for format detection and processing
enum ProcessorState {
    /// Haven't detected format yet
    Detecting,
    /// Plain FASTA (uncompressed)
    Plain(FastaProcessor),
    /// Gzipped FASTA - uses write::GzDecoder for true streaming decompression
    Gzipped(flate2::write::GzDecoder<FastaProcessor>),
}

/// A streaming FASTA hasher that processes data chunk-by-chunk.
///
/// This is designed for WASM environments where files are fetched in chunks.
/// Memory usage is constant (~100KB) regardless of file size:
/// - Internal state: ~200 bytes (hasher state, counters)
/// - Line buffer: ~8KB (handles long lines)
/// - Gzip decoder state: ~32KB if compressed
/// - Results: grows only with number of sequences (not sequence length)
///
/// # Example
/// ```
/// use gtars_refget::digest::stream::FastaStreamHasher;
///
/// let mut hasher = FastaStreamHasher::new();
///
/// // Process first chunk
/// hasher.update(b">chr1\nACGT").expect("update");
///
/// // Process second chunk
/// hasher.update(b"TGCA\n>chr2\nGGGG\n").expect("update");
///
/// // Finalize and get results
/// let collection = hasher.finish().expect("finish");
/// assert_eq!(collection.sequences.len(), 2);
/// ```
pub struct FastaStreamHasher {
    state: ProcessorState,
}

impl FastaStreamHasher {
    /// Create a new streaming FASTA hasher.
    pub fn new() -> Self {
        Self {
            state: ProcessorState::Detecting,
        }
    }

    /// Process a chunk of FASTA data.
    ///
    /// This method can be called multiple times with successive chunks of data.
    /// Handles both plain text and gzip-compressed FASTA with true streaming
    /// decompression (constant memory usage).
    ///
    /// # Arguments
    /// * `chunk` - A slice of bytes from the FASTA file
    ///
    /// # Returns
    /// Ok(()) on success, Err on parsing error
    pub fn update(&mut self, chunk: &[u8]) -> anyhow::Result<()> {
        if chunk.is_empty() {
            return Ok(());
        }

        // Detect format on first non-empty chunk
        if matches!(self.state, ProcessorState::Detecting) {
            let is_gzipped = chunk.len() >= 2 && chunk[0] == 0x1f && chunk[1] == 0x8b;

            if is_gzipped {
                // Create GzDecoder wrapping a FastaProcessor
                // Decompressed data flows directly into the processor
                let processor = FastaProcessor::new();
                let decoder = flate2::write::GzDecoder::new(processor);
                self.state = ProcessorState::Gzipped(decoder);
            } else {
                self.state = ProcessorState::Plain(FastaProcessor::new());
            }
        }

        // Process the chunk
        match &mut self.state {
            ProcessorState::Detecting => unreachable!(),
            ProcessorState::Plain(processor) => {
                processor.write_all(chunk)?;
            }
            ProcessorState::Gzipped(decoder) => {
                // Write compressed data to decoder
                // Decoder decompresses and writes to inner FastaProcessor
                decoder.write_all(chunk)?;
            }
        }

        Ok(())
    }

    /// Finalize processing and return the SequenceCollection.
    ///
    /// This must be called after all chunks have been processed via `update()`.
    pub fn finish(self) -> anyhow::Result<SequenceCollection> {
        match self.state {
            ProcessorState::Detecting => {
                // No data was ever provided - return empty collection
                let metadata = SequenceCollectionMetadata::from_sequences(&[], None);
                Ok(SequenceCollection {
                    metadata,
                    sequences: Vec::new(),
                })
            }
            ProcessorState::Plain(processor) => {
                processor.finish()
            }
            ProcessorState::Gzipped(decoder) => {
                // Finish decompression and get the inner processor
                let processor = decoder.finish()
                    .map_err(|e| anyhow::anyhow!("Gzip decompression error: {}", e))?;
                processor.finish()
            }
        }
    }

    /// Get the current number of completed sequences.
    pub fn sequence_count(&self) -> usize {
        match &self.state {
            ProcessorState::Detecting => 0,
            ProcessorState::Plain(p) => p.sequences.len(),
            ProcessorState::Gzipped(d) => d.get_ref().sequences.len(),
        }
    }

    /// Check if currently processing a sequence.
    pub fn in_sequence(&self) -> bool {
        match &self.state {
            ProcessorState::Detecting => false,
            ProcessorState::Plain(p) => p.current_name.is_some(),
            ProcessorState::Gzipped(d) => d.get_ref().current_name.is_some(),
        }
    }

    /// Get the name of the sequence currently being processed (if any).
    pub fn current_sequence_name(&self) -> Option<&str> {
        match &self.state {
            ProcessorState::Detecting => None,
            ProcessorState::Plain(p) => p.current_name.as_deref(),
            ProcessorState::Gzipped(d) => d.get_ref().current_name.as_deref(),
        }
    }

    /// Get the current length of the sequence being processed.
    pub fn current_sequence_length(&self) -> usize {
        match &self.state {
            ProcessorState::Detecting => 0,
            ProcessorState::Plain(p) => p.current_length,
            ProcessorState::Gzipped(d) => d.get_ref().current_length,
        }
    }
}

impl Default for FastaStreamHasher {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::digest::alphabet::AlphabetType;
    use crate::digest::fasta::digest_fasta_bytes;

    #[test]
    fn test_streaming_basic() {
        let mut hasher = FastaStreamHasher::new();
        hasher.update(b">chr1\nACGT\n>chr2\nTGCA\n").expect("update");
        let collection = hasher.finish().expect("finish");

        assert_eq!(collection.sequences.len(), 2);
        assert_eq!(collection.sequences[0].metadata().name, "chr1");
        assert_eq!(collection.sequences[0].metadata().length, 4);
        assert_eq!(collection.sequences[1].metadata().name, "chr2");
        assert_eq!(collection.sequences[1].metadata().length, 4);
    }

    #[test]
    fn test_streaming_chunked() {
        let mut hasher = FastaStreamHasher::new();

        // Split data across multiple chunks
        hasher.update(b">chr1\nAC").expect("chunk 1");
        hasher.update(b"GT\n>chr2\n").expect("chunk 2");
        hasher.update(b"TGCA\n").expect("chunk 3");

        let collection = hasher.finish().expect("finish");

        assert_eq!(collection.sequences.len(), 2);
        assert_eq!(collection.sequences[0].metadata().name, "chr1");
        assert_eq!(collection.sequences[0].metadata().length, 4);
        assert_eq!(collection.sequences[1].metadata().name, "chr2");
        assert_eq!(collection.sequences[1].metadata().length, 4);
    }

    #[test]
    fn test_streaming_split_header() {
        let mut hasher = FastaStreamHasher::new();

        // Header split across chunks
        hasher.update(b">ch").expect("chunk 1");
        hasher.update(b"r1 description\nACGT\n").expect("chunk 2");

        let collection = hasher.finish().expect("finish");

        assert_eq!(collection.sequences.len(), 1);
        assert_eq!(collection.sequences[0].metadata().name, "chr1");
        assert_eq!(
            collection.sequences[0].metadata().description,
            Some("description".to_string())
        );
    }

    #[test]
    fn test_streaming_matches_batch() {
        // Test that streaming produces same results as batch processing
        let fasta_data = b">chrX\nTTGGGGAA\n>chr1\nGGAA\n>chr2\nGCGC\n";

        // Batch processing
        let batch_result = digest_fasta_bytes(fasta_data).expect("batch");

        // Streaming processing (single chunk)
        let mut hasher = FastaStreamHasher::new();
        hasher.update(fasta_data).expect("streaming");
        let stream_result = hasher.finish().expect("finish");

        // Compare results
        assert_eq!(batch_result.metadata.digest, stream_result.metadata.digest);
        assert_eq!(
            batch_result.metadata.names_digest,
            stream_result.metadata.names_digest
        );
        assert_eq!(
            batch_result.metadata.sequences_digest,
            stream_result.metadata.sequences_digest
        );
        assert_eq!(
            batch_result.metadata.lengths_digest,
            stream_result.metadata.lengths_digest
        );

        for (batch_seq, stream_seq) in batch_result
            .sequences
            .iter()
            .zip(stream_result.sequences.iter())
        {
            assert_eq!(batch_seq.metadata().name, stream_seq.metadata().name);
            assert_eq!(batch_seq.metadata().length, stream_seq.metadata().length);
            assert_eq!(
                batch_seq.metadata().sha512t24u,
                stream_seq.metadata().sha512t24u
            );
            assert_eq!(batch_seq.metadata().md5, stream_seq.metadata().md5);
            assert_eq!(batch_seq.metadata().alphabet, stream_seq.metadata().alphabet);
        }
    }

    #[test]
    fn test_streaming_multiline_sequence() {
        let mut hasher = FastaStreamHasher::new();
        hasher
            .update(b">chr1\nACGT\nTGCA\nAAAA\n")
            .expect("update");
        let collection = hasher.finish().expect("finish");

        assert_eq!(collection.sequences.len(), 1);
        assert_eq!(collection.sequences[0].metadata().length, 12);
    }

    #[test]
    fn test_streaming_empty() {
        let hasher = FastaStreamHasher::new();
        let collection = hasher.finish().expect("finish");
        assert_eq!(collection.sequences.len(), 0);
    }

    #[test]
    fn test_streaming_known_digest() {
        let mut hasher = FastaStreamHasher::new();
        hasher.update(b">chrX\nTTGGGGAA\n").expect("update");
        let collection = hasher.finish().expect("finish");

        assert_eq!(
            collection.sequences[0].metadata().sha512t24u,
            "iYtREV555dUFKg2_agSJW6suquUyPpMw"
        );
        assert_eq!(
            collection.sequences[0].metadata().md5,
            "5f63cfaa3ef61f88c9635fb9d18ec945"
        );
        assert_eq!(
            collection.sequences[0].metadata().alphabet,
            AlphabetType::Dna2bit
        );
    }

    #[test]
    fn test_streaming_gzipped() {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let fasta = b">chr1\nACGT\n";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(fasta).expect("compress");
        let compressed = encoder.finish().expect("finish compression");

        let mut hasher = FastaStreamHasher::new();
        hasher.update(&compressed).expect("update");
        let collection = hasher.finish().expect("finish");

        assert_eq!(collection.sequences.len(), 1);
        assert_eq!(collection.sequences[0].metadata().name, "chr1");
        assert_eq!(collection.sequences[0].metadata().length, 4);
    }

    #[test]
    fn test_streaming_gzipped_chunked() {
        // Test that gzipped data can be split across chunks
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let fasta = b">chr1\nACGTTGCA\n>chr2\nGGGGAAAA\n";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(fasta).expect("compress");
        let compressed = encoder.finish().expect("finish compression");

        // Split compressed data into small chunks
        let mut hasher = FastaStreamHasher::new();
        for chunk in compressed.chunks(5) {
            hasher.update(chunk).expect("update chunk");
        }
        let collection = hasher.finish().expect("finish");

        assert_eq!(collection.sequences.len(), 2);
        assert_eq!(collection.sequences[0].metadata().name, "chr1");
        assert_eq!(collection.sequences[0].metadata().length, 8);
        assert_eq!(collection.sequences[1].metadata().name, "chr2");
        assert_eq!(collection.sequences[1].metadata().length, 8);
    }

    #[test]
    fn test_streaming_progress() {
        let mut hasher = FastaStreamHasher::new();

        assert_eq!(hasher.sequence_count(), 0);
        assert!(!hasher.in_sequence());
        assert!(hasher.current_sequence_name().is_none());

        hasher.update(b">chr1\n").expect("header");
        assert!(hasher.in_sequence());
        assert_eq!(hasher.current_sequence_name(), Some("chr1"));
        assert_eq!(hasher.current_sequence_length(), 0);

        hasher.update(b"ACGT\n").expect("sequence");
        assert_eq!(hasher.current_sequence_length(), 4);

        hasher.update(b">chr2\n").expect("next header");
        assert_eq!(hasher.sequence_count(), 1); // chr1 is finalized
        assert_eq!(hasher.current_sequence_name(), Some("chr2"));
    }

    #[test]
    fn test_streaming_chunked_matches_batch() {
        // Test that chunked streaming produces same results as batch
        let fasta_data = b">chrX\nTTGGGGAA\n>chr1\nGGAA\n>chr2\nGCGC\n";

        // Batch processing
        let batch_result = digest_fasta_bytes(fasta_data).expect("batch");

        // Streaming with various chunk sizes
        for chunk_size in [1, 2, 3, 5, 7, 11, 13, 17] {
            let mut hasher = FastaStreamHasher::new();
            for chunk in fasta_data.chunks(chunk_size) {
                hasher.update(chunk).expect("streaming chunk");
            }
            let stream_result = hasher.finish().expect("finish");

            assert_eq!(
                batch_result.metadata.digest,
                stream_result.metadata.digest,
                "Mismatch with chunk size {}",
                chunk_size
            );
        }
    }
}
