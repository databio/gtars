//! FASTA import pipeline for RefgetStore.
//!
//! Contains the multithreaded FASTA import pipeline and cached metadata fast path.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

use anyhow::{anyhow, Result};

use crate::collection::SequenceCollectionExt;
use crate::digest::{
    SequenceCollection, SequenceCollectionMetadata,
    SequenceEncoder, SequenceMetadata, SequenceRecord,
};
use crate::hashkeyable::HashKeyable;

// ============================================================================
// ReadonlyRefgetStore import methods
// ============================================================================

impl ReadonlyRefgetStore {
    /// Import a FASTA file into the store using a multithreaded pipeline.
    ///
    /// After the pipeline finishes, computes collection metadata, registers the
    /// collection, and inserts all sequences.
    pub fn add_sequence_collection_from_fasta<P: AsRef<Path>>(
        &mut self,
        file_path: P,
        opts: FastaImportOptions<'_>,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
        use crossbeam_channel::bounded;
        use md5::Md5;
        use sha2::{Digest, Sha512};

        if !self.quiet {
            println!("Processing {}...", file_path.as_ref().display());
        }

        let pipeline_start = Instant::now();

        // Check for RGSI cache
        let use_cache = self.local_path.is_some();
        use crate::utils::PathExtension;
        let rgsi_path = file_path.as_ref().replace_exts_with("rgsi");
        let have_rgsi = use_cache && rgsi_path.exists();

        if have_rgsi {
            match self.add_fasta_with_cached_metadata(&file_path, &rgsi_path, opts, pipeline_start) {
                Ok(result) => return Ok(result),
                Err(_) => {
                    // Cache was stale/empty, fall through to pipeline
                }
            }
        }

        // --- Channel message types ---

        /// Threshold above which Thread 1 processes sequences inline
        const LARGE_SEQ_THRESHOLD: usize = 500 * 1024 * 1024; // 500 MB

        struct DecompressedSequence {
            name: String,
            description: Option<String>,
            raw_header: String,
            raw_bytes: Vec<u8>,
        }

        struct DigestedSequence {
            metadata: SequenceMetadata,
            raw_bytes: Vec<u8>,
            aliases: Vec<(String, String)>,
        }

        struct ReadySequence {
            metadata: SequenceMetadata,
            sequence_data: Vec<u8>,
            aliases: Vec<(String, String)>,
        }

        enum ToDigest {
            NeedsWork(DecompressedSequence),
            AlreadyDone(ReadySequence),
        }

        enum ToEncode {
            NeedsWork(DigestedSequence),
            AlreadyDone(ReadySequence),
        }

        // --- Set up channels ---
        let (decompress_tx, decompress_rx) = bounded::<ToDigest>(1);
        let (digest_tx, digest_rx) = bounded::<ToEncode>(1);

        let file_path_buf = file_path.as_ref().to_path_buf();
        let namespaces: Vec<String> = opts.namespaces.iter().map(|s| s.to_string()).collect();
        let ns_for_digest = namespaces.clone();

        // --- Thread 1: Read FASTA ---
        let quiet = self.quiet;
        let mode_for_t1 = self.mode;
        let decompress_handle = std::thread::spawn(move || -> Result<()> {
            use sha2::{Digest, Sha512};
            use md5::Md5;

            let mut fasta_reader = crate::fasta::FastaReader::from_path(&file_path_buf)?;

            while let Some(record) = fasta_reader.next_record()? {
                if record.raw_bytes.len() > LARGE_SEQ_THRESHOLD {
                    if !quiet {
                        println!(
                            "  Large sequence '{}' ({} MB) -- processing inline to reduce memory",
                            record.name,
                            record.raw_bytes.len() / (1024 * 1024),
                        );
                    }
                    let crate::fasta::FastaRecord { name, description, raw_header, raw_bytes } = record;

                    let mut sha512_hasher = Sha512::new();
                    sha512_hasher.update(&raw_bytes);
                    let sha512 = base64_url::encode(&sha512_hasher.finalize()[0..24]);

                    let mut md5_hasher = Md5::new();
                    md5_hasher.update(&raw_bytes);
                    let md5 = format!("{:x}", md5_hasher.finalize());

                    let mut guesser = crate::digest::AlphabetGuesser::new();
                    guesser.update(&raw_bytes);
                    let alphabet = guesser.guess();

                    let length = raw_bytes.len();
                    let sequence_data = match mode_for_t1 {
                        StorageMode::Encoded => {
                            let mut encoder = SequenceEncoder::new(alphabet, length);
                            encoder.update(&raw_bytes);
                            drop(raw_bytes);
                            encoder.finalize()
                        }
                        StorageMode::Raw => raw_bytes,
                    };

                    let metadata = SequenceMetadata {
                        name,
                        description,
                        length,
                        sha512t24u: sha512,
                        md5,
                        alphabet,
                        fai: None,
                    };

                    let aliases = if !namespaces.is_empty() {
                        let ns_refs: Vec<&str> = namespaces.iter().map(|s| s.as_str()).collect();
                        crate::digest::fasta::extract_aliases_from_header(&raw_header, &ns_refs)
                    } else {
                        vec![]
                    };

                    decompress_tx.send(ToDigest::AlreadyDone(ReadySequence {
                        metadata,
                        sequence_data,
                        aliases,
                    })).map_err(|_| anyhow!("Digest thread stopped receiving"))?;
                } else {
                    decompress_tx.send(ToDigest::NeedsWork(DecompressedSequence {
                        name: record.name,
                        description: record.description,
                        raw_header: record.raw_header,
                        raw_bytes: record.raw_bytes,
                    })).map_err(|_| anyhow!("Digest thread stopped receiving"))?;
                }
            }
            Ok(())
        });

        // --- Thread 2: Digest ---
        let digest_handle = std::thread::spawn(move || -> Result<()> {
            let ns_refs: Vec<&str> = ns_for_digest.iter().map(|s| s.as_str()).collect();

            for msg in decompress_rx {
                match msg {
                    ToDigest::NeedsWork(seq) => {
                        let mut sha512_hasher = Sha512::new();
                        sha512_hasher.update(&seq.raw_bytes);
                        let sha512 = base64_url::encode(&sha512_hasher.finalize()[0..24]);

                        let mut md5_hasher = Md5::new();
                        md5_hasher.update(&seq.raw_bytes);
                        let md5 = format!("{:x}", md5_hasher.finalize());

                        let mut guesser = crate::digest::AlphabetGuesser::new();
                        guesser.update(&seq.raw_bytes);
                        let alphabet = guesser.guess();

                        let metadata = SequenceMetadata {
                            name: seq.name,
                            description: seq.description,
                            length: seq.raw_bytes.len(),
                            sha512t24u: sha512,
                            md5,
                            alphabet,
                            fai: None,
                        };

                        let aliases = if !ns_refs.is_empty() {
                            crate::digest::fasta::extract_aliases_from_header(&seq.raw_header, &ns_refs)
                        } else {
                            vec![]
                        };

                        digest_tx.send(ToEncode::NeedsWork(DigestedSequence {
                            metadata,
                            raw_bytes: seq.raw_bytes,
                            aliases,
                        })).map_err(|_| anyhow!("Encode thread stopped receiving"))?;
                    }
                    ToDigest::AlreadyDone(ready) => {
                        digest_tx.send(ToEncode::AlreadyDone(ready))
                            .map_err(|_| anyhow!("Encode thread stopped receiving"))?;
                    }
                }
            }
            drop(digest_tx);
            Ok(())
        });

        // --- Thread 3 (main thread): Encode + stream to disk ---
        let mode = self.mode;
        let mut sequence_metadata: Vec<SequenceMetadata> = Vec::new();
        let mut all_aliases: Vec<(String, String, String)> = Vec::new();

        for msg in digest_rx {
            match msg {
                ToEncode::NeedsWork(digested) => {
                    let metadata = digested.metadata;
                    let aliases = digested.aliases;
                    let raw_bytes = digested.raw_bytes;

                    let sequence_data = match mode {
                        StorageMode::Encoded => {
                            let mut encoder = SequenceEncoder::new(metadata.alphabet, metadata.length);
                            encoder.update(&raw_bytes);
                            drop(raw_bytes);
                            encoder.finalize()
                        }
                        StorageMode::Raw => raw_bytes,
                    };

                    for (ns, alias_value) in &aliases {
                        all_aliases.push((ns.clone(), alias_value.clone(), metadata.sha512t24u.clone()));
                    }

                    self.add_sequence_record(
                        SequenceRecord::Full {
                            metadata: metadata.clone(),
                            sequence: std::sync::Arc::new(sequence_data),
                        },
                        true,
                    )?;

                    sequence_metadata.push(metadata);
                }
                ToEncode::AlreadyDone(ready) => {
                    for (ns, alias_value) in &ready.aliases {
                        all_aliases.push((ns.clone(), alias_value.clone(), ready.metadata.sha512t24u.clone()));
                    }

                    self.add_sequence_record(
                        SequenceRecord::Full {
                            metadata: ready.metadata.clone(),
                            sequence: std::sync::Arc::new(ready.sequence_data),
                        },
                        true,
                    )?;

                    sequence_metadata.push(ready.metadata);
                }
            }
        }

        // Join threads and propagate errors
        decompress_handle.join()
            .map_err(|e| anyhow!("Decompress thread panicked: {:?}", e))??;
        digest_handle.join()
            .map_err(|e| anyhow!("Digest thread panicked: {:?}", e))??;

        let pipeline_elapsed = pipeline_start.elapsed();

        // --- Post-pipeline: compute collection metadata, register ---
        let metadata_start = std::time::Instant::now();
        let seq_count = sequence_metadata.len();

        let stub_records: Vec<SequenceRecord> = sequence_metadata
            .iter()
            .map(|meta| SequenceRecord::Stub(meta.clone()))
            .collect();

        let mut seqcol_metadata = SequenceCollectionMetadata::from_sequences(
            &stub_records,
            Some(file_path.as_ref().to_path_buf()),
        );
        if self.ancillary_digests {
            seqcol_metadata.compute_ancillary_digests(&stub_records);
        }

        let coll_key = seqcol_metadata.digest.to_key();
        let coll_digest_display = seqcol_metadata.digest.clone();
        let metadata = seqcol_metadata.clone();
        let metadata_elapsed = metadata_start.elapsed();

        if !opts.force && self.collections.contains_key(&coll_key) {
            if !self.quiet {
                println!("Skipped {} (already exists)", coll_digest_display);
            }
            return Ok((metadata, false));
        }

        let index_start = std::time::Instant::now();

        // Write RGSI cache for next time
        if use_cache {
            let seqcol_for_cache = SequenceCollection {
                metadata: seqcol_metadata.clone(),
                sequences: stub_records.clone(),
            };
            let _ = seqcol_for_cache.write_collection_rgsi(&rgsi_path);
        }

        // Register the collection
        let seqcol = SequenceCollection {
            metadata: seqcol_metadata,
            sequences: stub_records,
        };
        self.add_sequence_collection_internal(seqcol, opts.force)?;

        // Register aliases
        for (ns, alias_value, sha512t24u) in &all_aliases {
            self.add_sequence_alias(ns, alias_value, sha512t24u)?;
        }

        // Populate name_lookup for already-written sequences
        for meta in &sequence_metadata {
            self.name_lookup
                .entry(coll_key)
                .or_default()
                .insert(meta.name.clone(), meta.sha512t24u.to_key());
        }
        let index_elapsed = index_start.elapsed();

        if !self.quiet {
            let total = pipeline_elapsed + metadata_elapsed + index_elapsed;
            println!(
                "Added {} {} seqs in {:.1}s ({:.1} proc | {:.1} meta | {:.1} index)",
                coll_digest_display,
                seq_count,
                total.as_secs_f64(),
                pipeline_elapsed.as_secs_f64(),
                metadata_elapsed.as_secs_f64(),
                index_elapsed.as_secs_f64(),
            );
        }

        Ok((metadata, true))
    }

    /// Fast path when RGSI cache exists: skip digesting, only decompress for encoding.
    fn add_fasta_with_cached_metadata<P: AsRef<Path>>(
        &mut self,
        file_path: P,
        rgsi_path: &Path,
        opts: FastaImportOptions<'_>,
        start_time: Instant,
    ) -> Result<(SequenceCollectionMetadata, bool)> {
        let mut seqcol = crate::collection::read_rgsi_file(rgsi_path)?;

        if seqcol.sequences.is_empty() {
            let _ = std::fs::remove_file(rgsi_path);
            return Err(anyhow!("Empty RGSI cache"));
        }
        let coll_key = seqcol.metadata.digest.to_key();
        let coll_digest_display = seqcol.metadata.digest.clone();

        if !opts.force && self.collections.contains_key(&coll_key) {
            if !self.quiet {
                println!("Skipped {} (already exists)", coll_digest_display);
            }
            return Ok((seqcol.metadata.clone(), false));
        }

        if self.ancillary_digests {
            seqcol.metadata.compute_ancillary_digests(&seqcol.sequences);
        }
        let metadata = seqcol.metadata.clone();

        let seqmeta_hashmap: HashMap<String, SequenceMetadata> = seqcol
            .sequences
            .iter()
            .map(|r| {
                let meta = r.metadata().clone();
                (meta.name.clone(), meta)
            })
            .collect();

        self.add_sequence_collection_internal(seqcol, opts.force)?;

        // Decompress and encode sequences
        let mut fasta_reader = crate::fasta::FastaReader::from_path(file_path.as_ref())?;
        let mut seq_count = 0;

        while let Some(record) = fasta_reader.next_record()? {
            let (name, _) = crate::fasta::parse_fasta_header(&record.raw_header);

            if !opts.namespaces.is_empty() {
                let aliases = crate::digest::fasta::extract_aliases_from_header(&record.raw_header, opts.namespaces);
                for (ns, alias_value) in aliases {
                    if let Some(meta) = seqmeta_hashmap.get(&name) {
                        self.add_sequence_alias(&ns, &alias_value, &meta.sha512t24u)?;
                    }
                }
            }

            let dr = seqmeta_hashmap
                .get(&name)
                .ok_or_else(|| anyhow!("Sequence '{}' not found in cached metadata", name))?
                .clone();

            seq_count += 1;

            let sequence_data = match self.mode {
                StorageMode::Encoded => {
                    let mut encoder = SequenceEncoder::new(dr.alphabet, dr.length);
                    encoder.update(&record.raw_bytes);
                    encoder.finalize()
                }
                StorageMode::Raw => record.raw_bytes,
            };

            self.add_sequence(
                SequenceRecord::Full { metadata: dr, sequence: std::sync::Arc::new(sequence_data) },
                coll_key,
                true,
            )?;
        }

        let elapsed = start_time.elapsed();
        if !self.quiet {
            println!(
                "Added {} ({} seqs) in {:.1}s [cached metadata]",
                coll_digest_display, seq_count, elapsed.as_secs_f64()
            );
        }

        Ok((metadata, true))
    }
}
