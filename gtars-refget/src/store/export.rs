//! FASTA export and BED region extraction for RefgetStore.

use super::*;
use super::readonly::ReadonlyRefgetStore;

use std::ffi::OsStr;

use indexmap::IndexMap;
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
/// description), and sequence output.
///
/// `line_width` controls sequence wrapping: a value `> 0` wraps the sequence body at
/// that many bases per line; a value of `0` disables wrapping entirely, emitting the
/// whole sequence on a single line (one-sequence-per-line, as GGCAT/SSHash expect).
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

    if line_width == 0 {
        // Unwrapped: emit the entire sequence on a single line. `chunks(0)` would
        // panic, so this branch is required, not merely an optimization.
        writer.write_all(decoded_sequence.as_bytes())?;
        writer.write_all(b"\n")?;
    } else {
        for chunk in decoded_sequence.as_bytes().chunks(line_width) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }

    Ok(())
}

// ============================================================================
// Shared per-region substring retrieval
// ============================================================================

/// Retrieve the substring for a single region, shared by the file-based and
/// vectors-based iterators.
///
/// Mutates the caller-owned `previous_parsed_chr` / `current_seq_digest` cache
/// so the consecutive-same-chrom optimization works across either path.
/// `region_label` (e.g. `"Line 5"` or `"Region 5"`) is interpolated into error
/// messages.
#[allow(clippy::too_many_arguments)]
fn retrieve_substring_for_region<K: AsRef<[u8]>>(
    store: &ReadonlyRefgetStore,
    collection_digest: &K,
    parsed_chr: String,
    parsed_start: i64,
    parsed_end: i64,
    previous_parsed_chr: &mut String,
    current_seq_digest: &mut String,
    region_label: &str,
) -> Result<RetrievedSequence> {
    if parsed_start < 0 || parsed_end < 0 {
        return Err(anyhow!(
            "{} has invalid start or end coordinates: start={}, end={}",
            region_label,
            parsed_start,
            parsed_end
        ));
    }

    if *previous_parsed_chr != parsed_chr {
        *previous_parsed_chr = parsed_chr.clone();

        let result = store
            .get_sequence_by_name(collection_digest, &parsed_chr)
            .map_err(|e| {
                anyhow!(
                    "{}: sequence '{}' not found in collection '{}': {}",
                    region_label,
                    parsed_chr,
                    String::from_utf8_lossy(collection_digest.as_ref()),
                    e
                )
            })?;

        *current_seq_digest = result.metadata().sha512t24u.clone();
    }

    let retrieved_substring = store
        .get_substring(&*current_seq_digest, parsed_start as usize, parsed_end as usize)
        .map_err(|e| {
            anyhow!(
                "{}: failed to get substring for digest '{}' from {} to {}: {}",
                region_label,
                current_seq_digest,
                parsed_start,
                parsed_end,
                e
            )
        })?;

    Ok(RetrievedSequence {
        sequence: retrieved_substring,
        chrom_name: parsed_chr,
        start: parsed_start as u32,
        end: parsed_end as u32,
    })
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

        Some(retrieve_substring_for_region(
            self.store,
            &self.collection_digest,
            parsed_chr,
            parsed_start as i64,
            parsed_end as i64,
            &mut self.previous_parsed_chr,
            &mut self.current_seq_digest,
            &format!("Line {}", self.line_num + 1),
        ))
    }
}

// ============================================================================
// SubstringsFromRegionVectors Iterator impl
// ============================================================================

impl<K> Iterator for SubstringsFromRegionVectors<'_, K>
where
    K: AsRef<[u8]>,
{
    type Item = Result<RetrievedSequence, anyhow::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.chroms.len() {
            return None;
        }
        let i = self.index;
        self.index += 1;

        let parsed_chr = self.chroms[i].clone();
        let parsed_start = self.starts[i] as i64;
        let parsed_end = self.ends[i] as i64;

        Some(retrieve_substring_for_region(
            self.store,
            &self.collection_digest,
            parsed_chr,
            parsed_start,
            parsed_end,
            &mut self.previous_parsed_chr,
            &mut self.current_seq_digest,
            &format!("Region {}", i + 1),
        ))
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

    /// Get an iterator over substrings defined by in-memory region vectors.
    ///
    /// `chroms`, `starts`, and `ends` are parallel slices; element `i` defines
    /// region `chroms[i]:starts[i]-ends[i]`. Yields the same `RetrievedSequence`
    /// results as `substrings_from_regions` without any temp-file round-trip.
    pub fn substrings_from_region_vectors<'a, K, S>(
        &'a self,
        collection_digest: K,
        chroms: &[S],
        starts: &[u32],
        ends: &[u32],
    ) -> Result<SubstringsFromRegionVectors<'a, K>>
    where
        K: AsRef<[u8]>,
        S: AsRef<str>,
    {
        if chroms.len() != starts.len() || chroms.len() != ends.len() {
            return Err(anyhow!(
                "Mismatched region vector lengths: chroms={}, starts={}, ends={}",
                chroms.len(),
                starts.len(),
                ends.len()
            ));
        }

        Ok(SubstringsFromRegionVectors {
            store: self,
            collection_digest,
            chroms: chroms.iter().map(|c| c.as_ref().to_string()).collect(),
            starts: starts.to_vec(),
            ends: ends.to_vec(),
            index: 0,
            previous_parsed_chr: String::new(),
            current_seq_digest: String::new(),
        })
    }

    /// Export sequences from BED file regions to a FASTA file.
    ///
    /// Region export is inherently unwrapped: each region's sequence is written on a
    /// single line under a `>chrom:start-end` header, regardless of any line width.
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

        let seq_iter = self.substrings_from_regions(&collection_digest, bed_file_path)?;

        for rs in seq_iter.into_iter() {
            let rs = rs?;

            // Write one header per region with coordinates
            let header = format!(">{}:{}-{}\n", rs.chrom_name, rs.start, rs.end);
            writer.write_all(header.as_bytes())?;
            writer.write_all(rs.sequence.as_bytes())?;
            writer.write_all(b"\n")?;
        }

        writer.flush()?;

        Ok(())
    }

    /// Export sequences from a collection to a FASTA file.
    ///
    /// `line_width` controls sequence wrapping via the sentinel convention:
    /// `None` defaults to 80 bases per line; `Some(n)` with `n > 0` wraps at `n`;
    /// `Some(0)` disables wrapping, emitting one sequence per line (unwrapped) as
    /// GGCAT/SSHash and similar k-mer tooling expect.
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

        // Headers come from the collection's own per-sequence records so a
        // sequence shared across collections is exported under this collection's
        // name and description, not the label of whichever import came first.
        let records = self.collection_sequence_metadata(&collection_key).map_err(|e| {
            anyhow!(
                "Collection not found: {:?}: {}",
                String::from_utf8_lossy(collection_digest.as_ref()),
                e
            )
        })?;
        let name_to_meta: IndexMap<&str, &SequenceMetadata> = records
            .iter()
            .map(|r| (r.metadata().name.as_str(), r.metadata()))
            .collect();

        let names_to_export: Vec<&str> = if let Some(names) = sequence_names {
            names
        } else {
            name_to_meta.keys().copied().collect()
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
            let meta = *name_to_meta
                .get(seq_name)
                .ok_or_else(|| anyhow!("Sequence '{}' not found in collection", seq_name))?;
            let seq_digest = meta.sha512t24u.to_key();

            let record = self
                .sequence_store
                .get(&seq_digest)
                .ok_or_else(|| anyhow!("Sequence record not found for digest: {}", meta.sha512t24u))?;

            // Only the bytes come from the global record. `alphabet` and `length`
            // are content-derived and identical in the per-collection stub.
            let sequence_data: &[u8] = match record {
                SequenceRecord::Stub(_) => {
                    return Err(anyhow!("Sequence data not loaded for '{}'. Call load_sequence() or load_all_sequences() first.", seq_name));
                }
                SequenceRecord::Full { sequence, .. } => sequence.as_slice(),
            };

            write_fasta_record(&mut *writer, meta, sequence_data, self.mode, line_width)?;
        }

        writer.flush()?;

        Ok(())
    }

    /// Export sequences by their sequence digests to a FASTA file.
    ///
    /// There is no collection context here, so each header uses the store-wide
    /// record's name and description (first import wins for a sequence shared by
    /// several collections). Use [`Self::export_fasta`] for collection-specific names.
    ///
    /// `line_width` follows the same sentinel convention as [`Self::export_fasta`]:
    /// `None` defaults to 80; `Some(n)` with `n > 0` wraps at `n`; `Some(0)` disables
    /// wrapping, emitting one sequence per line (unwrapped).
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

            let (metadata, sequence_data): (&SequenceMetadata, &[u8]) = match record {
                SequenceRecord::Stub(_) => {
                    return Err(anyhow!(
                        "Sequence data not loaded for digest: {}",
                        digest_str
                    ));
                }
                SequenceRecord::Full { metadata, sequence } => (metadata, sequence.as_slice()),
            };

            write_fasta_record(&mut *writer, metadata, sequence_data, self.mode, line_width)?;
        }

        writer.flush()?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::store::{FastaImportOptions, RefgetStore};
    use std::io::Write;
    use tempfile::{NamedTempFile, TempDir};

    /// Build a small in-memory store from a temp FASTA and return it plus the
    /// digest of the single collection created.
    fn build_store() -> (RefgetStore, String) {
        let mut fasta = NamedTempFile::new().expect("create temp fasta");
        // A sequence longer than 80 bases so the wrapped export spans >1 line.
        writeln!(fasta, ">seq1").unwrap();
        writeln!(fasta, "{}", "ACGT".repeat(30)).unwrap(); // 120 bases
        writeln!(fasta, ">seq2").unwrap();
        writeln!(fasta, "TTGGCCAA").unwrap();
        fasta.flush().unwrap();

        let mut store = RefgetStore::in_memory();
        store
            .add_sequence_collection_from_fasta(fasta.path(), FastaImportOptions::new())
            .expect("import fasta");

        let collections = store.list_collections(0, usize::MAX, &[]).unwrap();
        let digest = collections.results[0].digest.clone();
        (store, digest)
    }

    /// Two single-record FASTAs with byte-identical sequences but different
    /// headers: `>chr1 ucsc desc` and `>1 flybase desc`. Mirrors the UCSC vs
    /// FlyBase situation from issue #270.
    fn shared_sequence_fastas() -> (NamedTempFile, NamedTempFile) {
        let seq = "ACGT".repeat(10);
        let mut a = NamedTempFile::new().expect("create temp fasta A");
        writeln!(a, ">chr1 ucsc desc").unwrap();
        writeln!(a, "{}", seq).unwrap();
        a.flush().unwrap();

        let mut b = NamedTempFile::new().expect("create temp fasta B");
        writeln!(b, ">1 flybase desc").unwrap();
        writeln!(b, "{}", seq).unwrap();
        b.flush().unwrap();
        (a, b)
    }

    /// Import `first` then `second` into `store`; return their collection digests.
    fn import_pair(store: &mut RefgetStore, first: &NamedTempFile, second: &NamedTempFile) -> (String, String) {
        let (meta_first, _) = store
            .add_sequence_collection_from_fasta(first.path(), FastaImportOptions::new())
            .expect("import first fasta");
        let (meta_second, _) = store
            .add_sequence_collection_from_fasta(second.path(), FastaImportOptions::new())
            .expect("import second fasta");
        (meta_first.digest, meta_second.digest)
    }

    /// Export a collection unwrapped and return the headers in order.
    fn export_headers(store: &crate::store::ReadonlyRefgetStore, digest: &str, dir: &TempDir) -> Vec<String> {
        let out = dir.path().join(format!("{}.fa", digest));
        store
            .export_fasta(digest, &out, None, Some(0))
            .expect("export collection");
        let content = std::fs::read_to_string(&out).unwrap();
        records(&content).into_iter().map(|(h, _)| h).collect()
    }

    /// Split a FASTA into (header, body-lines) pairs.
    fn records(fasta: &str) -> Vec<(String, Vec<String>)> {
        let mut out: Vec<(String, Vec<String>)> = Vec::new();
        for line in fasta.lines() {
            if let Some(header) = line.strip_prefix('>') {
                out.push((header.to_string(), Vec::new()));
            } else if let Some(last) = out.last_mut() {
                last.1.push(line.to_string());
            }
        }
        out
    }

    #[test]
    fn export_unwrapped_puts_each_sequence_on_one_line() {
        let (store, digest) = build_store();
        let readonly = store.into_readonly();

        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unwrapped.fa");
        readonly
            .export_fasta(&digest, &out, None, Some(0))
            .expect("export unwrapped");

        let content = std::fs::read_to_string(&out).unwrap();
        let recs = records(&content);
        assert_eq!(recs.len(), 2, "expected two records");
        for (header, body) in &recs {
            assert_eq!(
                body.len(),
                1,
                "record '{}' should have exactly one body line when unwrapped, got {:?}",
                header,
                body
            );
        }
        // seq1 body is the full 120 bases on a single line.
        assert_eq!(recs[0].1[0].len(), 120);
    }

    #[test]
    fn export_wrapped_at_80_splits_long_sequences() {
        let (store, digest) = build_store();
        let readonly = store.into_readonly();

        let dir = TempDir::new().unwrap();
        let out = dir.path().join("wrapped.fa");
        readonly
            .export_fasta(&digest, &out, None, Some(80))
            .expect("export wrapped");

        let content = std::fs::read_to_string(&out).unwrap();
        let recs = records(&content);
        // seq1 (120 bases) wraps to two lines at width 80; seq2 (8 bases) stays one.
        assert_eq!(recs[0].1.len(), 2, "120-base seq should wrap to 2 lines");
        assert_eq!(recs[0].1[0].len(), 80);
        assert_eq!(recs[0].1[1].len(), 40);
        assert_eq!(recs[1].1.len(), 1);
    }

    #[test]
    fn shared_sequence_exports_with_each_collections_own_header() {
        let (fasta_a, fasta_b) = shared_sequence_fastas();

        // Both import orders: the first-imported label must never leak into
        // the other collection's export.
        for swap in [false, true] {
            let mut store = RefgetStore::in_memory();
            let (digest_a, digest_b) = if swap {
                let (b, a) = import_pair(&mut store, &fasta_b, &fasta_a);
                (a, b)
            } else {
                import_pair(&mut store, &fasta_a, &fasta_b)
            };
            let readonly = store.into_readonly();
            let dir = TempDir::new().unwrap();

            assert_eq!(
                export_headers(&readonly, &digest_a, &dir),
                vec!["chr1 ucsc desc".to_string()],
                "swap={}",
                swap
            );
            assert_eq!(
                export_headers(&readonly, &digest_b, &dir),
                vec!["1 flybase desc".to_string()],
                "swap={}",
                swap
            );
        }
    }

    #[test]
    fn get_collection_returns_per_collection_description() {
        let (fasta_a, fasta_b) = shared_sequence_fastas();
        let mut store = RefgetStore::in_memory();
        let (digest_a, digest_b) = import_pair(&mut store, &fasta_a, &fasta_b);
        let readonly = store.into_readonly();

        let coll_b = readonly.get_collection(&digest_b).expect("get collection B");
        assert_eq!(coll_b.sequences.len(), 1);
        assert_eq!(coll_b.sequences[0].metadata().name, "1");
        assert_eq!(
            coll_b.sequences[0].metadata().description.as_deref(),
            Some("flybase desc")
        );

        let coll_a = readonly.get_collection(&digest_a).expect("get collection A");
        assert_eq!(coll_a.sequences[0].metadata().name, "chr1");
        assert_eq!(
            coll_a.sequences[0].metadata().description.as_deref(),
            Some("ucsc desc")
        );
    }

    #[test]
    fn shared_sequence_disk_roundtrip_exports_own_header() {
        let (fasta_a, fasta_b) = shared_sequence_fastas();
        let store_dir = TempDir::new().unwrap();

        let (digest_a, digest_b) = {
            let mut store = RefgetStore::on_disk(store_dir.path()).expect("create disk store");
            import_pair(&mut store, &fasta_a, &fasta_b)
        };

        // Reopen so collection stubs come from collections/<digest>.rgsi via
        // ensure_collection_loaded, not from the import-time in-memory records.
        let mut store = RefgetStore::open_local(store_dir.path()).expect("reopen store");
        store.load_collection(&digest_b).expect("load collection B");
        store.load_collection(&digest_a).expect("load collection A");
        store.load_all_sequences().expect("load sequences");
        let readonly = store.into_readonly();

        let dir = TempDir::new().unwrap();
        assert_eq!(
            export_headers(&readonly, &digest_b, &dir),
            vec!["1 flybase desc".to_string()]
        );
        assert_eq!(
            export_headers(&readonly, &digest_a, &dir),
            vec!["chr1 ucsc desc".to_string()]
        );

        let coll_b = readonly.get_collection(&digest_b).unwrap();
        assert_eq!(
            coll_b.sequences[0].metadata().description.as_deref(),
            Some("flybase desc")
        );
    }
}
