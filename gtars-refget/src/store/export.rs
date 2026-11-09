//! FASTA export and BED region extraction for RefgetStore.
//!
//! Exports stream records in order: a prefetch reads each record's stored
//! bytes (from RAM, a local `.seq` file, the shared cache, or the remote store,
//! many at once) ahead of the writer, within a fixed memory budget (1 GiB),
//! and the writer decodes them in chunks. An export never loads sequences into the store, and
//! writes downloaded `.seq` files to the shared cache only with
//! [`ExportOptions::write_cache`].

use super::*;
use super::label::{label_err, AliasMaps, LabelObject, LabelTemplate};
use super::readonly::ReadonlyRefgetStore;
use super::fetch::FetchCounters;
use super::prefetch::{PrefetchOptions, PrefetchOrder, DEFAULT_MEMORY_BUDGET};

use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{anyhow, Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use flate2::read::MultiGzDecoder;
use gtars_core::utils::get_file_info;

use crate::digest::{
    AlphabetType, SequenceMetadata, decode_substring_from_bytes, lookup_alphabet, strip_sq_prefix,
};
use crate::hashkeyable::{DigestKey, HashKeyable, key_to_digest_string};

/// Bases decoded per chunk when streaming a record.
pub(crate) const EXPORT_CHUNK_BASES: usize = 1 << 20;

/// Records between progress lines on stderr.
const PROGRESS_EVERY: usize = 1000;

// ============================================================================
// Public export types
// ============================================================================

/// Which records an export writes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExportSelection {
    /// One or more collections, optionally limited to some sequence names in each.
    Collections { digests: Vec<String>, names: Option<Vec<String>> },
    /// Specific sequences by digest. These records have no collection.
    SeqDigests(Vec<String>),
    /// Every collection, plus any sequence that is in no collection.
    All,
}

impl ExportSelection {
    /// Build a selection from binding-style arguments: exactly one of
    /// `collections`, `seq_digests` or `all`, and `names` only with collections.
    pub fn from_parts(
        collections: Option<Vec<String>>,
        names: Option<Vec<String>>,
        seq_digests: Option<Vec<String>>,
        all: bool,
    ) -> Result<Self> {
        if collections.is_some() as u8 + seq_digests.is_some() as u8 + all as u8 != 1 {
            return Err(label_err!("Set exactly one of collections, seq_digests or all"));
        }
        Ok(match (collections, seq_digests) {
            (Some(digests), _) => ExportSelection::Collections { digests, names },
            _ if names.is_some() => {
                return Err(label_err!("names only applies with collections"))
            }
            (_, Some(digests)) => ExportSelection::SeqDigests(digests),
            _ => ExportSelection::All,
        })
    }
}

/// Molecule filter: BLAST and similar tools need one molecule type per database.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Molecule {
    /// `Dna2bit`, `Dna3bit` and `DnaIupac` sequences.
    Dna,
    /// `Protein` sequences.
    Protein,
}

impl Molecule {
    /// Parse an optional `"dna"` / `"protein"` argument.
    pub fn from_opt(s: Option<&str>) -> Result<Option<Self>> {
        s.map(str::parse).transpose()
    }

    fn keeps(&self, alphabet: &AlphabetType) -> bool {
        match self {
            Molecule::Dna => matches!(
                alphabet,
                AlphabetType::Dna2bit | AlphabetType::Dna3bit | AlphabetType::DnaIupac
            ),
            Molecule::Protein => matches!(alphabet, AlphabetType::Protein),
        }
    }
}

impl std::str::FromStr for Molecule {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self> {
        match s.to_ascii_lowercase().as_str() {
            "dna" => Ok(Molecule::Dna),
            "protein" => Ok(Molecule::Protein),
            other => Err(label_err!("Invalid molecule '{}': expected 'dna' or 'protein'", other)),
        }
    }
}

/// How an export labels and formats its records.
#[derive(Debug, Clone, Default)]
pub struct ExportOptions {
    /// Label template; default `{seq}`.
    pub label: LabelTemplate,
    /// Keep only sequences of this molecule type.
    pub molecule: Option<Molecule>,
    /// Bases per line: `None` is 80, `Some(0)` writes each sequence on one line.
    pub line_width: Option<usize>,
    /// Publish sequences downloaded from a remote store to the shared cache,
    /// so an interrupted export can resume without downloading them again.
    pub write_cache: bool,
}

/// Counts reported by an export.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ExportSummary {
    /// Records written.
    pub n_records: usize,
    /// Bases written.
    pub n_bases: u64,
    /// Records whose label fell back to a digest (missing alias, or no collection).
    pub n_fallbacks: usize,
    /// Records skipped because the same label and sequence were already written.
    pub n_deduplicated: usize,
    /// What reading the sequences took: network, cache and local reads.
    pub fetch: FetchCounters,
}

// ============================================================================
// Output and record writing
// ============================================================================

/// Open an export destination. `"-"` is stdout; a `.gz` path is gzipped.
pub(crate) fn open_output(path: &Path) -> Result<Box<dyn Write>> {
    if path == Path::new("-") {
        return Ok(Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout().lock())));
    }
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent)?;
    }
    let file = File::create(path)
        .with_context(|| format!("Failed to create output file: {}", path.display()))?;
    Ok(if path.extension() == Some(OsStr::new("gz")) {
        Box::new(GzEncoder::new(BufWriter::new(file), Compression::default()))
    } else {
        Box::new(BufWriter::new(file))
    })
}

/// Write `bases` as lines of `line_width` (`0`: no wrapping, no newline).
/// Callers keep every chunk but the last a multiple of `line_width`.
fn write_lines(writer: &mut dyn Write, bases: &[u8], line_width: usize) -> Result<()> {
    if line_width == 0 {
        writer.write_all(bases)?;
    } else {
        for line in bases.chunks(line_width) {
            writer.write_all(line)?;
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

/// Read until `buf` is full or the reader ends; return the bytes read.
fn read_full(reader: &mut impl std::io::Read, buf: &mut [u8]) -> Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match reader.read(&mut buf[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e.into()),
        }
    }
    Ok(filled)
}

/// Write one FASTA record, decoding the body in chunks of about `chunk_bases`.
///
/// `header` is the label, written as `>{header}`. `data` is the stored bytes
/// (packed, raw or zstd per `mode`). `line_width > 0` wraps the body;
/// `0` writes it on a single line.
pub(crate) fn write_fasta_record(
    writer: &mut dyn Write,
    header: &str,
    meta: &SequenceMetadata,
    data: &[u8],
    mode: StorageMode,
    line_width: usize,
    chunk_bases: usize,
) -> Result<()> {
    // A multiple of line_width, so line breaks never straddle a chunk.
    let chunk = if line_width == 0 {
        chunk_bases.max(1)
    } else {
        (chunk_bases / line_width).max(1) * line_width
    };

    writer.write_all(b">")?;
    writer.write_all(header.as_bytes())?;
    writer.write_all(b"\n")?;

    match mode {
        StorageMode::Encoded => {
            let alphabet = lookup_alphabet(&meta.alphabet);
            let mut start = 0;
            while start < meta.length {
                let end = (start + chunk).min(meta.length);
                let decoded = decode_substring_from_bytes(data, start, end, alphabet);
                write_lines(writer, &decoded, line_width)?;
                start = end;
            }
        }
        StorageMode::Raw => {
            for piece in data.chunks(chunk) {
                write_lines(writer, piece, line_width)?;
            }
        }
        StorageMode::Zstd => {
            let mut decoder = zstd::stream::read::Decoder::with_buffer(data)
                .context("Failed to open zstd sequence record")?;
            let mut buf = vec![0u8; chunk];
            loop {
                let n = read_full(&mut decoder, &mut buf)
                    .context("Failed to zstd-decompress sequence record")?;
                if n == 0 {
                    break;
                }
                write_lines(writer, &buf[..n], line_width)?;
                if n < chunk {
                    break;
                }
            }
        }
    }

    if line_width == 0 {
        writer.write_all(b"\n")?;
    }
    Ok(())
}

/// One record to write: its label and the sequence behind it.
struct PlannedRecord {
    label: String,
    key: DigestKey,
}

/// A record before labeling: the sequence plus how it was reached.
struct Candidate<'a> {
    key: DigestKey,
    meta: &'a SequenceMetadata,
    /// `(collection digest, sequence name)` when reached through a collection.
    collection: Option<(String, &'a str)>,
}

fn source_of(collection: &Option<(String, &str)>) -> String {
    match collection {
        Some((col, _)) => format!("collection {}", col),
        None => "the sequence list".to_string(),
    }
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

    /// Digest -> alias maps for every namespace `label` uses, after checking
    /// each namespace exists in this store.
    pub(crate) fn label_alias_maps(&self, label: &LabelTemplate) -> Result<AliasMaps> {
        let mut maps = AliasMaps::default();
        let seq_namespaces = self.aliases.sequence_namespaces();
        for ns in label.namespaces(LabelObject::Seq) {
            if !seq_namespaces.contains(&ns) {
                return Err(label_err!(
                    "Unknown sequence-alias namespace '{}' in label '{}'. This store's sequence-alias namespaces: {}",
                    ns,
                    label,
                    list_or_none(&seq_namespaces)
                ));
            }
            maps.add_sequence_namespace(&ns, self.aliases.namespace_map(AliasKind::Sequence, &ns));
        }
        let col_namespaces = self.aliases.collection_namespaces();
        for ns in label.namespaces(LabelObject::Col) {
            if !col_namespaces.contains(&ns) {
                return Err(label_err!(
                    "Unknown collection-alias namespace '{}' in label '{}'. This store's collection-alias namespaces: {}",
                    ns,
                    label,
                    list_or_none(&col_namespaces)
                ));
            }
            maps.add_collection_namespace(&ns, self.aliases.namespace_map(AliasKind::Collection, &ns));
        }
        Ok(maps)
    }

    /// Export sequences from BED file regions to a FASTA file.
    ///
    /// Writes one record per BED line, each on a single line, labeled by
    /// `label` with `:start-end` appended (as samtools does). Labels repeat
    /// when regions repeat; no uniqueness rule applies.
    pub fn export_fasta_from_regions<K: AsRef<[u8]>>(
        &self,
        collection_digest: K,
        bed_file_path: &str,
        output_file_path: &str,
        label: &LabelTemplate,
    ) -> Result<ExportSummary> {
        let collection_key = collection_digest.as_ref().to_key();
        let collection_str = key_to_digest_string(&collection_key);
        let maps = self.label_alias_maps(label)?;

        // Label each name once, from this collection's own records. A repeated
        // name keeps its last record, as the region lookup does.
        let mut name_labels: HashMap<&str, (String, bool)> = HashMap::new();
        for rec in self.collection_sequence_metadata(&collection_key)? {
            let meta = rec.metadata();
            let rendered = label.render(&meta.sha512t24u, Some((&collection_str, &meta.name)), &maps);
            check_label_whitespace(&rendered.0)?;
            name_labels.insert(meta.name.as_str(), rendered);
        }

        let seq_iter = self.substrings_from_regions(&collection_digest, bed_file_path)?;
        let mut writer = open_output(Path::new(output_file_path))?;
        let mut summary = ExportSummary::default();

        for rs in seq_iter {
            let rs = rs?;
            let (name_label, fallback) = name_labels
                .get(rs.chrom_name.as_str())
                .ok_or_else(|| anyhow!("Sequence '{}' not found in collection", rs.chrom_name))?;
            writeln!(writer, ">{}:{}-{}", name_label, rs.start, rs.end)?;
            writer.write_all(rs.sequence.as_bytes())?;
            writer.write_all(b"\n")?;
            summary.n_records += 1;
            summary.n_bases += rs.sequence.len() as u64;
            summary.n_fallbacks += *fallback as usize;
        }

        writer.flush()?;
        Ok(summary)
    }

    /// Export sequences to a FASTA file, or to stdout when `output_path` is `"-"`.
    ///
    /// `selection` chooses which records to write; `opts.label` chooses how
    /// each is labeled, the same way whatever the selection. The whole
    /// `(label, sequence)` list is planned from metadata before any output is
    /// written, so label errors never leave half a file. A repeated label for
    /// the same sequence is written once (counted in `n_deduplicated`); a
    /// repeated label for different sequences is an error.
    ///
    /// Record order: collections in the order given (digest order for `All`),
    /// each in its own order; then, for `All`, sequences in no collection, in
    /// digest order. When the label uses no collection field (`{seq}`,
    /// `{seq:NS}`), `All` lists each sequence once in digest order, without
    /// reading any collection.
    ///
    /// Exports stream one sequence at a time and leave the store's loaded
    /// state unchanged. Collections named by `selection` (all of them for
    /// `All` with a collection field) must be loaded; [`RefgetStore::export_fasta`]
    /// loads them.
    pub fn export_fasta<P: AsRef<Path>>(
        &self,
        selection: &ExportSelection,
        output_path: P,
        opts: &ExportOptions,
    ) -> Result<ExportSummary> {
        let line_width = opts.line_width.unwrap_or(80);
        let maps = self.label_alias_maps(&opts.label)?;
        let candidates = self.export_candidates(selection, &opts.label)?;

        // Plan: label every record, apply the uniqueness rule.
        let mut summary = ExportSummary::default();
        let mut seen: HashMap<String, (DigestKey, String)> = HashMap::new();
        let mut plan: Vec<PlannedRecord> = Vec::new();
        for c in candidates {
            if let Some(m) = opts.molecule {
                if !m.keeps(&c.meta.alphabet) {
                    continue;
                }
            }
            let coll = c.collection.as_ref().map(|(col, name)| (col.as_str(), *name));
            let (label, fallback) = opts.label.render(&c.meta.sha512t24u, coll, &maps);
            check_label_whitespace(&label)?;
            match seen.get(&label) {
                Some((key, _)) if *key == c.key => {
                    summary.n_deduplicated += 1;
                    continue;
                }
                Some((key, first_source)) => {
                    return Err(label_err!(
                        "Label '{}' would name two different sequences: {} from {} and {} from {}. \
                         Add {{col:NS}} or {{col}} to the label (for example '{{col}}/{{name}}') to make labels unique.",
                        label,
                        key_to_digest_string(key),
                        first_source,
                        c.meta.sha512t24u,
                        source_of(&c.collection)
                    ));
                }
                None => {}
            }
            summary.n_fallbacks += fallback as usize;
            seen.insert(label.clone(), (c.key, source_of(&c.collection)));
            plan.push(PlannedRecord { label, key: c.key });
        }
        drop(seen);

        // Write: stream each record's bytes, read ahead in parallel.
        let remote = self.remote_source.is_some();
        let mut writer = open_output(output_path.as_ref())?;
        let keys: Vec<DigestKey> = plan.iter().map(|r| r.key).collect();
        let popts = PrefetchOptions {
            order: PrefetchOrder::Input,
            write_cache: opts.write_cache,
            progress: remote && !self.quiet,
        };
        let mut fetch = self.prefetch_keys(&keys, popts, DEFAULT_MEMORY_BUDGET)?;
        let total = plan.len();
        for (i, rec) in plan.iter().enumerate() {
            let item = fetch
                .next()
                .ok_or_else(|| anyhow!("Sequence data for record {} was not delivered", rec.label))??;
            debug_assert_eq!(item.digest.to_key(), rec.key);
            let meta = self
                .sequence_store
                .get(&rec.key)
                .map(|r| r.metadata())
                .ok_or_else(|| anyhow!("Sequence record not found: {}", key_to_digest_string(&rec.key)))?;
            write_fasta_record(
                &mut *writer,
                &rec.label,
                meta,
                item.bytes(),
                self.mode,
                line_width,
                EXPORT_CHUNK_BASES,
            )?;
            drop(item);
            summary.n_records += 1;
            summary.n_bases += meta.length as u64;
            if !self.quiet && ((i + 1) % PROGRESS_EVERY == 0 || i + 1 == total) {
                let fetched = if remote {
                    let c = fetch.counters();
                    format!(", fetched {:.2} GB at {:.1} MB/s", c.bytes_downloaded as f64 / 1e9, c.rate_mb_s())
                } else {
                    String::new()
                };
                eprintln!(
                    "exported {}/{} records ({:.2} GB){}",
                    i + 1,
                    total,
                    summary.n_bases as f64 / 1e9,
                    fetched
                );
            }
        }
        writer.flush()?;
        summary.fetch = fetch.counters();
        drop(fetch);
        if remote && !self.quiet {
            eprintln!("{}", summary.fetch.network_line());
        }
        Ok(summary)
    }

    /// Every record a selection reaches, in output order, before labeling.
    fn export_candidates(
        &self,
        selection: &ExportSelection,
        label: &LabelTemplate,
    ) -> Result<Vec<Candidate<'_>>> {
        let global = |key: &DigestKey| -> Result<&SequenceMetadata> {
            self.sequence_store
                .get(key)
                .map(|r| r.metadata())
                .ok_or_else(|| anyhow!("Sequence record not found for digest: {}", key_to_digest_string(key)))
        };
        let mut out: Vec<Candidate<'_>> = Vec::new();

        match selection {
            ExportSelection::Collections { digests, names } => {
                if digests.is_empty() {
                    return Err(anyhow!("No collections selected for export"));
                }
                for digest in digests {
                    self.collection_candidates(digest, names.as_deref(), &global, &mut out)?;
                }
            }
            ExportSelection::SeqDigests(digests) => {
                for d in digests {
                    let meta = self.get_sequence(strip_sq_prefix(d)).map(|r| r.metadata())?;
                    let key = meta.sha512t24u.to_key();
                    out.push(Candidate { key, meta: global(&key)?, collection: None });
                }
            }
            ExportSelection::All => {
                let mut in_collection: HashSet<DigestKey> = HashSet::new();
                if label.needs_collection() {
                    let mut coll_digests: Vec<String> = self
                        .collections
                        .values()
                        .map(|r| r.metadata().digest.clone())
                        .collect();
                    coll_digests.sort();
                    for digest in &coll_digests {
                        let start = out.len();
                        self.collection_candidates(digest, None, &global, &mut out)?;
                        in_collection.extend(out[start..].iter().map(|c| c.key));
                    }
                }
                let mut rest: Vec<(&String, DigestKey)> = self
                    .sequence_store
                    .iter()
                    .filter(|(k, _)| !in_collection.contains(*k))
                    .map(|(k, r)| (&r.metadata().sha512t24u, *k))
                    .collect();
                rest.sort();
                for (_, key) in rest {
                    out.push(Candidate { key, meta: global(&key)?, collection: None });
                }
            }
        }
        Ok(out)
    }

    /// Append one collection's records (optionally only `names`) to `out`.
    fn collection_candidates<'a>(
        &'a self,
        digest: &str,
        names: Option<&[String]>,
        global: &dyn Fn(&DigestKey) -> Result<&'a SequenceMetadata>,
        out: &mut Vec<Candidate<'a>>,
    ) -> Result<()> {
        let key = digest.to_key();
        let digest_str = key_to_digest_string(&key);
        let records = self.collection_sequence_metadata(&key)?;
        let wanted: Option<HashSet<&str>> = names.map(|ns| ns.iter().map(String::as_str).collect());
        if let Some(names) = names {
            let present: HashSet<&str> = records.iter().map(|r| r.metadata().name.as_str()).collect();
            if let Some(n) = names.iter().find(|n| !present.contains(n.as_str())) {
                return Err(anyhow!("Sequence '{}' not found in collection {}", n, digest_str));
            }
        }
        for rec in records {
            let m = rec.metadata();
            if wanted.as_ref().is_some_and(|w| !w.contains(m.name.as_str())) {
                continue;
            }
            let seq_key = m.sha512t24u.to_key();
            out.push(Candidate {
                key: seq_key,
                meta: global(&seq_key)?,
                collection: Some((digest_str.clone(), m.name.as_str())),
            });
        }
        Ok(())
    }
}

fn list_or_none(items: &[String]) -> String {
    if items.is_empty() {
        "(none)".to_string()
    } else {
        items.join(", ")
    }
}

fn check_label_whitespace(label: &str) -> Result<()> {
    if label.is_empty() {
        return Err(label_err!("A record's label is empty; FASTA IDs need at least one character"));
    }
    if label.chars().any(char::is_whitespace) {
        return Err(label_err!(
            "Label '{}' contains whitespace; FASTA IDs end at the first space. Change the label template or the alias.",
            label
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::store::{ExportOptions, ExportSelection, FastaImportOptions, LabelTemplate, RefgetStore};
    use std::io::Write;

    fn coll(digest: &str) -> ExportSelection {
        ExportSelection::Collections { digests: vec![digest.to_string()], names: None }
    }

    fn name_opts(line_width: usize) -> ExportOptions {
        ExportOptions {
            label: LabelTemplate::parse("name").unwrap(),
            molecule: None,
            line_width: Some(line_width),
            ..Default::default()
        }
    }
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
            .export_fasta(&coll(digest), &out, &name_opts(0))
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
            .export_fasta(&coll(&digest), &out, &name_opts(0))
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
            .export_fasta(&coll(&digest), &out, &name_opts(80))
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
    fn shared_sequence_exports_with_each_collections_own_name() {
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
                vec!["chr1".to_string()],
                "swap={}",
                swap
            );
            assert_eq!(
                export_headers(&readonly, &digest_b, &dir),
                vec!["1".to_string()],
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
    fn shared_sequence_disk_roundtrip_exports_own_name() {
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
            vec!["1".to_string()]
        );
        assert_eq!(
            export_headers(&readonly, &digest_a, &dir),
            vec!["chr1".to_string()]
        );

        let coll_b = readonly.get_collection(&digest_b).unwrap();
        assert_eq!(
            coll_b.sequences[0].metadata().description.as_deref(),
            Some("flybase desc")
        );
    }

    #[test]
    fn get_sequence_by_name_returns_per_collection_metadata() {
        let (fasta_a, fasta_b) = shared_sequence_fastas();

        // Both import orders: the first-imported label must never leak into
        // the other collection's get_sequence_by_name() result.
        for swap in [false, true] {
            let mut store = RefgetStore::in_memory();
            let (digest_a, digest_b) = if swap {
                let (b, a) = import_pair(&mut store, &fasta_b, &fasta_a);
                (a, b)
            } else {
                import_pair(&mut store, &fasta_a, &fasta_b)
            };
            let readonly = store.into_readonly();

            let rec_b = readonly
                .get_sequence_by_name(&digest_b, "1")
                .expect("get sequence '1' from collection B");
            assert_eq!(rec_b.metadata().name, "1", "swap={}", swap);
            assert_eq!(
                rec_b.metadata().description.as_deref(),
                Some("flybase desc"),
                "swap={}",
                swap
            );

            let rec_a = readonly
                .get_sequence_by_name(&digest_a, "chr1")
                .expect("get sequence 'chr1' from collection A");
            assert_eq!(rec_a.metadata().name, "chr1", "swap={}", swap);
            assert_eq!(
                rec_a.metadata().description.as_deref(),
                Some("ucsc desc"),
                "swap={}",
                swap
            );

            // Same underlying bytes, shared across both collections.
            assert_eq!(
                rec_a.metadata().sha512t24u,
                rec_b.metadata().sha512t24u,
                "swap={}",
                swap
            );

            // Names do not leak across collections.
            assert!(
                readonly.get_sequence_by_name(&digest_b, "chr1").is_err(),
                "swap={}",
                swap
            );
        }
    }

    #[test]
    fn get_sequence_by_name_disk_roundtrip() {
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

        let rec_b = readonly
            .get_sequence_by_name(&digest_b, "1")
            .expect("get sequence '1' from collection B");
        assert_eq!(rec_b.metadata().name, "1");
        assert_eq!(rec_b.metadata().description.as_deref(), Some("flybase desc"));
        assert!(rec_b.is_loaded());
        assert_eq!(rec_b.decode().unwrap(), "ACGT".repeat(10));

        let rec_a = readonly
            .get_sequence_by_name(&digest_a, "chr1")
            .expect("get sequence 'chr1' from collection A");
        assert_eq!(rec_a.metadata().name, "chr1");
        assert_eq!(rec_a.metadata().description.as_deref(), Some("ucsc desc"));
        assert!(rec_a.is_loaded());
        assert_eq!(rec_a.decode().unwrap(), "ACGT".repeat(10));
    }
}
