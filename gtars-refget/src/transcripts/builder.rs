//! Builder for creating transcript stores from cdot JSON (native-only).
//!
//! Gated behind `filesystem`: writes a `.reftx` file via `std::fs`.

use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
use serde::Deserialize;

use crate::transcripts::models::{Exon, ManeStatus, Strand, Transcript};
use crate::transcripts::store::build_reftx_bytes;

const READ_BUFFER_SIZE: usize = 256 * 1024;

/// cdot JSON file structure (subset of fields we need).
///
/// Real cdot files (0.2.x) keep per-build coordinates one level down:
/// `transcripts.<id>.genome_builds.<build>.{contig, strand, cds_start, cds_end, exons}`.
/// The top-level `genome_builds` lists every build present in the file.
#[derive(Deserialize)]
struct CdotFile {
    #[serde(default)]
    genome_builds: Vec<String>,
    transcripts: HashMap<String, CdotTranscript>,
}

#[derive(Deserialize)]
struct CdotTranscript {
    id: String,
    gene_name: Option<String>,
    #[serde(default)]
    genome_builds: HashMap<String, CdotBuild>,
}

/// One transcript's alignment to one genome build.
#[derive(Deserialize)]
struct CdotBuild {
    contig: String,
    /// `"+"` or `"-"`.
    strand: String,
    cds_start: Option<u32>,
    cds_end: Option<u32>,
    exons: Vec<CdotExon>,
}

/// A cdot exon: `[alt_start, alt_end, exon_ordinal, tx_start, tx_end, gap]`.
/// Genomic coordinates are 0-based half-open; we keep only those two.
struct CdotExon {
    start: u32,
    end: u32,
}

impl<'de> Deserialize<'de> for CdotExon {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct ExonVisitor;

        impl<'de> serde::de::Visitor<'de> for ExonVisitor {
            type Value = CdotExon;

            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("a cdot exon array starting with [start, end, ...]")
            }

            fn visit_seq<A: serde::de::SeqAccess<'de>>(
                self,
                mut seq: A,
            ) -> Result<CdotExon, A::Error> {
                use serde::de::Error;
                let start = seq
                    .next_element()?
                    .ok_or_else(|| A::Error::invalid_length(0, &self))?;
                let end = seq
                    .next_element()?
                    .ok_or_else(|| A::Error::invalid_length(1, &self))?;
                while seq.next_element::<serde::de::IgnoredAny>()?.is_some() {}
                Ok(CdotExon { start, end })
            }
        }

        deserializer.deserialize_seq(ExonVisitor)
    }
}

/// Builder for creating transcript stores.
pub struct TxStoreBuilder {
    /// Staged transcripts. Public for direct insertion in tests/benches.
    pub transcripts: Vec<Transcript>,
    chrom_to_digest: HashMap<String, [u8; 24]>,
    /// Map from transcript accession (versioned or base) to MANE flags byte.
    mane_flags: HashMap<String, u8>,
}

impl TxStoreBuilder {
    /// Create a new builder.
    pub fn new() -> Self {
        Self {
            transcripts: Vec::new(),
            chrom_to_digest: HashMap::new(),
            mane_flags: HashMap::new(),
        }
    }

    /// Load MANE flags from an NCBI MANE summary TSV (optionally `.gz`).
    ///
    /// Download from:
    /// `https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.vX.summary.txt.gz`
    ///
    /// Indexes both RefSeq (`RefSeq_nuc`) and Ensembl (`Ensembl_nuc`) accessions
    /// so cdot-imported transcripts from either source pick up the flags.
    ///
    /// Returns the number of accession→flags entries inserted.
    pub fn load_mane_summary<P: AsRef<Path>>(&mut self, path: P) -> Result<usize> {
        let file = File::open(path.as_ref())?;
        let reader: Box<dyn BufRead> = if path.as_ref().extension().is_some_and(|e| e == "gz") {
            Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut lines = reader.lines();
        let header = loop {
            let line = lines
                .next()
                .ok_or_else(|| anyhow!("MANE summary is empty"))??;
            if !line.is_empty() {
                break line;
            }
        };

        let columns: Vec<&str> = header.trim_start_matches('#').split('\t').collect();
        let refseq_idx = columns
            .iter()
            .position(|c| c.trim() == "RefSeq_nuc")
            .ok_or_else(|| anyhow!("RefSeq_nuc column not found in MANE summary"))?;
        let ensembl_idx = columns
            .iter()
            .position(|c| c.trim() == "Ensembl_nuc")
            .ok_or_else(|| anyhow!("Ensembl_nuc column not found in MANE summary"))?;
        let status_idx = columns
            .iter()
            .position(|c| c.trim() == "MANE_status")
            .ok_or_else(|| anyhow!("MANE_status column not found in MANE summary"))?;

        let mut count = 0;
        for line in lines {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let max_idx = status_idx.max(refseq_idx).max(ensembl_idx);
            if fields.len() <= max_idx {
                continue;
            }

            let flags = match fields[status_idx].trim() {
                "MANE Select" => 0x01,
                "MANE Plus Clinical" => 0x03,
                _ => continue,
            };

            for &raw in &[fields[refseq_idx], fields[ensembl_idx]] {
                let acc = raw.trim();
                if acc.is_empty() {
                    continue;
                }
                self.mane_flags.insert(acc.to_string(), flags);
                // Also key by base accession (without version) so cdot records
                // with a different version pick up MANE flags.
                if let Some((base, _)) = acc.split_once('.') {
                    self.mane_flags.insert(base.to_string(), flags);
                }
                count += 1;
            }
        }
        Ok(count)
    }

    /// True if any MANE flags have been loaded.
    #[inline]
    pub fn has_mane_flags(&self) -> bool {
        !self.mane_flags.is_empty()
    }

    /// Register chromosome name to refget digest mapping.
    pub fn add_chrom_mapping(&mut self, name: &str, digest: [u8; 24]) {
        self.chrom_to_digest.insert(name.to_string(), digest);
    }

    /// Ingest a cdot JSON file (optionally `.gz`).
    ///
    /// `genome_build` picks which build's coordinates to use (e.g. `"GRCh38"`).
    /// If `None`, the file must contain exactly one build. Transcripts without
    /// an alignment to that build, or on contigs not in the chrom_to_digest
    /// mapping, are skipped. Note that real cdot files name contigs by RefSeq
    /// accession (e.g. `NC_000007.14`), so mappings must use those names.
    pub fn ingest_cdot<P: AsRef<Path>>(
        &mut self,
        path: P,
        genome_build: Option<&str>,
    ) -> Result<usize> {
        let file = File::open(path.as_ref())?;

        let reader: Box<dyn std::io::Read> =
            if path.as_ref().extension().is_some_and(|e| e == "gz") {
                Box::new(flate2::read::GzDecoder::new(BufReader::with_capacity(
                    READ_BUFFER_SIZE,
                    file,
                )))
            } else {
                Box::new(BufReader::with_capacity(READ_BUFFER_SIZE, file))
            };

        let cdot: CdotFile = serde_json::from_reader(reader)?;

        let build = match genome_build {
            Some(b) => {
                if !cdot.genome_builds.iter().any(|x| x == b) {
                    return Err(anyhow!(
                        "genome build {:?} not in cdot file (available: {:?})",
                        b,
                        cdot.genome_builds
                    ));
                }
                b.to_string()
            }
            None => match cdot.genome_builds.as_slice() {
                [only] => only.clone(),
                builds => {
                    return Err(anyhow!(
                        "cdot file has {} genome builds {:?}; pass genome_build to pick one",
                        builds.len(),
                        builds
                    ))
                }
            },
        };

        let mut count = 0;
        for (_, mut tx) in cdot.transcripts {
            let Some(b) = tx.genome_builds.remove(&build) else {
                continue;
            };

            let chrom_digest = match self.chrom_to_digest.get(&b.contig) {
                Some(d) => *d,
                None => continue,
            };

            let strand = match b.strand.as_str() {
                "+" => Strand::Forward,
                "-" => Strand::Reverse,
                _ => continue,
            };

            let mut exons: Vec<Exon> = b
                .exons
                .into_iter()
                .map(|e| Exon { start: e.start, end: e.end })
                .collect();

            if exons.is_empty() {
                continue;
            }
            exons.sort_by_key(|e| e.start);

            let mane_flag_byte = self
                .mane_flags
                .get(&tx.id)
                .copied()
                .or_else(|| {
                    let base = tx.id.split('.').next().unwrap_or(&tx.id);
                    self.mane_flags.get(base).copied()
                })
                .unwrap_or(0);

            self.transcripts.push(Transcript {
                accession: tx.id,
                gene: tx.gene_name.unwrap_or_default(),
                chrom_digest,
                strand,
                cds_start: b.cds_start,
                cds_end: b.cds_end,
                exons,
                mane: ManeStatus::from_flags_byte(mane_flag_byte),
            });
            count += 1;
        }

        Ok(count)
    }

    /// Number of transcripts currently staged.
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    /// Returns true if no transcripts staged.
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }

    /// Build and write the binary `.reftx` store to disk atomically.
    ///
    /// Any reader (especially the mmap backend, whose `unsafe { Mmap::map }`
    /// requires the file to be immutable) must only ever observe a complete file:
    ///
    /// 1. Acquire an advisory build lock on a `<dest>.lock` sidecar (serializes
    ///    BUILDERS only — readers never need it since they open only a fully
    ///    published, immutable file).
    /// 2. Build the full `.reftx` byte image in memory via the shared
    ///    [`build_reftx_bytes`] encoder (the ONE encoder; header already
    ///    populated, no write-then-seek-back patch needed).
    /// 3. Publish it through [`crate::store::atomic::atomic_write`], the ONE
    ///    implementation of the temp-file → `fsync` → `rename(2)` → `fsync` dir
    ///    discipline, shared with the store index/manifest/alias writes.
    pub fn build<P: AsRef<Path>>(&mut self, output: P) -> Result<()> {
        if self.transcripts.is_empty() {
            return Err(anyhow!("No transcripts to write"));
        }

        let dest = output.as_ref();

        // (1) Advisory build lock (serializes concurrent builders only).
        let _lock = BuildLock::acquire(dest)?;

        // (2) Assemble the complete byte image with the shared encoder.
        let bytes = build_reftx_bytes(&self.transcripts)?;

        // (3) Atomic publish.
        crate::store::atomic::atomic_write_bytes(dest, &bytes)
            .with_context(|| format!("failed to atomically publish {:?}", dest))
    }
}

/// Advisory build lock backed by a create-exclusive `<dest>.lock` sidecar.
///
/// Portable across platforms (no `flock`/`fs2` dependency): the lock is held
/// for the lifetime of the guard and the sidecar is removed on drop. Serializes
/// only concurrent BUILDERS — readers never acquire it because they open only a
/// fully-published, immutable file.
struct BuildLock {
    path: PathBuf,
}

impl BuildLock {
    fn acquire(dest: &Path) -> Result<Self> {
        let mut path = dest.as_os_str().to_os_string();
        path.push(".lock");
        let path = PathBuf::from(path);
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .with_context(|| {
                format!(
                    "acquiring build lock {:?} (another builder may be running)",
                    path
                )
            })?;
        Ok(Self { path })
    }
}

impl Drop for BuildLock {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

impl Default for TxStoreBuilder {
    fn default() -> Self {
        Self::new()
    }
}
