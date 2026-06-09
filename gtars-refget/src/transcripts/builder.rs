//! Builder for creating transcript stores from cdot JSON (native-only).
//!
//! Gated behind `filesystem`: writes a `.reftx` file via `std::fs`.

use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
use serde::Deserialize;

use crate::transcripts::models::{Exon, ManeStatus, Strand, Transcript};
use crate::transcripts::store::build_reftx_bytes;

const READ_BUFFER_SIZE: usize = 256 * 1024;

/// cdot JSON transcript format (subset of fields we need).
#[derive(Deserialize)]
struct CdotTranscript {
    id: String,
    gene_name: Option<String>,
    contig: String,
    strand: i8,
    cds_start: Option<u32>,
    cds_end: Option<u32>,
    exons: Vec<(u32, u32)>,
}

/// cdot JSON file structure.
#[derive(Deserialize)]
struct CdotFile {
    transcripts: HashMap<String, CdotTranscript>,
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

    /// Ingest a cdot JSON file.
    ///
    /// Uses 256KB read buffer. Skips transcripts on chromosomes not in
    /// the chrom_to_digest mapping.
    pub fn ingest_cdot<P: AsRef<Path>>(&mut self, path: P) -> Result<usize> {
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

        let mut count = 0;
        for (_, tx) in cdot.transcripts {
            let chrom_digest = match self.chrom_to_digest.get(&tx.contig) {
                Some(d) => *d,
                None => continue,
            };

            let strand = match tx.strand {
                1 => Strand::Forward,
                -1 => Strand::Reverse,
                _ => continue,
            };

            let exons: Vec<Exon> = tx
                .exons
                .into_iter()
                .map(|(s, e)| Exon { start: s, end: e })
                .collect();

            if exons.is_empty() {
                continue;
            }

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
                cds_start: tx.cds_start,
                cds_end: tx.cds_end,
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
    /// Uses refget's write-once + atomic-publish discipline so any reader
    /// (especially the mmap backend, whose `unsafe { Mmap::map }` requires the
    /// file to be immutable) only ever observes a complete file:
    ///
    /// 1. Acquire an advisory build lock on a `<dest>.lock` sidecar (serializes
    ///    BUILDERS only — readers never need it since they open only a fully
    ///    published, immutable file).
    /// 2. Build the full `.reftx` byte image in memory via the shared
    ///    [`build_reftx_bytes`] encoder (the ONE encoder; header already
    ///    populated, no write-then-seek-back patch needed) and write it to a
    ///    temp file in the SAME directory as the destination (so the final
    ///    rename is same-filesystem and atomic).
    /// 3. `sync_all()` to force bytes+metadata to disk BEFORE publishing.
    /// 4. `NamedTempFile::persist` to atomically `rename(2)` onto the
    ///    destination: the dest either is unchanged or flips wholesale to the
    ///    complete file — never a torn intermediate.
    pub fn build<P: AsRef<Path>>(&mut self, output: P) -> Result<()> {
        if self.transcripts.is_empty() {
            return Err(anyhow!("No transcripts to write"));
        }

        let dest = output.as_ref();
        let dest_parent = dest.parent().filter(|p| !p.as_os_str().is_empty());
        let dest_dir: PathBuf = match dest_parent {
            Some(p) => p.to_path_buf(),
            None => PathBuf::from("."),
        };

        // (1) Advisory build lock (serializes concurrent builders only).
        let _lock = BuildLock::acquire(dest)?;

        // (2) Assemble the complete byte image with the shared encoder (header
        // already populated) and write it into a temp file in the dest dir.
        let bytes = build_reftx_bytes(&self.transcripts)?;
        let mut tmp = tempfile::NamedTempFile::new_in(&dest_dir)
            .with_context(|| format!("creating temp file in {:?}", dest_dir))?;
        {
            let file = tmp.as_file_mut();
            file.write_all(&bytes)?;
            // (3) Durability: force bytes+metadata to disk BEFORE publishing.
            file.sync_all()?;
        }

        // (4) Atomically publish via rename(2) onto the destination.
        tmp.persist(dest)
            .map_err(|e| anyhow!("failed to atomically publish {:?}: {}", dest, e.error))?;

        Ok(())
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
