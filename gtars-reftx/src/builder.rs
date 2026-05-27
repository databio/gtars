//! Builder for creating transcript stores from cdot JSON.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{anyhow, Result};
use serde::Deserialize;

use crate::models::{Exon, ManeStatus, Strand, Transcript};
use crate::store::{fnv1a_64, HEADER_SIZE, MAGIC, NONE_SENTINEL, VERSION};

const WRITE_BUFFER_SIZE: usize = 256 * 1024;
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
    record_buf: Vec<u8>,
}

impl TxStoreBuilder {
    /// Create a new builder.
    pub fn new() -> Self {
        Self {
            transcripts: Vec::new(),
            chrom_to_digest: HashMap::new(),
            mane_flags: HashMap::new(),
            record_buf: Vec::with_capacity(512),
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
        let reader: Box<dyn BufRead> =
            if path.as_ref().extension().map_or(false, |e| e == "gz") {
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
            if path.as_ref().extension().map_or(false, |e| e == "gz") {
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

    /// Build and write the binary store to disk.
    pub fn build<P: AsRef<Path>>(&mut self, output: P) -> Result<()> {
        if self.transcripts.is_empty() {
            return Err(anyhow!("No transcripts to write"));
        }

        self.transcripts
            .sort_by_key(|t| fnv1a_64(t.accession.as_bytes()));

        let file = File::create(output)?;
        let mut writer = BufWriter::with_capacity(WRITE_BUFFER_SIZE, file);

        let header_placeholder = [0u8; HEADER_SIZE];
        writer.write_all(&header_placeholder)?;

        let records_start = HEADER_SIZE as u64;
        let mut index_entries: Vec<(u64, u64)> = Vec::with_capacity(self.transcripts.len());
        let mut mane_entries: Vec<(u64, u64)> = Vec::new();
        let mut current_offset = records_start;

        // Take the record buffer out of self so we can simultaneously iterate
        // self.transcripts immutably and mutate the buffer.
        let mut record_buf = std::mem::take(&mut self.record_buf);
        for tx in &self.transcripts {
            let hash = fnv1a_64(tx.accession.as_bytes());
            index_entries.push((hash, current_offset));

            if tx.mane.mane_select {
                let gene_hash = fnv1a_64(tx.gene.to_uppercase().as_bytes());
                mane_entries.push((gene_hash, current_offset));
            }

            serialize_record_into(&mut record_buf, tx);
            writer.write_all(&record_buf)?;
            current_offset += record_buf.len() as u64;
        }
        self.record_buf = record_buf;

        let index_offset = current_offset;
        for (hash, offset) in &index_entries {
            writer.write_all(&hash.to_le_bytes())?;
            writer.write_all(&offset.to_le_bytes())?;
            current_offset += 16;
        }

        // MANE gene index: 0 if empty (sentinel for "no MANE info").
        let mane_index_offset = if mane_entries.is_empty() {
            0u64
        } else {
            mane_entries.sort_by_key(|(h, _)| *h);
            let off = current_offset;
            writer.write_all(&(mane_entries.len() as u64).to_le_bytes())?;
            for (hash, offset) in &mane_entries {
                writer.write_all(&hash.to_le_bytes())?;
                writer.write_all(&offset.to_le_bytes())?;
            }
            off
        };

        writer.flush()?;

        let mut file = writer.into_inner()?;
        file.seek(SeekFrom::Start(0))?;

        let mut header = [0u8; HEADER_SIZE];
        header[0..4].copy_from_slice(MAGIC);
        header[4..8].copy_from_slice(&VERSION.to_le_bytes());
        header[8..16].copy_from_slice(&(self.transcripts.len() as u64).to_le_bytes());
        header[16..24].copy_from_slice(&index_offset.to_le_bytes());
        header[24..32].copy_from_slice(&mane_index_offset.to_le_bytes());
        // bytes [32..40] reserved (already zero)

        file.write_all(&header)?;

        Ok(())
    }

}

/// Serialize a transcript record into a reusable buffer.
fn serialize_record_into(buf: &mut Vec<u8>, tx: &Transcript) {
    buf.clear();

    buf.push(tx.accession.len() as u8);
    buf.extend_from_slice(tx.accession.as_bytes());

    buf.push(tx.gene.len() as u8);
    buf.extend_from_slice(tx.gene.as_bytes());

    buf.extend_from_slice(&tx.chrom_digest);

    buf.push(tx.strand.to_byte());

    // MANE flags byte (v2 format).
    buf.push(tx.mane.to_flags_byte());

    buf.extend_from_slice(&tx.cds_start.unwrap_or(NONE_SENTINEL).to_le_bytes());
    buf.extend_from_slice(&tx.cds_end.unwrap_or(NONE_SENTINEL).to_le_bytes());

    buf.extend_from_slice(&(tx.exons.len() as u16).to_le_bytes());

    for exon in &tx.exons {
        buf.extend_from_slice(&exon.start.to_le_bytes());
        buf.extend_from_slice(&exon.end.to_le_bytes());
    }
}

impl Default for TxStoreBuilder {
    fn default() -> Self {
        Self::new()
    }
}
