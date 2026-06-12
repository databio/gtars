//! Transcript data models.

/// Strand orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i8)]
pub enum Strand {
    Forward = 1,
    Reverse = -1,
}

impl Strand {
    /// Parse from signed integer (+1 or -1).
    #[inline]
    pub fn from_i8(val: i8) -> Option<Self> {
        match val {
            1 => Some(Strand::Forward),
            -1 => Some(Strand::Reverse),
            _ => None,
        }
    }

    /// Byte representation for binary serialization.
    #[inline]
    pub fn to_byte(self) -> u8 {
        self as i8 as u8
    }
}

/// MANE (Matched Annotation from NCBI and EMBL-EBI) status flags.
///
/// Packed into a single byte for compact binary storage:
/// - bit 0: MANE Select (the single canonical transcript per gene)
/// - bit 1: MANE Plus Clinical (additional disease-relevant transcript)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ManeStatus {
    /// True if this is the MANE Select transcript for its gene.
    pub mane_select: bool,
    /// True if this is MANE Plus Clinical.
    pub mane_clinical: bool,
}

impl ManeStatus {
    /// Decode from packed flag byte.
    #[inline]
    pub fn from_flags_byte(byte: u8) -> Self {
        Self {
            mane_select: byte & 0x01 != 0,
            mane_clinical: byte & 0x02 != 0,
        }
    }

    /// Encode to packed flag byte.
    #[inline]
    pub fn to_flags_byte(&self) -> u8 {
        let mut flags = 0u8;
        if self.mane_select {
            flags |= 0x01;
        }
        if self.mane_clinical {
            flags |= 0x02;
        }
        flags
    }

    /// True if either MANE flag is set.
    #[inline]
    pub fn is_mane(&self) -> bool {
        self.mane_select || self.mane_clinical
    }
}

/// A single exon with genomic coordinates (0-based, half-open).
///
/// # Coordinate range limit
///
/// `start`/`end` are `u32`, so genomic coordinates are capped at
/// `u32::MAX` (~4.29 Gb). This is safe for every supported assembly
/// (the longest GRCh38 contig, chr1, is ~248.96 Mb) but would silently
/// truncate on a hypothetical genome larger than ~4.3 Gb. Downstream
/// mapping outputs widen to `u64` (see [`crate::transcripts::mapper`]),
/// so the limit lives only on the stored exon-coordinate side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Exon {
    /// Genomic start (0-based, inclusive).
    pub start: u32,
    /// Genomic end (0-based, exclusive).
    pub end: u32,
}

impl Exon {
    /// Length in bases.
    #[inline]
    pub fn len(&self) -> u32 {
        self.end - self.start
    }

    /// Returns true if zero-length.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }
}

/// Transcript annotation record.
///
/// Exons are stored in genomic order (5' to 3' on chromosome), regardless of strand.
/// For reverse-strand transcripts, transcript coordinates run in reverse genomic order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Transcript {
    /// Accession with version (e.g., "NM_004333.6").
    pub accession: String,
    /// Gene symbol (e.g., "BRAF").
    pub gene: String,
    /// Chromosome refget digest (SQ.xxx, truncated to 24 bytes binary).
    pub chrom_digest: [u8; 24],
    /// Strand orientation.
    pub strand: Strand,
    /// CDS start in genomic coordinates (None if non-coding RNA).
    pub cds_start: Option<u32>,
    /// CDS end in genomic coordinates (None if non-coding RNA).
    pub cds_end: Option<u32>,
    /// Exons in genomic order (always 5'->3' on chromosome).
    pub exons: Vec<Exon>,
    /// MANE Select / Plus Clinical status (defaults to all-false).
    pub mane: ManeStatus,
}

impl Transcript {
    /// Total transcript length in bases (sum of exon lengths).
    pub fn transcript_length(&self) -> u32 {
        self.exons.iter().map(|e| e.len()).sum()
    }

    /// CDS length in bases (0 if non-coding).
    pub fn cds_length(&self) -> u32 {
        match (self.cds_start, self.cds_end) {
            (Some(cds_s), Some(cds_e)) => self
                .exons
                .iter()
                .filter_map(|exon| {
                    // Clip exon to CDS bounds
                    let s = exon.start.max(cds_s);
                    let e = exon.end.min(cds_e);
                    if s < e {
                        Some(e - s)
                    } else {
                        None
                    }
                })
                .sum(),
            _ => 0,
        }
    }

    /// Returns true if this is a coding transcript.
    #[inline]
    pub fn is_coding(&self) -> bool {
        self.cds_start.is_some() && self.cds_end.is_some()
    }

    /// Get the accession without version (e.g., "NM_004333").
    pub fn accession_base(&self) -> &str {
        self.accession.split('.').next().unwrap_or(&self.accession)
    }
}
