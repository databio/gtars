use crate::errors::GtarsGenomicDistError;
use bio::io::fasta;
use gtars_core::models::{CoordinateMode, Region, RegionSet};
use memmap2::Mmap;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// A `RegionSet` that is guaranteed to be sorted by (chr, start).
///
/// Created by moving a `RegionSet` into `SortedRegionSet::new()`, which
/// sorts in place (no clone). Functions that require sorted input can
/// accept this type instead of re-sorting every time.
pub struct SortedRegionSet(pub RegionSet);

impl SortedRegionSet {
    /// Sort a RegionSet in place and wrap it. This is a move, not a clone.
    pub fn new(mut rs: RegionSet) -> Self {
        rs.sort();
        Self(rs)
    }
}

/// Genomic strand orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Strand {
    Plus,
    Minus,
    Unstranded,
}

impl Strand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => Strand::Plus,
            '-' => Strand::Minus,
            _ => Strand::Unstranded,
        }
    }
}

/// A `RegionSet` paired with a parallel `Vec<Strand>`.
///
/// Follows the `SortedRegionSet` wrapper pattern: wraps a `RegionSet` and
/// adds strand information without modifying `Region` itself (which lives
/// in gtars-core). Strand-aware operations (promoters, reduce, setdiff)
/// are implemented as methods on this type.
#[derive(Clone, Serialize, Deserialize)]
pub struct StrandedRegionSet {
    pub inner: RegionSet,
    pub strands: Vec<Strand>,
}

impl StrandedRegionSet {
    /// Create a new StrandedRegionSet from a RegionSet and parallel strand vector.
    ///
    /// # Panics
    /// Panics if `strands.len() != rs.regions.len()`.
    pub fn new(rs: RegionSet, strands: Vec<Strand>) -> Self {
        assert_eq!(
            rs.regions.len(),
            strands.len(),
            "StrandedRegionSet: regions and strands must have the same length"
        );
        StrandedRegionSet {
            inner: rs,
            strands,
        }
    }

    /// Wrap a RegionSet with all-Unstranded. Preserves existing behavior.
    pub fn unstranded(rs: RegionSet) -> Self {
        let n = rs.regions.len();
        StrandedRegionSet {
            strands: vec![Strand::Unstranded; n],
            inner: rs,
        }
    }

    /// Consume into the inner RegionSet, dropping strand information.
    pub fn into_regionset(self) -> RegionSet {
        self.inner
    }

    pub fn len(&self) -> usize {
        self.inner.regions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.regions.is_empty()
    }
}

/// Statistics summary for regions on a single chromosome.
///
/// Contains counts, bounds, and descriptive statistics for region lengths.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromosomeStatistics {
    /// Chromosome name
    pub chromosome: String,
    /// Total number of regions on this chromosome
    pub number_of_regions: u32,
    /// Leftmost start position across all regions
    pub start_nucleotide_position: u32,
    /// Rightmost end position across all regions
    pub end_nucleotide_position: u32,
    /// Length of the shortest region
    pub minimum_region_length: u32,
    /// Length of the longest region
    pub maximum_region_length: u32,
    /// Average region length
    pub mean_region_length: f64,
    /// Median region length
    pub median_region_length: f64,
}

/// A genomic bin with a count of overlapping regions.
///
/// Used to represent distribution of regions across fixed-size windows.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegionBin {
    /// Chromosome name
    pub chr: String,
    /// Start position of the bin
    pub start: u32,
    /// End position of the bin
    pub end: u32,
    /// Number of regions overlapping this bin
    pub n: u32,
    /// Rid: needed for plot to have correct order
    pub rid: u32,
}

/// Summary statistics over the distribution of inter-region spacings.
///
/// "Spacing" here means the bp distance between the end of one region and the
/// start of the next region on the same chromosome, counting only positive
/// gaps (overlapping / abutting neighbors are excluded, matching
/// `calc_neighbor_distances`). Cross-chromosome pairs are never counted.
///
/// Empty / singleton inputs return `n_gaps = 0` and NaN for all float fields,
/// matching numpy's behavior on empty arrays.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpacingStats {
    /// Number of positive inter-region gaps (excludes overlaps and abutting).
    pub n_gaps: usize,
    /// Mean gap length in bp.
    pub mean: f64,
    /// Median gap length in bp.
    pub median: f64,
    /// Population standard deviation of gap lengths in bp.
    pub std: f64,
    /// Interquartile range (Q3 − Q1) of gap lengths in bp.
    pub iqr: f64,
    /// Mean of `log10(gap + 1)`. Gap distributions are heavy-tailed;
    /// the log-transformed mean is a more robust central tendency.
    pub log_mean: f64,
    /// Population standard deviation of `log10(gap + 1)`.
    pub log_std: f64,
}

/// Summary statistics over peak clusters at a given stitching radius.
///
/// A cluster is a maximal set of regions connected via single-linkage,
/// where two regions link if the bp distance between `prev.end` and
/// `next.start` is at most `radius_bp`. Clusters are chromosome-scoped
/// (two regions on different chromosomes can never link).
///
/// # The `min_cluster_size` filter applies uniformly
///
/// Every size-dependent field except `max_cluster_size` is restricted
/// to clusters with size ≥ the `min_cluster_size` parameter passed to
/// `calc_peak_clusters`. The same threshold drives `n_clusters`,
/// `n_clustered_peaks`, `mean_cluster_size`, and `fraction_clustered`
/// so they always answer the same question about the same subset of
/// clusters.
///
/// **Default is `min_cluster_size = 2`** — the default ClusterStats
/// answers "how clustered are my peaks, counting only groups of at
/// least 2?". This matches the scientifically meaningful use case
/// (enhancer clustering, super-enhancer stitching) and makes the
/// arithmetic identity `n_clusters * mean_cluster_size ==
/// n_clustered_peaks` hold at the default.
///
/// **Pass `min_cluster_size = 1`** to get the "every connected
/// component including singletons" view: `n_clusters` then counts all
/// clusters, `n_clustered_peaks == total_peaks`,
/// `fraction_clustered == 1.0` trivially, and `mean_cluster_size` is
/// the simple average `total_peaks / n_clusters`. Useful when you want
/// the simple mean, but most size-related fields become tautological
/// at this threshold.
///
/// `max_cluster_size` is always the largest cluster in the input
/// regardless of filter — the maximum is inherent and unaffected by
/// the threshold.
///
/// Empty inputs return zero counts and NaN for `mean_cluster_size` and
/// `fraction_clustered`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClusterStats {
    /// Stitching radius in bp (pass-through from the call).
    pub radius_bp: u32,
    /// Number of clusters with size ≥ `min_cluster_size`. With the
    /// default `min_cluster_size = 2`, this counts only multi-peak
    /// clusters and excludes singletons.
    pub n_clusters: usize,
    /// Number of peaks belonging to a cluster with size ≥
    /// `min_cluster_size`. With default `min_cluster_size = 2`, this is
    /// "peaks with at least one neighbor within `radius_bp`". With
    /// `min_cluster_size = 1`, this degenerates to `total_peaks`.
    pub n_clustered_peaks: usize,
    /// Mean size of clusters with size ≥ `min_cluster_size`. The
    /// identity `n_clusters * mean_cluster_size == n_clustered_peaks`
    /// holds exactly. NaN when no clusters meet the threshold (empty
    /// input, or `min_cluster_size > max_cluster_size`).
    pub mean_cluster_size: f64,
    /// Size of the largest cluster in the input, regardless of
    /// `min_cluster_size`. 0 for empty input.
    pub max_cluster_size: usize,
    /// `n_clustered_peaks / total_peaks`, where `total_peaks` is the
    /// **raw input count** (not filtered). With default
    /// `min_cluster_size = 2`, this is the fraction of peaks in
    /// multi-peak clusters. NaN if input is empty.
    pub fraction_clustered: f64,
}

/// Dense per-window peak-count vector and the binning that produced it.
///
/// Unlike `region_distribution_with_chrom_sizes`, which returns only bins
/// that contain ≥1 region, this struct carries the full zero-padded vector
/// with one entry per window on every chromosome in `chrom_sizes`.
///
/// `counts` is ordered by karyotypic chromosome order (`chrom_karyotype_key`)
/// and then by bin index within each chromosome. `chrom_offsets` records the
/// start index of each chromosome's slice in `counts`, so callers can
/// recover per-chromosome subvectors without re-binning.
///
/// # `n_bins` is a target, not the total
///
/// `n_bins` is the **target bin count for the longest chromosome in
/// `chrom_sizes`**, not the length of `counts`. Bin width is derived as
/// `max(chrom_sizes.values()) / n_bins` (floored, minimum 1 bp), and every
/// chromosome is tiled with windows of that width. The length of `counts`
/// is `sum(ceil(chrom_size / bin_width))` over all chromosomes in
/// `chrom_sizes`, which can substantially exceed `n_bins` when many
/// chromosomes are present. To target a specific bin width in bp, set
/// `n_bins` to `max_chrom_len / desired_bin_width_bp`.
///
/// # Per-chromosome bin width
///
/// The last bin on each chromosome is narrower than `bin_width` whenever
/// `chrom_size` is not an exact multiple of `bin_width`. Chromosomes
/// shorter than `bin_width` (common with UCSC alt / random / unplaced
/// contigs) reduce to a single bin whose effective width equals the
/// chromosome size rather than `bin_width`. Individual entries in
/// `counts` are therefore counts per bin, not counts per `bin_width` bp —
/// bins of different effective widths are not directly comparable as
/// densities when `chrom_sizes` contains contigs significantly shorter
/// than `bin_width`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DensityVector {
    /// Target bin count for the longest chromosome in `chrom_sizes`.
    /// Shorter chromosomes get proportionally fewer bins. **Not** the
    /// total length of `counts`.
    pub n_bins: u32,
    /// Bin width in bp, computed as `max(chrom_sizes) / n_bins` (floored,
    /// minimum 1 to avoid zero-sized bins). The last bin on each
    /// chromosome and all bins on chromosomes shorter than `bin_width`
    /// have narrower effective widths — see the struct-level docs.
    pub bin_width: u32,
    /// Dense zero-padded per-window peak counts. Length ==
    /// `sum(ceil(chrom_size / bin_width))` over all chromosomes in
    /// `chrom_sizes`.
    pub counts: Vec<u32>,
    /// `(chr, start_index)` per chromosome, in karyotypic order. Slice
    /// `counts[start_index .. next_chrom_start_index]` for per-chromosome
    /// vectors.
    pub chrom_offsets: Vec<(String, usize)>,
}

/// Summary of how evenly peaks are distributed across genome windows.
///
/// Derived from a dense per-window count vector (see `DensityVector`).
/// A Poisson-distributed peak set has `cv ≈ 1`; clustered sets have
/// `cv >> 1`; evenly-spread sets have `cv << 1`.
///
/// See `DensityVector` for the definition of `n_bins` (it is the target
/// bin count for the longest chromosome, not the total window count)
/// and for the treatment of chromosomes shorter than the derived
/// `bin_width`. Both affect the interpretation of the statistics below —
/// short contigs each contribute a narrow single-bin entry which dilutes
/// `mean_count`, inflates `n_windows`, and raises `gini`.
///
/// **Note on Gini:** the Gini coefficient has a known bias toward high
/// values for sparse count distributions (many zero-count windows). For
/// typical TF-binding BED files over hg38 at 1 Mb bins the bias is small,
/// but for very sparse inputs (hundreds of peaks over thousands of bins)
/// the reported Gini will exaggerate concentration. `n_nonzero_windows`
/// is reported alongside so callers can detect the regime.
///
/// Empty-input convention (no chrom_sizes or zero regions in a populated
/// `chrom_sizes`): `n_windows` may still be > 0 (empty RegionSet over
/// populated chrom_sizes means all windows are zero); `mean = 0`,
/// `variance = 0`, `cv = NaN`, `gini = 0`. If `chrom_sizes` itself is
/// empty, all fields are zero or NaN.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DensityHomogeneity {
    /// Bin width in bp used for this summary. The last bin on each
    /// chromosome (and all bins on chromosomes shorter than this value)
    /// is narrower; counts are per-bin, not per-bp.
    pub bin_width: u32,
    /// Total number of bins (including zero-count bins). Each chromosome
    /// contributes `ceil(chrom_size / bin_width)` bins, so this can be
    /// substantially larger than the `n_bins` parameter when
    /// `chrom_sizes` has many entries.
    pub n_windows: usize,
    /// Number of bins containing at least one peak.
    pub n_nonzero_windows: usize,
    /// Mean peak count per window.
    pub mean_count: f64,
    /// Population variance of per-window counts.
    pub variance: f64,
    /// Coefficient of variation, `sqrt(variance) / mean_count`.
    /// NaN if `mean_count == 0`.
    pub cv: f64,
    /// Gini coefficient of per-window counts, in `[0, 1]`.
    /// `0` = perfectly even, `1` = all peaks in a single bin.
    pub gini: f64,
}

/// Trait for types that provide sequence access to a reference genome.
///
/// Implemented by [`GenomeAssembly`] (in-memory HashMap) and
/// [`BinaryGenomeAssembly`] (mmap .fab binary). Functions like
/// `calc_gc_content` accept `&impl SequenceAccess` to work with either.
pub trait SequenceAccess {
    /// Get the sequence for a genomic region. Returns owned bytes.
    fn get_sequence(&self, coords: &Region) -> Result<Vec<u8>, GtarsGenomicDistError>;

    /// Check whether a chromosome exists in the assembly.
    fn contains_chr(&self, chr: &str) -> bool;
}

/// In-memory genome assembly backed by a HashMap of chromosome sequences.
///
/// Loads the entire FASTA into memory on construction. Slower to construct
/// (~2s for hg38) but provides zero-copy `&[u8]` access to sequences,
/// making per-region operations like GC content highly vectorizable.
/// No `.fai` index required.
pub struct GenomeAssembly {
    seq_map: HashMap<String, Vec<u8>>,
}

impl TryFrom<&str> for GenomeAssembly {
    type Error = GtarsGenomicDistError;
    fn try_from(value: &str) -> Result<Self, GtarsGenomicDistError> {
        GenomeAssembly::try_from(Path::new(value))
    }
}

impl TryFrom<String> for GenomeAssembly {
    type Error = GtarsGenomicDistError;
    fn try_from(value: String) -> Result<Self, GtarsGenomicDistError> {
        GenomeAssembly::try_from(Path::new(&value))
    }
}

impl TryFrom<&Path> for GenomeAssembly {
    type Error = GtarsGenomicDistError;

    fn try_from(value: &Path) -> Result<GenomeAssembly, GtarsGenomicDistError> {
        let file = File::open(value)?;
        let genome = fasta::Reader::new(file);
        let records = genome.records();

        let mut seq_map: HashMap<String, Vec<u8>> = HashMap::new();
        for record in records {
            match record {
                Ok(record) => {
                    seq_map.insert(record.id().to_string(), record.seq().to_owned());
                }
                Err(e) => {
                    return Err(GtarsGenomicDistError::CustomError(format!(
                        "Error reading genome file: {}",
                        e
                    )));
                }
            }
        }
        Ok(GenomeAssembly { seq_map })
    }
}

impl GenomeAssembly {
    pub fn seq_from_region(&self, coords: &Region) -> Result<&[u8], GtarsGenomicDistError> {
        let chr = &coords.chr;
        let start = coords.start as usize;
        let end = coords.end as usize;

        if let Some(seq) = self.seq_map.get(chr) {
            if end <= seq.len() && start <= end {
                Ok(&seq[start..end])
            } else {
                Err(GtarsGenomicDistError::CustomError(format!(
                    "Invalid range: start={}, end={} for chromosome {} with length {}",
                    start, end, chr, seq.len()
                )))
            }
        } else {
            Err(GtarsGenomicDistError::CustomError(format!(
                "Unknown chromosome found in region set: {}",
                chr
            )))
        }
    }

    pub fn contains_chr(&self, chr: &str) -> bool {
        self.seq_map.contains_key(chr)
    }
}

impl SequenceAccess for GenomeAssembly {
    fn get_sequence(&self, coords: &Region) -> Result<Vec<u8>, GtarsGenomicDistError> {
        self.seq_from_region(coords).map(|s| s.to_vec())
    }

    fn contains_chr(&self, chr: &str) -> bool {
        self.seq_map.contains_key(chr)
    }
}

// --- Binary FASTA (.fab) format ---

const FAB_MAGIC: &[u8; 4] = b"GFAB";
const FAB_VERSION: u8 = 1;

/// Memory-mapped genome assembly backed by a binary FASTA (.fab) file.
///
/// The .fab format stores sequences contiguously without line wrapping,
/// enabling mmap + zero-copy `&[u8]` access with both instant construction
/// and zero-copy per-region performance.
///
/// Create .fab files with [`BinaryGenomeAssembly::write_from_fasta`] or
/// via `gtars prep --fasta <file>`.
#[derive(Debug)]
pub struct BinaryGenomeAssembly {
    mmap: Mmap,
    /// name → (offset_in_file, sequence_length)
    index: HashMap<String, (usize, usize)>,
}

impl BinaryGenomeAssembly {
    /// Open a .fab binary FASTA file via mmap.
    pub fn from_file(path: &Path) -> Result<Self, GtarsGenomicDistError> {
        let file = File::open(path).map_err(|e| {
            GtarsGenomicDistError::CustomError(format!(
                "Failed to open .fab file '{}': {}",
                path.display(), e
            ))
        })?;
        let mmap = unsafe { Mmap::map(&file) }.map_err(|e| {
            GtarsGenomicDistError::CustomError(format!(
                "Failed to mmap .fab file '{}': {}",
                path.display(), e
            ))
        })?;

        // Parse header
        if mmap.len() < 9 {
            return Err(GtarsGenomicDistError::CustomError(
                "Invalid .fab file: too short".into(),
            ));
        }
        if &mmap[0..4] != FAB_MAGIC {
            return Err(GtarsGenomicDistError::CustomError(
                "Invalid .fab file: bad magic bytes".into(),
            ));
        }
        let version = mmap[4];
        if version != FAB_VERSION {
            return Err(GtarsGenomicDistError::CustomError(format!(
                "Unsupported .fab version: {} (expected {})",
                version, FAB_VERSION
            )));
        }
        let n_chroms = u32::from_le_bytes(mmap[5..9].try_into().unwrap()) as usize;

        // Parse index
        let mut pos = 9;
        let mut index = HashMap::with_capacity(n_chroms);
        for _ in 0..n_chroms {
            if pos + 2 > mmap.len() {
                return Err(GtarsGenomicDistError::CustomError(
                    "Invalid .fab file: truncated index".into(),
                ));
            }
            let name_len = u16::from_le_bytes(mmap[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            if pos + name_len + 16 > mmap.len() {
                return Err(GtarsGenomicDistError::CustomError(
                    "Invalid .fab file: truncated index entry".into(),
                ));
            }
            let name = std::str::from_utf8(&mmap[pos..pos + name_len])
                .map_err(|e| {
                    GtarsGenomicDistError::CustomError(format!(
                        "Invalid .fab file: non-UTF8 chromosome name: {}",
                        e
                    ))
                })?
                .to_string();
            pos += name_len;
            let offset =
                u64::from_le_bytes(mmap[pos..pos + 8].try_into().unwrap()) as usize;
            pos += 8;
            let length =
                u64::from_le_bytes(mmap[pos..pos + 8].try_into().unwrap()) as usize;
            pos += 8;
            index.insert(name, (offset, length));
        }

        Ok(BinaryGenomeAssembly { mmap, index })
    }

    /// Get the sequence for a region as a zero-copy `&[u8]` slice.
    pub fn seq_from_region(&self, coords: &Region) -> Result<&[u8], GtarsGenomicDistError> {
        let chr = &coords.chr;
        let start = coords.start as usize;
        let end = coords.end as usize;

        let &(offset, length) = self.index.get(chr).ok_or_else(|| {
            GtarsGenomicDistError::CustomError(format!(
                "Unknown chromosome found in region set: {}",
                chr
            ))
        })?;

        if end > length || start > end {
            return Err(GtarsGenomicDistError::CustomError(format!(
                "Invalid range: start={}, end={} for chromosome {} with length {}",
                start, end, chr, length
            )));
        }

        let file_start = offset + start;
        let file_end = offset + end;
        if file_end > self.mmap.len() {
            return Err(GtarsGenomicDistError::CustomError(format!(
                "Corrupted .fab file: sequence data for {} extends beyond file boundary",
                chr
            )));
        }

        Ok(&self.mmap[file_start..file_end])
    }

    pub fn contains_chr(&self, chr: &str) -> bool {
        self.index.contains_key(chr)
    }

    /// Convert a FASTA file to .fab binary format.
    pub fn write_from_fasta(
        fasta_path: &Path,
        output_path: &Path,
    ) -> Result<(), GtarsGenomicDistError> {
        // Read all sequences into memory (same as GenomeAssembly)
        let file = File::open(fasta_path)?;
        let reader = fasta::Reader::new(file);

        let mut chroms: Vec<(String, Vec<u8>)> = Vec::new();
        for record in reader.records() {
            let record = record.map_err(|e| {
                GtarsGenomicDistError::CustomError(format!(
                    "Error reading FASTA: {}", e
                ))
            })?;
            chroms.push((record.id().to_string(), record.seq().to_owned()));
        }

        // Compute index: header size first
        let mut header_size: usize = 4 + 1 + 4; // magic + version + n_chroms
        for (name, _) in &chroms {
            header_size += 2 + name.len() + 8 + 8; // name_len + name + offset + length
        }

        // Write
        let out = File::create(output_path).map_err(|e| {
            GtarsGenomicDistError::CustomError(format!(
                "Failed to create .fab file '{}': {}",
                output_path.display(), e
            ))
        })?;
        let mut w = BufWriter::new(out);

        // Header
        w.write_all(FAB_MAGIC)?;
        w.write_all(&[FAB_VERSION])?;
        w.write_all(&(chroms.len() as u32).to_le_bytes())?;

        // Index
        let mut offset = header_size;
        for (name, seq) in &chroms {
            w.write_all(&(name.len() as u16).to_le_bytes())?;
            w.write_all(name.as_bytes())?;
            w.write_all(&(offset as u64).to_le_bytes())?;
            w.write_all(&(seq.len() as u64).to_le_bytes())?;
            offset += seq.len();
        }

        // Sequences
        for (_, seq) in &chroms {
            w.write_all(seq)?;
        }

        w.flush()?;
        Ok(())
    }
}

impl TryFrom<&str> for BinaryGenomeAssembly {
    type Error = GtarsGenomicDistError;
    fn try_from(value: &str) -> Result<Self, GtarsGenomicDistError> {
        BinaryGenomeAssembly::from_file(Path::new(value))
    }
}

impl TryFrom<String> for BinaryGenomeAssembly {
    type Error = GtarsGenomicDistError;
    fn try_from(value: String) -> Result<Self, GtarsGenomicDistError> {
        BinaryGenomeAssembly::from_file(Path::new(&value))
    }
}

impl TryFrom<&Path> for BinaryGenomeAssembly {
    type Error = GtarsGenomicDistError;
    fn try_from(value: &Path) -> Result<Self, GtarsGenomicDistError> {
        BinaryGenomeAssembly::from_file(value)
    }
}

impl SequenceAccess for BinaryGenomeAssembly {
    fn get_sequence(&self, coords: &Region) -> Result<Vec<u8>, GtarsGenomicDistError> {
        self.seq_from_region(coords).map(|s| s.to_vec())
    }

    fn contains_chr(&self, chr: &str) -> bool {
        self.index.contains_key(chr)
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub enum Dinucleotide {
    Aa,
    Ac,
    Ag,
    At,
    Ca,
    Cc,
    Cg,
    Ct,
    Ga,
    Gc,
    Gg,
    Gt,
    Ta,
    Tc,
    Tg,
    Tt,
}

impl Dinucleotide {
    pub fn from_bytes(bytes: &[u8]) -> Option<Dinucleotide> {
        if bytes.len() != 2 {
            return None;
        }
        // Normalize to uppercase for case-insensitive matching
        let normalized = [bytes[0].to_ascii_uppercase(), bytes[1].to_ascii_uppercase()];
        match &normalized {
            b"AA" => Some(Dinucleotide::Aa),
            b"AC" => Some(Dinucleotide::Ac),
            b"AG" => Some(Dinucleotide::Ag),
            b"AT" => Some(Dinucleotide::At),
            b"CA" => Some(Dinucleotide::Ca),
            b"CC" => Some(Dinucleotide::Cc),
            b"CG" => Some(Dinucleotide::Cg),
            b"CT" => Some(Dinucleotide::Ct),
            b"GA" => Some(Dinucleotide::Ga),
            b"GC" => Some(Dinucleotide::Gc),
            b"GG" => Some(Dinucleotide::Gg),
            b"GT" => Some(Dinucleotide::Gt),
            b"TA" => Some(Dinucleotide::Ta),
            b"TC" => Some(Dinucleotide::Tc),
            b"TG" => Some(Dinucleotide::Tg),
            b"TT" => Some(Dinucleotide::Tt),
            _ => None,
        }
    }

    pub fn to_string(&self) -> Result<String, GtarsGenomicDistError> {
        match self {
            Dinucleotide::Aa => Ok("Aa".to_string()),
            Dinucleotide::Ac => Ok("Ac".to_string()),
            Dinucleotide::Ag => Ok("Ag".to_string()),
            Dinucleotide::At => Ok("At".to_string()),
            Dinucleotide::Ca => Ok("Ca".to_string()),
            Dinucleotide::Cc => Ok("Cc".to_string()),
            Dinucleotide::Cg => Ok("Cg".to_string()),
            Dinucleotide::Ct => Ok("Ct".to_string()),
            Dinucleotide::Ga => Ok("Ga".to_string()),
            Dinucleotide::Gc => Ok("Gc".to_string()),
            Dinucleotide::Gg => Ok("Gg".to_string()),
            Dinucleotide::Gt => Ok("Gt".to_string()),
            Dinucleotide::Ta => Ok("Ta".to_string()),
            Dinucleotide::Tc => Ok("Tc".to_string()),
            Dinucleotide::Tg => Ok("Tg".to_string()),
            Dinucleotide::Tt => Ok("Tt".to_string()),
        }
    }
}

///
/// Struct to hold Tss information (RegionSet with additionally indexing) that is initialized from
/// RegionSet or BED file that holds tss regions
///
pub struct TssIndex {
    pub region_set: RegionSet,
    pub mid_points: HashMap<String, Vec<u32>>,
}

impl TryFrom<RegionSet> for TssIndex {
    type Error = GtarsGenomicDistError;
    fn try_from(value: RegionSet) -> Result<Self, GtarsGenomicDistError> {
        TssIndex::from_region_set(value, CoordinateMode::Bed)
    }
}

impl TssIndex {
    /// Create a TssIndex from a RegionSet using the specified coordinate mode for midpoints.
    pub fn from_region_set(
        value: RegionSet,
        mode: CoordinateMode,
    ) -> Result<Self, GtarsGenomicDistError> {
        let mut mid_points = value.calc_mid_points_with_mode(mode);

        for points in mid_points.values_mut() {
            points.sort_unstable();
        }

        Ok(TssIndex {
            region_set: value,
            mid_points,
        })
    }
}

impl TryFrom<&Path> for TssIndex {
    type Error = GtarsGenomicDistError;
    fn try_from(value: &Path) -> Result<Self, GtarsGenomicDistError> {
        let region_set = match RegionSet::try_from(value) {
            Ok(region_set) => region_set,
            Err(_e) => {
                return Err(GtarsGenomicDistError::TSSContentError(String::from(
                    "Unable to open Tss file",
                )));
            }
        };
        TssIndex::try_from(region_set)
    }
}

impl TryFrom<&str> for TssIndex {
    type Error = GtarsGenomicDistError;
    fn try_from(value: &str) -> Result<Self, GtarsGenomicDistError> {
        let region_set = match RegionSet::try_from(value) {
            Ok(region_set) => region_set,
            Err(_e) => {
                return Err(GtarsGenomicDistError::TSSContentError(String::from(
                    "Unable to open Tss file",
                )));
            }
        };
        TssIndex::try_from(region_set)
    }
}

impl TryFrom<String> for TssIndex {
    type Error = GtarsGenomicDistError;
    fn try_from(value: String) -> Result<Self, GtarsGenomicDistError> {
        let region_set = match RegionSet::try_from(value) {
            Ok(region_set) => region_set,
            Err(_e) => {
                return Err(GtarsGenomicDistError::TSSContentError(String::from(
                    "Unable to open Tss file",
                )));
            }
        };
        TssIndex::try_from(region_set)
    }
}

impl TssIndex {
    ///
    /// Calculate the distance from each region to the nearest TSS mid-point.
    ///
    /// Uses binary search for O(R * log M) complexity instead of O(R * M),
    /// where R is the number of regions and M is the number of TSS midpoints.
    ///
    pub fn calc_tss_distances(
        &self,
        rs: &RegionSet,
        mode: CoordinateMode,
    ) -> Result<Vec<u32>, GtarsGenomicDistError> {
        let mut distances: Vec<u32> = Vec::with_capacity(rs.len());

        for chromosome in rs.iter_chroms() {
            if let Some(chr_midpoints) = self.mid_points.get(chromosome.as_str()) {
                for region in rs.iter_chr_regions(chromosome.as_str()) {
                    let target = region.mid_point_with_mode(mode);

                    let min_distance = match chr_midpoints.binary_search(&target) {
                        Ok(_) => 0,
                        Err(idx) => {
                            let left = idx
                                .checked_sub(1)
                                .map(|i| target.abs_diff(chr_midpoints[i]));
                            let right = chr_midpoints.get(idx).map(|&v| target.abs_diff(v));

                            match (left, right) {
                                (Some(l), Some(r)) => l.min(r),
                                (Some(l), None) => l,
                                (None, Some(r)) => r,
                                (None, None) => continue,
                            }
                        }
                    };
                    distances.push(min_distance);
                }
            } else {
                // No features on this chromosome — push u32::MAX for each region
                for _ in rs.iter_chr_regions(chromosome.as_str()) {
                    distances.push(u32::MAX);
                }
            }
        }
        Ok(distances)
    }

    /// Calculate signed distances from each region to its nearest feature.
    ///
    /// Like `calc_tss_distances` but returns signed distances where:
    /// - Positive: nearest feature is downstream (right) of the query
    /// - Negative: nearest feature is upstream (left) of the query
    ///
    /// Sign convention: `nearest_feature_midpoint - query_midpoint`
    /// (matches R GenomicDistributions `calcFeatureDist()`).
    pub fn calc_feature_distances(
        &self,
        rs: &RegionSet,
        mode: CoordinateMode,
    ) -> Result<Vec<i64>, GtarsGenomicDistError> {
        let mut distances: Vec<i64> = Vec::with_capacity(rs.len());

        for chromosome in rs.iter_chroms() {
            if let Some(chr_midpoints) = self.mid_points.get(chromosome.as_str()) {
                for region in rs.iter_chr_regions(chromosome.as_str()) {
                    let target = region.mid_point_with_mode(mode) as i64;

                    let distance = match chr_midpoints.binary_search(&(target as u32)) {
                        Ok(_) => 0i64,
                        Err(idx) => {
                            // distance = feature_mid - query_mid
                            let left = idx
                                .checked_sub(1)
                                .map(|i| chr_midpoints[i] as i64 - target);
                            let right =
                                chr_midpoints.get(idx).map(|&v| v as i64 - target);

                            match (left, right) {
                                (Some(l), Some(r)) => {
                                    if l.unsigned_abs() <= r.unsigned_abs() {
                                        l
                                    } else {
                                        r
                                    }
                                }
                                (Some(l), None) => l,
                                (None, Some(r)) => r,
                                (None, None) => continue,
                            }
                        }
                    };
                    distances.push(distance);
                }
            } else {
                // No features on this chromosome — push i64::MAX for each region
                for _ in rs.iter_chr_regions(chromosome.as_str()) {
                    distances.push(i64::MAX);
                }
            }
        }
        Ok(distances)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Error;
    use std::path::PathBuf;

    use pretty_assertions::assert_eq;
    use rstest::*;

    fn get_test_path(file_name: &str) -> Result<PathBuf, Error> {
        let file_path: PathBuf = std::env::current_dir()
            .unwrap()
            .join("../tests/data/regionset")
            .join(file_name);
        Ok(file_path)
    }

    fn get_fasta_path(file_name: &str) -> PathBuf {
        std::env::current_dir()
            .unwrap()
            .join("../tests/data/fasta")
            .join(file_name)
    }

    // --- Strand ---

    #[test]
    fn test_strand_from_char() {
        assert_eq!(Strand::from_char('+'), Strand::Plus);
        assert_eq!(Strand::from_char('-'), Strand::Minus);
        assert_eq!(Strand::from_char('.'), Strand::Unstranded);
        assert_eq!(Strand::from_char('?'), Strand::Unstranded);
    }

    // --- Dinucleotide ---

    #[test]
    fn test_dinucleotide_from_bytes_all_variants() {
        let pairs = [
            (b"AA", Dinucleotide::Aa), (b"AC", Dinucleotide::Ac),
            (b"AG", Dinucleotide::Ag), (b"AT", Dinucleotide::At),
            (b"CA", Dinucleotide::Ca), (b"CC", Dinucleotide::Cc),
            (b"CG", Dinucleotide::Cg), (b"CT", Dinucleotide::Ct),
            (b"GA", Dinucleotide::Ga), (b"GC", Dinucleotide::Gc),
            (b"GG", Dinucleotide::Gg), (b"GT", Dinucleotide::Gt),
            (b"TA", Dinucleotide::Ta), (b"TC", Dinucleotide::Tc),
            (b"TG", Dinucleotide::Tg), (b"TT", Dinucleotide::Tt),
        ];
        for (bytes, expected) in &pairs {
            assert_eq!(Dinucleotide::from_bytes(&bytes[..]), Some(*expected));
        }
    }

    #[test]
    fn test_dinucleotide_case_insensitive() {
        assert_eq!(Dinucleotide::from_bytes(b"aa"), Some(Dinucleotide::Aa));
        assert_eq!(Dinucleotide::from_bytes(b"cG"), Some(Dinucleotide::Cg));
        assert_eq!(Dinucleotide::from_bytes(b"Tc"), Some(Dinucleotide::Tc));
    }

    #[test]
    fn test_dinucleotide_invalid() {
        assert_eq!(Dinucleotide::from_bytes(b"AN"), None);
        assert_eq!(Dinucleotide::from_bytes(b"A"), None);  // too short
        assert_eq!(Dinucleotide::from_bytes(b"ACG"), None); // too long
    }

    #[test]
    fn test_dinucleotide_to_string_round_trip() {
        let all = [
            Dinucleotide::Aa, Dinucleotide::Ac, Dinucleotide::Ag, Dinucleotide::At,
            Dinucleotide::Ca, Dinucleotide::Cc, Dinucleotide::Cg, Dinucleotide::Ct,
            Dinucleotide::Ga, Dinucleotide::Gc, Dinucleotide::Gg, Dinucleotide::Gt,
            Dinucleotide::Ta, Dinucleotide::Tc, Dinucleotide::Tg, Dinucleotide::Tt,
        ];
        for d in &all {
            let s = d.to_string().unwrap();
            assert_eq!(s.len(), 2);
            let round_tripped = Dinucleotide::from_bytes(s.as_bytes()).unwrap();
            assert_eq!(*d, round_tripped);
        }
    }

    // --- SortedRegionSet ---

    #[test]
    fn test_sorted_regionset_sorts_in_place() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 100, end: 200, rest: None },
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
            Region { chr: "chr2".into(), start: 5, end: 15, rest: None },
        ];
        let sorted = SortedRegionSet::new(RegionSet::from(regions));
        let starts: Vec<u32> = sorted.0.regions.iter().map(|r| r.start).collect();
        // chr1 regions should come first (sorted by chr then start)
        assert_eq!(starts, vec![10, 100, 5]);
    }

    // --- StrandedRegionSet ---

    #[test]
    fn test_stranded_regionset_new() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
            Region { chr: "chr1".into(), start: 30, end: 40, rest: None },
        ];
        let strands = vec![Strand::Plus, Strand::Minus];
        let srs = StrandedRegionSet::new(RegionSet::from(regions), strands);
        assert_eq!(srs.len(), 2);
        assert!(!srs.is_empty());
        assert_eq!(srs.strands[0], Strand::Plus);
        assert_eq!(srs.strands[1], Strand::Minus);
    }

    #[test]
    fn test_stranded_regionset_unstranded() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
        ];
        let srs = StrandedRegionSet::unstranded(RegionSet::from(regions));
        assert_eq!(srs.strands, vec![Strand::Unstranded]);
    }

    #[test]
    #[should_panic(expected = "regions and strands must have the same length")]
    fn test_stranded_regionset_mismatched_lengths() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
        ];
        StrandedRegionSet::new(RegionSet::from(regions), vec![]);
    }

    #[test]
    fn test_stranded_regionset_into_regionset() {
        let regions = vec![
            Region { chr: "chr1".into(), start: 10, end: 20, rest: None },
        ];
        let srs = StrandedRegionSet::unstranded(RegionSet::from(regions));
        let rs = srs.into_regionset();
        assert_eq!(rs.regions.len(), 1);
    }

    // --- GenomeAssembly ---

    #[test]
    fn test_genome_assembly_from_fasta() {
        // base.fa: chrX=TTGGGGAA, chr1=GGAA, chr2=GCGC
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.as_path()).unwrap();
        assert!(ga.contains_chr("chr1"));
        assert!(ga.contains_chr("chr2"));
        assert!(ga.contains_chr("chrX"));
        assert!(!ga.contains_chr("chr3"));
    }

    #[test]
    fn test_genome_assembly_seq_from_region() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.as_path()).unwrap();

        let region = Region { chr: "chr1".into(), start: 0, end: 4, rest: None };
        let seq = ga.seq_from_region(&region).unwrap();
        assert_eq!(seq, b"GGAA");

        let region2 = Region { chr: "chrX".into(), start: 2, end: 6, rest: None };
        let seq2 = ga.seq_from_region(&region2).unwrap();
        assert_eq!(seq2, b"GGGG");
    }

    #[test]
    fn test_genome_assembly_seq_unknown_chrom() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.as_path()).unwrap();

        let region = Region { chr: "chr99".into(), start: 0, end: 1, rest: None };
        assert!(ga.seq_from_region(&region).is_err());
    }

    #[test]
    fn test_genome_assembly_seq_out_of_bounds() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.as_path()).unwrap();

        // chr1 is only 4bp, request 0-100
        let region = Region { chr: "chr1".into(), start: 0, end: 100, rest: None };
        assert!(ga.seq_from_region(&region).is_err());
    }

    #[test]
    fn test_genome_assembly_try_from_str() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap());
        assert!(ga.is_ok());
    }

    #[test]
    fn test_genome_assembly_try_from_string() {
        let path = get_fasta_path("base.fa");
        let ga = GenomeAssembly::try_from(path.to_str().unwrap().to_string());
        assert!(ga.is_ok());
        assert!(ga.unwrap().contains_chr("chr1"));
    }

    #[test]
    fn test_binary_genome_assembly_round_trip() {
        // Write .fab from base.fa, then read it back and verify sequences match
        let fasta_path = get_fasta_path("base.fa");
        let fab_path = fasta_path.with_extension("fa.test.fab");

        BinaryGenomeAssembly::write_from_fasta(&fasta_path, &fab_path).unwrap();
        let bga = BinaryGenomeAssembly::from_file(&fab_path).unwrap();

        // Verify chromosomes exist
        assert!(bga.contains_chr("chr1"));
        assert!(bga.contains_chr("chr2"));
        assert!(bga.contains_chr("chrX"));
        assert!(!bga.contains_chr("chr3"));

        // Verify sequences match HashMap GenomeAssembly
        let ga = GenomeAssembly::try_from(fasta_path.as_path()).unwrap();

        let region1 = Region { chr: "chr1".into(), start: 0, end: 4, rest: None };
        assert_eq!(bga.seq_from_region(&region1).unwrap(), ga.seq_from_region(&region1).unwrap());
        assert_eq!(bga.seq_from_region(&region1).unwrap(), b"GGAA");

        let region2 = Region { chr: "chrX".into(), start: 2, end: 6, rest: None };
        assert_eq!(bga.seq_from_region(&region2).unwrap(), ga.seq_from_region(&region2).unwrap());
        assert_eq!(bga.seq_from_region(&region2).unwrap(), b"GGGG");

        // Out-of-bounds error
        let bad_region = Region { chr: "chr1".into(), start: 0, end: 100, rest: None };
        assert!(bga.seq_from_region(&bad_region).is_err());

        // Unknown chromosome error
        let unk_region = Region { chr: "chr99".into(), start: 0, end: 1, rest: None };
        assert!(bga.seq_from_region(&unk_region).is_err());

        // Clean up
        std::fs::remove_file(&fab_path).ok();
    }

    #[test]
    fn test_binary_genome_assembly_bad_magic() {
        let dir = tempfile::tempdir().unwrap();
        let fab_path = dir.path().join("bad.fab");
        // Need at least 9 bytes to pass the "too short" check
        std::fs::write(&fab_path, b"XXXX\x01\x00\x00\x00\x00").unwrap();
        let result = BinaryGenomeAssembly::from_file(&fab_path);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("bad magic"));
    }

    #[test]
    fn test_binary_genome_assembly_gc_parity() {
        // Verify calc_gc_content produces identical results from .fab and HashMap
        use crate::statistics::calc_gc_content;

        let fasta_path = get_fasta_path("base.fa");
        let fab_path = fasta_path.with_extension("fa.test2.fab");
        BinaryGenomeAssembly::write_from_fasta(&fasta_path, &fab_path).unwrap();

        let ga = GenomeAssembly::try_from(fasta_path.as_path()).unwrap();
        let bga = BinaryGenomeAssembly::from_file(&fab_path).unwrap();

        let regions = vec![
            Region { chr: "chr1".into(), start: 0, end: 4, rest: None },
            Region { chr: "chr2".into(), start: 0, end: 4, rest: None },
        ];
        let rs = RegionSet::from(regions);

        let gc_hashmap = calc_gc_content(&rs, &ga, false).unwrap();
        let gc_fab = calc_gc_content(&rs, &bga, false).unwrap();
        assert_eq!(gc_hashmap, gc_fab);

        std::fs::remove_file(&fab_path).ok();
    }

    #[test]
    fn test_tss_index_try_from_path() {
        let path = get_test_path("dummy_tss.bed").unwrap();
        let tss = TssIndex::try_from(path.as_path());
        assert!(tss.is_ok());
    }

    #[test]
    fn test_tss_index_try_from_path_invalid() {
        let path = PathBuf::from("/nonexistent/file.bed");
        let tss = TssIndex::try_from(path.as_path());
        assert!(tss.is_err());
    }

    #[test]
    fn test_tss_index_try_from_string() {
        let path = get_test_path("dummy_tss.bed").unwrap();
        let tss = TssIndex::try_from(path.to_str().unwrap().to_string());
        assert!(tss.is_ok());
    }

    // --- TssIndex sentinel behavior ---

    #[test]
    fn test_tss_distances_sentinel_for_missing_chrom() {
        // TSS features only on chr1
        let tss_regions = vec![
            Region { chr: "chr1".into(), start: 50, end: 51, rest: None },
        ];
        let tss_index = TssIndex::try_from(RegionSet::from(tss_regions)).unwrap();

        // Query has regions on chr1 and chr2 (no TSS on chr2)
        let query = RegionSet::from(vec![
            Region { chr: "chr1".into(), start: 40, end: 45, rest: None },
            Region { chr: "chr2".into(), start: 10, end: 20, rest: None },
        ]);

        let distances = tss_index.calc_tss_distances(&query, CoordinateMode::Bed).unwrap();
        assert_eq!(distances.len(), 2); // one per input region
        // One should be a real distance, the other should be u32::MAX sentinel
        // (order depends on HashSet iteration of iter_chroms)
        assert_eq!(distances.iter().filter(|&&d| d == u32::MAX).count(), 1);
        assert_eq!(distances.iter().filter(|&&d| d < u32::MAX).count(), 1);
    }

    #[test]
    fn test_feature_distances_sentinel_for_missing_chrom() {
        let tss_regions = vec![
            Region { chr: "chr1".into(), start: 50, end: 51, rest: None },
        ];
        let tss_index = TssIndex::try_from(RegionSet::from(tss_regions)).unwrap();

        let query = RegionSet::from(vec![
            Region { chr: "chr1".into(), start: 40, end: 45, rest: None },
            Region { chr: "chr2".into(), start: 10, end: 20, rest: None },
        ]);

        let distances = tss_index.calc_feature_distances(&query, CoordinateMode::Bed).unwrap();
        assert_eq!(distances.len(), 2);
        // One real distance, one i64::MAX sentinel
        assert_eq!(distances.iter().filter(|&&d| d == i64::MAX).count(), 1);
        assert_eq!(distances.iter().filter(|&&d| d != i64::MAX).count(), 1);
    }

    // --- Existing tests ---

    #[rstest]
    fn test_calc_tss_distances() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let tss_path = get_test_path("dummy_tss.bed").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
        let tss_index = TssIndex::try_from(tss_path.to_str().unwrap()).unwrap();

        let distances = tss_index.calc_tss_distances(&region_set, CoordinateMode::Bed).unwrap();

        assert_eq!(distances.len(), 9);
        assert_eq!(distances.iter().min(), Some(&2));
    }

    #[rstest]
    fn test_calc_feature_distances() {
        let file_path = get_test_path("dummy.narrowPeak").unwrap();
        let tss_path = get_test_path("dummy_tss.bed").unwrap();
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();
        let tss_index = TssIndex::try_from(tss_path.to_str().unwrap()).unwrap();

        let signed_distances = tss_index.calc_feature_distances(&region_set, CoordinateMode::Bed).unwrap();
        let abs_distances = tss_index.calc_tss_distances(&region_set, CoordinateMode::Bed).unwrap();

        // same number of results
        assert_eq!(signed_distances.len(), abs_distances.len());
        // absolute values should match calc_tss_distances
        for (signed, abs) in signed_distances.iter().zip(abs_distances.iter()) {
            assert_eq!(signed.unsigned_abs() as u32, *abs);
        }
        // should have both positive and negative distances
        assert!(signed_distances.iter().any(|d| *d > 0));
        assert!(signed_distances.iter().any(|d| *d < 0));
    }
}
