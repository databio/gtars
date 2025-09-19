pub const SCATRS_CMD: &str = "scatrs";
pub const VERSION: &str = "0.1.0";

pub const DEFAULT_EXTEND_BP: u32 = 250;
pub const DEFAULT_BIN_SIZE: u32 = 5000;
pub const DEFAULT_MERGE_DISTANCE: u32 = 20;
pub const DEFAULT_SIGNAL_TO_NOISE: f64 = 0.9;
pub const DEFAULT_FRAGMENTS_PER_CELL: u32 = 8000;

pub const VALID_CHROMOSOMES: &[&str] = &[
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
    "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"
];

pub const DEFAULT_THREAD_COUNT: usize = 4;

pub fn get_thread_count(requested: Option<usize>) -> usize {
    requested.unwrap_or(DEFAULT_THREAD_COUNT)
}