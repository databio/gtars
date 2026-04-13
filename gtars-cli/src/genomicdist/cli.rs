use clap::{Arg, Command, arg};

pub const GENOMICDIST_CMD: &str = "genomicdist";

pub fn create_genomicdist_cli() -> Command {
    Command::new(GENOMICDIST_CMD)
        .about("Compute genomic distribution statistics for a BED file.")
        .arg(
            arg!(--bed <BED>)
                .required(true)
                .help("Path to input BED file"),
        )
        .arg(
            arg!(--gtf <GTF>)
                .required(false)
                .help("Path to GTF/GTF.gz gene model (enables partitions; also TSS distances if no --tss)"),
        )
        .arg(
            arg!(--tss <TSS>)
                .required(false)
                .help("Path to TSS BED file (overrides GTF-derived TSS for distance calculation)"),
        )
        .arg(
            Arg::new("chrom-sizes")
                .long("chrom-sizes")
                .required(false)
                .help("Path to chrom.sizes file. When provided, region distribution uses a per-chromosome bin size derived from the reference genome (stable across files). Also enables expected partitions and promoter trimming."),
        )
        .arg(
            arg!(--output <OUTPUT>)
                .required(false)
                .help("Output JSON path (default: stdout)"),
        )
        .arg(
            arg!(--bins <BINS>)
                .required(false)
                .default_value("250")
                .help("Target bin count for the longest chromosome in --chrom-sizes. Bin width is derived as max_chrom_len/bins; shorter chromosomes get proportionally fewer bins. Total windows returned by density_vector/density_homogeneity is typically larger than this value (sum of ceil(size/bin_width) across chromosomes). Pass bins = max_chrom_len / desired_bin_width_bp to target a specific bin width."),
        )
        .arg(
            Arg::new("signal-matrix")
                .long("signal-matrix")
                .required(false)
                .help("Path to open signal matrix TSV (enables cell-type open chromatin enrichment)"),
        )
        .arg(
            Arg::new("fasta")
                .long("fasta")
                .required(false)
                .help("Path to genome FASTA (.fa) or binary FASTA (.fab) file. Enables GC content; also enables dinucleotide frequencies when --dinucl-freq is set. Use .fab format (via gtars prep --fasta) for best performance."),
        )
        .arg(
            Arg::new("ignore-unk-chroms")
                .long("ignore-unk-chroms")
                .action(clap::ArgAction::SetTrue)
                .help("When computing GC content, skip regions on chromosomes not in the FASTA (default: error)"),
        )
        .arg(
            Arg::new("dinucl-freq")
                .long("dinucl-freq")
                .action(clap::ArgAction::SetTrue)
                .help("Compute per-region dinucleotide frequencies (expensive for wide regions; opt-in even when --fasta is provided)"),
        )
        .arg(
            Arg::new("dinucl-raw-counts")
                .long("dinucl-raw-counts")
                .action(clap::ArgAction::SetTrue)
                .help("Return raw per-region dinucleotide counts instead of percentages (matches R GenomicDistributions' rawCounts=TRUE)"),
        )
        .arg(
            Arg::new("promoter-upstream")
                .long("promoter-upstream")
                .required(false)
                .default_value("200")
                .help("Upstream distance (bp) from TSS to define promoter regions"),
        )
        .arg(
            Arg::new("promoter-downstream")
                .long("promoter-downstream")
                .required(false)
                .default_value("2000")
                .help("Downstream distance (bp) from TSS to define promoter regions"),
        )
        .arg(
            Arg::new("cluster-radii")
                .long("cluster-radii")
                .required(false)
                .default_value("500,5000,50000")
                .help("Comma-separated list of stitching radii (bp) for peak_clusters summary statistics. Default probes promoter (500), enhancer (5000), and domain (50000) scales."),
        )
        .arg(
            Arg::new("compact")
                .long("compact")
                .action(clap::ArgAction::SetTrue)
                .help("Compact JSON output (default: pretty-printed)"),
        )
}
