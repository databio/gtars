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
                .help("Path to chrom.sizes file (enables expected partitions and promoter trimming)"),
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
                .help("Number of bins for region distribution"),
        )
}
