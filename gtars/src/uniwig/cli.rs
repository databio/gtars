use clap::{Arg, ArgAction, Command};

use crate::uniwig::consts::UNIWIG_CMD;

/// Creates the uniwig CLI Command object
///
/// Example to run uiwig
/// `cargo run uniwig -f /sourcefiles/test.bed -t "bed" -c /sourcefiles/hg38.chrom.sizes -m 5 -t 1 -l /numpy_arrays_created_with_rust/ -y npy`
pub fn create_uniwig_cli() -> Command {
    Command::new(UNIWIG_CMD)
        .author("DRC")
        .about("Create accumulation files from a BED or BAM file")
        .arg(
            Arg::new("file")
                .long("file")
                .short('f')
                .help("Path to the combined bed file we want to transform or a sorted bam file")
                .required(true),
        )
        .arg(
            Arg::new("filetype")
                .long("filetype")
                .short('t')
                .help("Input file type, 'bed' 'bam' or 'narrowpeak'")
                .default_value("bed"),
        )
        .arg(
            Arg::new("chromref")
                .long("chromref")
                .short('c')
                .help("Path to chromreference")
                .required(false),
        )
        .arg(
            Arg::new("smoothsize")
                .long("smoothsize")
                .short('m')
                .value_parser(clap::value_parser!(i32))
                .help("Integer value for smoothing")
                .required(true),
        )
        .arg(
            Arg::new("stepsize")
                .long("stepsize")
                .short('s')
                .value_parser(clap::value_parser!(i32))
                .help("Integer value for stepsize")
                .required(true),
        )
        .arg(
            Arg::new("bamscale")
                .long("bamscale")
                .short('e')
                .default_value("1.0")
                .value_parser(clap::value_parser!(f32))
                .help("Integer for scaling bam read values, default is 1")
                .required(false),
        )
        .arg(
            Arg::new("fileheader")
                .long("fileheader")
                .short('l')
                .help("Name of the file")
                .required(true),
        )
        .arg(
            Arg::new("outputtype")
                .long("outputtype")
                .short('y')
                .help("Output as wiggle or npy")
                .required(true),
        )
        .arg(
            Arg::new("counttype")
                .long("counttype")
                .short('u')
                .default_value("all")
                .help("Select to only output start, end, or core. Select `shift` for bam workflows. Defaults to all.")
                .required(false),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .short('p')
                .default_value("6")
                .value_parser(clap::value_parser!(i32))
                .help("Number of rayon threads to use for parallel processing")
                .required(false),
        )
        .arg(
            Arg::new("score")
                .long("score")
                .short('o')
                .help("Count via score (narrowPeak only!)")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("no-bamshift")
                .long("no-bamshift")
                .short('a')
                .help("Set bam shift to False, i.e. uniwig will count raw reads without considering read direction.")
                .action(ArgAction::SetFalse),
        )
        .arg(
            Arg::new("zoom")
                .long("zoom")
                .short('z')
                .default_value("1")
                .value_parser(clap::value_parser!(i32))
                .help("Number of zoom levels (for bw file output only")
                .required(false),
        )
        .arg(
            Arg::new("debug")
                .long("debug")
                .short('d')
                .help("Print more verbose debug messages?")
                .action(ArgAction::SetTrue),
        )
}
