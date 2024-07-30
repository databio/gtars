use clap::{Arg, Command};

use crate::uniwig::consts::UNIWIG_CMD;

/// Creates the uniwig CLI Command object
pub fn create_uniwig_cli() -> Command {
    Command::new(UNIWIG_CMD)
        .author("DRC")
        .about("Given a set of bed files, we want to produce 2")
        .arg(
            Arg::new("bed")
                .long("bed")
                .short('b')
                .help("Path to the combined bed file we want to transform")
                .required(true),
        )
        .arg(
            Arg::new("chromref")
                .long("chromref")
                .short('c')
                .help("Path to chromreference")
                .required(true),
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
                .short('t')
                .value_parser(clap::value_parser!(i32))
                .help("Integer value for stepsize")
                .required(true),
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
}
