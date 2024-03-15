use clap::{Arg, Command};

use crate::uniwig::consts::UNIWIG_CMD;

pub fn create_uniwig_cli() -> Command {
    Command::new(UNIWIG_CMD)
        .author("DRC")
        .about("Given a set of bed files, we want to produce 2")
        .arg(
            Arg::new("sorted")
                .long("sorted")
                .short('s')
                .help("Specify if the provided bed file is already sorted by the chromosome number.")
                .required(false)
        )

}