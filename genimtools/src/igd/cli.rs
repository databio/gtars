
use super::*;
use clap::{arg, ArgMatches, Command};
use crate::vocab::consts;

pub fn create_igd_cli() -> Command {
    Command::new("igd")
        .author("DRC")
        .about("Create a integrated genome database (IGD)")
        .arg(arg!(--output <VALUE> "Path to the output.").required(true))
        .arg(
            arg!(--filelist <VALUE> "Path to the list of files. This should be a folder of bed files.")
                .required(true),
        )
}