
use clap::{arg, ArgMatches, Command};
use crate::igd::consts::IGD_CMD;

pub fn create_igd_cli() -> Command {
    Command::new(IGD_CMD)
        .author("DRC")
        .about("Create a integrated genome database (IGD)")
        .arg(arg!(--output <VALUE> "Path to the output.").required(true))
        .arg(
            arg!(--filelist <VALUE> "Path to the list of files. This should be a folder of bed files.")
                .required(true),
        )
}