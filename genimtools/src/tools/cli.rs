use super::*;
use clap::{arg, ArgMatches, Command};

pub fn make_tools_cli() -> Command {
    Command::new(consts::DATA_DIR_STAT_CMD)
        .author("Nathan LeRoy")
        .about("Collect data statistics on all bed files in a directory.")
        .arg(arg!(<VALUE> "Path to the data directory.").required(true))
        .arg(
            arg!(--out <VALUE> "Path to the output file.")
                .required(false),
        )
}

pub mod handlers {
    use super::*;

    pub fn data_dir_stat_handler(matches: &ArgMatches) {
        let path = matches
            .get_one::<String>("path")
            .expect("Path is required");

        let out = matches
            .get_one::<String>("out")
            .unwrap_or(&consts::DEFAULT_OUTPUT.to_string())
            .to_owned();

        // core logic/algorithm here
        data_dir_stat(path, out.as_str());
    }
}