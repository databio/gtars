use super::*;
use clap::{arg, ArgMatches, Command};

fn make_data_dir_stat_cli() -> Command {
    Command::new(consts::DATA_DIR_STAT_CMD)
        .author("Nathan LeRoy")
        .about("Collect data statistics on all bed files in a directory.")
        .arg(arg!(--out <VALUE> "Path to the output file.").required(false))
        // positional path
        .arg(arg!(<path> "Path to the data directory.").required(true))
}

fn make_pre_tokenization_cli() -> Command {
    Command::new(consts::PRE_TOKENIZATION_CMD)
        .about("Pre-tokenize a bed file or folder of bed files into a specific universe.")
        .arg(arg!(--universe <VALUE> "Path to the output folder or file.").required(true))
        .arg(arg!(--out <VALUUE> "Path to output the new data.").required(false))
        // positional path
        .arg(arg!(<path> "Path to the data directory.").required(true))
}

pub fn make_tools_cli() -> Command {
    Command::new(consts::TOOLS_CMD)
        .author("Nathan LeRoy")
        .about("Tools for working with genomic data.")
        .subcommand(make_data_dir_stat_cli())
        .subcommand(make_pre_tokenization_cli())
}

pub mod handlers {

    use std::path::Path;

    use super::*;

    pub fn tools_handler(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
        match matches.subcommand() {
            Some((consts::DATA_DIR_STAT_CMD, matches)) => {
                data_dir_stat_handler(matches)?;
            }
            Some((consts::PRE_TOKENIZATION_CMD, matches)) => {
                pre_tokenization_handler(matches)?;
            }
            _ => unreachable!("Subcommand not found"),
        }

        Ok(())
    }

    pub fn data_dir_stat_handler(matches: &ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
        let path = matches.get_one::<String>("path").expect("Path is required");

        let out = matches
            .get_one::<String>("out")
            .unwrap_or(&consts::DEFAULT_DATA_DIR_STAT_OUTPUT.to_string())
            .to_owned();

        // core logic/algorithm here
        data_dir_stat(path, out.as_str());

        Ok(())
    }

    pub fn pre_tokenization_handler(
        matches: &ArgMatches,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let path = matches
            .get_one::<String>("path")
            .expect("Path to either a data file or a directory with data is required");

        let universe = matches
            .get_one::<String>("universe")
            .expect("Path to the universe file is required");

        let binding = consts::DEFAULT_PRETOKENIZE_OUT.to_string();
        let outdir = matches.get_one::<String>("out").unwrap_or(&binding);

        // check if the path is a file or a directory
        let path_to_data = Path::new(&path);
        let universe = Path::new(&universe);

        // create the tokenizer
        let tokenizer = crate::tokenizers::TreeTokenizer::try_from(universe)?;

        pre_tokenize_data(path_to_data, outdir, &tokenizer)?;

        Ok(())
    }
}
