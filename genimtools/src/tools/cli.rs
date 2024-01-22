use super::*;
use clap::{arg, ArgMatches, Command};

fn make_data_dir_stat_cli() -> Command {
    Command::new(consts::DATA_DIR_STAT_CMD)
        .author("Nathan LeRoy")
        .about("Collect data statistics on all bed files in a directory.")
        .arg(
            arg!(--out <VALUE> "Path to the output file.")
                .required(false),
        )
        // positional path
        .arg(arg!(<path> "Path to the data directory.").required(true))
}

fn make_pre_tokenization_cli() -> Command {
    Command::new(consts::PRE_TOKENIZATION_CMD)
        .about("Pre-tokenize a bed file or folder of bed files into a specific universe.")
        .arg(
            arg!(--universe <VALUE> "Path to the output folder or file.")
                .required(true),
        )
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
    
    use std::path::{Path, PathBuf};

    use crate::{tokenizers::{self, Tokenizer}, common::models::RegionSet};

    use super::*;

    pub fn tools_handler(matches: &ArgMatches) {
        match matches.subcommand() {
            Some((consts::DATA_DIR_STAT_CMD, matches)) => {
                data_dir_stat_handler(matches);
            },
            Some((consts::PRE_TOKENIZATION_CMD, matches)) => {
                pre_tokenization_handler(matches);
            }
            _ => unreachable!("Subcommand not found"),
        }
    }

    pub fn data_dir_stat_handler(matches: &ArgMatches) {
        let path = matches
            .get_one::<String>("path")
            .expect("Path is required");

        let out = matches
            .get_one::<String>("out")
            .unwrap_or(&consts::DEFAULT_DATA_DIR_STAT_OUTPUT.to_string())
            .to_owned();

        // core logic/algorithm here
        data_dir_stat(path, out.as_str());
    }

    pub fn pre_tokenization_handler(matches: &ArgMatches) {
        let path = matches
            .get_one::<String>("path")
            .expect("Path to either a data file or a directory with data is required");

        let universe = matches
            .get_one::<String>("universe")
            .expect("Path to the universe file is required");

        // check if the path is a file or a directory
        let path_to_data = Path::new(&path);
        let universe = Path::new(&universe);

        if path_to_data.is_file() {
            let file_name = path_to_data.file_stem();
            match file_name {
                Some(file_name) => {
                    let new_file = format!("{}.{}", file_name.to_str().unwrap(), consts::PRE_TOKENIZATION_EXT);
                    let new_file = Path::new(&new_file);
                    let tokenizer = tokenizers::TreeTokenizer::from(universe);
                    let data = RegionSet::try_from(path_to_data);
                    match data {
                        Ok(data) => {
                            let result = tokenizer.tokenize_region_set(&data).expect("Data couldn't be tokenized.");
                            let _ = result.to_gtok_file(&PathBuf::from(new_file));
                        },
                        Err(e) => panic!("There was an error readig the data file: {}", e)
                    }
                },
                None => panic!("There was an issue extracting the name of the file.")
            }
            
        }
        
    }
}