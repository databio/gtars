// use std::path::Path;
// use std::io::{self, Write};

// use anyhow::Result;
// use clap::{arg, Command, ArgMatches};

// use super::consts::TOKENIZERS_CMD;
// use super::config::TokenizerType;
// use super::utils::create_tokenize_core_from_universe;
// use super::utils::special_tokens::SpecialTokens;
// use super::universe::Universe;
// use super::tokenizer::Tokenizer;

// use crate::common::models::RegionSet;

// pub fn create_tokenizer_cli() -> Command {
//     Command::new(TOKENIZERS_CMD)
//         .author("NJL")
//         .about("Tokenize data into a universe")
//         .arg_required_else_help(true)
//         .arg(arg!(-a <query> "The file you are tokenizing"))
//         .arg(arg!(-b <universe> "The universe you are tokenizing into"))
//         .arg(arg!(-e --backend <backend> "Which backend to use (ailist or bits)"))
//         .arg(arg!(-i --ids "Print out input ids instead of regions"))
// }

// pub mod handlers {

//     use super::*;

//     pub fn tokenize_query_into_universe(matches: &ArgMatches) -> Result<()> {
//         let query_file = matches
//             .get_one::<String>("query")
//             .expect("A path to a query file is required.");

//         let universe_file = matches
//             .get_one::<String>("universe")
//             .expect("A path to a universe file is required.");

//         let output_ids = matches.get_one::<bool>("ids").unwrap_or(&false);

//         let default_backend = "bits".to_string();
//         let backend_str = matches.get_one::<String>("backend")
//             .unwrap_or(&default_backend);

//         let backend_type = match backend_str.as_str() {
//             "bits" => TokenizerType::Bits,
//             "ailist" => TokenizerType::AiList,
//             _ => return Err(anyhow::anyhow!("Invalid backend type: {}. Valid options are 'bits' or 'ailist'", backend_str)),
//         };

//         let universe_path = Path::new(universe_file);
//         let universe = Universe::try_from(universe_path)?;

//         let core = create_tokenize_core_from_universe(&universe, backend_type);
//         let tokenizer = Tokenizer::new(core, universe, SpecialTokens::default());

//         let query_path = Path::new(query_file);
//         let rs = RegionSet::try_from(query_path)?;
//         let tokens = tokenizer.tokenize(&rs.regions)?;

//         let stdout = io::stdout();
//         let mut handle = stdout.lock();

//         match output_ids {
//             true => {
//                 // Output token IDs to stdout
//                 for token in tokens {
//                     let id = tokenizer.convert_token_to_id(&token).unwrap();
//                     writeln!(handle, "{}", id)?;
//                 }
//             }
//             false => {
//                 // Output tokens formatted as chr\tstart\tend
//                 for token in tokens {
//                     // Tokens are formatted as chr:start-end, need to print chr\tstart\tend
//                     if let Some((chr, range)) = token.split_once(':') {
//                         if let Some((start, end)) = range.split_once('-') {
//                             writeln!(handle, "{}\t{}\t{}", chr, start, end)?;
//                         }
//                     }
//                 }
//             }
//         }

//         Ok(())
//     }
// }
