use std::path::Path;
use std::{collections::HashMap, io::BufRead};

use crate::common::models::Region;
use crate::common::utils::get_dynamic_reader;

use super::{Tokenizer, TokenizerError};

pub fn tokenize_fragment_file<P: AsRef<Path>>(
    file: P,
    tokenizer: &Tokenizer,
) -> Result<HashMap<String, Vec<u32>>, TokenizerError> {
    let reader = get_dynamic_reader(file.as_ref())?;
    let mut res = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let parts = line.split_whitespace().collect::<Vec<&str>>();

        if parts.len() < 5 {
            return Err(TokenizerError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Invalid fragment file detected at line: {i}"),
            )));
        }

        let chr = parts[0].to_string();
        let start = parts[1].parse::<u32>().unwrap();
        let end = parts[2].parse::<u32>().unwrap();
        let barcode = parts[3].to_string();
        // we skip read support

        let tokens_for_cell = res.entry(barcode).or_insert(vec![]);
        let tokens = tokenizer.tokenize(&[Region {
            chr,
            start,
            end,
            rest: None,
        }])?;
        let tokens = tokens
            .iter()
            .map(|t| tokenizer.convert_token_to_id(t).unwrap());

        tokens_for_cell.extend(tokens);
    }

    Ok(res)
}
