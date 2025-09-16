use std::path::Path;
use std::{collections::HashMap, io::BufRead};

use gtars_core::models::Region;
use gtars_core::utils::get_dynamic_reader;

use super::super::{Tokenizer, TokenizerError};

pub fn tokenize_fragment_file<P>(
    file: P,
    tokenizer: &Tokenizer,
) -> Result<HashMap<String, Vec<u32>>, TokenizerError>
where
    P: AsRef<Path>,
{
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

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    fn test_tokenize_frament() {
        let tokenizer = Tokenizer::from_bed("../tests/data/consensus/consensus1.bed").unwrap();
        let result = tokenize_fragment_file(
            "../tests/data/fragments/region_scoring/fragments1.bed.gz",
            &tokenizer,
        );
        assert_eq!(result.is_ok(), true);

        let result = result.unwrap();
        assert_eq!(result.clone().keys().len(), 2);
        println!("{:?}", result);
    }
}
