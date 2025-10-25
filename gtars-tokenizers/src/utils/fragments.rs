use std::path::Path;
use std::{collections::HashMap, io::BufRead};

use gtars_core::models::Region;
use gtars_core::utils::get_dynamic_reader;

use super::super::{Tokenizer, TokenizerError};

/// Helper function: process fragments and call a callback for each barcode's tokens
///
/// This factors out the common file reading/parsing/tokenization logic.
/// The callback determines what to do with the tokens (collect vs count).
fn process_fragments_by_barcode<P, F>(
    file: P,
    tokenizer: &Tokenizer,
    mut callback: F,
) -> Result<(), TokenizerError>
where
    P: AsRef<Path>,
    F: FnMut(&str, &[u32]) -> Result<(), TokenizerError>,
{
    let reader = get_dynamic_reader(file.as_ref())?;

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
        let start = parts[1].parse::<u32>().map_err(|e| {
            TokenizerError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Failed to parse start position at line {i}: {e}"),
            ))
        })?;
        let end = parts[2].parse::<u32>().map_err(|e| {
            TokenizerError::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Failed to parse end position at line {i}: {e}"),
            ))
        })?;
        let barcode = parts[3];

        // Get overlapping peaks using tokenizer
        let tokens = tokenizer.tokenize(&[Region {
            chr,
            start,
            end,
            rest: None,
        }])?;

        // Convert to token IDs
        let token_ids: Vec<u32> = tokens
            .iter()
            .map(|t| tokenizer.convert_token_to_id(t).unwrap())
            .collect();

        // Call callback with barcode and token IDs
        callback(barcode, &token_ids)?;
    }

    Ok(())
}

/// Tokenize fragment file (EXISTING - refactored to use helper)
///
/// Collects all token IDs for each barcode (can have duplicates).
pub fn tokenize_fragment_file<P>(
    file: P,
    tokenizer: &Tokenizer,
) -> Result<HashMap<String, Vec<u32>>, TokenizerError>
where
    P: AsRef<Path>,
{
    let mut res: HashMap<String, Vec<u32>> = HashMap::new();

    process_fragments_by_barcode(file, tokenizer, |barcode, token_ids| {
        let tokens_for_cell = res.entry(barcode.to_string()).or_insert_with(Vec::new);
        tokens_for_cell.extend_from_slice(token_ids);
        Ok(())
    })?;

    Ok(res)
}

/// Count fragment overlaps by barcode (for single-cell analysis)
///
/// Returns a sparse map: barcode → (peak_id → count)
pub fn count_fragments_by_barcode<P>(
    file: P,
    tokenizer: &Tokenizer,
) -> Result<HashMap<String, HashMap<u32, u32>>, TokenizerError>
where
    P: AsRef<Path>,
{
    let mut res: HashMap<String, HashMap<u32, u32>> = HashMap::new();

    process_fragments_by_barcode(file, tokenizer, |barcode, token_ids| {
        let counts = res.entry(barcode.to_string()).or_insert_with(HashMap::new);

        // COUNT occurrences
        for &token_id in token_ids {
            *counts.entry(token_id).or_insert(0) += 1;
        }

        Ok(())
    })?;

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

    #[rstest]
    fn test_count_fragments_by_barcode() {
        let tokenizer = Tokenizer::from_bed("../tests/data/consensus/consensus1.bed").unwrap();
        let result = count_fragments_by_barcode(
            "../tests/data/fragments/region_scoring/fragments1.bed.gz",
            &tokenizer,
        );
        assert_eq!(result.is_ok(), true);

        let counts = result.unwrap();

        // Verify we got barcodes
        assert!(counts.len() > 0);

        // Verify counts are sparse HashMaps
        for (barcode, peak_counts) in counts.iter() {
            assert!(peak_counts.len() > 0);  // Each cell has some peaks
            for (&peak_id, &count) in peak_counts.iter() {
                assert!(count > 0);  // No zeros stored
            }
        }
    }
}
