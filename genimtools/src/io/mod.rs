use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

///
/// Writes a vector of tokens to a file in the `.gtok` format.
/// # Arguments
/// - filename: the file to save the tokens to
/// - tokens: tokens to save
///
pub fn write_tokens_to_gtok(filename: &str, tokens: &[u32]) -> std::io::Result<()> {
    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    for &token in tokens {
        writer.write_all(&token.to_le_bytes())?;
    }

    Ok(())
}

///
/// Read in a vector of tokens from a file in the `.gtok` format.
/// # Arguments
/// - filename: filename to read the tokens from
///
/// # Returns
/// - vector of tokens in u32 format
pub fn read_tokens_from_gtok(filename: &str) -> std::io::Result<Vec<u32>> {
    let file = File::open(filename)?;
    let mut reader = BufReader::new(file);
    let mut tokens = Vec::new();
    let mut buffer = [0; 4];

    while let Ok(()) = reader.read_exact(&mut buffer) {
        tokens.push(u32::from_le_bytes(buffer));
    }

    Ok(tokens)
}
