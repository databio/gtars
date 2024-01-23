use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

pub mod consts;

use consts::{GTOK_HEADER, GTOK_U16_FLAG, GTOK_U32_FLAG};

///
/// Writes a vector of tokens to a file in the `.gtok` format.
/// # Arguments
/// - filename: the file to save the tokens to
/// - tokens: tokens to save
///
pub fn write_tokens_to_gtok(filename: &str, tokens: &[u32]) -> std::io::Result<()> {
    // make sure the path exists
    let path = std::path::Path::new(filename);
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    } else {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            "Failed to create directory",
        ));
    }

    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    // write the header
    writer.write_all(GTOK_HEADER)?;
    
    // determine size of tokens
    let is_small = tokens.iter().all(|&x| x <= u16::MAX as u32);
    let flag = if is_small { GTOK_U16_FLAG } else { GTOK_U32_FLAG };
    writer.write_all(&flag.to_le_bytes())?;

    for &token in tokens {
        if is_small {
            writer.write_all(&(token as u16).to_le_bytes())?;
            continue;
        }
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

    // check the header
    let mut header = [0; 4];
    reader.read_exact(&mut header)?;

    if &header != GTOK_HEADER {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            "File doesn't appear to be a valid .gtok file.",
        ));
    }

    let mut size_flag = [0; 1];
    reader.read_exact(&mut size_flag)?;

    let mut tokens = Vec::new();    

    match size_flag {
        [GTOK_U16_FLAG] => {
            let mut buffer = [0; 2];
            while let Ok(()) = reader.read_exact(&mut buffer) {
                tokens.push(u16::from_le_bytes(buffer) as u32);
            }
        },
        [GTOK_U32_FLAG] => {
            let mut buffer = [0; 4];
            while let Ok(()) = reader.read_exact(&mut buffer) {
                tokens.push(u32::from_le_bytes(buffer));
            }
        },
        _ => return Err(std::io::Error::new(std::io::ErrorKind::Other, "Invalid data format flag")),
    }

    

    Ok(tokens)
}
