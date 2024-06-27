use std::fs::File;
use std::fs::OpenOptions;
use std::io::{BufReader, BufWriter, Read, Write};

use anyhow::{Context, Result};

use super::consts::{GTOK_HEADER, GTOK_U16_FLAG, GTOK_U32_FLAG};

///
/// Writes a vector of tokens to a file in the `.gtok` format.
/// # Arguments
/// - filename: the file to save the tokens to
/// - tokens: tokens to save
///
pub fn write_tokens_to_gtok(filename: &str, tokens: &[u32]) -> Result<()> {
    // make sure the path exists
    let path = std::path::Path::new(filename);

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    } else {
        anyhow::bail!("Failed to create parent directories for gtok file!")
    }

    let file = File::create(filename).with_context(|| "Failed to create gtok file!")?;
    let mut writer = BufWriter::new(file);

    // write the header
    writer
        .write_all(GTOK_HEADER)
        .with_context(|| "Failed to write GTOK header to file!")?;

    // determine size of tokens
    let is_small = tokens.iter().all(|&x| x <= u16::MAX as u32);
    let flag = if is_small {
        GTOK_U16_FLAG
    } else {
        GTOK_U32_FLAG
    };
    writer
        .write_all(&flag.to_le_bytes())
        .with_context(|| "Failed to write GTOK size flag to file!")?;

    for &token in tokens {
        if is_small {
            writer
                .write_all(&(token as u16).to_le_bytes())
                .with_context(|| "Failed to write bytes to file!")?;
            continue;
        }
        writer
            .write_all(&token.to_le_bytes())
            .with_context(|| "Failed to write bytes to file!")?;
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
pub fn read_tokens_from_gtok(filename: &str) -> Result<Vec<u32>> {
    let file = File::open(filename)?;
    let mut reader = BufReader::new(file);

    // check the header
    let mut header = [0; 4];
    reader.read_exact(&mut header)?;

    if &header != GTOK_HEADER {
        anyhow::bail!("File doesn't appear to be a valid .gtok file.")
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
        }
        [GTOK_U32_FLAG] => {
            let mut buffer = [0; 4];
            while let Ok(()) = reader.read_exact(&mut buffer) {
                tokens.push(u32::from_le_bytes(buffer));
            }
        }
        _ => {
            anyhow::bail!("Invalid data format flag found in gtok file")
        }
    }

    Ok(tokens)
}

///
/// Initialize a `.gtok` file with a header and size flag.
/// # Arguments
/// - filename: the file to initialize
///
/// # Returns
/// - Result<(), anyhow::Error>
pub fn init_gtok_file(filename: &str) -> Result<()> {
    // make sure the path exists
    let path = std::path::Path::new(filename);

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    } else {
        anyhow::bail!("Failed to create parent directories for gtok file!")
    }

    let file = File::create(filename).with_context(|| "Failed to create gtok file!")?;
    let mut writer = BufWriter::new(file);

    writer
        .write_all(GTOK_HEADER)
        .with_context(|| "Failed to write GTOK header to file!")?;

    // assume large and write u32 flag
    writer
        .write_all(&GTOK_U32_FLAG.to_le_bytes())
        .with_context(|| "Failed to write GTOK size flag to file!")?;

    Ok(())
}

///
/// Add tokens to the end of an existing `.gtok` file.
///
/// # Arguments
/// - filename: the file to append to
///
/// # Returns
/// - Result<(), anyhow::Error>
pub fn append_tokens_to_gtok_file(filename: &str, tokens: &[u32]) -> Result<()> {
    let file = File::open(filename).with_context(|| "Failed to open gtok file!")?;

    let mut reader = BufReader::new(file);

    // check the header
    let mut header = [0; 4];
    reader.read_exact(&mut header)?;

    if &header != GTOK_HEADER {
        anyhow::bail!("File doesn't appear to be a valid .gtok file.")
    }

    // detect the size flag
    let mut size_flag = [0; 1];
    reader.read_exact(&mut size_flag)?;

    // start appending to the open file
    // must reopen because `Bufreader` takes ownership of `file`.
    let file = OpenOptions::new()
        .append(true)
        .open(filename)
        .with_context(|| "Failed to open gtok file for appending")?;
    let mut writer = BufWriter::new(file);

    match size_flag {
        [GTOK_U16_FLAG] => {
            for token in tokens {
                writer
                    .write_all(&(*token as u16).to_le_bytes())
                    .with_context(|| "Failed to write bytes to file!")?;
            }
        }
        [GTOK_U32_FLAG] => {
            for token in tokens {
                writer
                    .write_all(&token.to_le_bytes())
                    .with_context(|| "Failed to write bytes to file!")?;
            }
        }
        _ => {
            anyhow::bail!("Invalid data format flag found in gtok file")
        }
    }

    Ok(())
}
