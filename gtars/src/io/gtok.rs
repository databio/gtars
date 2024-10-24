use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

use anyhow::{Context, Result};

use crate::common::models::tokenized_region::TokenizedRegionPointer;

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
/// Writes a vector of tokens to a file in the `.gtokp` format. The  `gtokp`
/// format is the same as the `gtok` format the addition of positional information.
/// # Arguments
/// - filename: the file to save the tokens to
/// - tokens: tokens to save
///
pub fn write_tokens_to_gtokp(filename: &str, pointers: &[TokenizedRegionPointer]) -> Result<()> {
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
    let is_small = pointers.iter().all(|&x| x.id <= u16::MAX as u32);
    let flag = if is_small {
        GTOK_U16_FLAG
    } else {
        GTOK_U32_FLAG
    };
    writer
        .write_all(&flag.to_le_bytes())
        .with_context(|| "Failed to write GTOK size flag to file!")?;

    for pointer in pointers {
        if is_small {
            writer
                .write_all(&(pointer.id as u16).to_le_bytes())
                .with_context(|| "Failed to write bytes to file!")?;
        } else {
            writer
                .write_all(&pointer.id.to_le_bytes())
                .with_context(|| "Failed to write bytes to file!")?;
        }

        writer
            .write_all(&pointer.chrom_id.to_le_bytes())
            .with_context(|| "Failed to write chrom_id bytes to file!")?;

        writer
            .write_all(&pointer.source_start.to_le_bytes())
            .with_context(|| "Failed to write source_start bytes to file!")?;

        writer
            .write_all(&pointer.source_end.to_le_bytes())
            .with_context(|| "Failed to write source_start bytes to file!")?;
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
/// Read in a vector of tokens from a file in the `.gtokp` format.
/// # Arguments
/// - filename: filename to read the tokens from
///
/// # Returns
/// - vector of TokenizedRegionPointer
pub fn read_tokens_from_gtokp(filename: &str) -> Result<Vec<TokenizedRegionPointer>> {
    let file = File::open(filename)?;
    let mut reader = BufReader::new(file);

    // check the header
    let mut header = [0; 4];
    reader.read_exact(&mut header)?;

    if &header != GTOK_HEADER {
        anyhow::bail!("File doesn't appear to be a valid .gtokp file.")
    }

    let mut size_flag = [0; 1];
    reader.read_exact(&mut size_flag)?;

    let mut pointers = Vec::new();

    match size_flag {
        [GTOK_U16_FLAG] => {
            let mut buffer = [0; 2 + 4 + 4 + 4]; // u16 + u32 + u32 + u32
            while let Ok(()) = reader.read_exact(&mut buffer) {
                let id = u16::from_le_bytes([buffer[0], buffer[1]]) as u32;
                let chrom_id =
                    u32::from_le_bytes([buffer[2], buffer[3], buffer[4], buffer[5]]) as u16;
                let source_start = u32::from_le_bytes([buffer[6], buffer[7], buffer[8], buffer[9]]);
                let source_end =
                    u32::from_le_bytes([buffer[10], buffer[11], buffer[12], buffer[13]]);
                pointers.push(TokenizedRegionPointer {
                    id,
                    chrom_id,
                    source_start,
                    source_end,
                });
            }
        }
        [GTOK_U32_FLAG] => {
            let mut buffer = [0; 4 + 4 + 4 + 4]; // u32 + u32 + u32 + u32
            while let Ok(()) = reader.read_exact(&mut buffer) {
                let id = u32::from_le_bytes([buffer[0], buffer[1], buffer[2], buffer[3]]);
                let chrom_id =
                    u32::from_le_bytes([buffer[4], buffer[5], buffer[6], buffer[7]]) as u16;
                let source_start =
                    u32::from_le_bytes([buffer[8], buffer[9], buffer[10], buffer[11]]);
                let source_end =
                    u32::from_le_bytes([buffer[12], buffer[13], buffer[14], buffer[15]]);
                pointers.push(TokenizedRegionPointer {
                    id,
                    chrom_id,
                    source_start,
                    source_end,
                });
            }
        }
        _ => {
            anyhow::bail!("Invalid data format flag found in gtokp file")
        }
    }

    Ok(pointers)
}
