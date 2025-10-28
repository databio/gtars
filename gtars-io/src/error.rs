use std::io;
use thiserror::Error;

/// Error type for gtars-io operations.
#[derive(Error, Debug)]
pub enum GtokError {
    /// IO error occurred during file operations.
    #[error("IO error: {0}")]
    Io(#[from] io::Error),

    /// Failed to write GTOK header to file.
    #[error("Failed to write GTOK header to file")]
    HeaderWrite,

    /// Failed to write GTOK size flag to file.
    #[error("Failed to write GTOK size flag to file")]
    SizeFlagWrite,

    /// Failed to write token bytes to file.
    #[error("Failed to write token bytes to file")]
    TokenWrite,

    /// Failed to create parent directories for file.
    #[error("Failed to create parent directories for file")]
    ParentDirectoryCreation,

    /// File is not a valid GTOK file.
    #[error("File doesn't appear to be a valid .gtok file")]
    InvalidGtokFile,

    /// Invalid data format flag found in GTOK file.
    #[error("Invalid data format flag found in gtok file: {0:#x}")]
    InvalidFormatFlag(u8),
}

/// Result type alias for gtars-io operations.
pub type Result<T> = std::result::Result<T, GtokError>;
