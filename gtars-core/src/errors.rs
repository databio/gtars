use thiserror::Error;

#[derive(Error, Debug)]
pub enum RegionSetError {
    #[error("Can't read file: {0}")]
    FileReadError(String),

    #[error("Can't get file from path or url: {0}")]
    InvalidPathOrUrl(String),

    #[error("BEDbase identifier is not valid UTF-8: {0}")]
    InvalidBedbaseIdentifier(String),

    #[error("Can't get file from BEDbase: {0}")]
    BedbaseFetchError(String),

    #[error("Error parsing region: {0}")]
    RegionParseError(String),

    #[error("Corrupted file. 0 regions found in the file: {0}")]
    EmptyRegionSet(String),

    #[error("File not found and HTTP feature not enabled: {0}")]
    HttpFeatureDisabled(String),

    #[error("BigBed error: {0}")]
    BigBedError(String),

    #[error(transparent)]
    Io(#[from] std::io::Error),
}
