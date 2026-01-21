use thiserror::Error;

#[derive(Error, Debug)]
pub enum BedClassifierError {
    #[error("Failed to convert RegionSet to DataFrame")]
    DataFrameConversionError,
    #[error("Invalid column data at index {0}")]
    InvalidColumnData(usize),
    #[error("Regex compilation error: {0}")]
    RegexError(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

#[derive(Error, Debug)]
pub enum GtarsGenomicDistError {
    #[error("Custom Error.")]
    CustomError(String),
    #[error("Error getting sequence for region.")]
    GCContentError(String, u32, u32, String),
    #[error("No TSS's found for region. Double-check your index!")]
    TSSContentError(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}
