use thiserror::Error;

#[derive(Error, Debug)]
pub enum TokenizerConfigError {
    #[error("Custom Error for gtars-genomicdist")]
    CustomError,
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

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
