use thiserror::Error;

#[derive(Error, Debug)]
pub enum LolaError {
    #[error("Universe is empty (0 regions)")]
    EmptyUniverse,

    #[error("No database sets indexed (IGD has 0 files)")]
    EmptyDatabase,

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}
