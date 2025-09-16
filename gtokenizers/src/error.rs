use thiserror::Error;

use super::config::TokenizerConfigError;

#[derive(Error, Debug)]
pub enum TokenizerError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("Invalid special token configuration")]
    InvalidSpecialTokenConfig,
    #[error(transparent)]
    Config(#[from] TokenizerConfigError),
    #[error("Universe error: {0}")]
    UniverseError(#[from] crate::universe::UniverseError),
    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),
}
