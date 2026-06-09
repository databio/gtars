//! Errors produced by HGVS parsing.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum HgvsError {
    #[error("HGVS parse error at position {pos}: {msg} (in: {input:?})")]
    Parse {
        input: String,
        pos: usize,
        msg: String,
    },
}
