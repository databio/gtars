use thiserror::Error;

#[derive(Error, Debug)]
pub enum LolaError {
    #[error("Universe is empty (0 regions)")]
    EmptyUniverse,

    #[error("No database sets indexed (IGD has 0 files)")]
    EmptyDatabase,

    #[error("User set {0} is empty (0 regions)")]
    EmptyUserSet(usize),

    #[error("Negative contingency table value: {field} = {value} for db_set {db_set}")]
    NegativeContingency {
        field: &'static str,
        value: i64,
        db_set: usize,
    },

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}
