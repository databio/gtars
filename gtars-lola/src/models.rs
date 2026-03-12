/// Direction of the enrichment test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Test for enrichment: P(X >= a), alternative="greater"
    Enrichment,
    /// Test for depletion: P(X <= a), alternative="less"
    Depletion,
}

impl Default for Direction {
    fn default() -> Self {
        Direction::Enrichment
    }
}

/// Configuration for a LOLA run.
#[derive(Debug, Clone)]
pub struct LolaConfig {
    /// Minimum base-pair overlap to count as overlapping (default 1).
    pub min_overlap: i32,
    /// Test direction: enrichment or depletion.
    pub direction: Direction,
}

impl Default for LolaConfig {
    fn default() -> Self {
        LolaConfig {
            min_overlap: 1,
            direction: Direction::Enrichment,
        }
    }
}

/// A 2x2 contingency table for a single (user set, DB set) pair.
///
/// ```text
///                     In DB set    Not in DB set
/// In user set            a              c
/// Not in user set        b              d
/// ```
#[derive(Debug, Clone)]
pub struct ContingencyTable {
    /// Overlap count between user set and DB set (support).
    pub a: u64,
    /// Universe-DB overlap minus user-DB overlap.
    pub b: u64,
    /// User set size minus user-DB overlap.
    pub c: u64,
    /// Remainder: universe_size - a - b - c.
    pub d: u64,
}

/// Result for a single (user set, DB set) pair.
#[derive(Debug, Clone)]
pub struct LolaResult {
    /// Index of the user set (0-based).
    pub user_set: usize,
    /// Index of the DB set (0-based).
    pub db_set: usize,
    /// -log10(p-value) from Fisher's exact test.
    pub p_value_log: f64,
    /// Odds ratio from Fisher's exact test.
    pub odds_ratio: f64,
    /// Overlap count (contingency table 'a').
    pub support: u64,
    /// Rank by p-value (1-based, ascending by significance).
    pub rnk_pv: usize,
    /// Rank by odds ratio (1-based, descending).
    pub rnk_or: usize,
    /// Rank by support (1-based, descending).
    pub rnk_sup: usize,
    /// max(rnk_pv, rnk_or, rnk_sup) — worst rank across metrics.
    pub max_rnk: usize,
    /// mean(rnk_pv, rnk_or, rnk_sup).
    pub mean_rnk: f64,
    /// Contingency table values.
    pub b: u64,
    pub c: u64,
    pub d: u64,
    /// FDR-adjusted p-value (Benjamini-Hochberg). None until computed.
    pub q_value: Option<f64>,
    /// DB set filename (from IGD file_info).
    pub filename: String,
}
