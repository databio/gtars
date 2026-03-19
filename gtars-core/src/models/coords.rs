/// Coordinate convention for genomic intervals.
///
/// Different ecosystems use different coordinate systems:
/// - BED/Python: 0-based half-open, floor midpoint
/// - GRanges/R: 1-based closed, banker's rounding midpoint
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CoordinateMode {
    #[default]
    Bed,
    GRanges,
}
