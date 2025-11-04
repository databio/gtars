use num_traits::{PrimInt, Unsigned};

pub use gtars_core::models::Interval;

/// A trait for data structures that efficiently find overlapping genomic intervals.
///
/// This trait defines the interface for interval overlap query data structures. Implementors
/// provide efficient algorithms for finding all intervals that overlap with a query range.
///
/// # Type Parameters
///
/// * `I` - The integer type used for interval coordinates (e.g., `u32`, `usize`). Must be an
///   unsigned primitive integer type that is thread-safe.
/// * `T` - The type of value associated with each interval (e.g., gene name, peak score).
///   Must be cloneable, comparable, and thread-safe.
///
/// # Implementations
/// The crate provides two main implementations:
///
/// * [`AIList`](crate::AIList) - Augmented Interval List, optimized for genomic data with
///   high-coverage regions. Uses a decomposition strategy for efficient queries.
/// * [`Bits`](crate::Bits) - Binary Interval Search, uses binary search for fast overlap
///   detection. Particularly efficient for sorted sequential queries.
pub trait Overlapper<I, T>: Send + Sync
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// Constructs a new overlap query data structure from a vector of intervals.
    ///
    /// This method takes ownership of the intervals and builds an optimized data structure
    /// for efficient overlap queries. The construction time varies by implementation but is
    /// typically O(n log n) where n is the number of intervals.
    ///
    /// # Arguments
    ///
    /// * `intervals` - A vector of intervals to index. The intervals may be modified
    ///   (e.g., sorted) during construction.
    ///
    /// # Returns
    ///
    /// A new instance of the implementing data structure, ready for overlap queries.
    fn build(intervals: Vec<Interval<I, T>>) -> Self
    where
        Self: Sized;

    /// Finds all intervals that overlap with the query range `[start, end)`.
    ///
    /// Returns a vector containing clones of all intervals that have any overlap with the
    /// specified range. The query range is half-open: `start` is inclusive, `end` is exclusive.
    ///
    /// Two intervals overlap if they share at least one position. Formally, interval `[a, b)`
    /// overlaps with query `[start, end)` if `a < end` and `b > start`.
    ///
    /// # Arguments
    ///
    /// * `start` - The start position of the query range (inclusive)
    /// * `end` - The end position of the query range (exclusive)
    ///
    /// # Returns
    ///
    /// A vector of intervals that overlap with the query range. The order of results is
    /// implementation-dependent.
    ///
    /// # Performance
    ///
    /// For large result sets, consider using [`find_iter`](Self::find_iter) to avoid
    /// allocating and cloning all results at once.
    /// TODO: this should be a default implementation that we collect
    /// the below find_iter()
    fn find(&self, start: I, end: I) -> Vec<Interval<I, T>>;

    /// Returns an iterator over all intervals that overlap with the query range `[start, end)`.
    ///
    /// This method is more memory-efficient than [`find`](Self::find) as it returns references
    /// to intervals rather than cloning them. Use this when processing large result sets or when
    /// you don't need owned copies of the intervals.
    ///
    /// # Arguments
    ///
    /// * `start` - The start position of the query range (inclusive)
    /// * `end` - The end position of the query range (exclusive)
    ///
    /// # Returns
    ///
    /// A boxed iterator that yields references to overlapping intervals. The iterator is
    /// lazy and computes results on demand.
    fn find_iter<'a>(
        &'a self,
        start: I,
        end: I,
    ) -> Box<dyn Iterator<Item = &'a Interval<I, T>> + 'a>;
}
