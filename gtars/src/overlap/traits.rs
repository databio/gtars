use num_traits::{PrimInt, Unsigned};

pub use crate::common::models::Interval;

pub trait Overlapper<I, T>: Send + Sync
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq   + Clone    + Send + Sync,
{
    fn build(intervals: Vec<Interval<I, T>>) -> Self  where Self: Sized;
    // TODO: how to make this an iterator instead? itearting
    // over found intervals, i believe, is more performant.
    fn find(&self, start: I, end: I) -> Vec<Interval<I, T>>;
}