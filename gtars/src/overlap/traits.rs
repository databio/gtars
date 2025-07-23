use num_traits::{PrimInt, Unsigned};

pub use crate::common::models::Interval;

pub trait Overlapper<I, D>: Send + Sync
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    D: Eq + Clone + Send + Sync,
{
    fn build(intervals: Vec<Interval<I, D>>) -> Self  where Self: Sized;
    // TODO: how to make this an iterator instead? itearting
    // over found intervals, i believe, is more performant.
    fn find(&self, start: I, end: I) -> Vec<Interval<I, D>>;
}