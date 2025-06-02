use num_traits::{PrimInt, Unsigned};
use crate::tokenizers::intersect::intervals::Interval;

/// anything that can answer “which stored intervals touch [start,end)?”
pub trait Intersect<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq   + Clone   + Send + Sync,
{
    type Iter<'a>: Iterator<Item = &'a Interval<I, T>> + 'a
    where
        Self: 'a;

    fn find_overlaps<'a>(&'a self, start: I, end: I) -> Self::Iter<'a>;
}

pub trait IntervalCount {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;
}