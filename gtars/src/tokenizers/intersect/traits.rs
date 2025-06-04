use num_traits::{PrimInt, Unsigned};

use crate::tokenizers::intersect::Interval;

/// anything that can answer “which stored intervals touch [start,end)?”
pub trait Intersect<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq   + Clone   + Send + Sync,
{
    type Iter<'a>: Iterator<Item = &'a Interval<I, T>> + 'a
    where
        I: 'a,
        T: 'a,
        Self: 'a;

    fn find_overlaps(&self, start: I, end: I) -> Self::Iter<'_>;
}

pub trait IntervalCount {
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;
}