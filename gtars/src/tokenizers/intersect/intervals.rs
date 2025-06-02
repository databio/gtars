use std::cmp::Ordering;
use num_traits::{
    identities::zero,
    PrimInt, Unsigned,
};

/// Represent a range from [start, stop)
/// Inclusive start, exclusive of stop
#[derive(Eq, Debug, Clone)]
pub struct Interval<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    pub start: I,
    pub stop: I,
    pub val: T,
}

impl<I, T> Interval<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// Compute the intsect between two intervals
    #[inline]
    pub fn intersect(&self, other: &Interval<I, T>) -> I {
        std::cmp::min(self.stop, other.stop)
            .checked_sub(std::cmp::max(&self.start, &other.start))
            .unwrap_or_else(zero::<I>)
    }

    /// Check if two intervals overlap
    #[inline]
    pub fn overlap(&self, start: I, stop: I) -> bool {
        self.start < stop && self.stop > start
    }
}

impl<I, T> Ord for Interval<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    #[inline]
    fn cmp(&self, other: &Interval<I, T>) -> Ordering {
        match self.start.cmp(&other.start) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => self.stop.cmp(&other.stop),
        }
    }
}

impl<I, T> PartialOrd for Interval<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<I, T> PartialEq for Interval<I, T>
where
    I: PrimInt + Unsigned + Ord + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    #[inline]
    fn eq(&self, other: &Interval<I, T>) -> bool {
        self.start == other.start && self.stop == other.stop
    }
}