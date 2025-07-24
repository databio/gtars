// https://github.com/sstadick/rust-lapper/blob/7e3904daed85181f1faa39b15f51935f13945976/src/lib.rs#L92
use num_traits::{identities::zero, PrimInt, Unsigned};
use std::cmp::Ordering::{self};

/// Represent a range from [start, end)
/// Inclusive start, exclusive of end
#[derive(Eq, Debug, Clone)]
pub struct Interval<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    pub start: I,
    pub end: I,
    pub val: T,
}

impl<I, T> Ord for Interval<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    #[inline]
    fn cmp(&self, other: &Interval<I, T>) -> Ordering {
        match self.start.cmp(&other.start) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => self.end.cmp(&other.end),
        }
    }
}

impl<I, T> Interval<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    /// Compute the intsect between two intervals
    #[inline]
    pub fn intersect(&self, other: &Interval<I, T>) -> I {
        std::cmp::min(self.end, other.end)
            .checked_sub(std::cmp::max(&self.start, &other.start))
            .unwrap_or_else(zero::<I>)
    }

    /// Check if two intervals overlap
    #[inline]
    pub fn overlap(&self, start: I, end: I) -> bool {
        self.start < end && self.end > start
    }
}

impl<I, T> PartialOrd for Interval<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<I, T> PartialEq for Interval<I, T>
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    #[inline]
    fn eq(&self, other: &Interval<I, T>) -> bool {
        self.start == other.start && self.end == other.end
    }
}
