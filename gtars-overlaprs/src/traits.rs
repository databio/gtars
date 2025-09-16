use num_traits::{PrimInt, Unsigned};

pub use gtars_core::models::Interval;

pub trait Overlapper<I, T>: Send + Sync
where
    I: PrimInt + Unsigned + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    fn build(intervals: Vec<Interval<I, T>>) -> Self
    where
        Self: Sized;

    fn find(&self, start: I, end: I) -> Vec<Interval<I, T>>;

    fn find_iter<'a>(
        &'a self,
        start: I,
        end: I,
    ) -> Box<dyn Iterator<Item = &'a Interval<I, T>> + 'a>;
}
