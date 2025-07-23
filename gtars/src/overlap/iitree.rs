///
/// WIP
///
use std::cmp::Ordering::{self};

use num_traits::{PrimInt, Unsigned};

use super::Overlapper;
use crate::common::models::Interval;

/// An internal wrapper that stores the extra `max` field required by the iitree.
#[derive(Clone)]
struct Node<I, T>
where
    I: PrimInt + Unsigned + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    start: I,
    end: I,
    max: I,
    val: T,
}

impl<I, T> Ord for Node<I, T>
where
    I: PrimInt + Unsigned + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        self.start.cmp(&other.start)
    }
}
impl<I, T> PartialOrd for Node<I, T>
where
    I: PrimInt + Unsigned + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl<I, T> Eq for Node<I, T>
where
    I: PrimInt + Unsigned + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
}
impl<I, T> PartialEq for Node<I, T>
where
    I: PrimInt + Unsigned + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.start == other.start && self.end == other.end
    }
}

///
///
/// The Implicit interval tree, described in: https://academic.oup.com/bioinformatics/article/37/9/1315/5910546
///
pub struct IITree<I, T>
where
    I: PrimInt + Unsigned + Clone + Send + Sync,
    T: Eq + Clone + Send + Sync,
{
    max_level: usize,
    nodes: Vec<Node<I, T>>,
}
