use crate::common::models::region_set::RegionSet;

pub struct BedSet {
    pub region_sets: Vec<RegionSet>,
}

impl IntoIterator for BedSet {
    type Item = RegionSet;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.region_sets.into_iter()
    }
}