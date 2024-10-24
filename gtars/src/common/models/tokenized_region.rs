use std::fmt::Display;

use crate::common::models::{Region, Universe};

///
/// The `TokenizedRegionPointer` is a wrapper/convenience struct that
/// lets us store pointers to the original universe and stores some other metadata
/// that will be useful for modeling -- the genomic positional encodings.
pub struct TokenizedRegionPointer {
    /// The id of the token in the universe
    pub id: u32,

    /// The id of the chromosome in the universe
    pub chrom_id: u16,

    /// The original start position of the region in the query
    pub source_start: u32,

    /// The original end position of the region in the query
    pub source_end: u32
}

#[derive(Eq, PartialEq, Clone)]
pub struct TokenizedRegion<'a> {
    pub universe: &'a Universe,
    pub id: u32,
    // these attributes are for positional encodings
    pub source_chr_id: u16,
    pub source_start: u32,
    pub source_end: u32
}

impl From<TokenizedRegion<'_>> for Region {
    fn from(val: TokenizedRegion<'_>) -> Self {
        val.universe.convert_id_to_region(val.id).unwrap()
    }
}

impl TokenizedRegion<'_> {
    pub fn into_region(self) -> Region {
        self.into()
    }

    pub fn chrom_id(&self) -> u16 {
        let r = self.universe.convert_id_to_region(self.id).unwrap();
        self.universe.convert_chrom_to_id(&r.chr).unwrap()
    }
}

impl Display for TokenizedRegion<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let r = self.universe.convert_id_to_region(self.id).unwrap();
        write!(f, "{}:{}-{}", r.chr, r.start, r.end)
    }
}
