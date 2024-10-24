use std::fmt::Display;

use crate::common::models::{Region, Universe};

#[derive(Eq, PartialEq, Clone)]
pub struct TokenizedRegion<'a> {
    pub universe: &'a Universe,
    pub id: u32,
    pub genomic_pos: u32
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
