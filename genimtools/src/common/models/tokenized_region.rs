use crate::common::models::Region;

#[derive(Eq, PartialEq, Clone)]
pub struct TokenizedRegion {
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub id: u32,
}

impl From<TokenizedRegion> for Region {
    fn from(val: TokenizedRegion) -> Self {
        Region {
            chr: val.chr,
            start: val.start,
            end: val.end,
        }
    }
}

impl From<&TokenizedRegion> for Region {
    fn from(value: &TokenizedRegion) -> Self {
        Region {
            chr: value.chr,
            start: value.start,
            end: value.end,
        }
    }
}
