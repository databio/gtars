use crate::common::models::Universe;

#[derive(Eq, PartialEq, Clone)]
pub struct TokenizedRegion<'a>  {
    pub universe: &'a Universe,
    pub id: u32,
}