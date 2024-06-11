#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub struct Region {
    pub chr: String,
    pub start: u32,
    pub end: u32,
}
