#[derive(Eq, PartialEq, Hash, Debug)]
pub struct Region {
    pub chr: String,
    pub start: u32,
    pub end: u32,
}

impl Clone for Region {
    fn clone(&self) -> Self {
        Region {
            chr: self.chr.to_owned(),
            start: self.start,
            end: self.end,
        }
    }
}
