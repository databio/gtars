use std::fmt::{self, Display};

#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub struct Region {
    pub chr: String,
    pub start: u32,
    pub end: u32,

    pub rest: String,
}

impl Region {
    pub fn width(self) -> u32 {
        self.end - self.start
    }
}

impl Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Region: {}\t{}\t{}\t{}",
            self.chr, self.start, self.end, self.rest
        )
    }
}
