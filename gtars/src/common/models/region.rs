use std::fmt::{self, Display};

///
/// Region struct, representation of one Region in RegionSet files
///
#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub struct Region {
    pub chr: String,
    pub start: u32,
    pub end: u32,

    pub rest: Option<String>,
}

impl Region {
    ///
    /// Get length of the file
    ///
    pub fn width(&self) -> u32 {
        self.end - self.start
    }

    ///
    /// Get file string of Region
    ///
    pub fn as_string(&self) -> String {
        format!(
            "{}\t{}\t{}{}",
            self.chr,
            self.start,
            self.end,
            self.rest
                .as_deref()
                .map_or(String::new(), |s| format!("\t{}", s)),
        )
    }
}

impl Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_string())
    }
}
