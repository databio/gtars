use md5::{Digest, Md5};
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

    ///
    /// Calculate digest for the Region
    ///
    pub fn digest(&self) -> String {
        let digest_string = format!("{},{},{}", self.chr, self.start, self.end);

        let mut hasher = Md5::new();
        hasher.update(digest_string);
        let chrom_hash = hasher.finalize();
        format!("{:x}", chrom_hash)
    }

    ///
    /// Calculate mid point
    ///
    pub fn mid_point(&self) -> u32 {
        self.start + self.width() / 2
    }
}

impl Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_string())
    }
}

// TODO:
// impl Display for ChromosomeStats {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//
//     }
// }
