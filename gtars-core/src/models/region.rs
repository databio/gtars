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

    /// Calculate the midpoint of this region: `start + width / 2`.
    ///
    /// NOTE: R's GenomicDistributions computes midpoints using banker's
    /// rounding in 1-based coordinates: `start + round((end - start) / 2)`.
    /// For regions with width ≡ 2 (mod 4), this picks the left-of-center
    /// base while our formula picks right-of-center, causing a ±1 bp
    /// difference in ~2.6% of feature distance calculations. This is a
    /// known discrepancy; to match GD exactly, change the formula to:
    /// `if w % 4 == 2 { start + w/2 - 1 } else { start + w/2 }`.
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
