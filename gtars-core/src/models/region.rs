use md5::{Digest, Md5};
use std::fmt::{self, Display};

use super::coords::CoordinateMode;

///
/// Region struct, representation of one Region in RegionSet files
///
#[derive(Eq, PartialEq, Hash, Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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

    /// Calculate midpoint using the specified coordinate convention.
    ///
    /// - `Bed` (default): floor division → `start + width / 2`
    /// - `GRanges`: banker's rounding in 1-based coords →
    ///   `if w % 4 == 2 { start + w/2 - 1 } else { start + w/2 }`
    pub fn mid_point_with_mode(&self, mode: CoordinateMode) -> u32 {
        match mode {
            CoordinateMode::Bed => self.start + self.width() / 2,
            CoordinateMode::GRanges => {
                let w = self.width();
                if w % 4 == 2 {
                    self.start + w / 2 - 1
                } else {
                    self.start + w / 2
                }
            }
        }
    }

    /// Gap distance between two regions.
    ///
    /// Returns 0 if the regions overlap, otherwise returns the positive
    /// gap (in bases) between the closer edges of the two regions.
    pub fn distance_to(&self, other: &Region) -> i64 {
        if self.start < other.end && other.start < self.end {
            0i64
        } else if other.end <= self.start {
            (self.start as i64) - (other.end as i64)
        } else {
            (other.start as i64) - (self.end as i64)
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn region(start: u32, end: u32) -> Region {
        Region {
            chr: "chr1".to_string(),
            start,
            end,
            rest: None,
        }
    }

    #[rstest]
    #[case::overlap((100, 200), (150, 250), 0)]
    #[case::contained((100, 300), (150, 200), 0)]
    // half-open touching intervals are distance 0
    #[case::adjacent((100, 200), (200, 300), 0)]
    #[case::gap((100, 200), (250, 300), 50)]
    fn test_distance_to(#[case] a: (u32, u32), #[case] b: (u32, u32), #[case] expected: i64) {
        let (a, b) = (region(a.0, a.1), region(b.0, b.1));
        assert_eq!(a.distance_to(&b), expected);
        assert_eq!(b.distance_to(&a), expected);
    }
}
