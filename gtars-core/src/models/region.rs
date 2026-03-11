use md5::{Digest, Md5};
use std::fmt::{self, Display};

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

    fn region(chr: &str, start: u32, end: u32) -> Region {
        Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    #[test]
    fn test_distance_to_overlap() {
        let a = region("chr1", 100, 200);
        let b = region("chr1", 150, 250);
        assert_eq!(a.distance_to(&b), 0);
        assert_eq!(b.distance_to(&a), 0);
    }

    #[test]
    fn test_distance_to_contained() {
        let a = region("chr1", 100, 300);
        let b = region("chr1", 150, 200);
        assert_eq!(a.distance_to(&b), 0);
        assert_eq!(b.distance_to(&a), 0);
    }

    #[test]
    fn test_distance_to_adjacent() {
        let a = region("chr1", 100, 200);
        let b = region("chr1", 200, 300);
        assert_eq!(a.distance_to(&b), 0);
        assert_eq!(b.distance_to(&a), 0);
    }

    #[test]
    fn test_distance_to_gap_right() {
        let a = region("chr1", 100, 200);
        let b = region("chr1", 250, 300);
        assert_eq!(a.distance_to(&b), 50);
    }

    #[test]
    fn test_distance_to_gap_left() {
        let a = region("chr1", 250, 300);
        let b = region("chr1", 100, 200);
        assert_eq!(a.distance_to(&b), 50);
    }

    #[test]
    fn test_distance_to_symmetric() {
        let a = region("chr1", 100, 200);
        let b = region("chr1", 300, 400);
        assert_eq!(a.distance_to(&b), b.distance_to(&a));
    }
}
