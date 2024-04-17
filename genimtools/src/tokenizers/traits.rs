use crate::common::models::region::Region;
use crate::common::models::region_set::RegionSet;
// use crate::common::models::bed_set::BedSet;
use crate::common::consts::{PAD_CHR, PAD_END, PAD_START, MASK_CHR, MASK_START, MASK_END, UNKNOWN_CHR, UNKNOWN_END, UNKNOWN_START};
use crate::common::models::tokenized_regionset::TokenizedRegionSet;

pub trait Tokenizer {
    /// Tokenize a region into the vocabulary of the tokenizer
    ///
    /// # Arguments
    /// - `region` - the region to be tokenized
    ///
    /// # Returns
    /// A new region that corresponds to a region in the tokenizers vocab (or universe).
    ///
    fn tokenize_region(&self, region: &Region) -> TokenizedRegionSet;

    /// Tokenize a region set into the vocabulary of the tokenizer
    ///
    /// # Arguments
    /// - `region_set` - the region set to be tokenized
    ///
    /// # Returns
    /// A vector of regions that correspond to regions in the tokenizers vocab (or universe).
    ///
    fn tokenize_region_set(&self, region_set: &RegionSet) -> Option<TokenizedRegionSet>;

    fn unknown_token(&self) -> Region {
        Region {
            chr: UNKNOWN_CHR.to_string(),
            start: UNKNOWN_START as u32,
            end: UNKNOWN_END as u32,
        }
    }

    fn padding_token(&self) -> Region {
        Region {
            chr: PAD_CHR.to_string(),
            start: PAD_START as u32,
            end: PAD_END as u32,
        }
    }

    fn mask_token(&self) -> Region {
        Region {
            chr: MASK_CHR.to_string(),
            start: MASK_START as u32,
            end: MASK_END as u32,
        }
    }


}
