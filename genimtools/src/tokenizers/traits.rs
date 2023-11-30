use crate::common::models::region::Region;
use crate::common::models::region_set::RegionSet;
// use crate::common::models::bed_set::BedSet;
use crate::common::consts::{PAD_CHR, PAD_END, PAD_START, UNKNOWN_CHR, UNKNOWN_END, UNKNOWN_START};
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

    /// Tokenize a bed set into the vocabulary of the tokenizer
    ///
    /// # Arguments
    /// - `bed_set` - the bed set to be tokenized
    ///
    /// # Returns
    /// A vector of vectors of regions that correspond to regions in the tokenizers vocab (or universe).
    ///
    // fn tokenize_bed_set(&self, bed_set: &BedSet) -> Option<Vec<TokenizedRegionSet>> {
    //     let mut tokenized_region_sets = Vec::new();
    //     for region_set in bed_set.into_iter() {
    //         let tokenized_region_set = self.tokenize_region_set(&region_set)?;
    //         tokenized_region_sets.push(tokenized_region_set);
    //     }
    //     Some(tokenized_region_sets)
    // }

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
}
