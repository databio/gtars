use std::collections::HashMap;

use crate::common::models::region::Region;
use crate::common::models::region_set::RegionSet;
use crate::common::models::tokenized_regionset::TokenizedRegionSet;
use crate::tokenizers::special_tokens::SpecialToken;

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
    fn tokenize_region_set(&self, region_set: &RegionSet) -> TokenizedRegionSet;
}

pub trait SpecialTokens {
    fn unknown_token(&self) -> Region;
    fn padding_token(&self) -> Region;
    fn mask_token(&self) -> Region;
    fn cls_token(&self) -> Region;
    fn bos_token(&self) -> Region;
    fn eos_token(&self) -> Region;
    fn sep_token(&self) -> Region;
    fn special_tokens_map(&self) -> HashMap<SpecialToken, Region> {
        let mut map: HashMap<SpecialToken, Region> = HashMap::new();

        map.insert(SpecialToken::Unk, self.unknown_token());
        map.insert(SpecialToken::Pad, self.padding_token());
        map.insert(SpecialToken::Mask, self.mask_token());
        map.insert(SpecialToken::Cls, self.cls_token());
        map.insert(SpecialToken::Bos, self.bos_token());
        map.insert(SpecialToken::Eos, self.eos_token());
        map.insert(SpecialToken::Sep, self.sep_token());

        map
    }
}

pub trait AtttentionMask: SpecialTokens {
    fn attention_mask(&self) -> Vec<u32>;
}
