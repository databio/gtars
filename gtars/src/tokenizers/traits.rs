use std::collections::HashMap;
use std::path::Path;

use anyhow::Result;

use crate::common::models::region::Region;
use crate::common::models::region_set::RegionSet;
use crate::common::models::tokenized_regionset::TokenizedRegionSet;
use crate::common::models::universe::Universe;
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

    ///
    /// Get the vocabulary size of the tokenizer
    ///
    /// # Returns
    /// The size of the vocabulary as usize
    fn vocab_size(&self) -> usize;

    ///
    /// Get the universe of the tokenizer
    ///
    /// # Returns
    /// A reference to the universe of the tokenizer
    fn get_universe(&self) -> &Universe;

    ///
    /// Export the tokenizer to a toml file
    ///
    /// # Arguments
    /// - `path` - the path to the toml file
    ///
    /// # Returns
    /// A Result
    fn export(&self, path: &Path) -> Result<()>;
}

pub trait SingleCellTokenizer {
    ///
    /// Tokenize an AnnData object, this is single-cell data
    ///
    /// # Arguments
    /// - `anndata` - the path to the AnnData object
    ///
    /// # Returns
    /// A vector of TokenizedRegionSets
    fn tokenize_anndata(&self, anndata: &Path) -> Result<Vec<TokenizedRegionSet>>;
}

pub trait SpecialTokens {
    fn unknown_token(&self) -> Region;
    fn padding_token(&self) -> Region;
    fn mask_token(&self) -> Region;
    fn cls_token(&self) -> Region;
    fn bos_token(&self) -> Region;
    fn eos_token(&self) -> Region;
    fn sep_token(&self) -> Region;
    fn unknown_token_id(&self) -> u32;
    fn padding_token_id(&self) -> u32;
    fn mask_token_id(&self) -> u32;
    fn cls_token_id(&self) -> u32;
    fn bos_token_id(&self) -> u32;
    fn eos_token_id(&self) -> u32;
    fn sep_token_id(&self) -> u32;
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
    fn attention_mask(&self, tokens: &TokenizedRegionSet) -> Vec<u8> {
        let mut mask: Vec<u8> = Vec::with_capacity(tokens.len());
        let pad_token = self.padding_token();

        for token in tokens {
            if token.into_region() == pad_token {
                mask.push(1)
            } else {
                mask.push(0)
            }
        }

        mask
    }
}

pub trait Pad: SpecialTokens {
    fn pad(&self, tokens_list: &mut Vec<TokenizedRegionSet>) {
        let pad_token = self.padding_token_id();
        let longest = tokens_list.iter().map(|t| t.len()).max().unwrap();

        for token in tokens_list.iter_mut() {
            while token.len() < longest {
                token.ids.push(pad_token);
            }
        }
    }
}

pub trait FromPretrained: Tokenizer {
    fn from_pretrained(models: &str) -> Result<Self>
    where
        Self: Sized;
}
