use std::path::PathBuf;

use gtars_overlaprs::{OverlapperType, multi_chrom_overlapper::{IntoMultiChromOverlapper, MultiChromOverlapper}};
use gtars_core::models::region_set::RegionSet;
use gtars_core::models::region::Region;

use crate::sparse_vector::SparseVector;

pub struct BM25 {
    avg_doc_length: f32,
    b: f32,
    k: f32,
    tokenizer: MultiChromOverlapper<u32, Option<String>>,
}

pub struct BM25Builder {
    b: f32,
    k: f32,
    avg_doc_length: f32,
    tokenizer: MultiChromOverlapper<u32, Option<String>>,
}

impl BM25Builder {
    pub fn build(self) -> BM25 {

        BM25 {
            avg_doc_length: self.avg_doc_length,
            b: self.b,
            k: self.k,
            tokenizer: self.tokenizer,
        }
    }

    pub fn with_k(mut self, k: f32) -> Self {
        self.k = k;
        self
    }

    pub fn with_b(mut self, b: f32) -> Self {
        self.b = b;
        self
    }

    pub fn with_avg_doc_length(mut self, avg_doc_length: f32) -> Self {
        self.avg_doc_length = avg_doc_length;
        self
    }

    pub fn with_tokenizer(mut self, tokenizer: MultiChromOverlapper<u32, Option<String>>) -> Self {
        self.tokenizer = tokenizer;
        self
    }

    pub fn with_vocab(mut self, vocab: &PathBuf) -> Self {
        // vocab should be one of three things:
        // 1. a bed file on disk (canbe gzipped)
        // 2. a bed file at a remote URL (can be gzipped)
        // 3. a registry path to a huggingface repo containing `universe.bed` (can be gzipped)
        // We should instantiate a RegionSet from the bed file, and then convert it into a MultiChromOverlapper using the `into_multi_chrom_overlapper` method, with the `OverlapperType::AIList` option.
        self
    }
}

impl Default for BM25Builder {
    fn default() -> Self {
        Self {
            b: 0.75,
            k: 1.0,
            avg_doc_length: 1_000.0,
            tokenizer: RegionSet::from(vec![]).into_multi_chrom_overlapper(OverlapperType::AIList),
        }
    }
}

impl BM25 {
    pub fn encode(&self, query: &[Region]) -> SparseVector {
            // Tokenize the query using the tokenizer
            let mut indices = Vec::new();
            let mut values = Vec::new();
    
            for region in query {
                let tokens = self.tokenizer.find_overlaps(region);
                for token in tokens {
                    indices.push(token);
                    values.push(1.0); // Placeholder value, you can compute actual BM25 scores here
                }
            }
    
            SparseVector { indices, values }
    }
}
