use std::path::Path;
use std::collections::HashMap;

use gtars_core::models::region::Region;
use gtars_tokenizers::Tokenizer;

use crate::sparse_vector::SparseVector;

pub struct BM25 {
    avg_doc_length: f32,
    b: f32,
    k: f32,
    tokenizer: Tokenizer,
}

pub struct BM25Builder {
    b: f32,
    k: f32,
    avg_doc_length: f32,
    tokenizer: Option<Tokenizer>,
}

impl BM25Builder {
    /// Build the BM25 model from the builder.
    ///
    /// # Panics
    /// Panics if no tokenizer/vocabulary has been provided.
    pub fn build(self) -> BM25 {
        let tokenizer = self
            .tokenizer
            .expect("A tokenizer or vocabulary must be provided via with_tokenizer() or with_vocab()");

        BM25 {
            avg_doc_length: self.avg_doc_length,
            b: self.b,
            k: self.k,
            tokenizer,
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

    /// Set the tokenizer directly from an existing `Tokenizer` instance.
    pub fn with_tokenizer(mut self, tokenizer: Tokenizer) -> Self {
        self.tokenizer = Some(tokenizer);
        self
    }

    /// Load a vocabulary from a file path.
    ///
    /// The path can be:
    /// 1. A BED file on disk (optionally gzipped)
    /// 2. A TOML config file pointing to a universe
    ///
    /// The file type is auto-detected.
    pub fn with_vocab<P: AsRef<Path>>(mut self, vocab: P) -> Self {
        let tokenizer = Tokenizer::from_auto(vocab)
            .expect("Failed to load vocabulary. Ensure the path points to a valid .bed, .bed.gz, or .toml file.");
        self.tokenizer = Some(tokenizer);
        self
    }
}

impl Default for BM25Builder {
    fn default() -> Self {
        Self {
            b: 0.75,
            k: 1.0,
            avg_doc_length: 1_000.0,
            tokenizer: None,
        }
    }
}

impl BM25 {
    /// Create a new BM25Builder.
    pub fn builder() -> BM25Builder {
        BM25Builder::default()
    }

    /// Tokenize a set of regions into token IDs.
    ///
    /// This performs the overlap query against the vocabulary and returns
    /// the token IDs of all overlapping vocabulary regions.
    pub fn tokenize(&self, regions: &[Region]) -> Vec<u32> {
        let unk_id = self.tokenizer.get_unk_token_id();
        self.tokenizer
            .encode(regions)
            .unwrap_or_default()
            // filter OOV tokens (those that don't overlap any vocab region) by removing the unk_id
            .into_iter().filter(|&id| id != unk_id)
            .collect()
    }

    /// Returns a reference to the internal tokenizer.
    pub fn tokenizer(&self) -> &Tokenizer {
        &self.tokenizer
    }

    /// Returns the vocabulary size (number of regions in the vocabulary).
    pub fn vocab_size(&self) -> usize {
        self.tokenizer.get_vocab_size()
    }

    /// Returns the `k` parameter (term frequency saturation).
    pub fn k(&self) -> f32 {
        self.k
    }

    /// Returns the `b` parameter (document length normalization).
    pub fn b(&self) -> f32 {
        self.b
    }

    /// Returns the assumed average document length.
    pub fn avg_doc_length(&self) -> f32 {
        self.avg_doc_length
    }

    /// Encode a set of regions into a BM25 sparse vector.
    ///
    /// This tokenizes the regions and computes BM25-like term frequency
    /// scores for each token. The indices are token IDs and the values
    /// are the BM25 term-frequency component scores.
    /// Formula from: https://en.wikipedia.org/wiki/Okapi_BM25
    pub fn embed(&self, regions: &[Region]) -> SparseVector {
        // tokenize to get token IDs
        let token_ids = self.tokenize(regions);

        if token_ids.is_empty() {
            return SparseVector::empty();
        }

        // count term frequencies
        let mut tf_map: HashMap<u32, u32> = HashMap::new();
        for id in &token_ids {
            *tf_map.entry(*id).or_insert(0) += 1;
        }

        let doc_length = token_ids.len() as f32;

        // compute bm25 term-frequency scores
        let mut indices: Vec<u32> = Vec::with_capacity(tf_map.len());
        let mut values: Vec<f32> = Vec::with_capacity(tf_map.len());

        for (token_id, raw_tf) in tf_map {
            let tf = raw_tf as f32;
            let tf_score = (tf * (self.k + 1.0))
                / (tf + self.k * (1.0 - self.b + self.b * (doc_length / self.avg_doc_length)));

            indices.push(token_id);
            values.push(tf_score);
        }

        SparseVector::new(indices, values)
    }
}