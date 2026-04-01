use std::path::Path;
use std::collections::HashMap;

use gtars_core::models::region::Region;
use gtars_tokenizers::Tokenizer;

use crate::sparse_vector::SparseVector;

pub struct Bm25 {
    avg_doc_length: f32,
    b: f32,
    k: f32,
    tokenizer: Tokenizer,
}

pub struct Bm25Builder {
    b: f32,
    k: f32,
    avg_doc_length: f32,
    tokenizer: Option<Tokenizer>,
}

impl Bm25Builder {
    /// Build the BM25 model from the builder.
    ///
    /// # Panics
    /// Panics if no tokenizer/vocabulary has been provided.
    pub fn build(self) -> Bm25 {
        let tokenizer = self
            .tokenizer
            .expect("A tokenizer or vocabulary must be provided via with_tokenizer() or with_vocab()");

        Bm25 {
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

impl Default for Bm25Builder {
    fn default() -> Self {
        Self {
            b: 0.75,
            k: 1.0,
            avg_doc_length: 1_000.0,
            tokenizer: None,
        }
    }
}

impl Bm25 {
    /// Create a new Bm25Builder.
    pub fn builder() -> Bm25Builder {
        Bm25Builder::default()
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

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[fixture]
    fn peaks_path() -> String {
        "../tests/data/tokenizers/peaks.bed".to_string()
    }

    #[fixture]
    fn bm25(peaks_path: String) -> Bm25 {
        Bm25::builder()
            .with_vocab(&peaks_path)
            .with_k(1.5)
            .with_b(0.75)
            .with_avg_doc_length(1_000.0)
            .build()
    }

    #[rstest]
    fn test_builder_defaults() {
        let builder = Bm25Builder::default();
        assert_eq!(builder.k, 1.0);
        assert_eq!(builder.b, 0.75);
        assert_eq!(builder.avg_doc_length, 1_000.0);
    }

    #[rstest]
    fn test_builder_custom_params(peaks_path: String) {
        let model = Bm25::builder()
            .with_vocab(&peaks_path)
            .with_k(2.0)
            .with_b(0.5)
            .with_avg_doc_length(500.0)
            .build();

        assert_eq!(model.k(), 2.0);
        assert_eq!(model.b(), 0.5);
        assert_eq!(model.avg_doc_length(), 500.0);
    }

    #[rstest]
    fn test_vocab_size(bm25: Bm25) {
        // peaks.bed has 25 regions + 7 special tokens
        assert!(bm25.vocab_size() > 25);
    }

    #[rstest]
    fn test_tokenize_overlapping_regions(bm25: Bm25) {
        // query a region that overlaps the first entry in peaks.bed: chr17 7915738 7915777
        let regions = vec![Region {
            chr: "chr17".to_string(),
            start: 7915700,
            end: 7915800,
            rest: None,
        }];

        let token_ids = bm25.tokenize(&regions);
        assert!(!token_ids.is_empty(), "should find overlapping tokens");
    }

    #[rstest]
    fn test_tokenize_no_overlap(bm25: Bm25) {
        // query a region on a chromosome with no vocab entries
        let regions = vec![Region {
            chr: "chrZ".to_string(),
            start: 0,
            end: 100,
            rest: None,
        }];

        let token_ids = bm25.tokenize(&regions);
        assert!(token_ids.is_empty(), "should find no tokens for non-existent chrom");
    }

    #[rstest]
    fn test_embed_produces_sparse_vector(bm25: Bm25) {
        let regions = vec![Region {
            chr: "chr17".to_string(),
            start: 7915700,
            end: 7915800,
            rest: None,
        }];

        let sv = bm25.embed(&regions);
        assert!(!sv.is_empty());
        assert_eq!(sv.indices.len(), sv.values.len());
    }

    #[rstest]
    fn test_embed_empty_input(bm25: Bm25) {
        let regions: Vec<Region> = vec![];
        let sv = bm25.embed(&regions);
        assert!(sv.is_empty());
    }

    #[rstest]
    fn test_embed_no_overlap_returns_empty(bm25: Bm25) {
        let regions = vec![Region {
            chr: "chrZ".to_string(),
            start: 0,
            end: 100,
            rest: None,
        }];

        let sv = bm25.embed(&regions);
        assert!(sv.is_empty());
    }

    #[rstest]
    fn test_embed_values_are_positive(bm25: Bm25) {
        // query multiple overlapping regions
        let regions = vec![
            Region { chr: "chr17".to_string(), start: 7915700, end: 7915800, rest: None },
            Region { chr: "chr6".to_string(), start: 157381091, end: 157381200, rest: None },
        ];

        let sv = bm25.embed(&regions);
        for val in &sv.values {
            assert!(*val > 0.0, "Bm25 TF scores should be positive");
        }
    }

    #[rstest]
    fn test_embed_repeated_regions_increase_tf(bm25: Bm25) {
        let region = Region {
            chr: "chr17".to_string(),
            start: 7915700,
            end: 7915800,
            rest: None,
        };

        let sv_single = bm25.embed(&[region.clone()]);
        let sv_repeated = bm25.embed(&[region.clone(), region.clone(), region]);

        // with repeated regions, the tf score should be higher (but sublinear due to saturation)
        assert!(!sv_single.is_empty());
        assert!(!sv_repeated.is_empty());

        let val_single = sv_single.values[0];
        let val_repeated = sv_repeated.values[0];
        assert!(val_repeated > val_single, "repeated terms should have higher TF score");
    }

    #[rstest]
    fn test_with_tokenizer(peaks_path: String) {
        let tokenizer = Tokenizer::from_bed(&peaks_path).unwrap();
        let model = Bm25::builder()
            .with_tokenizer(tokenizer)
            .build();

        assert!(model.vocab_size() > 0);
    }

    #[rstest]
    #[should_panic(expected = "A tokenizer or vocabulary must be provided")]
    fn test_builder_panics_without_tokenizer() {
        Bm25Builder::default().build();
    }
}