/// A sparse vector representation for BM25 embeddings.
///
/// This is designed to be compatible with Qdrant's sparse vector format.
/// Each non-zero dimension is represented by an index-value pair.
pub struct SparseVector {
    pub indices: Vec<u32>,
    pub values: Vec<f32>,
}

impl SparseVector {
    /// Create a new sparse vector from indices and values.
    ///
    /// # Panics
    /// Panics if `indices` and `values` have different lengths.
    pub fn new(indices: Vec<u32>, values: Vec<f32>) -> Self {
        assert_eq!(
            indices.len(),
            values.len(),
            "indices and values must have the same length"
        );
        SparseVector { indices, values }
    }

    /// Create an empty sparse vector.
    pub fn empty() -> Self {
        SparseVector {
            indices: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Returns the number of non-zero entries.
    pub fn len(&self) -> usize {
        self.indices.len()
    }

    /// Returns true if the vector has no entries.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
}