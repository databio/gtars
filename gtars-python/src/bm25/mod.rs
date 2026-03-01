use pyo3::prelude::*;

use anyhow::Result;
use gtars_bm25::{BM25Builder, SparseVector};

use crate::utils::extract_regions_from_py_any;

#[pyclass(name = "SparseVector", module = "gtars.bm25")]
pub struct PySparseVector {
    inner: SparseVector,
}

#[pymethods]
impl PySparseVector {
    #[getter]
    fn indices(&self) -> Vec<u32> {
        self.inner.indices.clone()
    }

    #[getter]
    fn values(&self) -> Vec<f32> {
        self.inner.values.clone()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "SparseVector(len={}, indices={:?}, values={:?})",
            self.inner.len(),
            self.inner.indices,
            self.inner.values
        )
    }
}

#[pyclass(name = "Bm25", module = "gtars.bm25", subclass)]
pub struct PyBm25 {
    inner: gtars_bm25::BM25,
}

#[pymethods]
impl PyBm25 {
    #[new]
    #[pyo3(signature = (tokenizer, k=1.0, b=0.75, avg_doc_length=1000.0))]
    fn new(tokenizer: &Bound<'_, PyAny>, k: f32, b: f32, avg_doc_length: f32) -> Result<Self> {
        let mut builder = BM25Builder::default()
            .with_k(k)
            .with_b(b)
            .with_avg_doc_length(avg_doc_length);

        // Accept either a string path or a Tokenizer object
        if let Ok(path) = tokenizer.extract::<String>() {
            builder = builder.with_vocab(&path);
        } else {
            // Try to use an existing Tokenizer's vocab path
            // The Tokenizer is opaque from Python, so we need to accept a path
            return Err(anyhow::anyhow!(
                "tokenizer must be a string path to a BED file, BED.GZ file, or TOML config"
            ));
        }

        Ok(PyBm25 {
            inner: builder.build(),
        })
    }

    /// Embed a set of regions into a BM25 sparse vector.
    fn embed(&self, regions: &Bound<'_, PyAny>) -> Result<PySparseVector> {
        let rs = extract_regions_from_py_any(regions)?;
        let sv = self.inner.embed(&rs.regions);
        Ok(PySparseVector { inner: sv })
    }

    /// Tokenize regions into token IDs without computing BM25 scores.
    fn tokenize(&self, regions: &Bound<'_, PyAny>) -> Result<Vec<u32>> {
        let rs = extract_regions_from_py_any(regions)?;
        Ok(self.inner.tokenize(&rs.regions))
    }

    #[getter]
    fn vocab_size(&self) -> usize {
        self.inner.vocab_size()
    }

    #[getter]
    fn k(&self) -> f32 {
        self.inner.k()
    }

    #[getter]
    fn b(&self) -> f32 {
        self.inner.b()
    }

    #[getter]
    fn avg_doc_length(&self) -> f32 {
        self.inner.avg_doc_length()
    }

    fn __repr__(&self) -> String {
        format!(
            "Bm25(vocab_size={}, k={}, b={}, avg_doc_length={})",
            self.inner.vocab_size(),
            self.inner.k(),
            self.inner.b(),
            self.inner.avg_doc_length()
        )
    }
}

#[pymodule]
pub fn bm25(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBm25>()?;
    m.add_class::<PySparseVector>()?;
    Ok(())
}