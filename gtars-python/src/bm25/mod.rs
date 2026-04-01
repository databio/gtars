use pyo3::prelude::*;

use anyhow::Result;
use gtars_bm25::{Bm25Builder, SparseVector};
use gtars_tokenizers::{Tokenizer, create_tokenize_core_from_universe, config::TokenizerType};

use crate::tokenizers::py_tokenizers::PyTokenizer;
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
    inner: gtars_bm25::Bm25,
}

#[pymethods]
impl PyBm25 {
    #[new]
    #[pyo3(signature = (tokenizer, k=1.0, b=0.75, avg_doc_length=1000.0))]
    fn new(tokenizer: &Bound<'_, PyAny>, k: f32, b: f32, avg_doc_length: f32) -> Result<Self> {
        let mut builder = Bm25Builder::default()
            .with_k(k)
            .with_b(b)
            .with_avg_doc_length(avg_doc_length);

        if let Ok(path) = tokenizer.extract::<String>() {
            // Accept a string path to a BED/BED.GZ/TOML file
            builder = builder.with_vocab(&path);
        } else if let Ok(py_tok) = tokenizer.cast::<PyTokenizer>() {
            // Accept an existing Tokenizer object — rebuild from its universe
            let borrowed = py_tok.borrow();
            let inner_tok = borrowed.inner();
            let universe = inner_tok.get_universe().clone();
            let special_tokens = inner_tok.get_special_tokens().clone();
            let core = create_tokenize_core_from_universe(&universe, TokenizerType::AIList);
            let tok = Tokenizer::new(core, universe, special_tokens);
            builder = builder.with_tokenizer(tok);
        } else {
            return Err(anyhow::anyhow!(
                "tokenizer must be a string path (BED, BED.GZ, TOML) or a Tokenizer object"
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