// This is intended to provide minimal Python bindings to functions in the `digests` module of the `gtars` crate.

use pyo3::prelude::*;
use gtars::digests::{sha512t24u, md5, DigestResult};

#[pyfunction]
pub fn sha512t24u_digest(readable: &str) -> String {
    return sha512t24u(readable);
}

#[pyfunction]
pub fn md5_digest(readable: &str) -> String {
    return md5(readable);
}

#[pyfunction]
pub fn digest_fasta(fasta: &str) -> PyResult<Vec<PyDigestResult>> {
    match gtars::digests::digest_fasta(fasta) {
        Ok(digest_results) => {
            let py_digest_results: Vec<PyDigestResult> = digest_results.into_iter().map(PyDigestResult::from).collect();
            Ok(py_digest_results)
        },
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error processing FASTA file: {}", e))),
    }
}

#[pyclass]
#[pyo3(name="DigestResult")]
pub struct PyDigestResult {
    #[pyo3(get,set)]
    pub id: String,
    #[pyo3(get,set)]
    pub length: usize,
    #[pyo3(get,set)]
    pub sha512t24u: String,
    #[pyo3(get,set)]
    pub md5: String
}

#[pymethods]
impl PyDigestResult {
    fn __repr__(&self) -> String {
        format!("<DigestResult for {}>", self.id)
    }

    fn __str__(&self) -> PyResult<String> {
        Ok(format!("DigestResult for sequence {}\n  length: {}\n  sha512t24u: {}\n  md5: {}", self.id, self.length, self.sha512t24u, self.md5))
    }
}

impl From<DigestResult> for PyDigestResult {
    fn from(value: DigestResult) -> Self {
        PyDigestResult {
            id: value.id,
            length: value.length,
            sha512t24u: value.sha512t24u,
            md5: value.md5
        }
    }
}

// This represents the Python module to be created
#[pymodule]
pub fn digests(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sha512t24u_digest, m)?)?;
    m.add_function(wrap_pyfunction!(md5_digest, m)?)?;
    m.add_function(wrap_pyfunction!(digest_fasta, m)?)?;
    m.add_class::<PyDigestResult>()?;
    Ok(())
}

