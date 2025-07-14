// This mode provides Python bindings to the `refget rust module of gtars.
// It will allow computing ga4gh digests, creating sequence store objects,
// and sequence collection objects from Python.

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};
use pyo3::exceptions::PyTypeError;


// use gtars::refget::{GlobalRefgetStore, SequenceCollection};
use gtars::refget::digest::{md5, sha512t24u};
use gtars::refget::collection::{SequenceCollection, SequenceMetadata, SequenceRecord, SeqColDigestLvl1};
use gtars::refget::alphabet::AlphabetType;


#[pyfunction]
pub fn sha512t24u_digest(readable: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(s) = readable.downcast::<PyString>() {
        Ok(sha512t24u(s.encode_utf8()?.as_bytes())) // Borrowed, no copying
    } else if let Ok(b) = readable.downcast::<PyBytes>() {
        Ok(sha512t24u(b.as_bytes())) // Borrowed, no copying
    } else {
        Err(PyTypeError::new_err("Expected str or bytes"))
    }
}

#[pyfunction]
pub fn md5_digest(readable: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(s) = readable.downcast::<PyString>() {
        Ok(md5(s.encode_utf8()?.as_bytes())) // Borrowed, no copying
    } else if let Ok(b) = readable.downcast::<PyBytes>() {
        Ok(md5(b.as_bytes())) // Borrowed, no copying
    } else {
        Err(PyTypeError::new_err("Expected str or bytes"))
    }
}

// This can take either a PosixPath or a string
// The `&Bound<'_, PyAny>` references any Python object, bound to the Python runtime.
#[pyfunction]
pub fn digest_fasta(fasta: &Bound<'_, PyAny>) -> PyResult<PySequenceCollection> {
    let fasta = fasta.to_string();
    match gtars::refget::fasta::digest_fasta(&fasta) {
        Ok(sequence_collection) => {
            Ok(PySequenceCollection::from(sequence_collection))
        }
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
            "Error processing FASTA file: {}",
            e
        ))),
    }
}

#[pyclass]
#[derive(Clone)]
#[pyo3(name = "AlphabetType")]
pub enum PyAlphabetType {
    Dna2bit,
    Dna3bit,
    DnaIupac,
    Protein,
    Ascii,
    Unknown,
}

#[pyclass]
#[derive(Clone)]
#[pyo3(name = "SequenceMetadata")]
pub struct PySequenceMetadata {
    #[pyo3(get, set)]
    pub name: String,
    #[pyo3(get, set)]
    pub length: usize,
    #[pyo3(get, set)]
    pub sha512t24u: String,
    #[pyo3(get, set)]
    pub md5: String,
    #[pyo3(get, set)]
    pub alphabet: PyAlphabetType,
}

#[pyclass]
#[derive(Clone)]
#[pyo3(name = "SequenceRecord")]
pub struct PySequenceRecord {
    #[pyo3(get, set)]
    pub metadata: PySequenceMetadata,
    #[pyo3(get, set)]
    pub data: Option<Vec<u8>>,
}

#[pyclass]
#[derive(Clone)]
#[pyo3(name = "SeqColDigestLvl1")]
pub struct PySeqColDigestLvl1 {
    #[pyo3(get, set)]
    pub sequences_digest: String,
    #[pyo3(get, set)]
    pub names_digest: String,
    #[pyo3(get, set)]
    pub lengths_digest: String,
}

#[pyclass]
#[derive(Clone)]
#[pyo3(name = "SequenceCollection")]
pub struct PySequenceCollection {
    #[pyo3(get, set)]
    pub sequences: Vec<PySequenceRecord>,
    #[pyo3(get, set)]
    pub digest: String,
    #[pyo3(get, set)]
    pub lvl1: PySeqColDigestLvl1,
    #[pyo3(get, set)]
    pub file_path: Option<String>,
    #[pyo3(get, set)]
    pub has_data: bool,
}

#[pymethods]
impl PyAlphabetType {
    fn __str__(&self) -> String {
        match self {
            PyAlphabetType::Dna2bit => "dna2bit".to_string(),
            PyAlphabetType::Dna3bit => "dna3bit".to_string(),
            PyAlphabetType::DnaIupac => "dnaio".to_string(),
            PyAlphabetType::Protein => "protein".to_string(),
            PyAlphabetType::Ascii => "ASCII".to_string(),
            PyAlphabetType::Unknown => "Unknown".to_string(),
        }
    }
}

#[pymethods]
impl PySequenceMetadata {
    fn __repr__(&self) -> String {
        format!("<SequenceMetadata for {}>", self.name)
    }

    fn __str__(&self) -> PyResult<String> {
        Ok(format!(
            "SequenceMetadata for sequence {}\n  length: {}\n  sha512t24u: {}\n  md5: {}\n  alphabet: {}",
            self.name, self.length, self.sha512t24u, self.md5, self.alphabet.__str__()
        ))
    }
}

#[pymethods]
impl PySequenceRecord {
    fn __repr__(&self) -> String {
        format!("<SequenceRecord for {}>", self.metadata.name)
    }

    fn __str__(&self) -> String {
        format!("SequenceRecord for {}", self.metadata.name)
    }
}

#[pymethods]
impl PySeqColDigestLvl1 {
    fn __repr__(&self) -> String {
        "<SeqColDigestLvl1>".to_string()
    }

    fn __str__(&self) -> String {
        format!(
            "SeqColDigestLvl1:\n  sequences_digest: {}\n  names_digest: {}\n  lengths_digest: {}",
            self.sequences_digest, self.names_digest, self.lengths_digest
        )
    }
}

#[pymethods]
impl PySequenceCollection {
    fn __repr__(&self) -> String {
        format!("<SequenceCollection with {} sequences>", self.sequences.len())
    }

    fn __str__(&self) -> String {
        format!(
            "SequenceCollection with {} sequences, digest: {}",
            self.sequences.len(), self.digest
        )
    }
}

// Conversion from Rust AlphabetType to Python PyAlphabetType
impl From<AlphabetType> for PyAlphabetType {
    fn from(value: AlphabetType) -> Self {
        match value {
            AlphabetType::Dna2bit => PyAlphabetType::Dna2bit,
            AlphabetType::Dna3bit => PyAlphabetType::Dna3bit,
            AlphabetType::DnaIupac => PyAlphabetType::DnaIupac,
            AlphabetType::Protein => PyAlphabetType::Protein,
            AlphabetType::Ascii => PyAlphabetType::Ascii,
            AlphabetType::Unknown => PyAlphabetType::Unknown,
        }
    }
}

// Conversion from Rust SequenceMetadata to Python PySequenceMetadata
impl From<SequenceMetadata> for PySequenceMetadata {
    fn from(value: SequenceMetadata) -> Self {
        PySequenceMetadata {
            name: value.name,
            length: value.length,
            sha512t24u: value.sha512t24u,
            md5: value.md5,
            alphabet: PyAlphabetType::from(value.alphabet),
        }
    }
}

// Conversion from Rust SequenceRecord to Python PySequenceRecord
impl From<SequenceRecord> for PySequenceRecord {
    fn from(value: SequenceRecord) -> Self {
        PySequenceRecord {
            metadata: PySequenceMetadata::from(value.metadata),
            data: value.data,
        }
    }
}

// Conversion from Rust SeqColDigestLvl1 to Python PySeqColDigestLvl1
impl From<SeqColDigestLvl1> for PySeqColDigestLvl1 {
    fn from(value: SeqColDigestLvl1) -> Self {
        PySeqColDigestLvl1 {
            sequences_digest: value.sequences_digest,
            names_digest: value.names_digest,
            lengths_digest: value.lengths_digest,
        }
    }
}

// Conversion from Rust SequenceCollection to Python PySequenceCollection
impl From<SequenceCollection> for PySequenceCollection {
    fn from(value: SequenceCollection) -> Self {
        PySequenceCollection {
            sequences: value.sequences.into_iter().map(PySequenceRecord::from).collect(),
            digest: value.digest,
            lvl1: PySeqColDigestLvl1::from(value.lvl1),
            file_path: value.file_path.map(|p| p.to_string_lossy().to_string()),
            has_data: value.has_data,
        }
    }
}

// This represents the Python module to be created
#[pymodule]
pub fn refget(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sha512t24u_digest, m)?)?;
    m.add_function(wrap_pyfunction!(md5_digest, m)?)?;
    m.add_function(wrap_pyfunction!(digest_fasta, m)?)?;
    m.add_class::<PyAlphabetType>()?;
    m.add_class::<PySequenceMetadata>()?;
    m.add_class::<PySequenceRecord>()?;
    m.add_class::<PySeqColDigestLvl1>()?;
    m.add_class::<PySequenceCollection>()?;
    Ok(())
}
