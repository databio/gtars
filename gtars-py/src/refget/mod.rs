// This mode provides Python bindings to the `refget rust module of gtars.
// It will allow computing ga4gh digests, creating sequence store objects,
// and sequence collection objects from Python.

use pyo3::exceptions::{PyIndexError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString, PyType};

use gtars_refget::alphabet::AlphabetType;
use gtars_refget::collection::{
    SeqColDigestLvl1, SequenceCollection, SequenceMetadata, SequenceRecord,
};
use gtars_refget::digest::{md5, sha512t24u};
use gtars_refget::store::GlobalRefgetStore;
use gtars_refget::store::StorageMode;
// use gtars::refget::store::RetrievedSequence; // This is the Rust-native struct

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
    match gtars_refget::fasta::digest_fasta(&fasta) {
        Ok(sequence_collection) => Ok(PySequenceCollection::from(sequence_collection)),
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
            "Error processing FASTA file: {}",
            e
        ))),
    }
}

#[pyclass(name = "AlphabetType")]
#[derive(Clone)]
pub enum PyAlphabetType {
    Dna2bit,
    Dna3bit,
    DnaIupac,
    Protein,
    Ascii,
    Unknown,
}

#[pyclass(name = "SequenceMetadata")]
#[derive(Clone)]
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

#[pyclass(name = "SequenceRecord")]
#[derive(Clone)]
pub struct PySequenceRecord {
    #[pyo3(get, set)]
    pub metadata: PySequenceMetadata,
    #[pyo3(get, set)]
    pub data: Option<Vec<u8>>,
}

#[pyclass(name = "SeqColDigestLvl1")]
#[derive(Clone)]
pub struct PySeqColDigestLvl1 {
    #[pyo3(get, set)]
    pub sequences_digest: String,
    #[pyo3(get, set)]
    pub names_digest: String,
    #[pyo3(get, set)]
    pub lengths_digest: String,
}

#[derive(Clone)]
#[pyclass(name = "SequenceCollection")]
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

#[pyclass(name = "RetrievedSequence")]
#[derive(Debug, Clone, PartialEq)]
pub struct PyRetrievedSequence {
    #[pyo3(get, set)]
    pub sequence: String,
    #[pyo3(get, set)]
    pub chrom_name: String,
    #[pyo3(get, set)]
    pub start: u32,
    #[pyo3(get, set)]
    pub end: u32,
}

// This `From` implementation converts the Rust-native `RetrievedSequence`
// into the Python-exposed `PyRetrievedSequence`.
impl From<gtars_refget::store::RetrievedSequence> for PyRetrievedSequence {
    fn from(value: gtars_refget::store::RetrievedSequence) -> Self {
        PyRetrievedSequence {
            sequence: value.sequence,
            chrom_name: value.chrom_name,
            start: value.start,
            end: value.end,
        }
    }
}

#[pymethods]
impl PyRetrievedSequence {
    #[new]
    fn new(sequence: String, chrom_name: String, start: u32, end: u32) -> Self {
        PyRetrievedSequence {
            sequence,
            chrom_name,
            start,
            end,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "RetrievedSequence(chrom_name='{}', start={}, end={}, sequence='{}')",
            self.chrom_name, self.start, self.end, self.sequence
        )
    }

    fn __str__(&self) -> String {
        format!(
            "{}|{}-{}: {}",
            self.chrom_name, self.start, self.end, self.sequence
        )
    }
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
        format!(
            "<SequenceCollection with {} sequences>",
            self.sequences.len()
        )
    }

    fn __str__(&self) -> String {
        format!(
            "SequenceCollection with {} sequences, digest: {}",
            self.sequences.len(),
            self.digest
        )
    }
    fn __len__(&self) -> PyResult<usize> {
        Ok(self.sequences.len())
    }
    fn __getitem__(&self, idx: isize, py: Python) -> PyResult<PyObject> {
        let len = self.sequences.len() as isize;

        // Handle negative indexing like Python lists do
        let index = if idx < 0 { len + idx } else { idx };

        if index >= 0 && (index as usize) < self.sequences.len() {
            // Convert the PySequenceRecord to a PyObject before returning
            Ok(self.sequences[index as usize].clone().into_py(py))
        } else {
            Err(PyIndexError::new_err(
                "SequenceCollection index out of range",
            ))
        }
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
            sequences: value
                .sequences
                .into_iter()
                .map(PySequenceRecord::from)
                .collect(),
            digest: value.digest,
            lvl1: PySeqColDigestLvl1::from(value.lvl1),
            file_path: value.file_path.map(|p| p.to_string_lossy().to_string()),
            has_data: value.has_data,
        }
    }
}

#[pyclass(name = "StorageMode")]
#[derive(Clone)]
pub enum PyStorageMode {
    Raw,
    Encoded,
}

impl From<StorageMode> for PyStorageMode {
    fn from(mode: StorageMode) -> Self {
        match mode {
            StorageMode::Raw => PyStorageMode::Raw,
            StorageMode::Encoded => PyStorageMode::Encoded,
        }
    }
}

impl From<PyStorageMode> for StorageMode {
    fn from(mode: PyStorageMode) -> Self {
        match mode {
            PyStorageMode::Raw => StorageMode::Raw,
            PyStorageMode::Encoded => StorageMode::Encoded,
        }
    }
}

#[pyclass(name = "GlobalRefgetStore")]
pub struct PyGlobalRefgetStore {
    inner: GlobalRefgetStore,
}

#[pymethods]
impl PyGlobalRefgetStore {
    #[new]
    fn new(mode: PyStorageMode) -> Self {
        Self {
            inner: GlobalRefgetStore::new(mode.into()),
        }
    }

    fn import_fasta(&mut self, file_path: &str) -> PyResult<()> {
        self.inner.import_fasta(file_path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error importing FASTA: {}", e))
        })
    }

    fn get_sequence_by_id(&self, digest: &str) -> PyResult<Option<PySequenceRecord>> {
        // Try as SHA512t24u first (32 bytes)
        let result = self
            .inner
            .get_sequence_by_id(digest.as_bytes())
            .map(|record| PySequenceRecord::from(record.clone()));

        // If not found and input looks like MD5 (32 hex chars), try MD5 lookup
        if result.is_none() && digest.len() == 32 {
            return Ok(self
                .inner
                .get_sequence_by_md5(digest.as_bytes())
                .map(|record| PySequenceRecord::from(record.clone())));
        }

        Ok(result)
    }

    fn get_sequence_by_collection_and_name(
        &self,
        collection_digest: &str,
        sequence_name: &str,
    ) -> Option<PySequenceRecord> {
        self.inner
            .get_sequence_by_collection_and_name(collection_digest, sequence_name)
            .map(|record| PySequenceRecord::from(record.clone()))
    }

    fn get_substring(&self, seq_digest: &str, start: usize, end: usize) -> Option<String> {
        self.inner.get_substring(seq_digest, start, end)
    }

    fn write_store_to_directory(
        &self,
        root_path: &str,
        seqdata_path_template: &str,
    ) -> PyResult<()> {
        self.inner
            .write_store_to_directory(root_path, seqdata_path_template)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error writing store: {}", e))
            })
    }

    #[classmethod]
    fn load_from_directory(_cls: &Bound<'_, PyType>, root_path: &str) -> PyResult<Self> {
        let store = GlobalRefgetStore::load_from_directory(root_path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error loading store: {}", e))
        })?;
        Ok(Self { inner: store })
    }

    fn get_seqs_bed_file(
        &self,
        collection_digest: &str,
        bed_file_path: &str,
        output_file_path: &str,
    ) -> PyResult<()> {
        // Rust function expects K: AsRef<[u8]>, &str works directly
        self.inner
            .get_seqs_bed_file(collection_digest, bed_file_path, output_file_path)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                    "Error retrieving sequences and writing to file: {}",
                    e
                ))
            })
    }

    fn get_seqs_bed_file_to_vec(
        &self,
        collection_digest: &str,
        bed_file_path: &str,
    ) -> PyResult<Vec<PyRetrievedSequence>> {
        // Corrected: use `let ... = ...?;` to bind the result
        let rust_results = self
            .inner
            .get_seqs_bed_file_to_vec(collection_digest, bed_file_path)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                    "Error retrieving sequences to list: {}",
                    e
                ))
            })?;

        // Now `rust_results` is available
        let py_results: Vec<PyRetrievedSequence> = rust_results
            .into_iter()
            .map(PyRetrievedSequence::from)
            .collect();

        Ok(py_results)
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __repr__(&self) -> String {
        format!(
            "<GlobalRefgetStore with {} sequences>",
            self.inner.list_sequence_digests().len()
        )
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
    m.add_class::<PyStorageMode>()?;
    m.add_class::<PyGlobalRefgetStore>()?;
    m.add_class::<PyRetrievedSequence>()?;
    Ok(())
}
