// This mode provides Python bindings to the `refget rust module of gtars.
// It will allow computing ga4gh digests, creating sequence store objects,
// and sequence collection objects from Python.

use pyo3::exceptions::{PyIndexError, PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString, PyType};

use gtars_refget::alphabet::AlphabetType;
use gtars_refget::collection::{
    FaiMetadata, SeqColDigestLvl1, SequenceCollection, SequenceCollectionMetadata, SequenceMetadata, SequenceRecord,
};
use gtars_refget::digest::{md5, sha512t24u};
use gtars_refget::fasta::FaiRecord;
use gtars_refget::store::RefgetStore;
use gtars_refget::store::StorageMode;
// use gtars::refget::store::RetrievedSequence; // This is the Rust-native struct

/// Compute the GA4GH SHA-512/24u digest for a sequence.
///
/// Accepts either a string or bytes and returns the truncated SHA-512 digest
/// as a 32-character base64url string.
#[pyfunction]
pub fn sha512t24u_digest(readable: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(s) = readable.cast::<PyString>() {
        Ok(sha512t24u(s.encode_utf8()?.as_bytes())) // Borrowed, no copying
    } else if let Ok(b) = readable.cast::<PyBytes>() {
        Ok(sha512t24u(b.as_bytes())) // Borrowed, no copying
    } else {
        Err(PyTypeError::new_err("Expected str or bytes"))
    }
}

/// Compute the MD5 digest for a sequence.
///
/// Accepts either a string or bytes and returns the MD5 hash as a 32-character
/// hexadecimal string. Supported for backward compatibility with legacy systems.
#[pyfunction]
pub fn md5_digest(readable: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(s) = readable.cast::<PyString>() {
        Ok(md5(s.encode_utf8()?.as_bytes())) // Borrowed, no copying
    } else if let Ok(b) = readable.cast::<PyBytes>() {
        Ok(md5(b.as_bytes())) // Borrowed, no copying
    } else {
        Err(PyTypeError::new_err("Expected str or bytes"))
    }
}

/// Read a FASTA file and compute GA4GH digests for all sequences.
///
/// Returns a SequenceCollection with computed sequence-level and collection-level
/// digests following the GA4GH refget specification.
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

/// Compute FASTA index (FAI) metadata for all sequences in a FASTA file.
///
/// Returns a list of FAI records compatible with samtools faidx format.
/// Only works with uncompressed FASTA files.
#[pyfunction]
pub fn compute_fai(fasta: &Bound<'_, PyAny>) -> PyResult<Vec<PyFaiRecord>> {
    let fasta = fasta.to_string();
    match gtars_refget::fasta::compute_fai(&fasta) {
        Ok(fai_records) => Ok(fai_records.into_iter().map(PyFaiRecord::from).collect()),
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
            "Error computing FAI: {}",
            e
        ))),
    }
}

/// Load a FASTA file with sequence data into a SequenceCollection.
///
/// Unlike digest_fasta(), this loads the actual sequence data into memory,
/// not just metadata and digests.
#[pyfunction]
pub fn load_fasta(fasta: &Bound<'_, PyAny>) -> PyResult<PySequenceCollection> {
    let fasta = fasta.to_string();
    match gtars_refget::fasta::load_fasta(&fasta) {
        Ok(sequence_collection) => Ok(PySequenceCollection::from(sequence_collection)),
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
            "Error loading FASTA file: {}",
            e
        ))),
    }
}

/// Create a SequenceRecord from raw data, computing all metadata.
///
/// This is the sequence-level parallel to `digest_fasta()` for collections.
/// It computes the GA4GH sha512t24u digest, MD5 digest, detects the alphabet,
/// and returns a SequenceRecord with the computed metadata and original data.
///
/// Args:
///     name: The sequence name (e.g., "chr1")
///     data: The raw sequence bytes (e.g., b"ACGTACGT")
///     description: Optional description text
///
/// Returns:
///     A SequenceRecord with computed metadata and the original data (uppercased)
///
/// Example:
///     >>> from gtars.refget import digest_sequence
///     >>> seq = digest_sequence("chr1", b"ACGTACGT")
///     >>> print(seq.name, seq.length)
///     chr1 8
#[pyfunction]
#[pyo3(signature = (name, data, description=None))]
pub fn digest_sequence(name: &str, data: &[u8], description: Option<&str>) -> PySequenceRecord {
    let seq_record = match description {
        Some(desc) => gtars_refget::collection::digest_sequence_with_description(name, Some(desc), data),
        None => gtars_refget::collection::digest_sequence(name, data),
    };
    PySequenceRecord::from(seq_record)
}

/// The type of alphabet for a biological sequence.
///
/// Used to determine encoding strategy and validation rules.
///
/// Variants:
///     Dna2bit: Standard DNA (A, C, G, T only) - can use 2-bit encoding.
///     Dna3bit: DNA with N (A, C, G, T, N) - requires 3-bit encoding.
///     DnaIupac: Full IUPAC DNA alphabet with ambiguity codes.
///     Protein: Amino acid sequences.
///     Ascii: Generic ASCII text.
///     Unknown: Alphabet could not be determined.
#[pyclass(name = "AlphabetType", module = "gtars.refget")]
#[derive(Clone)]
pub enum PyAlphabetType {
    Dna2bit,
    Dna3bit,
    DnaIupac,
    Protein,
    Ascii,
    Unknown,
}

/// Metadata for a biological sequence.
///
/// Contains identifying information and computed digests for a sequence,
/// without the actual sequence data.
///
/// Attributes:
///     name (str): Sequence name (first word of FASTA header).
///     description (str | None): Description from FASTA header (text after first whitespace).
///     length (int): Length of the sequence in bases.
///     sha512t24u (str): GA4GH SHA-512/24u digest (32-char base64url).
///     md5 (str): MD5 digest (32-char hex string).
///     alphabet (AlphabetType): Detected alphabet type (DNA, protein, etc.).
///     fai (FaiMetadata | None): FASTA index metadata if available.
#[pyclass(name = "SequenceMetadata", module = "gtars.refget")]
#[derive(Clone)]
pub struct PySequenceMetadata {
    #[pyo3(get, set)]
    pub name: String,
    #[pyo3(get, set)]
    pub description: Option<String>,
    #[pyo3(get, set)]
    pub length: usize,
    #[pyo3(get, set)]
    pub sha512t24u: String,
    #[pyo3(get, set)]
    pub md5: String,
    #[pyo3(get, set)]
    pub alphabet: PyAlphabetType,
    #[pyo3(get, set)]
    pub fai: Option<PyFaiMetadata>,
}

/// FASTA index (FAI) metadata for a sequence.
///
/// Contains the information needed to quickly seek to a sequence
/// in a FASTA file, compatible with samtools faidx format.
///
/// Attributes:
///     offset (int): Byte offset of the first base in the FASTA file.
///     line_bases (int): Number of bases per line.
///     line_bytes (int): Number of bytes per line (including newline).
#[pyclass(name = "FaiMetadata", module = "gtars.refget")]
#[derive(Clone)]
pub struct PyFaiMetadata {
    #[pyo3(get, set)]
    pub offset: u64,
    #[pyo3(get, set)]
    pub line_bases: u32,
    #[pyo3(get, set)]
    pub line_bytes: u32,
}

/// A FASTA index record for a single sequence.
///
/// Represents one line of a .fai index file with sequence name,
/// length, and FAI metadata for random access.
///
/// Attributes:
///     name (str): Sequence name.
///     length (int): Sequence length in bases.
///     fai (FaiMetadata | None): FAI metadata (None for gzipped files).
#[pyclass(name = "FaiRecord", module = "gtars.refget")]
#[derive(Clone)]
pub struct PyFaiRecord {
    #[pyo3(get, set)]
    pub name: String,
    #[pyo3(get, set)]
    pub length: usize,
    #[pyo3(get, set)]
    pub fai: Option<PyFaiMetadata>,
}

/// A record representing a biological sequence with metadata and optional data.
///
/// SequenceRecord can be either a "stub" (metadata only) or "full" (metadata + data).
/// Stubs are used for lazy-loading where sequence data is fetched on demand.
///
/// Attributes:
///     metadata (SequenceMetadata): Sequence metadata (name, length, digests).
///     sequence (bytes | None): Raw sequence data if loaded, None for stubs.
#[pyclass(name = "SequenceRecord", module = "gtars.refget")]
#[derive(Clone)]
pub struct PySequenceRecord {
    #[pyo3(get, set)]
    pub metadata: PySequenceMetadata,
    #[pyo3(get, set)]
    pub sequence: Option<Vec<u8>>,
}

/// Level 1 digests for a sequence collection.
///
/// These are intermediate digests computed over the arrays of sequence
/// properties, used in the GA4GH seqcol specification.
///
/// Attributes:
///     sequences_digest (str): Digest of the array of sequence digests.
///     names_digest (str): Digest of the array of sequence names.
///     lengths_digest (str): Digest of the array of sequence lengths.
#[pyclass(name = "SeqColDigestLvl1", module = "gtars.refget")]
#[derive(Clone)]
pub struct PySeqColDigestLvl1 {
    #[pyo3(get, set)]
    pub sequences_digest: String,
    #[pyo3(get, set)]
    pub names_digest: String,
    #[pyo3(get, set)]
    pub lengths_digest: String,
}

/// Metadata for a sequence collection.
///
/// Contains the collection digest and level 1 digests for names, sequences, and lengths.
/// This is a lightweight representation of a collection without the actual sequence list.
///
/// Attributes:
///     digest (str): The collection's SHA-512/24u digest.
///     n_sequences (int): Number of sequences in the collection.
///     names_digest (str): Level 1 digest of the names array.
///     sequences_digest (str): Level 1 digest of the sequences array.
///     lengths_digest (str): Level 1 digest of the lengths array.
#[pyclass(name = "SequenceCollectionMetadata", module = "gtars.refget")]
#[derive(Clone)]
pub struct PySequenceCollectionMetadata {
    #[pyo3(get, set)]
    pub digest: String,
    #[pyo3(get, set)]
    pub n_sequences: usize,
    #[pyo3(get, set)]
    pub names_digest: String,
    #[pyo3(get, set)]
    pub sequences_digest: String,
    #[pyo3(get, set)]
    pub lengths_digest: String,
}

#[pymethods]
impl PySequenceCollectionMetadata {
    fn __repr__(&self) -> String {
        format!(
            "SequenceCollectionMetadata(digest='{}', n_sequences={})",
            self.digest, self.n_sequences
        )
    }

    fn __str__(&self) -> String {
        format!(
            "SequenceCollectionMetadata:\n  digest: {}\n  n_sequences: {}\n  names_digest: {}\n  sequences_digest: {}\n  lengths_digest: {}",
            self.digest, self.n_sequences, self.names_digest, self.sequences_digest, self.lengths_digest
        )
    }
}

impl From<SequenceCollectionMetadata> for PySequenceCollectionMetadata {
    fn from(value: SequenceCollectionMetadata) -> Self {
        PySequenceCollectionMetadata {
            digest: value.digest,
            n_sequences: value.n_sequences,
            names_digest: value.names_digest,
            sequences_digest: value.sequences_digest,
            lengths_digest: value.lengths_digest,
        }
    }
}

/// A collection of biological sequences (e.g., a genome assembly).
///
/// SequenceCollection represents a set of sequences with collection-level
/// digests following the GA4GH seqcol specification. Supports iteration
/// and indexing.
///
/// Attributes:
///     sequences (list[SequenceRecord]): List of sequence records.
///     digest (str): Collection-level SHA-512/24u digest (Level 2).
///     lvl1 (SeqColDigestLvl1): Level 1 digests for names, lengths, sequences.
///     file_path (str | None): Source file path if loaded from FASTA.
///
/// Examples:
///     Iterate over sequences::
///
///         for seq in collection:
///             print(f"{seq.metadata.name}: {seq.metadata.length} bp")
///
///     Access by index::
///
///         first_seq = collection[0]
///         last_seq = collection[-1]
#[derive(Clone)]
#[pyclass(name = "SequenceCollection", module = "gtars.refget")]
pub struct PySequenceCollection {
    #[pyo3(get, set)]
    pub sequences: Vec<PySequenceRecord>,
    #[pyo3(get, set)]
    pub digest: String,
    #[pyo3(get, set)]
    pub lvl1: PySeqColDigestLvl1,
    #[pyo3(get, set)]
    pub file_path: Option<String>,
}

/// A retrieved sequence segment with its genomic coordinates.
///
/// Returned by methods that extract subsequences from specific regions,
/// such as substrings_from_regions().
///
/// Attributes:
///     sequence (str): The extracted sequence string.
///     chrom_name (str): Chromosome/sequence name (e.g., "chr1").
///     start (int): Start position (0-based, inclusive).
///     end (int): End position (0-based, exclusive).
#[pyclass(name = "RetrievedSequence", module = "gtars.refget")]
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
        format!(
            "SequenceMetadata(name='{}', length={}, sha512t24u='{}', md5='{}', alphabet={})",
            self.name, self.length, self.sha512t24u, self.md5, self.alphabet.__str__()
        )
    }

    fn __str__(&self) -> PyResult<String> {
        Ok(format!(
            "SequenceMetadata for sequence {}\n  length: {}\n  sha512t24u: {}\n  md5: {}\n  alphabet: {}",
            self.name, self.length, self.sha512t24u, self.md5, self.alphabet.__str__()
        ))
    }
}

#[pymethods]
impl PyFaiMetadata {
    fn __repr__(&self) -> String {
        format!("<FaiMetadata offset={} line_bases={} line_bytes={}>",
                self.offset, self.line_bases, self.line_bytes)
    }

    fn __str__(&self) -> String {
        format!(
            "FaiMetadata:\n  offset: {}\n  line_bases: {}\n  line_bytes: {}",
            self.offset, self.line_bases, self.line_bytes
        )
    }
}

#[pymethods]
impl PyFaiRecord {
    fn __repr__(&self) -> String {
        format!("<FaiRecord name='{}' length={}>", self.name, self.length)
    }

    fn __str__(&self) -> String {
        let fai_str = if let Some(ref fai) = self.fai {
            format!("\n  FAI offset: {}\n  FAI line_bases: {}\n  FAI line_bytes: {}",
                    fai.offset, fai.line_bases, fai.line_bytes)
        } else {
            "\n  FAI: None (gzipped file)".to_string()
        };
        format!("FaiRecord:\n  name: {}\n  length: {}{}",
                self.name, self.length, fai_str)
    }
}

#[pymethods]
impl PySequenceRecord {
    fn __repr__(&self) -> String {
        format!(
            "SequenceRecord(name='{}', length={}, sha512t24u='{}', is_loaded={})",
            self.metadata.name,
            self.metadata.length,
            self.metadata.sha512t24u,
            self.sequence.is_some()
        )
    }

    fn __str__(&self) -> String {
        format!("SequenceRecord for {}", self.metadata.name)
    }

    /// Whether sequence data is loaded (true) or just metadata (false).
    #[getter]
    fn is_loaded(&self) -> bool {
        self.sequence.is_some()
    }

    /// Decode the sequence data to a string.
    ///
    /// This method decodes the sequence data stored in this record. It handles
    /// both raw (uncompressed UTF-8) and encoded (bit-packed) data automatically
    /// based on the alphabet type.
    ///
    /// Returns:
    ///     Optional[str]: The decoded sequence string if data is loaded, None otherwise.
    ///
    /// Example:
    ///     >>> record = store.get_sequence_by_collection_and_name(digest, "chr1")
    ///     >>> sequence = record.decode()
    ///     >>> if sequence:
    ///     ...     print(f"First 50 bases: {sequence[:50]}")
    pub fn decode(&self) -> Option<String> {
        // Convert PySequenceRecord to Rust SequenceRecord
        let metadata = SequenceMetadata {
            name: self.metadata.name.clone(),
            description: self.metadata.description.clone(),
            length: self.metadata.length,
            sha512t24u: self.metadata.sha512t24u.clone(),
            md5: self.metadata.md5.clone(),
            alphabet: match self.metadata.alphabet {
                PyAlphabetType::Dna2bit => AlphabetType::Dna2bit,
                PyAlphabetType::Dna3bit => AlphabetType::Dna3bit,
                PyAlphabetType::DnaIupac => AlphabetType::DnaIupac,
                PyAlphabetType::Protein => AlphabetType::Protein,
                PyAlphabetType::Ascii => AlphabetType::Ascii,
                PyAlphabetType::Unknown => AlphabetType::Unknown,
            },
            fai: self.metadata.fai.as_ref().map(|fai| FaiMetadata {
                offset: fai.offset,
                line_bases: fai.line_bases,
                line_bytes: fai.line_bytes,
            }),
        };

        let rust_record = match &self.sequence {
            None => SequenceRecord::Stub(metadata),
            Some(seq) => SequenceRecord::Full {
                metadata,
                sequence: seq.clone(),
            },
        };

        // Call the Rust decode method
        rust_record.decode()
    }
}

#[pymethods]
impl PySeqColDigestLvl1 {
    fn __repr__(&self) -> String {
        format!(
            "SeqColDigestLvl1(sequences='{}', names='{}', lengths='{}')",
            self.sequences_digest, self.names_digest, self.lengths_digest
        )
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
            "SequenceCollection(n_sequences={}, digest='{}')",
            self.sequences.len(),
            self.digest
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
    fn __getitem__(&self, idx: isize) -> PyResult<PySequenceRecord> {
        let len = self.sequences.len() as isize;

        // Handle negative indexing like Python lists do
        let index = if idx < 0 { len + idx } else { idx };

        if index >= 0 && (index as usize) < self.sequences.len() {
            // Convert the PySequenceRecord to a PyObject before returning
            let record = self.sequences[index as usize].clone();
            Ok(record)
        } else {
            Err(PyIndexError::new_err(
                "SequenceCollection index out of range",
            ))
        }
    }

    /// Write the collection to a FASTA file.
    ///
    /// Args:
    ///     file_path (str): Path to the output FASTA file
    ///     line_width (int, optional): Number of bases per line (default: 70)
    ///
    /// Raises:
    ///     IOError: If any sequence doesn't have data loaded
    ///
    /// Example:
    ///     >>> collection = load_fasta("genome.fa")
    ///     >>> collection.write_fasta("output.fa")
    ///     >>> collection.write_fasta("output.fa", line_width=60)
    fn write_fasta(&self, file_path: &str, line_width: Option<usize>) -> PyResult<()> {
        // Convert Python sequences back to Rust SequenceCollection
        let rust_collection = SequenceCollection {
            sequences: self.sequences.iter().map(|py_rec| {
                let metadata = SequenceMetadata {
                    name: py_rec.metadata.name.clone(),
                    description: py_rec.metadata.description.clone(),
                    length: py_rec.metadata.length,
                    sha512t24u: py_rec.metadata.sha512t24u.clone(),
                    md5: py_rec.metadata.md5.clone(),
                    alphabet: match &py_rec.metadata.alphabet {
                        PyAlphabetType::Dna2bit => AlphabetType::Dna2bit,
                        PyAlphabetType::Dna3bit => AlphabetType::Dna3bit,
                        PyAlphabetType::DnaIupac => AlphabetType::DnaIupac,
                        PyAlphabetType::Protein => AlphabetType::Protein,
                        PyAlphabetType::Ascii => AlphabetType::Ascii,
                        PyAlphabetType::Unknown => AlphabetType::Unknown,
                    },
                    fai: py_rec.metadata.fai.as_ref().map(|f| FaiMetadata {
                        offset: f.offset,
                        line_bases: f.line_bases,
                        line_bytes: f.line_bytes,
                    }),
                };

                match &py_rec.sequence {
                    None => SequenceRecord::Stub(metadata),
                    Some(seq) => SequenceRecord::Full {
                        metadata,
                        sequence: seq.clone(),
                    },
                }
            }).collect(),
            metadata: SequenceCollectionMetadata {
                digest: self.digest.clone(),
                n_sequences: self.sequences.len(),
                sequences_digest: self.lvl1.sequences_digest.clone(),
                names_digest: self.lvl1.names_digest.clone(),
                lengths_digest: self.lvl1.lengths_digest.clone(),
                file_path: None,
            },
        };

        rust_collection
            .write_fasta(file_path, line_width)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to write FASTA: {}",
                e
            )))
    }

    /// Iterate over sequences in the collection.
    ///
    /// Allows Pythonic iteration: `for seq in collection:`
    ///
    /// Yields:
    ///     SequenceRecord: Each sequence record in the collection.
    ///
    /// Example:
    ///     >>> collection = digest_fasta("genome.fa")
    ///     >>> for seq in collection:
    ///     ...     print(f"{seq.metadata.name}: {seq.metadata.length} bp")
    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<PySequenceCollectionIterator>> {
        let py = slf.py();
        Py::new(py, PySequenceCollectionIterator {
            iter: slf.sequences.clone().into_iter(),
        })
    }
}

/// Iterator for SequenceCollection - simple wrapper around Vec iterator.
#[pyclass]
pub struct PySequenceCollectionIterator {
    iter: std::vec::IntoIter<PySequenceRecord>,
}

#[pymethods]
impl PySequenceCollectionIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PySequenceRecord> {
        slf.iter.next()
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

// Conversion from Rust FaiMetadata to Python PyFaiMetadata
impl From<FaiMetadata> for PyFaiMetadata {
    fn from(value: FaiMetadata) -> Self {
        PyFaiMetadata {
            offset: value.offset,
            line_bases: value.line_bases,
            line_bytes: value.line_bytes,
        }
    }
}

// Conversion from Rust FaiRecord to Python PyFaiRecord
impl From<FaiRecord> for PyFaiRecord {
    fn from(value: FaiRecord) -> Self {
        PyFaiRecord {
            name: value.name,
            length: value.length,
            fai: value.fai.map(PyFaiMetadata::from),
        }
    }
}

// Conversion from Rust SequenceMetadata to Python PySequenceMetadata
impl From<SequenceMetadata> for PySequenceMetadata {
    fn from(value: SequenceMetadata) -> Self {
        PySequenceMetadata {
            name: value.name,
            description: value.description,
            length: value.length,
            sha512t24u: value.sha512t24u,
            md5: value.md5,
            alphabet: PyAlphabetType::from(value.alphabet),
            fai: value.fai.map(PyFaiMetadata::from),
        }
    }
}

// Conversion from Rust SequenceRecord to Python PySequenceRecord
impl From<SequenceRecord> for PySequenceRecord {
    fn from(value: SequenceRecord) -> Self {
        match value {
            SequenceRecord::Stub(metadata) => PySequenceRecord {
                metadata: PySequenceMetadata::from(metadata),
                sequence: None,
            },
            SequenceRecord::Full { metadata, sequence } => PySequenceRecord {
                metadata: PySequenceMetadata::from(metadata),
                sequence: Some(sequence),
            },
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
        let metadata = value.metadata;
        PySequenceCollection {
            sequences: value
                .sequences
                .into_iter()
                .map(PySequenceRecord::from)
                .collect(),
            digest: metadata.digest,
            lvl1: PySeqColDigestLvl1 {
                names_digest: metadata.names_digest,
                sequences_digest: metadata.sequences_digest,
                lengths_digest: metadata.lengths_digest,
            },
            file_path: metadata.file_path.map(|p| p.to_string_lossy().to_string()),
        }
    }
}

/// Defines how sequence data is stored in the RefgetStore.
///
/// StorageMode.Encoded uses 2-bit packing for DNA sequences (A=00, C=01, G=10, T=11),
/// reducing memory usage by 75% compared to raw storage. StorageMode.Raw stores
/// sequences as plain bytes.
///
/// Variants:
///     Raw: Store sequences as raw bytes (1 byte per base).
///     Encoded: Store sequences with 2-bit encoding (4 bases per byte).
#[pyclass(name = "StorageMode", module = "gtars.refget")]
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

/// A global store for GA4GH refget sequences with lazy-loading support.
///
/// RefgetStore provides content-addressable storage for reference genome
/// sequences following the GA4GH refget specification. Supports both local and
/// remote stores with on-demand sequence loading.
///
/// Attributes:
///     cache_path (str | None): Local directory path where the store is located or cached.
///         None for in-memory stores.
///     remote_url (str | None): Remote URL of the store if loaded remotely, None otherwise.
///
/// Examples:
///     Create a new in-memory store and import sequences::
///
///         from gtars.refget import RefgetStore
///         store = RefgetStore.in_memory()
///         store.add_sequence_collection_from_fasta("genome.fa")
///
///     Load an existing local store::
///
///         store = RefgetStore.load_local("/data/hg38")
///         seq = store.get_substring("chr1_digest", 0, 1000)
///
///     Load a remote store with caching::
///
///         store = RefgetStore.load_remote(
///             "/local/cache",
///             "https://example.com/hg38"
///         )
#[pyclass(name = "RefgetStore", module = "gtars.refget")]
pub struct PyRefgetStore {
    inner: RefgetStore,
}

#[pymethods]
impl PyRefgetStore {
    /// Create a disk-backed RefgetStore.
    ///
    /// Sequences are written to disk immediately and loaded on-demand (lazy loading).
    /// Only metadata is kept in memory.
    ///
    /// Args:
    ///     cache_path (str or Path): Directory for storing sequences and metadata
    ///     mode: Storage mode (StorageMode.Raw or StorageMode.Encoded)
    ///
    /// Returns:
    ///     RefgetStore: A configured disk-backed store
    ///
    /// Example:
    ///     >>> from gtars.refget import RefgetStore
    ///     >>> store = RefgetStore.on_disk("/data/store")
    ///     >>> store.add_sequence_collection_from_fasta("genome.fa")
    #[classmethod]
    fn on_disk(_cls: &Bound<'_, PyType>, cache_path: &Bound<'_, PyAny>) -> PyResult<Self> {
        let cache_path = cache_path.to_string();
        let store = RefgetStore::on_disk(cache_path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error with disk-backed store: {}", e))
        })?;
        Ok(Self { inner: store })
    }

    /// Create an in-memory RefgetStore.
    ///
    /// All sequences kept in RAM for fast access.
    /// Defaults to Encoded storage mode (2-bit packing for space efficiency).
    /// Use set_encoding_mode() to change storage mode after creation.
    ///
    /// Returns:
    ///     RefgetStore: A new in-memory store
    ///
    /// Example:
    ///     >>> from gtars.refget import RefgetStore
    ///     >>> store = RefgetStore.in_memory()
    ///     >>> store.add_sequence_collection_from_fasta("genome.fa")
    #[classmethod]
    fn in_memory(_cls: &Bound<'_, PyType>) -> Self {
        Self {
            inner: RefgetStore::in_memory(),
        }
    }

    /// Change the storage mode, re-encoding/decoding existing sequences as needed.
    ///
    /// When switching from Raw to Encoded:
    /// - All Full sequences in memory are encoded (2-bit packed)
    ///
    /// When switching from Encoded to Raw:
    /// - All Full sequences in memory are decoded back to raw bytes
    ///
    /// Args:
    ///     mode: The storage mode to switch to (StorageMode.Raw or StorageMode.Encoded)
    ///
    /// Example:
    ///     >>> from gtars.refget import RefgetStore, StorageMode
    ///     >>> store = RefgetStore.in_memory()
    ///     >>> store.set_encoding_mode(StorageMode.Raw)
    fn set_encoding_mode(&mut self, mode: PyStorageMode) {
        self.inner.set_encoding_mode(mode.into());
    }

    /// Enable 2-bit encoding for space efficiency.
    /// Re-encodes any existing Raw sequences in memory.
    ///
    /// Example:
    ///     >>> store = RefgetStore.in_memory()
    ///     >>> store.disable_encoding()  # Switch to Raw
    ///     >>> store.enable_encoding()   # Back to Encoded
    fn enable_encoding(&mut self) {
        self.inner.enable_encoding();
    }

    /// Disable encoding, use raw byte storage.
    /// Decodes any existing Encoded sequences in memory.
    ///
    /// Example:
    ///     >>> store = RefgetStore.in_memory()
    ///     >>> store.disable_encoding()  # Switch to Raw mode
    fn disable_encoding(&mut self) {
        self.inner.disable_encoding();
    }

    /// Set whether to suppress progress output.
    ///
    /// When quiet is True, operations like add_sequence_collection_from_fasta
    /// will not print progress messages.
    ///
    /// Args:
    ///     quiet (bool): Whether to suppress progress output.
    ///
    /// Example:
    ///     >>> store = RefgetStore.in_memory()
    ///     >>> store.set_quiet(True)  # Suppress output
    ///     >>> store.add_sequence_collection_from_fasta("genome.fa")  # No output
    fn set_quiet(&mut self, quiet: bool) {
        self.inner.set_quiet(quiet);
    }

    /// Returns whether the store is in quiet mode.
    ///
    /// Returns:
    ///     bool: True if quiet mode is enabled.
    #[getter]
    fn quiet(&self) -> bool {
        self.inner.is_quiet()
    }

    /// Enable disk persistence for this store.
    ///
    /// Sets up the store to write sequences to disk. Any in-memory Full sequences
    /// are flushed to disk and converted to Stubs.
    ///
    /// Args:
    ///     path (str or Path): Directory for storing sequences and metadata
    ///
    /// Raises:
    ///     IOError: If the directory cannot be created or written to.
    ///
    /// Example:
    ///     >>> store = RefgetStore.in_memory()
    ///     >>> store.add_sequence_collection_from_fasta("genome.fa")
    ///     >>> store.enable_persistence("/data/store")  # Flush to disk
    #[pyo3(signature = (path))]
    fn enable_persistence(&mut self, path: &Bound<'_, PyAny>) -> PyResult<()> {
        let path = path.to_string();
        self.inner.enable_persistence(path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error enabling persistence: {}", e))
        })
    }

    /// Disable disk persistence for this store.
    ///
    /// New sequences will be kept in memory only. Existing Stub sequences
    /// can still be loaded from disk if local_path is set.
    ///
    /// Example:
    ///     >>> store = RefgetStore.load_remote("/cache", "https://example.com")
    ///     >>> store.disable_persistence()  # Stop caching new sequences
    fn disable_persistence(&mut self) {
        self.inner.disable_persistence();
    }

    /// Add a sequence collection from a FASTA file.
    ///
    /// Reads a FASTA file, digests the sequences, creates a SequenceCollection,
    /// and adds it to the store along with all its sequences.
    ///
    /// Args:
    ///     file_path (str or Path): Path to the FASTA file to import.
    ///     force (bool, optional): If True, overwrite existing collections/sequences.
    ///                            If False (default), skip duplicates.
    ///
    /// Returns:
    ///     tuple[SequenceCollectionMetadata, bool]: A tuple containing:
    ///         - SequenceCollectionMetadata: Metadata for the collection (digest, n_sequences, level 1 digests)
    ///         - bool: True if the collection was newly added, False if it already existed
    ///
    /// Raises:
    ///     IOError: If the file cannot be read or processed.
    ///
    /// Example:
    ///     >>> store = RefgetStore.in_memory()
    ///     >>> metadata, was_new = store.add_sequence_collection_from_fasta("genome.fa")
    ///     >>> print(f"{'Added' if was_new else 'Skipped'}: {metadata.digest} ({metadata.n_sequences} seqs)")
    #[pyo3(signature = (file_path, force=false))]
    fn add_sequence_collection_from_fasta(&mut self, file_path: &Bound<'_, PyAny>, force: bool) -> PyResult<(PySequenceCollectionMetadata, bool)> {
        let file_path = file_path.to_string();
        let result = if force {
            self.inner.add_sequence_collection_from_fasta_force(file_path)
        } else {
            self.inner.add_sequence_collection_from_fasta(file_path)
        };
        result
            .map(|(metadata, was_new)| (PySequenceCollectionMetadata::from(metadata), was_new))
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error importing FASTA: {}", e))
            })
    }

    /// Retrieve a sequence record by its digest (SHA-512/24u or MD5).
    ///
    /// Searches for a sequence by its GA4GH SHA-512/24u digest. If not found
    /// and the input looks like an MD5 digest (32 hex characters), tries MD5 lookup.
    ///
    /// Args:
    ///     digest: Sequence digest (SHA-512/24u).
    ///
    /// Returns:
    ///     Optional[SequenceRecord]: The sequence record with data if found, None otherwise.
    ///
    /// Example:
    ///     >>> record = store.get_sequence("aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2")
    ///     >>> if record:
    ///     ...     print(f"Found: {record.metadata.name}")
    fn get_sequence(&mut self, digest: &str) -> Option<PySequenceRecord> {
        self.inner
            .get_sequence(digest.as_bytes())
            .ok()
            .map(|record| PySequenceRecord::from(record.clone()))
    }

    /// Retrieve a sequence by collection digest and sequence name.
    ///
    /// Looks up a sequence within a specific collection using its name
    /// (e.g., "chr1", "chrM"). Loads sequence data if needed.
    ///
    /// Args:
    ///     collection_digest: The collection's SHA-512/24u digest.
    ///     sequence_name: Name of the sequence within that collection.
    ///
    /// Returns:
    ///     Optional[SequenceRecord]: The sequence record with data if found, None otherwise.
    ///
    /// Example:
    ///     >>> record = store.get_sequence_by_name(
    ///     ...     "uC_UorBNf3YUu1YIDainBhI94CedlNeH",
    ///     ...     "chr1"
    ///     ... )
    fn get_sequence_by_name(
        &mut self,
        collection_digest: &str,
        sequence_name: &str,
    ) -> Option<PySequenceRecord> {
        self.inner
            .get_sequence_by_name(collection_digest, sequence_name)
            .ok()
            .map(|record| PySequenceRecord::from(record.clone()))
    }

    /// Get metadata for a single sequence by digest (no sequence data).
    ///
    /// Use this for lightweight lookups when you don't need the actual sequence.
    ///
    /// Args:
    ///     digest: Sequence digest (SHA-512/24u).
    ///
    /// Returns:
    ///     Optional[SequenceMetadata]: Sequence metadata if found, None otherwise.
    fn get_sequence_metadata(&self, digest: &str) -> Option<PySequenceMetadata> {
        self.inner
            .get_sequence_metadata(digest.as_bytes())
            .map(|meta| PySequenceMetadata::from(meta.clone()))
    }

    /// Extract a substring from a sequence.
    ///
    /// Retrieves a specific region from a sequence using 0-based, half-open
    /// coordinates [start, end). Automatically loads sequence data if not
    /// already cached (for lazy-loaded stores).
    ///
    /// Args:
    ///     seq_digest: Sequence digest (SHA-512/24u).
    ///     start: Start position (0-based, inclusive).
    ///     end: End position (0-based, exclusive).
    ///
    /// Returns:
    ///     Optional[str]: The substring sequence if found, None otherwise.
    ///
    /// Example:
    ///     >>> # Get first 1000 bases of chr1
    ///     >>> seq = store.get_substring("chr1_digest", 0, 1000)
    ///     >>> print(f"First 50bp: {seq[:50]}")
    fn get_substring(&mut self, seq_digest: &str, start: usize, end: usize) -> Option<String> {
        self.inner.get_substring(seq_digest, start, end).ok()
    }

    #[getter]
    fn cache_path(&self) -> Option<String> {
        self.inner.local_path().map(|p| p.display().to_string())
    }

    #[getter]
    fn remote_url(&self) -> Option<String> {
        self.inner.remote_source().map(|s| s.to_string())
    }

    #[getter]
    fn storage_mode(&self) -> PyStorageMode {
        self.inner.storage_mode().into()
    }

    /// List all sequences in the store (metadata only, no sequence data).
    ///
    /// Returns metadata for all sequences without loading sequence data.
    /// Use this for browsing/inventory operations.
    ///
    /// Returns:
    ///     list[SequenceMetadata]: List of sequence metadata objects.
    ///
    /// Example:
    ///     >>> for meta in store.list_sequences():
    ///     ...     print(f"{meta.name}: {meta.length} bp")
    fn list_sequences(&self) -> Vec<PySequenceMetadata> {
        self.inner
            .list_sequences()
            .into_iter()
            .map(|meta| PySequenceMetadata::from(meta))
            .collect()
    }

    /// List all collections in the store (metadata only, no sequence data).
    ///
    /// Returns metadata for all collections without loading sequence data.
    /// Use this for browsing/inventory operations.
    ///
    /// Returns:
    ///     list[SequenceCollectionMetadata]: List of collection metadata objects.
    ///
    /// Example:
    ///     >>> for meta in store.list_collections():
    ///     ...     print(f"{meta.digest}: {meta.n_sequences} sequences")
    fn list_collections(&self) -> Vec<PySequenceCollectionMetadata> {
        self.inner
            .list_collections()
            .into_iter()
            .map(|meta| PySequenceCollectionMetadata::from(meta))
            .collect()
    }

    /// Get a collection with all its sequences loaded.
    ///
    /// This loads the collection metadata and all sequence data, returning
    /// a complete `SequenceCollection` ready for use.
    ///
    /// Args:
    ///     collection_digest: The collection's SHA-512/24u digest.
    ///
    /// Returns:
    ///     SequenceCollection: The loaded collection with sequence data.
    ///
    /// Example:
    ///     >>> collection = store.get_collection("abc123")
    ///     >>> for seq in collection.sequences:
    ///     ...     print(f"{seq.metadata.name}: {seq.decode()[:20]}...")
    fn get_collection(&mut self, collection_digest: &str) -> PyResult<PySequenceCollection> {
        let collection = self.inner.get_collection(collection_digest).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error loading collection: {}", e))
        })?;
        Ok(PySequenceCollection::from(collection))
    }

    /// Get metadata for a collection by digest.
    ///
    /// Returns lightweight metadata without loading the full collection.
    /// Use this for quick lookups of collection information.
    ///
    /// Args:
    ///     collection_digest: The collection's SHA-512/24u digest.
    ///
    /// Returns:
    ///     Optional[SequenceCollectionMetadata]: Collection metadata if found, None otherwise.
    ///
    /// Example:
    ///     >>> meta = store.get_collection_metadata("uC_UorBNf3YUu1YIDainBhI94CedlNeH")
    ///     >>> if meta:
    ///     ...     print(f"Collection has {meta.n_sequences} sequences")
    fn get_collection_metadata(&self, collection_digest: &str) -> Option<PySequenceCollectionMetadata> {
        self.inner
            .get_collection_metadata(collection_digest)
            .map(|meta| PySequenceCollectionMetadata::from(meta.clone()))
    }

    /// Check if a collection is fully loaded.
    ///
    /// Returns True if the collection's sequence list is loaded in memory,
    /// False if it's only metadata (stub).
    ///
    /// Args:
    ///     collection_digest: The collection's SHA-512/24u digest.
    ///
    /// Returns:
    ///     bool: True if loaded, False otherwise.
    fn is_collection_loaded(&self, collection_digest: &str) -> bool {
        self.inner.is_collection_loaded(collection_digest)
    }

    /// Iterate over all collections with their sequences loaded.
    ///
    /// This loads all collection data upfront and returns a list of
    /// SequenceCollection objects with full sequence data.
    ///
    /// For browsing without loading data, use list_collections() instead.
    ///
    /// Returns:
    ///     list[SequenceCollection]: All collections with loaded sequences.
    ///
    /// Example:
    ///     >>> for coll in store.iter_collections():
    ///     ...     print(f"{coll.digest}: {len(coll.sequences)} sequences")
    fn iter_collections(&mut self) -> Vec<PySequenceCollection> {
        self.inner
            .iter_collections()
            .map(|coll| PySequenceCollection::from(coll))
            .collect()
    }

    /// Iterate over all sequences with their data loaded.
    ///
    /// This ensures all sequence data is loaded and returns a list of
    /// SequenceRecord objects with full sequence data.
    ///
    /// For browsing without loading data, use list_sequences() instead.
    ///
    /// Returns:
    ///     list[SequenceRecord]: All sequences with loaded data.
    ///
    /// Example:
    ///     >>> for seq in store.iter_sequences():
    ///     ...     print(f"{seq.metadata.name}: {seq.decode()[:20]}...")
    fn iter_sequences(&mut self) -> Vec<PySequenceRecord> {
        self.inner
            .iter_sequences()
            .map(|rec| PySequenceRecord::from(rec))
            .collect()
    }

    /// Returns statistics about the store.
    ///
    /// Returns:
    ///     dict: Dictionary with keys 'n_sequences', 'n_collections', 'n_collections_loaded', 'storage_mode'
    ///
    /// Note:
    ///     n_collections is the total number of collections (both loaded and stubs).
    ///     n_collections_loaded only reflects collections fully loaded in memory.
    ///     For remote stores, collections are loaded on-demand when accessed.
    ///
    /// Example:
    ///     >>> stats = store.stats()
    ///     >>> print(f"Store has {stats['n_sequences']} sequences")
    ///     >>> print(f"Collections: {stats['n_collections']} total, {stats['n_collections_loaded']} loaded")
    fn stats(&self) -> std::collections::HashMap<String, String> {
        let extended_stats = self.inner.stats_extended();
        let mut stats = std::collections::HashMap::new();
        stats.insert("n_sequences".to_string(), extended_stats.n_sequences.to_string());
        stats.insert("n_sequences_loaded".to_string(), extended_stats.n_sequences_loaded.to_string());
        stats.insert("n_collections".to_string(), extended_stats.n_collections.to_string());
        stats.insert("n_collections_loaded".to_string(), extended_stats.n_collections_loaded.to_string());
        stats.insert("storage_mode".to_string(), extended_stats.storage_mode);
        stats.insert("total_disk_size".to_string(), extended_stats.total_disk_size.to_string());
        stats
    }

    /// Write the store using its configured paths.
    ///
    /// Convenience method for disk-backed stores.
    /// Uses the store's own local_path and seqdata_path_template.
    fn write(&self) -> PyResult<()> {
        self.inner.write().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error writing store: {}", e))
        })
    }

    #[pyo3(signature = (root_path, seqdata_path_template=None))]
    fn write_store_to_dir(
        &self,
        root_path: &Bound<'_, PyAny>,
        seqdata_path_template: Option<&str>,
    ) -> PyResult<()> {
        let root_path = root_path.to_string();
        self.inner
            .write_store_to_dir(root_path, seqdata_path_template)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error writing store: {}", e))
            })
    }

    /// Load a local RefgetStore from a directory.
    ///
    /// Loads metadata from the local store immediately; sequence data is loaded
    /// on-demand when first accessed. This is efficient for large genomes where
    /// you may only need specific sequences.
    ///
    /// Args:
    ///     path (str or Path): Local directory containing the refget store (must have
    ///         rgstore.json and sequences.rgsi files).
    ///
    /// Returns:
    ///         RefgetStore: Store with metadata loaded, sequences lazy-loaded.
    ///
    /// Raises:
    ///     IOError: If the store directory or index files cannot be read.
    ///
    /// Example:
    ///     >>> from gtars.refget import RefgetStore
    ///     >>> store = RefgetStore.open_local("/data/hg38_store")
    ///     >>> seq = store.get_substring("chr1_digest", 0, 1000)
    #[classmethod]
    fn open_local(_cls: &Bound<'_, PyType>, path: &Bound<'_, PyAny>) -> PyResult<Self> {
        let path = path.to_string();
        let store = RefgetStore::open_local(path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening local store: {}", e))
        })?;
        Ok(Self { inner: store })
    }

    /// Open a remote RefgetStore with local caching.
    ///
    /// Fetches metadata (rgstore.json, sequences.rgsi) from a remote URL immediately.
    /// Sequence data (.seq files) are downloaded on-demand when first accessed and
    /// cached locally. This is ideal for working with large remote genomes where
    /// you only need specific sequences.
    ///
    /// By default, persistence is enabled (sequences are cached to disk).
    /// Call `disable_persistence()` after loading to keep only in memory.
    ///
    /// Args:
    ///     cache_path (str or Path): Local directory to cache downloaded metadata and sequences.
    ///         Created if it doesn't exist.
    ///     remote_url (str): Base URL of the remote refget store (e.g.,
    ///         "https://example.com/hg38" or "s3://bucket/hg38").
    ///
    /// Returns:
    ///     RefgetStore: Store with metadata loaded, sequences fetched on-demand.
    ///
    /// Raises:
    ///     IOError: If remote metadata cannot be fetched or cache cannot be written.
    ///
    /// Example:
    ///     >>> from gtars.refget import RefgetStore
    ///     >>> # With disk caching (default)
    ///     >>> store = RefgetStore.open_remote(
    ///     ...     "/data/cache/hg38",
    ///     ...     "https://refget-server.com/hg38"
    ///     ... )
    ///     >>> # Memory-only mode (no sequence data caching to disk)
    ///     >>> store = RefgetStore.open_remote(
    ///     ...     "/tmp/cache",
    ///     ...     "https://refget-server.com/hg38"
    ///     ... )
    ///     >>> store.disable_persistence()
    #[classmethod]
    fn open_remote(_cls: &Bound<'_, PyType>, cache_path: &Bound<'_, PyAny>, remote_url: &Bound<'_, PyAny>) -> PyResult<Self> {
        let cache_path_str = cache_path.to_string();
        let remote_url = remote_url.to_string();

        // Validate cache_path is not empty
        if cache_path_str.trim().is_empty() {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "cache_path cannot be empty"
            ));
        }

        // Expand tilde in path (Python users often pass paths with ~)
        let cache_path_expanded = if cache_path_str.starts_with("~/") {
            if let Some(home) = std::env::var_os("HOME") {
                let home_str = home.to_string_lossy();
                cache_path_str.replacen("~", &home_str, 1)
            } else {
                cache_path_str
            }
        } else if cache_path_str == "~" {
            if let Some(home) = std::env::var_os("HOME") {
                home.to_string_lossy().to_string()
            } else {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "Cannot expand '~': HOME environment variable not set"
                ));
            }
        } else {
            cache_path_str
        };

        let store = RefgetStore::open_remote(&cache_path_expanded, remote_url).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening remote store: {}", e))
        })?;
        Ok(Self { inner: store })
    }

    /// Export sequences from BED file regions to a FASTA file.
    ///
    /// Reads a BED file defining genomic regions and exports the sequences
    /// for those regions to a FASTA file.
    ///
    /// Args:
    ///     collection_digest: The collection's SHA-512/24u digest.
    ///     bed_file_path: Path to BED file defining regions.
    ///     output_file_path: Path to write the output FASTA file.
    ///
    /// Raises:
    ///     IOError: If files cannot be read/written or sequences not found.
    ///
    /// Example:
    ///     >>> store.export_fasta_from_regions(
    ///     ...     "uC_UorBNf3YUu1YIDainBhI94CedlNeH",
    ///     ...     "regions.bed",
    ///     ...     "output.fa"
    ///     ... )
    fn export_fasta_from_regions(
        &mut self,
        collection_digest: &str,
        bed_file_path: &Bound<'_, PyAny>,
        output_file_path: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let bed_file_path = bed_file_path.to_string();
        let output_file_path = output_file_path.to_string();
        self.inner
            .export_fasta_from_regions(collection_digest, &bed_file_path, &output_file_path)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                    "Error exporting FASTA from regions: {}",
                    e
                ))
            })
    }

    /// Get substrings for BED file regions as a list.
    ///
    /// Reads a BED file and returns a list of sequences for each region.
    /// This is a convenience wrapper around the iterator for cases where
    /// you want all results in memory at once.
    ///
    /// Note: For large BED files, consider using the Rust API's iterator
    /// directly to avoid loading all results into memory.
    ///
    /// Args:
    ///     collection_digest: The collection's SHA-512/24u digest.
    ///     bed_file_path: Path to BED file defining regions.
    ///
    /// Returns:
    ///     list[RetrievedSequence]: List of sequences with region metadata.
    ///
    /// Raises:
    ///     IOError: If files cannot be read or sequences not found.
    ///
    /// Example:
    ///     >>> sequences = store.substrings_from_regions(
    ///     ...     "uC_UorBNf3YUu1YIDainBhI94CedlNeH",
    ///     ...     "regions.bed"
    ///     ... )
    ///     >>> for seq in sequences:
    ///     ...     print(f"{seq.chrom_name}:{seq.start}-{seq.end}")
    fn substrings_from_regions(
        &mut self,
        collection_digest: &str,
        bed_file_path: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<PyRetrievedSequence>> {
        let bed_file_path = bed_file_path.to_string();
        // Get iterator and collect results
        let iter = self
            .inner
            .substrings_from_regions(collection_digest, &bed_file_path)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                    "Error getting substrings from regions: {}",
                    e
                ))
            })?;

        // Collect results, logging errors to stderr (matches Rust export_fasta_from_regions behavior)
        let mut py_results: Vec<PyRetrievedSequence> = Vec::new();
        for item in iter {
            match item {
                Ok(retrieved) => py_results.push(PyRetrievedSequence::from(retrieved)),
                Err(e) => eprintln!("{}", e),
            }
        }

        Ok(py_results)
    }

    /// Export sequences from a collection to a FASTA file.
    ///
    /// This method exports sequences from a specific sequence collection (genome assembly)
    /// using the collection digest. You can optionally filter which sequences to export
    /// by providing sequence names.
    ///
    /// **Use this method when:** You have a collection digest and want to export sequences
    /// from that specific collection, optionally filtering by chromosome/sequence names.
    ///
    /// **Contrast with export_fasta_by_digests:** That method exports sequences by their
    /// individual sequence digests (not collection digest) and bypasses collection
    /// information entirely.
    ///
    /// Args:
    ///     collection_digest: The collection's SHA-512/24u digest (identifies the genome
    ///         assembly/sequence collection).
    ///     output_path: Path to write the output FASTA file.
    ///     sequence_names: Optional list of sequence names to export (e.g., ["chr1", "chr2"]).
    ///         If None, exports all sequences in the collection.
    ///     line_width: Number of bases per line in output FASTA (default: 80).
    ///
    /// Raises:
    ///     IOError: If file cannot be written or collection/sequences not found.
    ///
    /// Example:
    ///     >>> # Export all sequences from a collection
    ///     >>> store.export_fasta("uC_UorBNf3YUu1YIDainBhI94CedlNeH", "output.fa", None, None)
    ///     >>> # Export only chr1 and chr2
    ///     >>> store.export_fasta("uC_UorBNf3YUu1YIDainBhI94CedlNeH", "output.fa", ["chr1", "chr2"], None)
    fn export_fasta(
        &mut self,
        collection_digest: &str,
        output_path: &Bound<'_, PyAny>,
        sequence_names: Option<Vec<String>>,
        line_width: Option<usize>,
    ) -> PyResult<()> {
        let output_path = output_path.to_string();
        let sequence_names_refs = sequence_names.as_ref().map(|names| {
            names.iter().map(|s| s.as_str()).collect::<Vec<&str>>()
        });
        self.inner
            .export_fasta(collection_digest, &output_path, sequence_names_refs, line_width)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                    "Error exporting FASTA: {}",
                    e
                ))
            })
    }

    /// Export sequences by their individual sequence digests to a FASTA file.
    ///
    /// This method exports sequences directly by their SHA-512/24u sequence digests,
    /// bypassing collection information. Useful when you have specific sequence
    /// identifiers and want to retrieve them regardless of which collection(s) they
    /// belong to.
    ///
    /// **Use this method when:** You have a list of sequence digests (SHA-512/24u)
    /// and want to export those specific sequences from the global sequence store.
    ///
    /// **Contrast with export_fasta:** That method uses a collection digest and
    /// optionally filters by sequence names within that collection.
    ///
    /// Args:
    ///     seq_digests: List of SHA-512/24u sequence digests (not collection digests)
    ///         to export. Each digest identifies an individual sequence.
    ///     output_path: Path to write the output FASTA file.
    ///     line_width: Number of bases per line in output FASTA (default: 80).
    ///
    /// Raises:
    ///     IOError: If file cannot be written or sequences not found.
    ///
    /// Example:
    ///     >>> # Export specific sequences by their digests
    ///     >>> seq_digests = [
    ///     ...     "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2",
    ///     ...     "bXE123dAxcJAqme6QYQ7EZ07-fiw8Kw2"
    ///     ... ]
    ///     >>> store.export_fasta_by_digests(seq_digests, "output.fa", None)
    fn export_fasta_by_digests(
        &mut self,
        seq_digests: Vec<String>,
        output_path: &Bound<'_, PyAny>,
        line_width: Option<usize>,
    ) -> PyResult<()> {
        let output_path = output_path.to_string();
        let digests_refs: Vec<&str> = seq_digests.iter().map(|s| s.as_str()).collect();
        self.inner
            .export_fasta_by_digests(digests_refs, &output_path, line_width)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                    "Error exporting FASTA by digests: {}",
                    e
                ))
            })
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __repr__(&self) -> String {
        let (n_sequences, n_collections_loaded, mode_str) = self.inner.stats();
        let persist_str = if self.inner.is_persisting() { "persist=on" } else { "persist=off" };
        let quiet_str = if self.inner.is_quiet() { "quiet=on" } else { "quiet=off" };

        let location = if let Some(remote) = self.inner.remote_source() {
            let cache = self.inner.local_path()
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| "None".to_string());
            format!("cache='{}', remote='{}'", cache, remote)
        } else if let Some(path) = self.inner.local_path() {
            format!("cache='{}'", path.display())
        } else {
            "memory-only".to_string()
        };

        format!(
            "RefgetStore(n_sequences={}, n_collections_loaded={}, mode={}, {}, {}, {})",
            n_sequences, n_collections_loaded, mode_str, persist_str, quiet_str, location
        )
    }

    /// Return the number of sequences in the store.
    ///
    /// Allows using `len(store)` in Python to get the total number of sequences
    /// across all collections.
    ///
    /// Returns:
    ///     int: Total number of sequences in the store.
    ///
    /// Example:
    ///     >>> store = RefgetStore(StorageMode.Encoded)
    ///     >>> store.import_fasta("genome.fa")
    ///     >>> print(f"Store contains {len(store)} sequences")
    fn __len__(&self) -> usize {
        self.inner.sequence_digests().count()
    }

    /// Iterate over all sequences in the store.
    ///
    /// Allows using `for sequence in store:` in Python to iterate over all
    /// sequence metadata in the store.
    ///
    /// Yields:
    ///     SequenceMetadata: Metadata for each sequence in the store.
    ///
    /// Example:
    ///     >>> store = RefgetStore(StorageMode.Encoded)
    ///     >>> store.import_fasta("genome.fa")
    ///     >>> for seq_meta in store:
    ///     ...     print(f"{seq_meta.name}: {seq_meta.length} bp")
    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<PyRefgetStoreIterator> {
        let sequences = slf.inner.list_sequences()
            .into_iter()
            .map(|meta| PySequenceMetadata::from(meta))
            .collect();
        Ok(PyRefgetStoreIterator {
            sequences,
            index: 0,
        })
    }
}

/// Iterator for RefgetStore that yields SequenceMetadata.
#[pyclass]
pub struct PyRefgetStoreIterator {
    sequences: Vec<PySequenceMetadata>,
    index: usize,
}

#[pymethods]
impl PyRefgetStoreIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PySequenceMetadata> {
        if slf.index < slf.sequences.len() {
            let item = slf.sequences[slf.index].clone();
            slf.index += 1;
            Some(item)
        } else {
            None
        }
    }
}

// This represents the Python module to be created
#[pymodule]
pub fn refget(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sha512t24u_digest, m)?)?;
    m.add_function(wrap_pyfunction!(md5_digest, m)?)?;
    m.add_function(wrap_pyfunction!(digest_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(compute_fai, m)?)?;
    m.add_function(wrap_pyfunction!(load_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(digest_sequence, m)?)?;
    m.add_class::<PyAlphabetType>()?;
    m.add_class::<PySequenceMetadata>()?;
    m.add_class::<PyFaiMetadata>()?;
    m.add_class::<PyFaiRecord>()?;
    m.add_class::<PySequenceRecord>()?;
    m.add_class::<PySeqColDigestLvl1>()?;
    m.add_class::<PySequenceCollectionMetadata>()?;
    m.add_class::<PySequenceCollection>()?;
    m.add_class::<PyStorageMode>()?;
    m.add_class::<PyRefgetStore>()?;
    m.add_class::<PyRetrievedSequence>()?;
    Ok(())
}
