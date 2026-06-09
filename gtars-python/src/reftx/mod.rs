//! Python bindings for the `gtars-refget` transcript store
//! (`gtars_refget::transcripts`): transcript store, MANE index, and
//! coordinate mapping for HGVS-to-genomic projection.
//!
//! Exposed as the `gtars.reftx` submodule.

use std::fmt::Display;
use std::sync::Arc;

use pyo3::create_exception;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyType};

use gtars_refget::transcripts::builder::TxStoreBuilder as RsTxStoreBuilder;
use gtars_refget::transcripts::mapper::{
    CoordinateMapper as RsCoordinateMapper, MappingError as RsMappingError,
};
use gtars_refget::transcripts::mmap::TxStore as RsTxStore;
use gtars_refget::transcripts::models::{
    Exon as RsExon, ManeStatus as RsManeStatus, Strand as RsStrand,
};
use gtars_refget::transcripts::store::{ReadonlyTxStore as RsReadonlyTxStore, TranscriptRef};

#[cfg(feature = "vrs")]
use gtars_vrs::provider::TxProvider as RsReftxProvider;

// ---------------------------------------------------------------------------
// Exception classes
// ---------------------------------------------------------------------------

create_exception!(gtars.reftx, MappingError, pyo3::exceptions::PyException);
create_exception!(gtars.reftx, TxStoreError, pyo3::exceptions::PyException);

#[inline]
fn map_mapping_err<E: Display>(e: E) -> PyErr {
    MappingError::new_err(format!("{}", e))
}

#[inline]
fn map_store_err<E: Display>(e: E) -> PyErr {
    TxStoreError::new_err(format!("{}", e))
}

// ---------------------------------------------------------------------------
// Strand
// ---------------------------------------------------------------------------

#[pyclass(name = "Strand", module = "gtars.reftx", eq, eq_int)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    fn from_rs(s: RsStrand) -> Self {
        match s {
            RsStrand::Forward => Strand::Plus,
            RsStrand::Reverse => Strand::Minus,
        }
    }
}

#[pymethods]
impl Strand {
    fn __repr__(&self) -> String {
        match self {
            Strand::Plus => "Strand.Plus".to_string(),
            Strand::Minus => "Strand.Minus".to_string(),
        }
    }

    fn __str__(&self) -> String {
        match self {
            Strand::Plus => "+".to_string(),
            Strand::Minus => "-".to_string(),
        }
    }

    #[classmethod]
    fn from_str(_cls: &Bound<'_, PyType>, s: &str) -> PyResult<Self> {
        match s {
            "+" | "Plus" | "forward" | "Forward" | "1" | "+1" => Ok(Strand::Plus),
            "-" | "Minus" | "reverse" | "Reverse" | "-1" => Ok(Strand::Minus),
            other => Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Unrecognized strand: {:?}",
                other
            ))),
        }
    }
}

// ---------------------------------------------------------------------------
// ManeStatus
// ---------------------------------------------------------------------------

#[pyclass(name = "ManeStatus", module = "gtars.reftx")]
#[derive(Debug, Clone, Copy)]
pub struct ManeStatusPy {
    #[pyo3(get, set)]
    pub select: bool,
    #[pyo3(get, set)]
    pub plus_clinical: bool,
}

#[pymethods]
impl ManeStatusPy {
    #[new]
    #[pyo3(signature = (select=false, plus_clinical=false))]
    fn new(select: bool, plus_clinical: bool) -> Self {
        Self {
            select,
            plus_clinical,
        }
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("select", self.select)?;
        d.set_item("plus_clinical", self.plus_clinical)?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "ManeStatus(select={}, plus_clinical={})",
            if self.select { "True" } else { "False" },
            if self.plus_clinical { "True" } else { "False" }
        )
    }
}

impl ManeStatusPy {
    fn from_rs(m: RsManeStatus) -> Self {
        Self {
            select: m.mane_select,
            plus_clinical: m.mane_clinical,
        }
    }
}

// ---------------------------------------------------------------------------
// Exon
// ---------------------------------------------------------------------------

#[pyclass(name = "Exon", module = "gtars.reftx")]
#[derive(Debug, Clone, Copy)]
pub struct ExonPy {
    #[pyo3(get, set)]
    pub start: u32,
    #[pyo3(get, set)]
    pub end: u32,
}

#[pymethods]
impl ExonPy {
    #[new]
    #[pyo3(signature = (start, end))]
    fn new(start: u32, end: u32) -> Self {
        Self { start, end }
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("start", self.start)?;
        d.set_item("end", self.end)?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!("Exon(start={}, end={})", self.start, self.end)
    }
}

impl ExonPy {
    fn from_rs(e: RsExon) -> Self {
        Self {
            start: e.start,
            end: e.end,
        }
    }
}

// ---------------------------------------------------------------------------
// Transcript
// ---------------------------------------------------------------------------

#[pyclass(name = "Transcript", module = "gtars.reftx")]
#[derive(Debug, Clone)]
pub struct TranscriptPy {
    #[pyo3(get)]
    pub accession: String,
    #[pyo3(get)]
    pub gene: Option<String>,
    /// Chromosome refget digest as a `SQ.<base64url>` accession.
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub strand: Strand,
    #[pyo3(get)]
    pub cds_start: Option<u32>,
    #[pyo3(get)]
    pub cds_end: Option<u32>,
    #[pyo3(get)]
    pub exons: Vec<ExonPy>,
    #[pyo3(get)]
    pub mane: Option<ManeStatusPy>,
}

#[pymethods]
impl TranscriptPy {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("accession", &self.accession)?;
        d.set_item("gene", &self.gene)?;
        d.set_item("chrom", &self.chrom)?;
        d.set_item(
            "strand",
            match self.strand {
                Strand::Plus => "+",
                Strand::Minus => "-",
            },
        )?;
        d.set_item("cds_start", self.cds_start)?;
        d.set_item("cds_end", self.cds_end)?;
        let exons: Vec<Bound<'_, PyDict>> = self
            .exons
            .iter()
            .map(|e| e.to_dict(py))
            .collect::<PyResult<_>>()?;
        d.set_item("exons", exons)?;
        d.set_item(
            "mane",
            self.mane.map(|m| m.to_dict(py)).transpose()?,
        )?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "Transcript(accession={:?}, gene={:?}, chrom={:?}, strand={:?}, exons={})",
            self.accession,
            self.gene,
            self.chrom,
            self.strand,
            self.exons.len()
        )
    }
}

impl TranscriptPy {
    fn from_rs(tx: gtars_refget::transcripts::models::Transcript) -> Self {
        let chrom = format!("SQ.{}", base64_url::encode(&tx.chrom_digest));
        let mane = if tx.mane.mane_select || tx.mane.mane_clinical {
            Some(ManeStatusPy::from_rs(tx.mane))
        } else {
            None
        };
        Self {
            accession: tx.accession,
            gene: if tx.gene.is_empty() {
                None
            } else {
                Some(tx.gene)
            },
            chrom,
            strand: Strand::from_rs(tx.strand),
            cds_start: tx.cds_start,
            cds_end: tx.cds_end,
            exons: tx.exons.into_iter().map(ExonPy::from_rs).collect(),
            mane,
        }
    }

    fn from_ref(tx: TranscriptRef<'_>) -> Self {
        Self::from_rs((*tx).clone())
    }
}

// Decode a 24-byte chromosome digest from a refget-style "SQ.<base64url>"
// accession or a plain base64url string.
fn decode_chrom_digest(s: &str) -> PyResult<[u8; 24]> {
    let body = s.strip_prefix("SQ.").unwrap_or(s);
    let bytes = base64_url::decode(body).map_err(|e| {
        TxStoreError::new_err(format!(
            "Invalid base64url chrom accession {:?}: {}",
            s, e
        ))
    })?;
    if bytes.len() != 24 {
        return Err(TxStoreError::new_err(format!(
            "Chrom accession must decode to 24 bytes, got {} (input: {:?})",
            bytes.len(),
            s
        )));
    }
    let mut out = [0u8; 24];
    out.copy_from_slice(&bytes);
    Ok(out)
}

// ---------------------------------------------------------------------------
// TxStoreBuilder
// ---------------------------------------------------------------------------

#[pyclass(name = "TxStoreBuilder", module = "gtars.reftx", unsendable)]
pub struct TxStoreBuilderPy {
    inner: RsTxStoreBuilder,
}

#[pymethods]
impl TxStoreBuilderPy {
    #[new]
    fn new() -> Self {
        Self {
            inner: RsTxStoreBuilder::new(),
        }
    }

    /// Add a transcript directly. Accepts either a `Transcript` instance or a
    /// dict with the same keys as `Transcript.to_dict()`.
    fn add_transcript(&mut self, py: Python<'_>, value: PyObject) -> PyResult<()> {
        let tx = if let Ok(t) = value.extract::<TranscriptPy>(py) {
            t
        } else if let Ok(d) = value.downcast_bound::<PyDict>(py) {
            transcript_from_dict(d)?
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "add_transcript expects a Transcript or dict",
            ));
        };

        let chrom_digest = decode_chrom_digest(&tx.chrom)?;
        let strand = match tx.strand {
            Strand::Plus => RsStrand::Forward,
            Strand::Minus => RsStrand::Reverse,
        };
        let exons = tx
            .exons
            .iter()
            .map(|e| RsExon {
                start: e.start,
                end: e.end,
            })
            .collect();
        let mane = match tx.mane {
            Some(m) => RsManeStatus {
                mane_select: m.select,
                mane_clinical: m.plus_clinical,
            },
            None => RsManeStatus::default(),
        };
        self.inner.transcripts.push(gtars_refget::transcripts::models::Transcript {
            accession: tx.accession,
            gene: tx.gene.unwrap_or_default(),
            chrom_digest,
            strand,
            cds_start: tx.cds_start,
            cds_end: tx.cds_end,
            exons,
            mane,
        });
        Ok(())
    }

    /// Register a chromosome name to its `SQ.<base64url>` refget accession.
    fn add_chrom_mapping(&mut self, name: &str, accession: &str) -> PyResult<()> {
        let digest = decode_chrom_digest(accession)?;
        self.inner.add_chrom_mapping(name, digest);
        Ok(())
    }

    /// Ingest a cdot JSON file. Auto-detects `.gz`.
    fn load_cdot(&mut self, py: Python<'_>, path: &str) -> PyResult<usize> {
        let path = path.to_string();
        py.allow_threads(|| self.inner.ingest_cdot(&path))
            .map_err(map_store_err)
    }

    /// Alias for `load_cdot`; the auto-detected gzip path also works through it.
    fn load_cdot_gz(&mut self, py: Python<'_>, path: &str) -> PyResult<usize> {
        self.load_cdot(py, path)
    }

    /// Load MANE flags from an NCBI MANE summary TSV (auto-detects `.gz`).
    fn load_mane_summary(&mut self, py: Python<'_>, path: &str) -> PyResult<usize> {
        let path = path.to_string();
        py.allow_threads(|| self.inner.load_mane_summary(&path))
            .map_err(map_store_err)
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Build and write the binary store at `out_path`.
    fn build(&mut self, py: Python<'_>, out_path: &str) -> PyResult<()> {
        let out_path = out_path.to_string();
        py.allow_threads(|| self.inner.build(&out_path))
            .map_err(map_store_err)
    }
}

fn transcript_from_dict(d: &Bound<'_, PyDict>) -> PyResult<TranscriptPy> {
    fn get<'a, 'py, T>(d: &'a Bound<'py, PyDict>, key: &str) -> PyResult<T>
    where
        T: for<'b, 'c> FromPyObject<'b, 'c>,
    {
        let item = d
            .get_item(key)?
            .ok_or_else(|| {
                pyo3::exceptions::PyKeyError::new_err(format!("Missing field {:?}", key))
            })?;
        item.extract::<T>().map_err(Into::into)
    }
    fn get_opt<'a, 'py, T>(d: &'a Bound<'py, PyDict>, key: &str) -> PyResult<Option<T>>
    where
        T: for<'b, 'c> FromPyObject<'b, 'c>,
    {
        match d.get_item(key)? {
            None => Ok(None),
            Some(v) if v.is_none() => Ok(None),
            Some(v) => v.extract::<T>().map(Some).map_err(Into::into),
        }
    }

    let accession: String = get(d, "accession")?;
    let gene: Option<String> = get_opt(d, "gene")?;
    let chrom: String = get(d, "chrom")?;
    let strand_raw = d.get_item("strand")?.ok_or_else(|| {
        pyo3::exceptions::PyKeyError::new_err("Missing field \"strand\"")
    })?;
    let strand: Strand = if let Ok(s) = strand_raw.extract::<Strand>() {
        s
    } else {
        let s: String = strand_raw.extract()?;
        match s.as_str() {
            "+" | "Plus" | "forward" | "Forward" | "1" | "+1" => Strand::Plus,
            "-" | "Minus" | "reverse" | "Reverse" | "-1" => Strand::Minus,
            other => {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "Unrecognized strand: {:?}",
                    other
                )))
            }
        }
    };
    let cds_start: Option<u32> = get_opt(d, "cds_start")?;
    let cds_end: Option<u32> = get_opt(d, "cds_end")?;

    let exons_obj = d
        .get_item("exons")?
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err("Missing field \"exons\""))?;
    let mut exons = Vec::new();
    for item in exons_obj.try_iter()? {
        let item = item?;
        let exon = if let Ok(e) = item.extract::<ExonPy>() {
            e
        } else if let Ok(d) = item.downcast::<PyDict>() {
            let start: u32 = get(d, "start")?;
            let end: u32 = get(d, "end")?;
            ExonPy { start, end }
        } else if let Ok((s, e)) = item.extract::<(u32, u32)>() {
            ExonPy { start: s, end: e }
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "exons entries must be Exon, dict, or (start, end) tuple",
            ));
        };
        exons.push(exon);
    }

    let mane: Option<ManeStatusPy> = match d.get_item("mane")? {
        None => None,
        Some(v) if v.is_none() => None,
        Some(v) => {
            if let Ok(m) = v.extract::<ManeStatusPy>() {
                Some(m)
            } else if let Ok(dd) = v.downcast::<PyDict>() {
                let select: bool = get_opt(dd, "select")?.unwrap_or(false);
                let plus_clinical: bool = get_opt(dd, "plus_clinical")?.unwrap_or(false);
                Some(ManeStatusPy {
                    select,
                    plus_clinical,
                })
            } else {
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "mane must be ManeStatus or dict",
                ));
            }
        }
    };

    Ok(TranscriptPy {
        accession,
        gene,
        chrom,
        strand,
        cds_start,
        cds_end,
        exons,
        mane,
    })
}

// ---------------------------------------------------------------------------
// ReadonlyTxStore
// ---------------------------------------------------------------------------

#[pyclass(name = "ReadonlyTxStore", module = "gtars.reftx")]
#[derive(Clone)]
pub struct ReadonlyTxStorePy {
    pub(crate) inner: Arc<RsReadonlyTxStore>,
}

#[pymethods]
impl ReadonlyTxStorePy {
    /// Open a `.reftx` store from disk (mmap-backed).
    #[classmethod]
    fn open(_cls: &Bound<'_, PyType>, py: Python<'_>, path: &str) -> PyResult<Self> {
        let path = path.to_string();
        let store = py
            .allow_threads(|| RsTxStore::open(&path))
            .map_err(map_store_err)?;
        Ok(Self {
            inner: Arc::new(store.into_readonly()),
        })
    }

    fn lookup(&self, py: Python<'_>, accession: &str) -> Option<TranscriptPy> {
        py.allow_threads(|| {
            self.inner
                .lookup(accession)
                .map(TranscriptPy::from_ref)
        })
    }

    fn lookup_mane(&self, py: Python<'_>, gene: &str) -> Option<TranscriptPy> {
        py.allow_threads(|| self.inner.lookup_mane(gene).map(TranscriptPy::from_rs))
    }

    fn has_mane_index(&self) -> bool {
        self.inner.has_mane_index()
    }

    /// Returns the `SQ.<base64url>` accession of the chromosome for the given
    /// transcript, or `None` if the accession is not present.
    fn get_chrom_accession(&self, py: Python<'_>, accession: &str) -> Option<String> {
        py.allow_threads(|| {
            self.inner
                .lookup(accession)
                .map(|tx| format!("SQ.{}", base64_url::encode(&tx.chrom_digest)))
        })
    }

    fn get_strand(&self, py: Python<'_>, accession: &str) -> Option<Strand> {
        py.allow_threads(|| self.inner.lookup(accession).map(|tx| Strand::from_rs(tx.strand)))
    }

    fn __len__(&self) -> usize {
        self.inner.len() as usize
    }
}

// ---------------------------------------------------------------------------
// CoordinateMapper
// ---------------------------------------------------------------------------

#[pyclass(name = "CoordinateMapper", module = "gtars.reftx")]
pub struct CoordinateMapperPy {
    store: Arc<RsReadonlyTxStore>,
}

#[pymethods]
impl CoordinateMapperPy {
    #[new]
    fn new(store: &ReadonlyTxStorePy) -> Self {
        Self {
            store: store.inner.clone(),
        }
    }

    /// Map a coding (`c.`) coordinate to a 0-based genomic position.
    #[pyo3(signature = (accession, c_pos, datum=None))]
    fn c_to_g(
        &self,
        py: Python<'_>,
        accession: &str,
        c_pos: i64,
        datum: Option<i32>,
    ) -> PyResult<u64> {
        let is_cds_end = matches!(datum, Some(1));
        let acc = accession.to_string();
        py.allow_threads(|| {
            let mapper = RsCoordinateMapper::new(&self.store);
            mapper
                .c_to_g_full(&acc, c_pos, 0, is_cds_end)
                .map(|r| r.position)
        })
        .map_err(map_mapping_err)
    }

    /// Map a non-coding (`n.`) coordinate to a 0-based genomic position.
    fn n_to_g(&self, py: Python<'_>, accession: &str, n_pos: u32) -> PyResult<u64> {
        let acc = accession.to_string();
        py.allow_threads(|| {
            let mapper = RsCoordinateMapper::new(&self.store);
            mapper.n_to_g(&acc, n_pos as u64).map(|r| r.position)
        })
        .map_err(map_mapping_err)
    }

    /// Map a `c.` coordinate and return full metadata.
    ///
    /// Returns `{"chrom": str, "chrom_accession": str, "genomic_pos": int, "strand": Strand}`.
    /// Note: `chrom` and `chrom_accession` are the same `SQ.<base64url>` value;
    /// chromosome name lookup is not available from the binary store.
    #[pyo3(signature = (accession, c_pos, datum=None))]
    fn c_to_g_full<'py>(
        &self,
        py: Python<'py>,
        accession: &str,
        c_pos: i64,
        datum: Option<i32>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let is_cds_end = matches!(datum, Some(1));
        let acc = accession.to_string();
        let store = self.store.clone();
        let (result, strand) = py
            .allow_threads(|| {
                let mapper = RsCoordinateMapper::new(&store);
                let r = mapper.c_to_g_full(&acc, c_pos, 0, is_cds_end)?;
                let strand = store.lookup(&acc).map(|tx| Strand::from_rs(tx.strand));
                Ok::<_, RsMappingError>((r, strand))
            })
            .map_err(map_mapping_err)?;
        build_full_dict(py, &result, strand)
    }

    /// Map a `n.` coordinate and return full metadata.
    fn n_to_g_full<'py>(
        &self,
        py: Python<'py>,
        accession: &str,
        n_pos: u32,
    ) -> PyResult<Bound<'py, PyDict>> {
        let acc = accession.to_string();
        let store = self.store.clone();
        let (result, strand) = py
            .allow_threads(|| {
                let mapper = RsCoordinateMapper::new(&store);
                let r = mapper.n_to_g(&acc, n_pos as u64)?;
                let strand = store.lookup(&acc).map(|tx| Strand::from_rs(tx.strand));
                Ok::<_, RsMappingError>((r, strand))
            })
            .map_err(map_mapping_err)?;
        build_full_dict(py, &result, strand)
    }

    /// Map a `c.` coordinate via a gene's MANE Select transcript.
    ///
    /// Returns `{"accession": str, "chrom": str, "chrom_accession": str,
    /// "genomic_pos": int, "strand": Strand}`.
    #[pyo3(signature = (gene, c_pos, datum=None))]
    fn c_to_g_by_gene<'py>(
        &self,
        py: Python<'py>,
        gene: &str,
        c_pos: i64,
        datum: Option<i32>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let is_cds_end = matches!(datum, Some(1));
        let g = gene.to_string();
        let store = self.store.clone();
        let (acc, result, strand) = py
            .allow_threads(|| {
                let mapper = RsCoordinateMapper::new(&store);
                let (acc, r) = mapper.c_to_g_by_gene(&g, c_pos, 0, is_cds_end)?;
                let strand = store.lookup(&acc).map(|tx| Strand::from_rs(tx.strand));
                Ok::<_, RsMappingError>((acc, r, strand))
            })
            .map_err(map_mapping_err)?;
        let d = build_full_dict(py, &result, strand)?;
        d.set_item("accession", acc)?;
        Ok(d)
    }
}

fn build_full_dict<'py>(
    py: Python<'py>,
    result: &gtars_refget::transcripts::mapper::MappingResult,
    strand: Option<Strand>,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    let chrom_acc = format!("SQ.{}", base64_url::encode(&result.chrom_digest));
    // `chrom` is exposed for parity with hgvs APIs that key off chromosome
    // labels; here we only have the digest, so both keys carry it.
    d.set_item("chrom", &chrom_acc)?;
    d.set_item("chrom_accession", &chrom_acc)?;
    d.set_item("genomic_pos", result.position)?;
    d.set_item("strand", strand)?;
    Ok(d)
}

// ---------------------------------------------------------------------------
// ReftxProvider (gated on the `vrs` feature)
// ---------------------------------------------------------------------------

#[cfg(feature = "vrs")]
#[pyclass(name = "ReftxProvider", module = "gtars.reftx")]
#[derive(Clone)]
pub struct ReftxProviderPy {
    pub(crate) inner: RsReftxProvider,
}

#[cfg(feature = "vrs")]
#[pymethods]
impl ReftxProviderPy {
    #[new]
    fn new(store: &ReadonlyTxStorePy) -> Self {
        Self {
            inner: RsReftxProvider::new(store.inner.clone()),
        }
    }
}

// ---------------------------------------------------------------------------
// Module registration
// ---------------------------------------------------------------------------

#[pymodule]
pub fn reftx(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Strand>()?;
    m.add_class::<ManeStatusPy>()?;
    m.add_class::<ExonPy>()?;
    m.add_class::<TranscriptPy>()?;
    m.add_class::<TxStoreBuilderPy>()?;
    m.add_class::<ReadonlyTxStorePy>()?;
    m.add_class::<CoordinateMapperPy>()?;

    #[cfg(feature = "vrs")]
    m.add_class::<ReftxProviderPy>()?;

    m.add("MappingError", py.get_type::<MappingError>())?;
    m.add("TxStoreError", py.get_type::<TxStoreError>())?;

    Ok(())
}
