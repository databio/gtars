use gtars::tokenizers::FragmentTokenizer;
use gtars::tokenizers::TreeTokenizer;
use pyo3::prelude::*;

use super::PyTokenizedRegionSet;
use super::PyUniverse;

#[pyclass(name = "FragmentTokenizer", module="gtars.tokenizers")]
pub struct PyFragmentTokenizer {
    pub tokenizer: gtars::tokenizers::FragmentTokenizer<TreeTokenizer>,
    pub universe: Py<PyUniverse>, // this is a Py-wrapped version self.tokenizer.universe for performance reasons
}

#[pymethods]
impl PyFragmentTokenizer {
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        Python::with_gil(|py| {
            let path = std::path::Path::new(&path);
            let tokenizer = gtars::tokenizers::TreeTokenizer::try_from(path)?;
            let frag_tokenizer = FragmentTokenizer::new(tokenizer);
            let py_universe: PyUniverse = frag_tokenizer.tokenizer.universe.to_owned().into();
            let py_universe_bound = Py::new(py, py_universe)?;
            Ok(PyFragmentTokenizer {
                tokenizer: frag_tokenizer,
                universe: py_universe_bound,
            })
        })
    }

    pub fn tokenize_fragments_to_gtoks(
        &self,
        file: String,
        out_path: Option<String>,
        filter: Option<Vec<String>>,
    ) -> PyResult<()> {
        let path = std::path::Path::new(&file);
        let out_path = out_path.unwrap_or("".to_string());
        let out_path = std::path::Path::new(&out_path);
        match filter {
            Some(filter) => self
                .tokenizer
                .tokenize_fragments_to_gtoks_with_filter(path, out_path, filter),
            None => self.tokenizer.tokenize_fragments_to_gtoks(path, out_path),
        }?;
        Ok(())
    }

    pub fn tokenize_fragments(
        &self,
        file: String,
        filter: Option<Vec<String>>,
    ) -> PyResult<Vec<PyTokenizedRegionSet>> {
        let path = std::path::Path::new(&file);
        match filter {
            Some(filter) => {
                let tokenized_region_sets = self
                    .tokenizer
                    .tokenize_fragments_with_filter(path, filter)?;
                Python::with_gil(|py| {
                    let py_tokenized_regions_sets = tokenized_region_sets
                        .into_iter()
                        .map(|trs| PyTokenizedRegionSet {
                            ids: trs.ids,
                            curr: 0,
                            universe: self.universe.clone_ref(py),
                        })
                        .collect();
                    Ok(py_tokenized_regions_sets)
                })
            }
            None => {
                let tokenized_region_sets = self.tokenizer.tokenize_fragments(path)?;
                Python::with_gil(|py| {
                    let py_tokenized_regions_sets = tokenized_region_sets
                        .into_iter()
                        .map(|trs| PyTokenizedRegionSet {
                            ids: trs.ids,
                            curr: 0,
                            universe: self.universe.clone_ref(py),
                        })
                        .collect();
                    Ok(py_tokenized_regions_sets)
                })
            }
        }
    }
}
