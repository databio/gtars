use pyo3::prelude::*;

#[pyclass(name="FragmentTokenizer")]
pub struct PyFragmentTokenizer {
    pub tokenizer: genimtools::tokenizers::FragmentTokenizer,
}

#[pymethods]
impl PyFragmentTokenizer {
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let path = std::path::Path::new(&path);
        let tokenizer = genimtools::tokenizers::FragmentTokenizer::try_from(path)?;
        Ok(PyFragmentTokenizer { tokenizer })
    }

    pub fn tokenize_fragments(&self, file: String, out_path: Option<String>, filter: Option<Vec<String>>) -> PyResult<()> {
        let path = std::path::Path::new(&file);
        let out_path = out_path.unwrap_or("".to_string());
        let out_path = std::path::Path::new(&out_path);
        match filter {
            Some(filter) => self.tokenizer.tokenize_fragments_with_filter(path, out_path, filter),
            None => self.tokenizer.tokenize_fragments(path, out_path),
        }?;
        Ok(())
    }
}