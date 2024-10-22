// TODO: stil a work in progress
use pyo3::prelude::*;

use anyhow::Result;

use std::path::Path;


use gtars::tokenizers::TokenizerConfig;

use super::{
    PyMetaTokenizer,
    PyTreeTokenizer
};

#[pyclass(name="TokenizerBuilder")]
pub struct PyTokenizerBuilder;

#[pymethods]
impl PyTokenizerBuilder {

    #[classmethod]
    pub fn from_toml(path: String) -> Result<PyObject> {
        let config = TokenizerConfig::new(Path::new(&path))?;
        
        match config.tokenizer_type {
            Some(tokenizer_type) => {
                match tokenizer_type.as_str() {
                    "tree" => {
                        let t = PyTreeTokenizer::new(path)?;
                        t.to_object()
                    },
                    "meta" => {
                        PyMetaTokenizer::new(path)
                    },
                    _ => {
                        anyhow::bail!("Tokenizer type {} not supported", tokenizer_type)
                    }
                }
            },
            None => {
                println!("No tokenizer type found in config file. Instantiating a default TreeTokenizer. Note that this may lead to unexpected behavior.");
                PyTreeTokenizer::new(path)
            }
        };
    
    }
}