use std::path::Path;
use pyo3::prelude::*;

use gtars::igd::search::igd_search;

#[pyclass(name="IGD")]
pub struct IGD;

#[pymethods]
impl IGD {

    #[classmethod]
    pub fn search(database_path: &String, query_file_path: &String) ->  Ok() {

        igd_search(database_path, query_file_path).unwrap()


    }
}