use pyo3::prelude::*;
use pyo3::types::PyDict;

use rayon::prelude::*;

use gtars::tokenizers::utils::r2v::{create_instances, Algorithm, Instance};

#[pyfunction(name = "create_instances")]
pub fn py_create_instances(
    sequences: &Bound<'_, PyAny>,
    window_size: usize,
    algorithm: &str,
) -> PyResult<Vec<Py<PyDict>>> {
    Python::with_gil(|py| {
        let algorithm = match algorithm {
            "cbow" => Algorithm::Cbow,
            "sg" => Algorithm::Sg,
            _ => return Err(pyo3::exceptions::PyValueError::new_err("Invalid algorithm")),
        };

        if let Ok(sequence) = sequences.extract::<Vec<u32>>() {
            let result = create_instances(&sequence, window_size, algorithm);
            let mapped_dicts = result
                .into_iter()
                .map(|instance| {
                    let dict = PyDict::new_bound(py);
                    match instance {
                        Instance::Cbow {
                            context_ids,
                            target_id,
                        } => {
                            dict.set_item("context_ids", context_ids).unwrap();
                            dict.set_item("target_id", target_id).unwrap();
                        }
                        Instance::Sg {
                            center_id,
                            context_ids,
                        } => {
                            dict.set_item("center_id", center_id).unwrap();
                            dict.set_item("context_ids", context_ids).unwrap();
                        }
                    }
                    dict.into()
                })
                .collect::<Vec<Py<PyDict>>>();
            Ok(mapped_dicts)
        } else if let Ok(sequences) = sequences.extract::<Vec<Vec<u32>>>() {
            let result: Vec<Vec<Instance>> = sequences
                .par_iter()
                .map(|sequence| create_instances(sequence, window_size, algorithm))
                .collect();

            let mapped_dicts = result
                .into_iter()
                .flat_map(|instances| {
                    instances.into_iter().map(|instance| {
                        let dict = PyDict::new_bound(py);
                        match instance {
                            Instance::Cbow {
                                context_ids,
                                target_id,
                            } => {
                                dict.set_item("context_ids", context_ids).unwrap();
                                dict.set_item("target_id", target_id).unwrap();
                            }
                            Instance::Sg {
                                center_id,
                                context_ids,
                            } => {
                                dict.set_item("center_id", center_id).unwrap();
                                dict.set_item("context_ids", context_ids).unwrap();
                            }
                        }
                        dict.into()
                    })
                })
                .collect::<Vec<Py<PyDict>>>();
            return Ok(mapped_dicts);
        } else {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "Invalid input type. Must be a sequence or list of sequences.",
            ));
        }
    })
}
