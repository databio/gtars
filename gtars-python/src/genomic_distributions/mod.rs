use pyo3::prelude::*;

mod tools;
pub use self::tools::{py_calc_dinucleotide_frequency, py_calc_gc_content};

#[pymodule]
pub fn genomic_distributions(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_calc_gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(py_calc_dinucleotide_frequency, m)?)?;
    Ok(())
}
