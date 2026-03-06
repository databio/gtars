use pyo3::prelude::*;

mod tools;
pub use self::tools::{
    py_calc_dinucleotide_frequency, py_calc_expected_partitions, py_calc_gc_content,
    py_calc_partitions, py_calc_summary_signal, py_consensus, py_median_abs_distance,
};

#[pymodule]
pub fn genomic_distributions(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_calc_gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(py_calc_dinucleotide_frequency, m)?)?;
    m.add_function(wrap_pyfunction!(py_calc_partitions, m)?)?;
    m.add_function(wrap_pyfunction!(py_calc_expected_partitions, m)?)?;
    m.add_function(wrap_pyfunction!(py_calc_summary_signal, m)?)?;
    m.add_function(wrap_pyfunction!(py_median_abs_distance, m)?)?;
    m.add_function(wrap_pyfunction!(py_consensus, m)?)?;
    Ok(())
}
