use std::path::Path;
use std::path::PathBuf;

use pyo3::exceptions::PyNotImplementedError;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyIterator};

use anyhow::Result;
use gtars::common::models::{Region, RegionSet};
use gtars::uniwig::utils::{
    read_bw_file, get_max_val_chr_bw
};

// this is for internal use only
pub fn extract_regions_from_py_any(regions: &Bound<'_, PyAny>) -> Result<RegionSet> {
    // is a string?
    if let Ok(regions) = regions.extract::<String>() {
        let regions = Path::new(&regions);

        if !regions.exists() {
            return Err(pyo3::exceptions::PyFileNotFoundError::new_err(format!(
                "The file {} does not exist.",
                regions.display()
            ))
            .into());
        }

        let regions = gtars::common::utils::extract_regions_from_bed_file(regions);
        match regions {
            Ok(regions) => return Ok(RegionSet::from(regions)),
            Err(e) => return Err(pyo3::exceptions::PyValueError::new_err(e.to_string()).into()),
        }
    }

    let regions = PyIterator::from_bound_object(regions)?;

    // attempt to map the list to a vector of regions
    let regions = regions
        .iter()?
        .map(|x| {
            let x = match x {
                Ok(x) => x,
                Err(e) => anyhow::bail!("There was an error iterating over regions: {}", e),
            };

            // extract chr, start, end
            // this lets us interface any python object with chr, start, end attributes
            let chr = x.getattr("chr").unwrap().extract::<String>().unwrap();
            let start = x.getattr("start").unwrap().extract::<u32>().unwrap();
            let end = x.getattr("end").unwrap().extract::<u32>().unwrap();

            Ok(Region { chr, start, end })
        })
        .collect::<Vec<_>>();

    let regions = regions.into_iter().collect::<Result<Vec<_>>>()?;

    Ok(RegionSet::from(regions))
}

#[pyclass(name = "Coverage", module="gtars.utils")]
#[derive(Clone, Debug)]
pub struct PyCoverage {
    bw_file: String
}

#[pymethods]
impl PyCoverage {
    #[new]
    pub fn new(path: String) -> Result<Self> {
        Ok(PyCoverage {
            bw_file: path
        })
    }

    pub fn stats(&self, chr: String, start: usize, end: usize, stat_type: String) -> Result<f64> {
        let bigwig = read_bw_file(&self.bw_file)?;
        match stat_type.as_str() {
            "max" => {
                if let Some(val) = get_max_val_chr_bw(bigwig, &chr, start as u32, end as u32) {
                    Ok(val as f64)
                } else {
                    Err(PyValueError::new_err("Could not get value for {chr}:{start}-{end}. Please check your chrom name and coordinates.").into())
                }
            },
            _ => {
                Err(PyNotImplementedError::new_err("Only max is supported at this time.").into())
            }
        }
    }

    pub fn __repr__(&self) -> String {
        format!(
            "Coverage({})",
            self.bw_file
        )
    }

    pub fn __str__(&self) -> String {
        format!(
            "Coverage({})",
            self.bw_file
        )
    }


}

#[pyfunction]
pub fn write_tokens_to_gtok(filename: &str, tokens: Vec<u32>) -> PyResult<()> {
    gtars::io::write_tokens_to_gtok(filename, &tokens)?;
    Ok(())
}

#[pyfunction]
pub fn read_tokens_from_gtok(filename: &str) -> PyResult<Vec<u32>> {
    let tokens = gtars::io::read_tokens_from_gtok(filename)?;
    Ok(tokens)
}

#[pyfunction]
pub fn read_tokens_from_gtok_as_strings(filename: &str) -> PyResult<Vec<String>> {
    let tokens = gtars::io::read_tokens_from_gtok(filename)?;
    let tokens = tokens.iter().map(|t| t.to_string()).collect();
    Ok(tokens)
}

#[pymodule]
pub fn utils(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(write_tokens_to_gtok))?;
    m.add_wrapped(wrap_pyfunction!(read_tokens_from_gtok))?;
    m.add_wrapped(wrap_pyfunction!(read_tokens_from_gtok_as_strings))?;
    Ok(())
}
