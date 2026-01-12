use pyo3::prelude::*;

use crate::models::PyRegionSet;
use gtars_genomicdist::statistics;
use std::collections::HashMap;
use std::path::Path;

use crate::models::PyGenomeAssembly;

#[pyfunction(name = "calc_gc_content")]
pub fn py_calc_gc_content(
    rs: &PyRegionSet,
    genome: &PyGenomeAssembly,
    ignore_unk_chroms: Option<bool>,
) -> anyhow::Result<Vec<f64>> {
    statistics::calc_gc_content(
        &rs.regionset,
        &genome.genome_assembly,
        ignore_unk_chroms.unwrap_or(false),
    )
}

#[pyfunction(name = "calc_dincleotide_frequency")]
pub fn py_calc_dinucleotide_frequency(
    rs: &PyRegionSet,
    genome: &PyGenomeAssembly,
) -> anyhow::Result<HashMap<String, f64>> {
    let frequencies = statistics::calc_dinucl_freq(&rs.regionset, &genome.genome_assembly)?;

    let mut freq_map: HashMap<String, f64> = HashMap::new();

    // Convert Dinucleotide to String and push to HashMap
    for (di, freq) in frequencies {
        freq_map.insert(di.to_string()?, freq);
    }

    Ok(freq_map)
}
