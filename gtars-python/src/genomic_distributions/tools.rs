use pyo3::prelude::*;

use crate::models::PyRegionSet;
use gtars_genomicdist::statistics;
use std::collections::HashMap;

use crate::models::PyGenomeAssembly;

#[pyfunction(name = "calc_gc_content")]
#[pyo3(signature = (rs, genome, ignore_unk_chroms = Some(false)))]
pub fn py_calc_gc_content(
    rs: &PyRegionSet,
    genome: &PyGenomeAssembly,
    ignore_unk_chroms: Option<bool>,
) -> anyhow::Result<Vec<f64>> {
    let result = statistics::calc_gc_content(
        &rs.regionset,
        &genome.genome_assembly,
        ignore_unk_chroms.unwrap_or(false),
    )?;
    Ok(result)
}

#[pyfunction(name = "calc_dinucleotide_frequency")]
pub fn py_calc_dinucleotide_frequency(
    rs: &PyRegionSet,
    genome: &PyGenomeAssembly,
) -> anyhow::Result<HashMap<String, u64>> {
    println!("Calculating dinucleotide_frequency...");
    let frequencies = statistics::calc_dinucl_freq(&rs.regionset, &genome.genome_assembly)?;
    let mut freq_map: HashMap<String, u64> = HashMap::new();
    // Convert Dinucleotide to String and push to HashMap
    for (di, freq) in frequencies {
        freq_map.insert(di.to_string()?, freq);
    }

    Ok(freq_map)
}
