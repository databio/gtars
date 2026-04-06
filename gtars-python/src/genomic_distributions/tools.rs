use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::models::{PyBinaryGenomeAssembly, PyGenomeAssembly, PyPartitionList, PyRegionSet, PySignalMatrix};
use gtars_genomicdist::{self, statistics};
use std::collections::HashMap;

#[pyfunction(name = "calc_gc_content")]
#[pyo3(signature = (rs, genome, ignore_unk_chroms = Some(false)))]
pub fn py_calc_gc_content(
    rs: &PyRegionSet,
    genome: &Bound<'_, PyAny>,
    ignore_unk_chroms: Option<bool>,
) -> anyhow::Result<Vec<f64>> {
    let ignore = ignore_unk_chroms.unwrap_or(false);
    if let Ok(g) = genome.cast::<PyBinaryGenomeAssembly>() {
        let g = g.borrow();
        Ok(statistics::calc_gc_content(&rs.regionset, &g.assembly, ignore)?)
    } else if let Ok(g) = genome.cast::<PyGenomeAssembly>() {
        let g = g.borrow();
        Ok(statistics::calc_gc_content(&rs.regionset, &g.genome_assembly, ignore)?)
    } else {
        Err(anyhow::anyhow!("genome must be a GenomeAssembly or BinaryGenomeAssembly"))
    }
}

/// Per-region dinucleotide frequencies, matching R GenomicDistributions `calcDinuclFreq`.
///
/// Arguments:
///   - ``rs``: RegionSet
///   - ``genome``: GenomeAssembly (reference)
///   - ``raw_counts``: if False (default), return percentages (0–100) per row;
///     if True, return raw integer-valued counts
///
/// Returns a dict:
///   - ``region_labels``: list of ``chr_start_end`` strings (one per region)
///   - ``dinucleotides``: list of 16 dinucleotide names in canonical order
///   - ``frequencies``: list of lists — one row per region, 16 values per row
///     matching ``dinucleotides`` order
///
/// For pooled global counts, sum each column of the raw-counts matrix.
#[pyfunction(name = "calc_dinucl_freq")]
#[pyo3(signature = (rs, genome, raw_counts = false, ignore_unk_chroms = false))]
pub fn py_calc_dinucl_freq<'py>(
    py: Python<'py>,
    rs: &PyRegionSet,
    genome: &Bound<'py, PyAny>,
    raw_counts: bool,
    ignore_unk_chroms: bool,
) -> anyhow::Result<Bound<'py, PyDict>> {
    let (labels, matrix) = if let Ok(g) = genome.cast::<PyBinaryGenomeAssembly>() {
        let g = g.borrow();
        statistics::calc_dinucl_freq(&rs.regionset, &g.assembly, raw_counts, ignore_unk_chroms)?
    } else if let Ok(g) = genome.cast::<PyGenomeAssembly>() {
        let g = g.borrow();
        statistics::calc_dinucl_freq(&rs.regionset, &g.genome_assembly, raw_counts, ignore_unk_chroms)?
    } else {
        return Err(anyhow::anyhow!("genome must be a GenomeAssembly or BinaryGenomeAssembly"));
    };
    let dinucl_names: Vec<String> = statistics::DINUCL_ORDER
        .iter()
        .map(|d| d.to_string().unwrap_or_default())
        .collect();
    let freqs_nested: Vec<Vec<f64>> = matrix.into_iter().map(|row| row.to_vec()).collect();
    let result = PyDict::new(py);
    result.set_item("region_labels", labels)?;
    result.set_item("dinucleotides", dinucl_names)?;
    result.set_item("frequencies", freqs_nested)?;
    Ok(result)
}

#[pyfunction(name = "calc_partitions")]
#[pyo3(signature = (rs, partition_list, bp_proportion = false))]
pub fn py_calc_partitions<'py>(
    py: Python<'py>,
    rs: &PyRegionSet,
    partition_list: &PyPartitionList,
    bp_proportion: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let result =
        gtars_genomicdist::calc_partitions(&rs.regionset, &partition_list.partition_list, bp_proportion);
    let dict = PyDict::new(py);
    let partitions: Vec<String> = result.counts.iter().map(|(name, _)| name.clone()).collect();
    let counts: Vec<u32> = result.counts.iter().map(|(_, count)| *count).collect();
    dict.set_item("partition", partitions)?;
    dict.set_item("count", counts)?;
    dict.set_item("total", result.total)?;
    Ok(dict)
}

#[pyfunction(name = "calc_expected_partitions")]
#[pyo3(signature = (rs, partition_list, chrom_sizes, bp_proportion = false))]
pub fn py_calc_expected_partitions<'py>(
    py: Python<'py>,
    rs: &PyRegionSet,
    partition_list: &PyPartitionList,
    chrom_sizes: HashMap<String, u32>,
    bp_proportion: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let result = gtars_genomicdist::calc_expected_partitions(
        &rs.regionset,
        &partition_list.partition_list,
        &chrom_sizes,
        bp_proportion,
    );
    let dict = PyDict::new(py);
    let partitions: Vec<String> = result.rows.iter().map(|r| r.partition.clone()).collect();
    let observed: Vec<f64> = result.rows.iter().map(|r| r.observed).collect();
    let expected: Vec<f64> = result.rows.iter().map(|r| r.expected).collect();
    let log10oe: Vec<f64> = result.rows.iter().map(|r| r.log10_oe).collect();
    let pvalue: Vec<f64> = result.rows.iter().map(|r| r.chi_sq_pval).collect();
    dict.set_item("partition", partitions)?;
    dict.set_item("observed", observed)?;
    dict.set_item("expected", expected)?;
    dict.set_item("log10OE", log10oe)?;
    dict.set_item("pvalue", pvalue)?;
    Ok(dict)
}

#[pyfunction(name = "calc_summary_signal")]
pub fn py_calc_summary_signal<'py>(
    py: Python<'py>,
    rs: &PyRegionSet,
    signal_matrix: &PySignalMatrix,
) -> anyhow::Result<Bound<'py, PyDict>> {
    let result =
        gtars_genomicdist::calc_summary_signal(&rs.regionset, &signal_matrix.signal_matrix, gtars_genomicdist::CoordinateMode::Bed)?;
    let dict = PyDict::new(py);
    dict.set_item("condition_names", &result.condition_names)?;
    let region_labels: Vec<&str> = result.signal_matrix.iter().map(|(label, _)| label.as_str()).collect();
    dict.set_item("region_labels", region_labels)?;

    // Signal matrix as list of lists (each inner list = per-condition max values for one region)
    let matrix: Vec<Vec<f64>> = result
        .signal_matrix
        .iter()
        .map(|(_, vals)| vals.clone())
        .collect();
    dict.set_item("signal_matrix", matrix)?;

    // Matrix stats as list of dicts
    let mut stats = Vec::with_capacity(result.matrix_stats.len());
    for cs in &result.matrix_stats {
        let s = PyDict::new(py);
        s.set_item("condition", &cs.condition)?;
        s.set_item("lower_whisker", cs.lower_whisker)?;
        s.set_item("lower_hinge", cs.lower_hinge)?;
        s.set_item("median", cs.median)?;
        s.set_item("upper_hinge", cs.upper_hinge)?;
        s.set_item("upper_whisker", cs.upper_whisker)?;
        stats.push(s);
    }
    dict.set_item("matrix_stats", stats)?;
    Ok(dict)
}

#[pyfunction(name = "median_abs_distance")]
pub fn py_median_abs_distance(distances: Vec<f64>) -> Option<f64> {
    let as_i64: Vec<i64> = distances
        .iter()
        .map(|&d| {
            if d.is_nan() || d.is_infinite() {
                i64::MAX
            } else {
                d as i64
            }
        })
        .collect();
    gtars_genomicdist::median_abs_distance(&as_i64)
}

#[pyfunction(name = "consensus")]
pub fn py_consensus<'py>(
    py: Python<'py>,
    region_sets: Vec<PyRef<'_, PyRegionSet>>,
) -> PyResult<Vec<Bound<'py, PyDict>>> {
    let sets: Vec<_> = region_sets.iter().map(|rs| rs.regionset.clone()).collect();
    let result = gtars_genomicdist::consensus(&sets);
    let mut out = Vec::with_capacity(result.len());
    for cr in result {
        let dict = PyDict::new(py);
        dict.set_item("chr", &cr.chr)?;
        dict.set_item("start", cr.start)?;
        dict.set_item("end", cr.end)?;
        dict.set_item("count", cr.count)?;
        out.push(dict);
    }
    Ok(out)
}
