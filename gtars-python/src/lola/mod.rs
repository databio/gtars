use std::path::Path;

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::models::PyRegionSetList;
use gtars_core::models::{Region, RegionSet};
use gtars_igd::igd::Igd;
use gtars_lola::database::{RegionDB, RegionSetAnno};
use gtars_lola::enrichment::run_lola;
use gtars_lola::models::{Direction, LolaConfig, LolaResult};
use gtars_lola::output::{annotate_results, apply_fdr_correction};
use gtars_lola::universe;

// =========================================================================
// PyRegionDB
// =========================================================================

/// A LOLA region database containing indexed genomic region sets.
#[pyclass(name = "RegionDB", module = "gtars.lola")]
pub struct PyRegionDB {
    inner: RegionDB,
}

#[pymethods]
impl PyRegionDB {
    /// Load a RegionDB from a LOLA-format folder.
    #[staticmethod]
    #[pyo3(signature = (db_path, collections=None, limit=None))]
    fn from_folder(
        db_path: &str,
        collections: Option<Vec<String>>,
        limit: Option<usize>,
    ) -> PyResult<Self> {
        let path = Path::new(db_path);
        let coll_refs: Option<Vec<&str>> = collections
            .as_ref()
            .map(|v| v.iter().map(|s| s.as_str()).collect());

        let db = RegionDB::from_lola_folder(path, coll_refs.as_deref(), limit)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to load RegionDB: {}", e)))?;

        Ok(PyRegionDB { inner: db })
    }

    /// Load a RegionDB from a list of BED file paths.
    #[staticmethod]
    #[pyo3(signature = (bed_files, filenames=None))]
    fn from_bed_files(
        bed_files: Vec<String>,
        filenames: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let names: Vec<String> = filenames.unwrap_or_else(|| {
            bed_files
                .iter()
                .map(|p| {
                    Path::new(p)
                        .file_name()
                        .map(|f| f.to_string_lossy().into_owned())
                        .unwrap_or_else(|| p.clone())
                })
                .collect()
        });

        let mut region_sets = Vec::new();
        let mut igd_inputs = Vec::new();
        let mut region_anno = Vec::new();

        for (i, bed_path) in bed_files.iter().enumerate() {
            let path = Path::new(bed_path);
            let rs = RegionSet::try_from(path)
                .map_err(|e| PyRuntimeError::new_err(format!("Failed to read {}: {}", bed_path, e)))?;

            let regions: Vec<(String, i32, i32)> = rs
                .regions
                .iter()
                .map(|r| (r.chr.clone(), r.start as i32, r.end as i32))
                .collect();

            let name = names.get(i).cloned().unwrap_or_default();
            igd_inputs.push((name.clone(), regions));
            region_sets.push(rs);
            region_anno.push(RegionSetAnno {
                filename: name,
                ..Default::default()
            });
        }

        let igd = Igd::from_region_sets(igd_inputs);
        let db = RegionDB::from_igd_with_regions(igd, region_sets, region_anno);

        Ok(PyRegionDB { inner: db })
    }

    /// Number of region sets in this database.
    #[getter]
    fn num_region_sets(&self) -> usize {
        self.inner.num_region_sets()
    }

    /// List region set filenames.
    #[pyo3(signature = (collections=None))]
    fn list_region_sets(&self, collections: Option<Vec<String>>) -> Vec<String> {
        let coll_refs: Option<Vec<&str>> = collections
            .as_ref()
            .map(|v| v.iter().map(|s| s.as_str()).collect());
        self.inner.list_region_sets(coll_refs.as_deref())
    }

    /// Extract region sets by 0-based index as a RegionSetList.
    ///
    /// If indices is None, returns all region sets.
    #[pyo3(signature = (indices=None))]
    fn get_region_sets(&self, indices: Option<Vec<usize>>) -> PyRegionSetList {
        let idx = indices.unwrap_or_else(|| (0..self.inner.num_region_sets()).collect());
        let rsl = self.inner.get_region_set_list(&idx);
        PyRegionSetList::from_inner(rsl)
    }

    /// Get per-file annotations as a list of dicts.
    #[getter]
    fn region_anno<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let mut result = Vec::new();
        for a in &self.inner.region_anno {
            let d = PyDict::new(py);
            d.set_item("filename", &a.filename)?;
            d.set_item("cellType", &a.cell_type)?;
            d.set_item("description", &a.description)?;
            d.set_item("tissue", &a.tissue)?;
            d.set_item("dataSource", &a.data_source)?;
            d.set_item("antibody", &a.antibody)?;
            d.set_item("treatment", &a.treatment)?;
            d.set_item("collection", &a.collection)?;
            result.push(d);
        }
        Ok(result)
    }

    /// Get collection-level annotations as a list of dicts.
    #[getter]
    fn collection_anno<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        let mut result = Vec::new();
        for a in &self.inner.collection_anno {
            let d = PyDict::new(py);
            d.set_item("collectionname", &a.collection_name)?;
            d.set_item("collector", &a.collector)?;
            d.set_item("date", &a.date)?;
            d.set_item("source", &a.source)?;
            d.set_item("description", &a.description)?;
            result.push(d);
        }
        Ok(result)
    }

    fn __repr__(&self) -> String {
        format!(
            "RegionDB({} region sets, {} contigs)",
            self.inner.num_region_sets(),
            self.inner.igd.num_contigs()
        )
    }
}

// =========================================================================
// run_lola
// =========================================================================

/// Run LOLA enrichment analysis.
///
/// Args:
///     user_sets: List of dicts with 'chr', 'start', 'end' keys (each a list),
///                or list of lists of (chr, start, end) tuples.
///     universe: Dict with 'chr', 'start', 'end' keys, or list of tuples.
///     region_db: A RegionDB object.
///     min_overlap: Minimum base-pair overlap (default 1).
///     direction: "enrichment" or "depletion".
///
/// Returns:
///     List of dicts, one per (user_set, db_set) pair.
#[pyfunction(name = "run_lola")]
#[pyo3(signature = (user_sets, universe, region_db, min_overlap=1, direction="enrichment"))]
fn py_run_lola<'py>(
    py: Python<'py>,
    user_sets: Vec<Vec<(String, u32, u32)>>,
    universe: Vec<(String, u32, u32)>,
    region_db: &PyRegionDB,
    min_overlap: i32,
    direction: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let rs_user: Vec<RegionSet> = user_sets
        .into_iter()
        .map(|regions| {
            RegionSet::from(
                regions
                    .into_iter()
                    .map(|(chr, start, end)| Region {
                        chr,
                        start,
                        end,
                        rest: None,
                    })
                    .collect::<Vec<Region>>(),
            )
        })
        .collect();

    let rs_universe = RegionSet::from(
        universe
            .into_iter()
            .map(|(chr, start, end)| Region {
                chr,
                start,
                end,
                rest: None,
            })
            .collect::<Vec<Region>>(),
    );

    let dir = match direction {
        "depletion" | "less" => Direction::Depletion,
        "enrichment" | "greater" => Direction::Enrichment,
        _ => return Err(PyValueError::new_err("direction must be 'enrichment' or 'depletion'")),
    };

    let config = LolaConfig {
        min_overlap,
        direction: dir,
    };

    let mut results = run_lola(&region_db.inner.igd, &rs_user, &rs_universe, &config)
        .map_err(|e| PyRuntimeError::new_err(format!("LOLA error: {}", e)))?;

    annotate_results(&mut results, &region_db.inner);
    apply_fdr_correction(&mut results);

    results_to_dict(py, &results)
}

/// Convert results to a column-oriented Python dict (DataFrame-friendly).
fn results_to_dict<'py>(
    py: Python<'py>,
    results: &[LolaResult],
) -> PyResult<Bound<'py, PyDict>> {
    use gtars_lola::output::results_to_columns;

    let c = results_to_columns(results);
    let dict = PyDict::new(py);
    dict.set_item("userSet", c.user_set)?;
    dict.set_item("dbSet", c.db_set)?;
    dict.set_item("collection", c.collection)?;
    dict.set_item("pValueLog", c.p_value_log)?;
    dict.set_item("oddsRatio", c.odds_ratio)?;
    dict.set_item("support", c.support)?;
    dict.set_item("rnkPV", c.rnk_pv)?;
    dict.set_item("rnkOR", c.rnk_or)?;
    dict.set_item("rnkSup", c.rnk_sup)?;
    dict.set_item("maxRnk", c.max_rnk)?;
    dict.set_item("meanRnk", c.mean_rnk)?;
    dict.set_item("b", c.b)?;
    dict.set_item("c", c.c)?;
    dict.set_item("d", c.d)?;
    dict.set_item("description", c.description)?;
    dict.set_item("cellType", c.cell_type)?;
    dict.set_item("tissue", c.tissue)?;
    dict.set_item("antibody", c.antibody)?;
    dict.set_item("treatment", c.treatment)?;
    dict.set_item("dataSource", c.data_source)?;
    dict.set_item("filename", c.filename)?;
    dict.set_item("qValue", c.q_value)?;
    dict.set_item("size", c.db_set_size)?;
    Ok(dict)
}

// =========================================================================
// Universe functions
// =========================================================================

/// Check universe appropriateness for user sets.
#[pyfunction(name = "check_universe")]
fn py_check_universe<'py>(
    py: Python<'py>,
    user_sets: Vec<Vec<(String, u32, u32)>>,
    universe_regions: Vec<(String, u32, u32)>,
) -> PyResult<Bound<'py, PyDict>> {
    let rs_user: Vec<RegionSet> = user_sets
        .into_iter()
        .map(|regions| tuples_to_regionset(regions))
        .collect();

    let rs_universe = tuples_to_regionset(universe_regions);
    let universe_igd = gtars_igd::igd::Igd::from_single_region_set(&rs_universe);

    let report = universe::check_universe_appropriateness(&rs_user, &universe_igd);

    let dict = PyDict::new(py);
    let mut us_index = Vec::new();
    let mut total = Vec::new();
    let mut in_univ = Vec::new();
    let mut coverage = Vec::new();
    let mut m2m = Vec::new();
    let mut warnings: Vec<String> = Vec::new();

    for ur in &report.user_set_reports {
        us_index.push(ur.user_set_index);
        total.push(ur.total_regions);
        in_univ.push(ur.regions_in_universe);
        coverage.push(ur.coverage);
        m2m.push(ur.many_to_many_count);
        warnings.extend(ur.warnings.clone());
    }

    dict.set_item("userSet", us_index)?;
    dict.set_item("totalRegions", total)?;
    dict.set_item("regionsInUniverse", in_univ)?;
    dict.set_item("coverage", coverage)?;
    dict.set_item("manyToMany", m2m)?;
    dict.set_item("warnings", warnings)?;

    Ok(dict)
}

/// Redefine user sets in terms of universe regions.
#[pyfunction(name = "redefine_user_sets")]
fn py_redefine_user_sets(
    user_sets: Vec<Vec<(String, u32, u32)>>,
    universe_regions: Vec<(String, u32, u32)>,
) -> Vec<Vec<(String, u32, u32)>> {
    let rs_user: Vec<RegionSet> = user_sets
        .into_iter()
        .map(|regions| tuples_to_regionset(regions))
        .collect();

    let rs_universe = tuples_to_regionset(universe_regions);
    let universe_igd = gtars_igd::igd::Igd::from_single_region_set(&rs_universe);
    let redefined = universe::redefine_user_sets(&rs_user, &rs_universe, &universe_igd);

    redefined
        .into_iter()
        .map(|rs| regionset_to_tuples(&rs))
        .collect()
}

/// Build a restricted universe from user sets.
#[pyfunction(name = "build_restricted_universe")]
fn py_build_restricted_universe(
    user_sets: Vec<Vec<(String, u32, u32)>>,
) -> Vec<(String, u32, u32)> {
    let rs_user: Vec<RegionSet> = user_sets
        .into_iter()
        .map(|regions| tuples_to_regionset(regions))
        .collect();

    let restricted = universe::build_restricted_universe(&rs_user);
    regionset_to_tuples(&restricted)
}

// =========================================================================
// Helpers
// =========================================================================

fn tuples_to_regionset(tuples: Vec<(String, u32, u32)>) -> RegionSet {
    RegionSet::from(
        tuples
            .into_iter()
            .map(|(chr, start, end)| Region {
                chr,
                start,
                end,
                rest: None,
            })
            .collect::<Vec<Region>>(),
    )
}

fn regionset_to_tuples(rs: &RegionSet) -> Vec<(String, u32, u32)> {
    rs.regions
        .iter()
        .map(|r| (r.chr.clone(), r.start, r.end))
        .collect()
}

// =========================================================================
// Module registration
// =========================================================================

#[pymodule]
pub fn lola(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRegionDB>()?;
    m.add_function(wrap_pyfunction!(py_run_lola, m)?)?;
    m.add_function(wrap_pyfunction!(py_check_universe, m)?)?;
    m.add_function(wrap_pyfunction!(py_redefine_user_sets, m)?)?;
    m.add_function(wrap_pyfunction!(py_build_restricted_universe, m)?)?;
    Ok(())
}
