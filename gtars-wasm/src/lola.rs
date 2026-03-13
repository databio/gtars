use wasm_bindgen::prelude::*;

use crate::regionset::JsRegionSetList;
use gtars_core::models::{Region, RegionSet};
use gtars_igd::igd::Igd;
use gtars_lola::database::{RegionDB, RegionSetAnno};
use gtars_lola::enrichment::run_lola;
use gtars_lola::models::{Direction, LolaConfig, LolaResult};
use gtars_lola::output::apply_fdr_correction;
use gtars_lola::universe;

// =========================================================================
// JsLolaRegionDB
// =========================================================================

/// A LOLA region database for enrichment testing.
///
/// Built in-memory from region data — no filesystem needed.
#[wasm_bindgen(js_name = "LolaRegionDB")]
pub struct JsLolaRegionDB {
    inner: RegionDB,
}

#[wasm_bindgen(js_class = "LolaRegionDB")]
impl JsLolaRegionDB {
    /// Create a LolaRegionDB from an array of {regions, name} objects.
    ///
    /// @param entries - Array of { regions: [["chr1", 100, 200], ...], name: "filename.bed" }
    #[wasm_bindgen(constructor)]
    pub fn new(entries: &JsValue) -> Result<JsLolaRegionDB, JsValue> {
        let entries: Vec<DbEntry> =
            serde_wasm_bindgen::from_value(entries.clone()).map_err(|e| {
                JsValue::from_str(&format!(
                    "Expected array of {{regions, name}} objects: {}",
                    e
                ))
            })?;

        let mut region_sets = Vec::with_capacity(entries.len());
        let mut igd_inputs = Vec::with_capacity(entries.len());
        let mut region_anno = Vec::with_capacity(entries.len());

        for entry in entries {
            let regions: Vec<Region> = entry
                .regions
                .into_iter()
                .map(|(chr, start, end)| Region {
                    chr,
                    start,
                    end,
                    rest: None,
                })
                .collect();

            let igd_regions: Vec<(String, i32, i32)> = regions
                .iter()
                .map(|r| (r.chr.clone(), r.start as i32, r.end as i32))
                .collect();

            igd_inputs.push((entry.name.clone(), igd_regions));
            region_sets.push(RegionSet::from(regions));
            region_anno.push(RegionSetAnno {
                filename: entry.name,
                ..Default::default()
            });
        }

        let igd = Igd::from_region_sets(igd_inputs);
        let db = RegionDB::from_igd_with_regions(igd, region_sets, region_anno);
        Ok(JsLolaRegionDB { inner: db })
    }

    /// Number of region sets in this database.
    #[wasm_bindgen(getter, js_name = "numRegionSets")]
    pub fn num_region_sets(&self) -> usize {
        self.inner.num_region_sets()
    }

    /// List region set filenames.
    #[wasm_bindgen(js_name = "listRegionSets")]
    pub fn list_region_sets(&self) -> Vec<String> {
        self.inner.list_region_sets(None)
    }

    /// Extract region sets by 0-based indices as a RegionSetList.
    ///
    /// @param indices - Optional array of 0-based indices. If omitted, returns all.
    #[wasm_bindgen(js_name = "getRegionSets")]
    pub fn get_region_sets(&self, indices: Option<Vec<usize>>) -> JsRegionSetList {
        let idx = indices.unwrap_or_else(|| (0..self.inner.num_region_sets()).collect());
        let rsl = self.inner.get_region_set_list(&idx);
        JsRegionSetList { inner: rsl }
    }
}

// =========================================================================
// runLOLA
// =========================================================================

/// Run LOLA enrichment analysis.
///
/// @param userSets - Array of arrays of [chr, start, end] tuples
/// @param universe - Array of [chr, start, end] tuples
/// @param regionDb - A LolaRegionDB object
/// @param minOverlap - Minimum bp overlap (default 1)
/// @param direction - "enrichment" or "depletion" (default "enrichment")
/// @returns Object with column arrays: userSet, dbSet, pValueLog, oddsRatio, support, etc.
#[wasm_bindgen(js_name = "runLOLA")]
pub fn js_run_lola(
    user_sets: &JsValue,
    universe_regions: &JsValue,
    region_db: &JsLolaRegionDB,
    min_overlap: Option<i32>,
    direction: Option<String>,
) -> Result<JsValue, JsValue> {
    let user_sets_data: Vec<Vec<(String, u32, u32)>> =
        serde_wasm_bindgen::from_value(user_sets.clone()).map_err(|e| {
            JsValue::from_str(&format!(
                "userSets: expected array of arrays of [chr, start, end]: {}",
                e
            ))
        })?;

    let universe_data: Vec<(String, u32, u32)> =
        serde_wasm_bindgen::from_value(universe_regions.clone()).map_err(|e| {
            JsValue::from_str(&format!(
                "universe: expected array of [chr, start, end]: {}",
                e
            ))
        })?;

    let rs_user: Vec<RegionSet> = user_sets_data
        .into_iter()
        .map(|regions| tuples_to_regionset(regions))
        .collect();

    let rs_universe = tuples_to_regionset(universe_data);

    let dir = match direction.as_deref() {
        Some("depletion") | Some("less") => Direction::Depletion,
        _ => Direction::Enrichment,
    };

    let config = LolaConfig {
        min_overlap: min_overlap.unwrap_or(1),
        direction: dir,
    };

    let mut results = run_lola(&region_db.inner.igd, &rs_user, &rs_universe, &config)
        .map_err(|e| JsValue::from_str(&format!("LOLA error: {}", e)))?;

    apply_fdr_correction(&mut results);

    results_to_js(&results)
}

// =========================================================================
// checkUniverseAppropriateness
// =========================================================================

/// Check whether the universe is appropriate for the given user sets.
///
/// @param userSets - Array of arrays of [chr, start, end] tuples
/// @param universe - Array of [chr, start, end] tuples
/// @returns Object with column arrays: userSet, totalRegions, regionsInUniverse, coverage, manyToMany, warnings
#[wasm_bindgen(js_name = "checkUniverseAppropriateness")]
pub fn js_check_universe(
    user_sets: &JsValue,
    universe_regions: &JsValue,
) -> Result<JsValue, JsValue> {
    let user_sets_data: Vec<Vec<(String, u32, u32)>> =
        serde_wasm_bindgen::from_value(user_sets.clone()).map_err(|e| {
            JsValue::from_str(&format!("userSets: {}", e))
        })?;

    let universe_data: Vec<(String, u32, u32)> =
        serde_wasm_bindgen::from_value(universe_regions.clone()).map_err(|e| {
            JsValue::from_str(&format!("universe: {}", e))
        })?;

    let rs_user: Vec<RegionSet> = user_sets_data
        .into_iter()
        .map(|regions| tuples_to_regionset(regions))
        .collect();

    let rs_universe = tuples_to_regionset(universe_data);
    let universe_igd = gtars_igd::igd::Igd::from_single_region_set(&rs_universe);

    let report = universe::check_universe_appropriateness(&rs_user, &universe_igd);

    let result = UniverseCheckResult {
        user_set: report.user_set_reports.iter().map(|r| r.user_set_index).collect(),
        total_regions: report.user_set_reports.iter().map(|r| r.total_regions).collect(),
        regions_in_universe: report.user_set_reports.iter().map(|r| r.regions_in_universe).collect(),
        coverage: report.user_set_reports.iter().map(|r| r.coverage).collect(),
        many_to_many: report.user_set_reports.iter().map(|r| r.many_to_many_count).collect(),
        warnings: report
            .user_set_reports
            .iter()
            .flat_map(|r| r.warnings.clone())
            .collect(),
    };

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&format!("{}", e)))
}

// =========================================================================
// Helpers
// =========================================================================

#[derive(serde::Deserialize)]
struct DbEntry {
    regions: Vec<(String, u32, u32)>,
    name: String,
}

#[derive(serde::Serialize)]
#[serde(rename_all = "camelCase")]
struct UniverseCheckResult {
    user_set: Vec<usize>,
    total_regions: Vec<usize>,
    regions_in_universe: Vec<usize>,
    coverage: Vec<f64>,
    many_to_many: Vec<usize>,
    warnings: Vec<String>,
}

#[derive(serde::Serialize)]
#[serde(rename_all = "camelCase")]
struct LolaResults {
    user_set: Vec<usize>,
    db_set: Vec<usize>,
    p_value_log: Vec<f64>,
    odds_ratio: Vec<f64>,
    support: Vec<u64>,
    rnk_pv: Vec<usize>,
    rnk_or: Vec<usize>,
    rnk_sup: Vec<usize>,
    max_rnk: Vec<usize>,
    mean_rnk: Vec<f64>,
    b: Vec<u64>,
    c: Vec<u64>,
    d: Vec<u64>,
    filename: Vec<String>,
    q_value: Vec<Option<f64>>,
}

fn results_to_js(results: &[LolaResult]) -> Result<JsValue, JsValue> {
    let out = LolaResults {
        user_set: results.iter().map(|r| r.user_set).collect(),
        db_set: results.iter().map(|r| r.db_set).collect(),
        p_value_log: results.iter().map(|r| r.p_value_log).collect(),
        odds_ratio: results.iter().map(|r| r.odds_ratio).collect(),
        support: results.iter().map(|r| r.support).collect(),
        rnk_pv: results.iter().map(|r| r.rnk_pv).collect(),
        rnk_or: results.iter().map(|r| r.rnk_or).collect(),
        rnk_sup: results.iter().map(|r| r.rnk_sup).collect(),
        max_rnk: results.iter().map(|r| r.max_rnk).collect(),
        mean_rnk: results.iter().map(|r| r.mean_rnk).collect(),
        b: results.iter().map(|r| r.b).collect(),
        c: results.iter().map(|r| r.c).collect(),
        d: results.iter().map(|r| r.d).collect(),
        filename: results.iter().map(|r| r.filename.clone()).collect(),
        q_value: results.iter().map(|r| r.q_value).collect(),
    };

    serde_wasm_bindgen::to_value(&out).map_err(|e| JsValue::from_str(&format!("{}", e)))
}

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
