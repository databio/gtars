use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::{calc_summary_signal, SignalMatrix};
use serde::Serialize;
use wasm_bindgen::prelude::*;

use crate::regionset::JsRegionSet;

// ── JsSignalMatrix ──────────────────────────────────────────────────

#[wasm_bindgen(js_name = "SignalMatrix")]
pub struct JsSignalMatrix {
    matrix: SignalMatrix,
}

#[wasm_bindgen(js_class = "SignalMatrix")]
impl JsSignalMatrix {
    /// Construct a SignalMatrix from JS data.
    ///
    /// - `regionIds`: array of strings in `"chr_start_end"` format
    /// - `conditionNames`: array of condition/cell-type names
    /// - `values`: flat row-major array of signal values
    /// - `nRegions`: number of rows (regions)
    /// - `nConditions`: number of columns (conditions)
    #[wasm_bindgen(constructor)]
    pub fn new(
        region_ids: &JsValue,
        condition_names: &JsValue,
        values: &[f64],
        n_regions: u32,
        n_conditions: u32,
    ) -> Result<JsSignalMatrix, JsValue> {
        let ids: Vec<String> = serde_wasm_bindgen::from_value(region_ids.clone())?;
        let cond_names: Vec<String> = serde_wasm_bindgen::from_value(condition_names.clone())?;

        let n_r = n_regions as usize;
        let n_c = n_conditions as usize;

        if ids.len() != n_r {
            return Err(JsValue::from_str(&format!(
                "regionIds length ({}) != nRegions ({})",
                ids.len(),
                n_r
            )));
        }
        if cond_names.len() != n_c {
            return Err(JsValue::from_str(&format!(
                "conditionNames length ({}) != nConditions ({})",
                cond_names.len(),
                n_c
            )));
        }
        if values.len() != n_r * n_c {
            return Err(JsValue::from_str(&format!(
                "values length ({}) != nRegions * nConditions ({})",
                values.len(),
                n_r * n_c
            )));
        }

        // Parse region IDs into regions
        let mut regions = Vec::with_capacity(n_r);
        for id in &ids {
            let parts: Vec<&str> = id.split('_').collect();
            if parts.len() != 3 {
                return Err(JsValue::from_str(&format!(
                    "Invalid region ID format '{}': expected 'chr_start_end'",
                    id
                )));
            }
            let chr = parts[0].to_string();
            let start: u32 = parts[1]
                .parse()
                .map_err(|_| JsValue::from_str(&format!("Invalid start in '{}'", id)))?;
            let end: u32 = parts[2]
                .parse()
                .map_err(|_| JsValue::from_str(&format!("Invalid end in '{}'", id)))?;
            regions.push(Region {
                chr,
                start,
                end,
                rest: None,
            });
        }

        // Reshape flat values into row-major Vec<Vec<f64>>
        let row_values: Vec<Vec<f64>> = (0..n_r)
            .map(|i| values[i * n_c..(i + 1) * n_c].to_vec())
            .collect();

        Ok(JsSignalMatrix {
            matrix: SignalMatrix {
                regions: RegionSet::from(regions),
                condition_names: cond_names,
                values: row_values,
            },
        })
    }
}

// ── Calc function ───────────────────────────────────────────────────

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct SignalSummaryResultJs {
    signal_matrix: Vec<SignalRowJs>,
    matrix_stats: Vec<ConditionStatsJs>,
    condition_names: Vec<String>,
}

#[derive(Serialize)]
struct SignalRowJs {
    region: String,
    values: Vec<f64>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ConditionStatsJs {
    condition: String,
    lower_whisker: f64,
    lower_hinge: f64,
    median: f64,
    upper_hinge: f64,
    upper_whisker: f64,
}

/// Compute summary signal: overlap query regions with a signal matrix,
/// take MAX per condition, and compute boxplot statistics.
#[wasm_bindgen(js_name = "calcSummarySignal")]
pub fn js_calc_summary_signal(
    region_set: &JsRegionSet,
    signal_matrix: &JsSignalMatrix,
) -> Result<JsValue, JsValue> {
    let result = calc_summary_signal(&region_set.region_set, &signal_matrix.matrix)
        .map_err(|e| JsValue::from_str(&e.to_string()))?;

    let js_result = SignalSummaryResultJs {
        signal_matrix: result
            .signal_matrix
            .into_iter()
            .map(|(region, values)| SignalRowJs { region, values })
            .collect(),
        matrix_stats: result
            .matrix_stats
            .into_iter()
            .map(|s| ConditionStatsJs {
                condition: s.condition,
                lower_whisker: s.lower_whisker,
                lower_hinge: s.lower_hinge,
                median: s.median,
                upper_hinge: s.upper_hinge,
                upper_whisker: s.upper_whisker,
            })
            .collect(),
        condition_names: result.condition_names,
    };

    serde_wasm_bindgen::to_value(&js_result).map_err(|e| e.into())
}
