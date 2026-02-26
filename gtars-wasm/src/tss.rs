use gtars_genomicdist::models::TssIndex;
use wasm_bindgen::prelude::*;

use crate::regionset::JsRegionSet;

// ── JsTssIndex ──────────────────────────────────────────────────────

#[wasm_bindgen(js_name = "TssIndex")]
pub struct JsTssIndex {
    index: TssIndex,
}

#[wasm_bindgen(js_class = "TssIndex")]
impl JsTssIndex {
    /// Build a TssIndex from a RegionSet of features (e.g. TSS sites).
    ///
    /// Clones the inner RegionSet and computes sorted midpoints per chromosome.
    #[wasm_bindgen(constructor)]
    pub fn new(features: &JsRegionSet) -> Result<JsTssIndex, JsValue> {
        let index = TssIndex::try_from(features.region_set.clone())
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        Ok(JsTssIndex { index })
    }

    /// Calculate unsigned distance from each query region to its nearest feature midpoint.
    ///
    /// Returns an array where `null` indicates no features on that chromosome
    /// (sentinel `u32::MAX` in Rust → `null` in JS).
    #[wasm_bindgen(js_name = "calcTssDistances")]
    pub fn calc_tss_distances(&self, query: &JsRegionSet) -> Result<JsValue, JsValue> {
        let distances = self
            .index
            .calc_tss_distances(&query.region_set)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;

        // Convert u32::MAX sentinel to None (JS null)
        let nullable: Vec<Option<f64>> = distances
            .into_iter()
            .map(|d| {
                if d == u32::MAX {
                    None
                } else {
                    Some(d as f64)
                }
            })
            .collect();

        serde_wasm_bindgen::to_value(&nullable).map_err(|e| e.into())
    }

    /// Calculate signed distance from each query region to its nearest feature.
    ///
    /// Positive = feature is downstream, negative = feature is upstream.
    /// Returns an array where `null` indicates no features on that chromosome
    /// (sentinel `i64::MAX` in Rust → `null` in JS).
    #[wasm_bindgen(js_name = "calcFeatureDistances")]
    pub fn calc_feature_distances(&self, query: &JsRegionSet) -> Result<JsValue, JsValue> {
        let distances = self
            .index
            .calc_feature_distances(&query.region_set)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;

        // Convert i64::MAX sentinel to None (JS null)
        let nullable: Vec<Option<f64>> = distances
            .into_iter()
            .map(|d| {
                if d == i64::MAX {
                    None
                } else {
                    Some(d as f64)
                }
            })
            .collect();

        serde_wasm_bindgen::to_value(&nullable).map_err(|e| e.into())
    }
}
