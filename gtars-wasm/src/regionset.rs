use crate::models::BedEntries;
use gtars_core::models::{Region, RegionSet};
use std::collections::HashMap;
use wasm_bindgen::prelude::*;

#[wasm_bindgen(js_name = "ChromosomeStats")]
#[derive(serde::Serialize)]
pub struct JsChromosomeStats {
    chromosome: String,
    count: u32,
    minimum: u32,
    maximum: u32,
    mean: f64,
    median: f64,
    start: u32,
    end: u32,
}

#[wasm_bindgen(js_class = "ChromosomeStats")]
impl JsChromosomeStats {
    #[wasm_bindgen(getter)]
    pub fn chromosome(&self) -> String {
        self.chromosome.clone()
    }

    #[wasm_bindgen(getter)]
    pub fn count(&self) -> u32 {
        self.count
    }

    #[wasm_bindgen(getter)]
    pub fn minimum(&self) -> u32 {
        self.minimum
    }

    #[wasm_bindgen(getter)]
    pub fn maximum(&self) -> u32 {
        self.maximum
    }

    #[wasm_bindgen(getter)]
    pub fn mean(&self) -> f64 {
        self.mean
    }

    #[wasm_bindgen(getter)]
    pub fn median(&self) -> f64 {
        self.median
    }
}

#[wasm_bindgen(js_name = "RegionSet")]
pub struct JsRegionSet {
    region_set: RegionSet,
}

#[wasm_bindgen(js_class = "RegionSet")]
impl JsRegionSet {
    // #[wasm_bindgen(constructor)]
    // pub fn new(word: &str) -> Result<JsRegionSet, JsValue> {
    //     println!("Phineas is ready");
    //     Ok(JsRegionSet {
    //         internal: String::from(word),
    //     })
    // }
    #[wasm_bindgen(constructor)]
    pub fn new(regions: &JsValue) -> Result<JsRegionSet, JsValue> {
        let regions: BedEntries = serde_wasm_bindgen::from_value(regions.clone())?;
        let regions: Vec<Region> = regions
            .0
            .into_iter()
            .map(|be| Region {
                chr: be.0,
                start: be.1,
                end: be.2,
                rest: None,
            })
            .collect::<Vec<Region>>();
        let mut region_set: RegionSet = RegionSet::from(regions);
        region_set.sort();
        Ok(JsRegionSet { region_set })
    }

    #[wasm_bindgen(getter, js_name = "number_of_regions")]
    pub fn get_region(&self) -> i32 {
        self.region_set.len() as i32
    }

    #[wasm_bindgen(getter, js_name = "mean_region_width")]
    pub fn get_mean_region_width(&self) -> f64 {
        self.region_set.mean_region_width()
    }

    #[wasm_bindgen(getter, js_name = "total_nucleotides")]
    pub fn get_total_nucleotides(&self) -> i32 {
        self.region_set.nucleotides_length() as i32
    }

    #[wasm_bindgen(getter, js_name = "digest")]
    pub fn get_digest(&self) -> String {
        self.region_set.identifier()
    }

    #[wasm_bindgen(js_name = "calculate_statistics")]
    pub fn calculate_statistics(&self) -> Result<JsValue, JsValue> {
        let stats = self.region_set.calculate_statistics();
        let mut result: HashMap<String, JsChromosomeStats> = HashMap::new();

        for (key, value) in stats {
            result.insert(
                key.clone(),
                JsChromosomeStats {
                    chromosome: value.chromosome.clone(),
                    count: value.count,
                    minimum: value.minimum,
                    maximum: value.maximum,
                    mean: value.mean,
                    median: value.median,
                    start: value.start,
                    end: value.end,
                },
            );
        }

        serde_wasm_bindgen::to_value(&result).map_err(|e| e.into())
    }
}
