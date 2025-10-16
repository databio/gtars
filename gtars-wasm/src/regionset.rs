use wasm_bindgen::prelude::*;
use gtars_core::models::{Region, RegionSet};
use crate::models::BedEntries;

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
        Ok(JsRegionSet { region_set: RegionSet::from(regions) })
    }


    #[wasm_bindgen(getter, js_name = "number_of_regions")]
    pub fn get_region(&self) -> i32 {
        self.region_set.len() as i32

    }

    #[wasm_bindgen(getter, js_name = "get_phi")]
    pub fn get_phi(&mut self) -> String {
        String::from("phi-neaaas")
    }
    #[wasm_bindgen(getter, js_name = "mean_region_width")]
    pub fn get_mean_region_width(&self) -> f64 {
        self.region_set.mean_region_width()
    }

    #[wasm_bindgen(getter, js_name = "total_nucleotides")]
    pub fn get_total_nucleotides(&self) -> i32 {
        self.region_set.nucleotides_length() as i32
    }
}
