use std::collections::HashMap;

use gtars_core::models::Region;
use wasm_bindgen::prelude::*;
use crate::models::BedEntriesFull;


#[wasm_bindgen(js_name = "RegionSet")]
pub struct RegionSet {
    internal: String,
}

#[wasm_bindgen(js_class = "RegionSet")]
impl RegionSet {

    #[wasm_bindgen(constructor)]
    pub fn new(word: &str) -> Result<RegionSet, JsValue> {
        println!("Phineas is ready");
        Ok(RegionSet {internal: String::from(word)})
    }

    #[wasm_bindgen(getter, js_name = "get_region")]
    pub fn get_region(&mut self) -> String{
        return self.internal.clone();
    }

    #[wasm_bindgen(getter, js_name = "get_phi")]
    pub fn get_phi(&mut self) -> String{
        return String::from("phi-neaaas");
    }
}