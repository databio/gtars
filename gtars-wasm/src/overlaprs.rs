use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use gtars_core::models::Interval;
use gtars_overlaprs::{AiList, Bits, Overlapper};

#[derive(Serialize, Deserialize)]
pub struct BedEntries(pub Vec<(String, u32, u32)>);

#[wasm_bindgen(js_name = "Overlapper")]
pub struct JsOverlapper {
    internal: HashMap<String, Box<dyn Overlapper<u32, u32>>>,
    backend: String,
}

#[wasm_bindgen]
impl JsOverlapper {
    #[wasm_bindgen(constructor)]
    pub fn new(universe: &JsValue, backend: String) -> Result<JsOverlapper, JsValue> {
        let mut chr_to_region_store: HashMap<String, Vec<Interval<u32, u32>>> = HashMap::new();
        let universe_regions: BedEntries = serde_wasm_bindgen::from_value(universe.to_owned())?;
        for region in universe_regions.0 {
            let chrom = region.0;
            let start = region.1;
            let end = region.2;

            let iv_store = chr_to_region_store.entry(chrom).or_default();
            iv_store.push(Interval { start, end, val: 0 })
        }
        let mut core: HashMap<String, Box<dyn Overlapper<u32, u32>>> = HashMap::new();

        match backend.as_str() {
            "ailist" => {
                for (chr, iv_list) in chr_to_region_store {
                    core.insert(chr, Box::new(AiList::build(iv_list)));
                }
                Ok(JsOverlapper {
                    internal: core,
                    backend,
                })
            }
            "bits" => {
                for (chr, iv_list) in chr_to_region_store {
                    core.insert(chr, Box::new(Bits::build(iv_list)));
                }
                Ok(JsOverlapper {
                    internal: core,
                    backend,
                })
            }
            _ => Err(JsValue::from_str(&format!(
                "Invalid backend specified: {}",
                backend
            ))),
        }
    }

    pub fn get_backend(&self) -> String {
        self.backend.to_string()
    }

    pub fn intersect(&self, regions: &JsValue) -> Result<JsValue, JsValue> {
        let regions: BedEntries = serde_wasm_bindgen::from_value(regions.to_owned())?;
        let tokens = regions
            .0
            .into_iter()
            .filter_map(|e| {
                let chrom = e.0;
                let start = e.1;
                let end = e.2;

                let olapper = &self.internal.get(&chrom);
                olapper.as_ref().map(|olapper| {
                    olapper
                        .find_iter(start, end)
                        .map(|iv| (chrom.to_string(), iv.start, iv.end))
                        .collect::<Vec<(String, u32, u32)>>()
                })
            })
            .flatten()
            .collect::<Vec<(String, u32, u32)>>();

        serde_wasm_bindgen::to_value(&tokens)
            .map_err(|e| JsValue::from_str(&format!("Error serializing tokens: {}", e)))
    }
}
