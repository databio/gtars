use wasm_bindgen::prelude::*;

#[wasm_bindgen(js_name = "RegionSet")]
pub struct JsRegionSet {
    internal: String,
}

#[wasm_bindgen(js_class = "RegionSet")]
impl JsRegionSet {

    #[wasm_bindgen(constructor)]
    pub fn new(word: &str) -> Result<JsRegionSet, JsValue> {
        println!("Phineas is ready");
        Ok(JsRegionSet {internal: String::from(word)})
    }

    #[wasm_bindgen(getter, js_name = "get_region")]
    pub fn get_region(&mut self) -> String{
        return self.internal.clone();
    }

    #[wasm_bindgen(getter, js_name = "get_phi")]
    pub fn get_phi(&mut self) -> String{
        String::from("phi-neaaas")
    }
}