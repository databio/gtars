mod models;
mod overlaprs;
mod tokenizers;
mod utils;
mod regionset;

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn greet(name: &str) {
    alert(&format!("Hello, {}!", name));
}

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
}
