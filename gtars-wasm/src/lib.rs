mod models;
mod overlaprs;
mod refget;
mod regionset;
mod tokenizers;
mod utils;

use wasm_bindgen::prelude::*;

// Re-export refget functions at the top level
pub use refget::*;

#[wasm_bindgen]
pub fn greet(name: &str) {
    alert(&format!("Hello, {}!", name));
}

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
}
