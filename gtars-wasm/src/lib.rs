mod bed_stream;
mod models;
mod overlaprs;
mod partitions;
mod refget;
mod regionset;
mod signal;
mod tokenizers;
mod tss;
mod utils;

use wasm_bindgen::prelude::*;

// Re-export functions at the top level
pub use bed_stream::*;
pub use refget::*;

#[wasm_bindgen]
pub fn greet(name: &str) {
    alert(&format!("Hello, {}!", name));
}

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
}
