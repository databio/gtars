use extendr_api::prelude::*;

pub mod io;

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod gtars;
    use io;
}
