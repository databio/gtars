// bindings/r/src/rust/src/lib.rs
use extendr_api::prelude::*;
pub mod igd;
pub mod io;

#[extendr]
fn __init__() {}

extendr_module! {
    mod gtars;
    use io;
    use igd;
    fn __init__;
}
