// bindings/r/src/rust/src/lib.rs
use extendr_api::prelude::*;
pub mod io;
pub mod igd;

#[extendr]
fn __init__() {}

extendr_module! {
    mod gtars;
    use io;
    use igd;
    fn __init__;
}
