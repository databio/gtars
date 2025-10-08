// bindings/r/src/rust/src/lib.rs
use extendr_api::prelude::*;
pub mod igd;
pub mod io;
pub mod refget;

#[extendr]
fn __init__() {}

extendr_module! {
    mod gtars;
    use io;
    use igd;
    use refget;
    fn __init__;
}
