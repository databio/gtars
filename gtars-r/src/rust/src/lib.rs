// bindings/r/src/rust/src/lib.rs
use extendr_api::prelude::*;
pub mod genomicdist;
pub mod igd;
pub mod io;
pub mod refget;

#[extendr]
fn __init__() {}

extendr_module! {
    mod gtars;
    use genomicdist;
    use io;
    use igd;
    use refget;
    fn __init__;
}
