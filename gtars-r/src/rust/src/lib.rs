// bindings/r/src/rust/src/lib.rs
use extendr_api::prelude::*;
pub mod genomicdist;
pub mod igd;
pub mod io;
pub mod refget;
pub mod sc;

#[extendr]
fn __init__() {}

extendr_module! {
    mod gtars;
    use genomicdist;
    use io;
    use igd;
    use refget;
    use sc;
    fn __init__;
}
