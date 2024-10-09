#![allow(nonstandard_style)]

pub mod cli;
pub mod create;
pub mod search;

pub mod consts {
    pub const IGD_CMD: &str = "igd";
    pub const IGD_CREATE: &str = "create";
    pub const IGD_SEARCH: &str = "search";
}
