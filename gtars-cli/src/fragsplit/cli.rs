use clap::{Arg, Command, arg};

pub const FRAGSPLIT_CMD: &str = "pb";
pub const DEFAULT_OUT: &str = "out/";

pub fn create_fragsplit_cli() -> Command {
    Command::new(FRAGSPLIT_CMD)
        .author("Nathan LeRoy")
        .about("Split fragment files into pseudobulks based on cluster labels.")
        .arg(Arg::new("fragments"))
        .arg(Arg::new("mapping"))
        .arg(arg!(--output <output>))
}
