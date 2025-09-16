use clap::{Arg, Command, arg};

pub use fragsplit::consts::*;

pub fn create_fragsplit_cli() -> Command {
    Command::new(FRAGSPLIT_CMD)
        .author("Nathan LeRoy")
        .about("Split fragment files into pseudobulks based on cluster labels.")
        .arg(Arg::new("fragments"))
        .arg(Arg::new("mapping"))
        .arg(arg!(--output <output>))
}
