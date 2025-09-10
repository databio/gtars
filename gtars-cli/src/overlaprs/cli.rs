use clap::{Command, arg};

pub const OVERLAP_CMD: &str = "overlaprs";

pub fn create_overlap_cli() -> Command {
    Command::new(OVERLAP_CMD)
        .author("NJL")
        .about("Tokenize data into a universe")
        .arg_required_else_help(true)
        .arg(arg!(-q <query> "The file you are tokenizing"))
        .arg(arg!(-u <universe> "The universe you are tokenizing into"))
        .arg(arg!(-e --backend <backend> "Which backend to use (ailist or bits)"))
        .arg(arg!(--streaming "Use streaming mode for very large universes (lower memory usage)"))
}
