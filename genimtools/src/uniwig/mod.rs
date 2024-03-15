use clap::ArgMatches;

pub mod cli;

pub fn run_uniwig(matches: &ArgMatches) {
    println!("Im running. Here are the arguments: {:?}", matches)
}

pub mod consts {
    pub const UNIWIG_CMD: &str = "uniwig";

}