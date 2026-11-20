use anyhow::Result;
use clap::{ArgMatches, FromArgMatches};

/// Run one `gtars blast` verb. The numeric exit code passes through unchanged: 2 for usage
/// errors and a missing engine, 130 when cancelled.
pub fn run_blast(matches: &ArgMatches) -> Result<()> {
    let code = gtars_rgblast::cli::BlastCommand::from_arg_matches(matches)?.run_code("gtars blast");
    if code != 0 {
        std::process::exit(code);
    }
    Ok(())
}
