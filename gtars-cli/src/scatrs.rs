// Re-export the scatrs library components
pub use gtars_scatrs::{cli, SCATRS_CMD};

pub mod handlers {
    use anyhow::Result;
    use clap::ArgMatches;

    pub fn run_scatrs(matches: &ArgMatches) -> Result<()> {
        // Delegate to the scatrs library's handler
        gtars_scatrs::cli::handle_scatrs_command(matches)
    }
}