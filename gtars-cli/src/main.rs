#[cfg(feature = "igd")]
mod igd;

#[cfg(feature = "overlaprs")]
mod overlaprs;

#[cfg(feature = "uniwig")]
mod uniwig;

#[cfg(feature = "bbcache")]
mod bbcache;

#[cfg(feature = "fragsplit")]
mod fragsplit;

#[cfg(feature = "scoring")]
mod scoring;

#[cfg(feature = "genomicdist")]
mod genomicdist;

#[cfg(feature = "genomicdist")]
mod ranges;

#[cfg(feature = "genomicdist")]
mod consensus;

#[cfg(feature = "genomicdist")]
mod prep;

use anyhow::Result;
use clap::Command;

pub mod consts {
    pub const VERSION: &str = env!("CARGO_PKG_VERSION");
    pub const PKG_NAME: &str = "gtars";
    pub const BIN_NAME: &str = "gtars";
}

fn build_parser() -> Command {
    let cmd = Command::new(consts::BIN_NAME)
        .bin_name(consts::BIN_NAME)
        .version(consts::VERSION)
        .author("Databio")
        .about("Performance critical tools for working with genomic interval data with an emphasis on preprocessing for machine learning pipelines.")
        .subcommand_required(true);

    #[cfg(feature = "igd")]
    let cmd = cmd.subcommand(igd::cli::create_igd_cli());

    #[cfg(feature = "overlaprs")]
    let cmd = cmd.subcommand(overlaprs::cli::create_overlap_cli());

    #[cfg(feature = "uniwig")]
    let cmd = cmd.subcommand(uniwig::cli::create_uniwig_cli());

    #[cfg(feature = "bbcache")]
    let cmd = cmd.subcommand(bbcache::cli::create_bbcache_cli());

    #[cfg(feature = "fragsplit")]
    let cmd = cmd.subcommand(fragsplit::cli::create_fragsplit_cli());

    #[cfg(feature = "scoring")]
    let cmd = cmd.subcommand(scoring::cli::create_scoring_cli());

    #[cfg(feature = "genomicdist")]
    let cmd = cmd.subcommand(genomicdist::cli::create_genomicdist_cli());

    #[cfg(feature = "genomicdist")]
    let cmd = cmd.subcommand(ranges::cli::create_ranges_cli());

    #[cfg(feature = "genomicdist")]
    let cmd = cmd.subcommand(consensus::cli::create_consensus_cli());

    #[cfg(feature = "genomicdist")]
    let cmd = cmd.subcommand(prep::cli::create_prep_cli());

    cmd
}

#[allow(unreachable_code)]
fn main() -> Result<()> {
    let app = build_parser();
    let matches = app.get_matches();

    match matches.subcommand() {
        //
        // UNIWIG
        //
        #[cfg(feature = "uniwig")]
        Some((uniwig::cli::UNIWIG_CMD, matches)) => {
            uniwig::handlers::run_uniwig(matches);
        }

        //
        // OVERLAPPES
        //
        #[cfg(feature = "overlaprs")]
        Some((overlaprs::cli::OVERLAP_CMD, matches)) => {
            let _ = overlaprs::handlers::run_overlaprs(matches);
        }

        //
        // BBCACHE
        //
        #[cfg(feature = "bbcache")]
        Some((bbcache::cli::BBCACHE_CMD, matches)) => {
            bbcache::handlers::run_bbcache(matches);
        }

        //
        // IGD
        //
        #[cfg(feature = "igd")]
        Some((igd::cli::IGD_CMD, matches)) => match matches.subcommand() {
            Some((igd::cli::IGD_CREATE, matches)) => {
                igd::handlers::igd_get_create_matches(matches);
            }
            Some((igd::cli::IGD_SEARCH, matches)) => {
                igd::handlers::igd_get_search_matches(matches);
            }
            _ => unreachable!("IGD Subcommand not found"),
        },

        //
        // FRAGMENT SPLITTING UTIL
        //
        #[cfg(feature = "fragsplit")]
        Some((fragsplit::cli::FRAGSPLIT_CMD, matches)) => {
            fragsplit::handlers::run_fragsplit(matches)?;
        }

        //
        // REGION SCORING MATRIX
        //
        #[cfg(feature = "scoring")]
        Some((scoring::cli::FSCORING_CMD, matches)) => {
            scoring::handlers::run_scoring(matches)?;
        }

        //
        // GENOMIC DISTRIBUTIONS
        //
        #[cfg(feature = "genomicdist")]
        Some((genomicdist::cli::GENOMICDIST_CMD, matches)) => {
            genomicdist::handlers::run_genomicdist(matches)?;
        }

        //
        // INTERVAL RANGES
        //
        #[cfg(feature = "genomicdist")]
        Some((ranges::cli::RANGES_CMD, matches)) => {
            ranges::handlers::run_ranges(matches)?;
        }

        //
        // CONSENSUS
        //
        #[cfg(feature = "genomicdist")]
        Some((consensus::cli::CONSENSUS_CMD, matches)) => {
            consensus::handlers::run_consensus(matches)?;
        }

        //
        // PREP (bincode serialization)
        //
        #[cfg(feature = "genomicdist")]
        Some((prep::cli::PREP_CMD, matches)) => {
            prep::handlers::run_prep(matches)?;
        }

        _ => unreachable!("Subcommand not found"),
    };

    Ok(())
}
