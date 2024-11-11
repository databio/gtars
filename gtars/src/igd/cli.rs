use crate::igd::consts::IGD_CMD;
use clap::{arg, Command};

pub fn create_igd_cli() -> Command {
    Command::new(IGD_CMD)
        .author("DRC")
        .about("Create or search an integrated genome database (IGD)")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("create")
                .about("Create igd database")
                .arg(arg!(--output <VALUE> "Path to the output.").required(true))
                .arg(
                    arg!(--filelist <VALUE> "Path to the list of files. This should be a folder of bed files.")
                        .required(true),
                )
                .arg(
                    arg!(--dbname <VALUE> "Database name")
                        .required(false).default_value("igd_database"),
                )
        )
        .subcommand(
            Command::new("search")
                .about("Search igd database")
                .arg(arg!(--database <VALUE> "Path to the igd database.").required(true).short('d'))
                .arg(
                    arg!(--query <VALUE> "Path to the query file (.bed or .bed.gz)")
                        .required(true).short('q'),
                )
                .arg(
                    arg!(--singlequery <VALUE> "chrN start end (a single query)")
                        .required(false).short('r'),
                )
                .arg(
                    arg!(--signalvalue <VALUE> "signal value 0-1000 (signal value > v)")
                        .required(false).short('v'),
                )
                .arg(
                    arg!(--output <VALUE> "output file path and name")
                        .required(false).short('o'),
                )
                .arg(
                    arg!(--outseqpare <VALUE> "output seqpare similarity")
                        .required(false).short('s'),
                )
                .arg(
                    arg!(--full <VALUE> "output full overlaps, for -q and -r only")
                        .required(false).short('f'),
                )
                .arg(
                    arg!(--hitsmap <VALUE> "hitsmap of igd datasets")
                        .required(false).short('m'),
                )
        )
}
