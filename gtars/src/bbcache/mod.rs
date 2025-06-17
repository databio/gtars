
use clap::ArgMatches;
use std::collections::HashMap;
use crate::bbclient::bbclient::BBClient;



pub mod utils;
pub mod bbclient;

pub fn run_bbclient(matches: &ArgMatches) {
    let cache_folder = matches.get_one::<String>("cache_folder").unwrap();
    let bedbase_api = matches.get_one::<String>("bedbase_api").unwrap();
    let bedset_id = matches.get_one::<String>("bedset_id").unwrap();

    let bbclient = bbclient::BBClient::new(cache_folder, bedbase_api);
    let bedset_data = bbclient.load_bedset(bedset_id);

    for entry in bedset_data {
        println!("{}", entry);
    }
}