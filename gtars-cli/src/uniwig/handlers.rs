use clap::ArgMatches;

use gtars_uniwig::uniwig_main;

/// Matches items from CLAP args before running uniwig_main
pub fn run_uniwig(matches: &ArgMatches) {
    //println!("I am running. Here are the arguments: {:?}", matches);

    let filepath = matches
        .get_one::<String>("file")
        .expect("file path is required");

    let filetype = matches
        .get_one::<String>("filetype")
        .expect("file type is required");

    let chromsizerefpath = matches
        .get_one::<String>("chromref")
        .cloned()
        .expect("chrom sizes is required!");

    let bwfileheader = matches
        .get_one::<String>("fileheader")
        .expect("fileheader is required");

    let smoothsize = matches
        .get_one::<i32>("smoothsize")
        .expect("smoothsize required");

    let output_type = matches
        .get_one::<String>("outputtype")
        .expect("output type is required");

    //let default_vec = &vec!["start", "end", "core"];
    let count_types = matches
        .get_one::<String>("counttype")
        .expect("output type is required");

    // let mut vec_count_type: Vec<&str> = Vec::new();
    let vec_count_type = match count_types.as_str() {
        "all" => {
            vec!["start", "end", "core"]
        }
        "start" => {
            vec!["start"]
        }
        "end" => {
            vec!["end"]
        }
        "core" => {
            vec!["core"]
        }
        "shift" => {
            vec!["shift"]
        }

        _ => {
            vec!["start", "end", "core"]
        }
    };

    //println!("FOUND count_type {:?}", vec_count_type);

    let num_threads = matches
        .get_one::<i32>("threads")
        .expect("requires integer value");

    let bam_scale = matches
        .get_one::<f32>("bamscale")
        .expect("requires int value");

    let score = matches.get_one::<bool>("score").unwrap_or(&false);
    let bam_shift = matches.get_one::<bool>("no-bamshift").unwrap_or(&true);

    let debug = matches.get_one::<bool>("debug").unwrap_or(&false);

    let stepsize = matches
        .get_one::<i32>("stepsize")
        .expect("requires integer value");

    let zoom = matches
        .get_one::<i32>("zoom")
        .expect("requires integer value");

    uniwig_main(
        vec_count_type,
        *smoothsize,
        filepath,
        chromsizerefpath.as_str(),
        bwfileheader,
        output_type,
        filetype,
        *num_threads,
        *score,
        *stepsize,
        *zoom,
        *debug,
        *bam_shift,
        *bam_scale,
    )
    .expect("Uniwig failed.");
}
