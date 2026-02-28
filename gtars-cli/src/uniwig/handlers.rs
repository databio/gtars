use clap::ArgMatches;

use gtars_uniwig::stream::{CountType, OutputFormat, read_chrom_sizes, uniwig_streaming};
use gtars_uniwig::uniwig_main;
use std::fs::File;
use std::io::{self, BufReader, Cursor, Read, Write};

/// Matches items from CLAP args before running uniwig_main
pub fn run_uniwig(matches: &ArgMatches) {
    let streaming = matches.get_one::<bool>("streaming").unwrap_or(&false);

    if *streaming {
        run_uniwig_streaming(matches);
        return;
    }

    // --- existing batch path (unchanged from here down) ---
    let filepath = matches
        .get_one::<String>("file")
        .expect("--file is required when not using --streaming mode");

    let filetype = matches
        .get_one::<String>("filetype")
        .expect("file type is required");

    let chromsizerefpath = matches
        .get_one::<String>("chromref")
        .cloned()
        .expect("chrom sizes is required!");

    let bwfileheader = matches
        .get_one::<String>("fileheader")
        .expect("--fileheader is required when not using --streaming mode");

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

    let wigstep = matches
        .get_one::<String>("wigstep")
        .expect("wigstep is required");

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
        wigstep,
    )
    .expect("Uniwig failed.");
}

fn run_uniwig_streaming(matches: &ArgMatches) {
    // 1. Parse chromsizes
    let chrom_sizes = match matches.get_one::<String>("chromref") {
        Some(path) => {
            let file = File::open(path).expect("Cannot open chrom.sizes file");
            read_chrom_sizes(BufReader::new(file)).expect("Cannot parse chrom.sizes")
        }
        None => std::collections::HashMap::new(),
    };

    // 2. Parse count types to run
    let count_type_str = matches
        .get_one::<String>("counttype")
        .map(|s| s.as_str())
        .unwrap_or("start");

    let count_types_to_run: Vec<(&str, CountType)> = match count_type_str {
        "all" => vec![
            ("start", CountType::Start),
            ("end", CountType::End),
            ("core", CountType::Core),
        ],
        "start" => vec![("start", CountType::Start)],
        "end" => vec![("end", CountType::End)],
        "core" => vec![("core", CountType::Core)],
        other => {
            eprintln!("Error: unknown count type '{}' for streaming mode", other);
            std::process::exit(1);
        }
    };

    // 3. Parse other params
    let smooth_size = *matches
        .get_one::<i32>("smoothsize")
        .expect("smoothsize required");
    let step_size = *matches
        .get_one::<i32>("stepsize")
        .expect("stepsize required");

    let output_format = match matches
        .get_one::<String>("outputtype")
        .map(|s| s.as_str())
        .unwrap_or("wig")
    {
        "wig" => OutputFormat::Wig,
        "bedgraph" | "bg" => OutputFormat::BedGraph,
        other => {
            eprintln!(
                "Error: output type '{}' not supported in streaming mode (use wig or bedgraph)",
                other
            );
            std::process::exit(1);
        }
    };

    let use_stdout = *matches.get_one::<bool>("stdout").unwrap_or(&false);
    let max_gap = *matches.get_one::<i64>("dense").unwrap_or(&0);
    let filepath = matches.get_one::<String>("file").map(|s| s.as_str());
    let is_stdin = filepath.is_none() || filepath == Some("-");

    // 4. For stdin with multiple count types, buffer input into memory
    let input_bytes: Option<Vec<u8>> = if is_stdin && count_types_to_run.len() > 1 {
        let mut buf = Vec::new();
        io::stdin()
            .read_to_end(&mut buf)
            .expect("Failed to read stdin");
        Some(buf)
    } else {
        None
    };

    // 5. Run each count type
    for (label, ct) in &count_types_to_run {
        let input: Box<dyn Read> = if let Some(ref bytes) = input_bytes {
            Box::new(Cursor::new(bytes.clone()))
        } else if is_stdin {
            Box::new(io::stdin())
        } else {
            Box::new(File::open(filepath.unwrap()).expect("Cannot open input file"))
        };

        let mut output: Box<dyn Write> = if use_stdout {
            if count_types_to_run.len() > 1 {
                // Write a separator comment when outputting multiple types to stdout
                let stdout = io::stdout();
                let mut handle = stdout.lock();
                writeln!(handle, "# count_type={}", label).expect("Failed to write separator");
            }
            Box::new(io::stdout().lock())
        } else {
            let header = matches
                .get_one::<String>("fileheader")
                .expect("--fileheader required for file output in streaming mode");
            let filename = format!("{}_{}.wig", header, label);
            Box::new(File::create(&filename).expect("Cannot create output file"))
        };

        uniwig_streaming(
            input,
            &mut output,
            chrom_sizes.clone(),
            smooth_size,
            step_size,
            *ct,
            output_format,
            max_gap,
        )
        .expect("Streaming uniwig failed");
    }
}
