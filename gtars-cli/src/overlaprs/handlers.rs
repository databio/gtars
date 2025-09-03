use std::fmt::Write as FmtWrite;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;

use anyhow::Result;
use clap::ArgMatches;
use fxhash::FxHashMap as HashMap;

use overlaprs::Overlapper;
use overlaprs::{AiList, Bits};

use gtars_core::models::Interval;
use gtars_core::utils::get_dynamic_reader;

enum BackendType {
    Bits,
    AiList,
}

type OverlapperMap = HashMap<String, Box<dyn Overlapper<u32, Option<u8>>>>;

pub fn run_overlaprs(matches: &ArgMatches) -> Result<()> {
    let query_file = matches
        .get_one::<String>("query")
        .expect("A path to a query file is required.");

    let universe_file = matches
        .get_one::<String>("universe")
        .expect("A path to a universe file is required.");

    let default_backend = "bits".to_string();
    let backend_str = matches
        .get_one::<String>("backend")
        .unwrap_or(&default_backend);

    let backend_type = match backend_str.as_str() {
        "bits" => BackendType::Bits,
        "ailist" => BackendType::AiList,
        _ => {
            return Err(anyhow::anyhow!(
                "Invalid backend type: {}. Valid options are 'bits' or 'ailist'",
                backend_str
            ));
        }
    };

    // Build overlap data structures directly from universe file
    let core = build_overlap_structures(universe_file, backend_type)?;

    // Process queries with buffered output
    process_queries(query_file, &core)?;

    Ok(())
}

fn build_overlap_structures(
    universe_file: &str,
    backend_type: BackendType,
) -> Result<OverlapperMap> {
    let universe_path = Path::new(universe_file);
    let universe_reader = get_dynamic_reader(universe_path)?;

    let mut intervals_by_chr: HashMap<String, Vec<Interval<u32, Option<u8>>>> = HashMap::default();

    // Single pass: collect intervals by chromosome
    for line in universe_reader.lines() {
        let line = line?;

        // Pre-allocate and reuse split iterator
        let mut fields = line.split('\t');

        let chr = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing chromosome field"))?;
        let start_str = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing start field"))?;
        let end_str = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing end field"))?;

        let start = start_str.parse::<u32>()?;
        let end = end_str.parse::<u32>()?;

        let interval = Interval {
            start,
            end,
            val: None,
        };
        intervals_by_chr
            .entry(chr.to_string())
            .or_default()
            .push(interval);
    }

    // Build overlap structures for each chromosome
    let mut core: OverlapperMap = HashMap::default();
    core.reserve(intervals_by_chr.len());

    for (chr, chr_intervals) in intervals_by_chr {
        if chr_intervals.is_empty() {
            continue;
        }

        let overlapper: Box<dyn Overlapper<u32, Option<u8>>> = match backend_type {
            BackendType::Bits => Box::new(Bits::build(chr_intervals)),
            BackendType::AiList => Box::new(AiList::build(chr_intervals)),
        };
        core.insert(chr, overlapper);
    }

    Ok(core)
}

fn process_queries(query_file: &str, core: &OverlapperMap) -> Result<()> {
    let query_path = Path::new(query_file);
    let query_reader = get_dynamic_reader(query_path)?;

    let stdout = io::stdout();
    let mut writer = BufWriter::new(stdout.lock());
    let mut output_buffer = String::with_capacity(1024); // Pre-allocate buffer

    for line in query_reader.lines() {
        let line = line?;

        let mut fields = line.split('\t');

        let chr = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing chromosome field"))?;
        let start_str = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing start field"))?;
        let end_str = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing end field"))?;

        let start = start_str.parse::<u32>()?;
        let end = end_str.parse::<u32>()?;

        // Skip if chromosome not in universe
        if let Some(overlapper) = core.get(chr) {
            let hits = overlapper.find_iter(start, end);

            // Batch output writing
            for hit in hits {
                output_buffer.clear();
                writeln!(&mut output_buffer, "{}\t{}\t{}", chr, hit.start, hit.end)?;
                writer.write_all(output_buffer.as_bytes())?;
            }
        }
    }

    writer.flush()?;

    Ok(())
}
