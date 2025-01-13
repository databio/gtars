use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use rstest::*;

#[fixture]
fn path_to_data() -> &'static str {
    "tests/data"
}

#[fixture]
fn path_to_bed_file() -> &'static str {
    "tests/data/peaks.bed"
}

#[fixture]
fn path_to_sorted_small_bed_file() -> &'static str {
    "tests/data/test_sorted_small.bed"
}

#[fixture]
fn path_to_small_bam_file() -> &'static str {
    "tests/data/test_chr22_small.bam"
    //"/home/drc/Downloads/bam files for rust test/test1_sort_dedup.bam"
}

#[fixture]
fn path_to_chrom_sizes_file() -> &'static str {
    "tests/hg38.chrom.sizes"
}

#[fixture]
fn path_to_bed_file_gzipped() -> &'static str {
    "tests/data/peaks.bed.gz"
}

#[fixture]
fn path_to_dummy_bed_file() -> &'static str {
    "tests/data/dummy.bed"
}

#[fixture]
fn path_to_dummy_chromsizes() -> &'static str {
    "tests/data/dummy.chrom.sizes"
}

#[fixture]
fn path_to_dummy_narrowpeak() -> &'static str {
    "tests/data/dummy.narrowPeak"
}

#[fixture]
fn path_to_start_wig_output() -> &'static str {
    "tests/data/out/_start.wig"
}

#[fixture]
fn path_to_core_wig_output() -> &'static str {
    "tests/data/out/_core.wig"
}

#[fixture]
fn path_to_start_bedgraph_output() -> &'static str {
    "tests/data/out/_start.bedGraph"
}

#[fixture]
fn path_to_core_bedgraph_output() -> &'static str {
    "tests/data/out/_core.bedGraph"
}

mod tests {
    use super::*;
    use gtars::igd::create::{
        create_igd_f, gdata_t, igd_add, igd_saveT, igd_save_db, igd_t, parse_bed,
    };
    use gtars::igd::search::{
        getOverlaps, get_file_info_tsv, get_igd_info, get_tsv_path, igd_search, igd_t_from_disk,
    };

    use gtars::uniwig::{uniwig_main, Chromosome};

    use gtars::uniwig::counting::{core_counts, start_end_counts};
    use gtars::uniwig::reading::{
        parse_bed_file, read_bam_header, read_bed_vec, read_chromosome_sizes, read_narrow_peak_vec,
    };

    use gtars::uniwig::writing::write_bw_files;

    use anyhow::Context;
    use byteorder::{LittleEndian, ReadBytesExt};
    use std::collections::HashMap;
    use std::collections::HashSet;
    use std::fs::OpenOptions;
    use std::io::{Seek, SeekFrom};
    // IGD TESTS

    #[rstest]
    fn test_igd_parse_bed_file() {
        // Given some random line from a  bed file...
        let bed_file_string =
            String::from("chr1	32481	32787	SRX4150706.05_peak_1	92	.	7.69231	13.22648	9.25988	155");

        //Placeholder start and end values
        let mut start = 0;
        let mut end = 0;
        let mut va = 0;

        let result = parse_bed(&bed_file_string, &mut start, &mut end, &mut va).unwrap(); // this will return

        let unwrapped_result = result.as_str();

        assert_eq!(unwrapped_result, "chr1");

        // Ensure start and end is modified via parse_bed
        assert_eq!(start, 32481);
        assert_eq!(end, 32787);
    }

    #[rstest]
    fn test_igd_create_short_long_regions() {
        // Depending on start and end coordinates which are divided by nbp=16384
        // the number of tiles per ctg are adjusted, this tests to ensure they are created appropriately
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());
        let db_path_unwrapped = path.into_os_string().into_string().unwrap();
        let db_output_path = db_path_unwrapped;

        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let testfilelists = format!("{}{}", path_to_crate, "/tests/data/igd_file_list_01/");

        let demo_name = String::from("demo");

        let igd = create_igd_f(&db_output_path, &testfilelists, &demo_name);
        assert_eq!(igd.ctg[0].name, "chr1");
        assert_eq!(igd.ctg[1].name, "chr2");
        assert_eq!(igd.ctg[2].name, "chr3");
        assert_eq!(igd.nctg, 3);

        assert_eq!(igd.ctg[0].mTiles, 4); // chr1 has 4 Tiles because of the 32768, and 49152 starts
        assert_eq!(igd.ctg[1].mTiles, 1); // chr only has 1 Tile due to the 200 start

        assert_eq!(igd.ctg[0].gTile[0].gList[0].start, 1); // look specific tile's start
        assert_eq!(
            igd.ctg[0].gTile[(igd.ctg[0].mTiles - 1) as usize].gList[0].start,
            49152
        ); // look specific tile's start

        assert_eq!(igd.ctg[0].gTile[0].nCnts, 2); // look at nCnts
        assert_eq!(igd.ctg[0].gTile[1].nCnts, 0); // look at nCnts
        assert_eq!(igd.ctg[0].gTile[2].nCnts, 1); // look at nCnts

        // Overall stats
        assert_eq!(igd.total_regions, 8);
        assert_eq!(igd.total_average, 998.0);
        assert_eq!(igd.average_length, 124.75);
    }

    #[rstest]
    fn test_igd_create_then_load_from_disk() {
        // Depending on start and end coordinates which are divided by nbp=16384
        // the number of tiles per ctg are adjusted, this tests to ensure they are created appropriately
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());
        let mut db_path_unwrapped = path.into_os_string().into_string().unwrap();
        db_path_unwrapped.push_str("/");
        let db_output_path = db_path_unwrapped.clone();

        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let testfilelists = format!("{}{}", path_to_crate, "/tests/data/igd_file_list_01/");

        let demo_name = String::from("demo");

        let igd_saved = create_igd_f(&db_output_path, &testfilelists, &demo_name);

        println!("dboutput_path {}", db_output_path);

        db_path_unwrapped.push_str("/demo.igd");

        let mut hash_table: HashMap<String, i32> = HashMap::new();

        // Create IGD Struct from database
        let mut igd_from_disk: igd_t_from_disk =
            get_igd_info(&db_path_unwrapped, &mut hash_table).expect("Could not open IGD");
        let tsv_path = get_tsv_path(db_path_unwrapped.as_str()).unwrap();
        get_file_info_tsv(tsv_path, &mut igd_from_disk).unwrap(); //sets igd.finfo

        assert_eq!(igd_saved.ctg.len(), igd_from_disk.nCtg as usize);

        assert_eq!(igd_from_disk.nFiles, 1);

        assert_eq!(
            igd_from_disk.nCnt[0].len(),
            igd_saved.ctg[0].mTiles as usize
        );
        assert_eq!(
            igd_from_disk.nCnt[1].len(),
            igd_saved.ctg[1].mTiles as usize
        );
        assert_eq!(
            igd_from_disk.nCnt[2].len(),
            igd_saved.ctg[2].mTiles as usize
        );

        assert_eq!(igd_from_disk.nCnt[0][0], igd_saved.ctg[0].gTile[0].nCnts);
        assert_eq!(igd_from_disk.nCnt[0][1], igd_saved.ctg[0].gTile[1].nCnts);
        assert_eq!(igd_from_disk.nCnt[0][2], igd_saved.ctg[0].gTile[2].nCnts);
        assert_eq!(igd_from_disk.nCnt[0][3], igd_saved.ctg[0].gTile[3].nCnts);

        // Check to see if the regions on disk are the same as the original igd (minus the unused zeros)
        let dbpath = std::path::Path::new(&db_path_unwrapped);
        let db_file = OpenOptions::new()
            .create(true)
            .append(true)
            .read(true)
            .open(dbpath)
            .unwrap();
        let mut db_reader = BufReader::new(db_file);

        for k in 0..3 {
            let nCnt_len = igd_from_disk.nCnt[k].len();

            for l in 0..nCnt_len {
                let mut a: HashSet<i32> = Default::default();
                let mut b: HashSet<i32> = Default::default();

                let tmpi = igd_from_disk.nCnt[k][l]; // number of gdata_t to read

                //println!("Here is k {}, l {}, and igd_from_disk.tIdx[k][l] {}",k,l, igd_from_disk.tIdx[k][l]);
                db_reader
                    .seek(SeekFrom::Start(igd_from_disk.tIdx[k][l] as u64)) // [k]contig [l] tile position
                    .unwrap();

                let mut gData: Vec<gdata_t> = Vec::new();

                //println!("Creating gData with tmpi {}", tmpi);
                for j in 0..tmpi {
                    gData.push(gdata_t::default())
                }

                for i in 0..tmpi {
                    // number of gdata_t to read
                    //println!("Iterating with i {} of tmpi {} ",i,tmpi);
                    let mut buf = [0u8; 16];

                    let n = db_reader.read(&mut buf).unwrap();

                    if n == 0 {
                        //println!("Breaking loop while reading tempfile");
                        break;
                    } else if n != 16 {
                        //panic!("Cannot read temp file.");
                        break;
                    }

                    let mut rdr = &buf[..] as &[u8];
                    let idx = rdr.read_i32::<LittleEndian>().unwrap();
                    let start = rdr.read_i32::<LittleEndian>().unwrap();
                    let end = rdr.read_i32::<LittleEndian>().unwrap();
                    let value = rdr.read_i32::<LittleEndian>().unwrap();

                    //println!("Looping through g_datat in temp files");
                    //println!("Chr_name: {} Filename: {}  start: {} end: {}", igd_from_disk.cName[k], igd_from_disk.file_info[idx as usize].fileName, start, end);

                    gData[i as usize] = gdata_t {
                        idx: idx,
                        start,
                        end,
                        value,
                    };
                }

                //println!("here is k {}, l {}",k,l);
                for g in gData.iter() {
                    //println!("Inserting {} from gData on Disk", g.start);
                    a.insert(g.start);
                }

                for g in igd_saved.ctg[k].gTile[l].gList.iter() {
                    //println!("Inserting {} from original gList ", g.start);
                    b.insert(g.start);
                }
                //println!("A: {:?}", a);
                //println!("B: {:?}", b);
                // There difference should at most be a 0 from unused tiles, therefore the difference length should at MOST be 1.
                let diff = b.difference(&a).collect::<Vec<&i32>>();
                //println!("Difference: {:?}", diff);
                assert!(diff.len() <= 1)
            }
        }
    }

    #[rstest]
    #[case(
        "/tests/data/igd_file_list_01/",
        "/tests/data/igd_query_files/query1.bed",
        8,
        8
    )]
    // #[case(
    //     "/tests/data/igd_file_list_02/",
    //     "/tests/data/igd_query_files/query2.bed",
    //     4,
    //     1
    // )]
    fn test_igd_create_then_search(
        #[case] input: &str,
        #[case] query_file: &str,
        #[case] expected_regions: u32,
        #[case] expected_hits: u32,
    ) {
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());
        let mut db_path_unwrapped = path.into_os_string().into_string().unwrap();
        db_path_unwrapped.push_str("/");
        let db_output_path = db_path_unwrapped.clone();

        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let testfilelists = format!("{}{}", path_to_crate, input);

        let demo_name = String::from("demo");

        let igd_saved = create_igd_f(&db_output_path, &testfilelists, &demo_name);

        println!("dboutput_path {}", db_output_path);

        db_path_unwrapped.push_str("/demo.igd");

        let queryfile = format!("{}{}", path_to_crate, query_file);
        let res = igd_search(&db_path_unwrapped, &queryfile).expect("Error during testing:");
        let mut res_iter = res[1].split('\t');

        // Skip the first two columns
        res_iter.next().unwrap();

        // Extract the third and fourth columns
        let second_column = res_iter.next().unwrap().to_string();
        let third_column = res_iter.next().unwrap().to_string();

        println!("Number of Regions: {}", second_column);
        println!("Number of Hits: {}", third_column);

        assert_eq!(second_column, expected_regions.to_string());
        assert_eq!(third_column, expected_hits.to_string());
    }

    #[rstest]
    fn test_igd_add() {
        // First create a new igd struct

        let mut igd = igd_t::new();
        // create hash table
        let mut hash_table: HashMap<String, i32> = HashMap::new();

        // Set values of struct
        igd.gType = 1;
        igd.nbp = 16384; // from og code tile_size = 16384;  -> this is the bin size (2^14) from the original paper
        igd.nctg = 0;
        igd.mctg = 32;
        igd.total = 0;

        // Given some random line from a bed file...
        let bed_file_string =
            String::from("chr1	32481	32787	SRX4150706.05_peak_1	92	.	7.69231	13.22648	9.25988	155");
        //Placeholder start and end values
        let mut start = 0;
        let mut end = 0;
        let mut va = 0;

        // We've now parsed to get the chromosome and the new start and end of the current contig.
        let result = parse_bed(&bed_file_string, &mut start, &mut end, &mut va).unwrap();
        let chromosome = result;

        // Add to the database (hash table)
        igd_add(&mut igd, &mut hash_table, chromosome, start, end, 0, 0);
    }

    #[rstest]
    fn test_igd_saving() {
        let mut igd = igd_t::new();
        // create hash table
        let mut hash_table: HashMap<String, i32> = HashMap::new();

        // Set values of struct
        igd.gType = 1;
        igd.nbp = 16384; // from og code tile_size = 16384;  -> this is the bin size (2^14) from the original paper
        igd.nctg = 0;
        igd.mctg = 32;
        igd.total = 0;

        // Given some random line from a bed file...
        let bed_file_string =
            String::from("chr1	32481	32787	SRX4150706.05_peak_1	92	.	7.69231	13.22648	9.25988	155");
        //Placeholder start and end values
        let mut start = 0;
        let mut end = 0;
        let mut va = 0;

        // We've now parsed to get the chromosome and the new start and end of the current contig.
        let result = parse_bed(&bed_file_string, &mut start, &mut end, &mut va).unwrap();
        let chromosome = result;

        // Add to the database (hash table)
        igd_add(&mut igd, &mut hash_table, chromosome, start, end, 0, 0);

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let db_path_unwrapped = path.into_os_string().into_string().unwrap();
        let db_output_path = &db_path_unwrapped;

        // First test igd_saveT
        igd_saveT(&mut igd, db_output_path);

        // then test saveing main databse

        igd_save_db(&mut igd, db_output_path, &String::from("randomname"));
    }

    // UNIWIG TESTS
    #[rstest]
    fn test_uniwig_parsed_bed_file(path_to_bed_file: &str) {
        let path = Path::new(path_to_bed_file);
        let file = File::open(path).unwrap();

        let mut reader = BufReader::new(file);
        let first_line = reader.by_ref().lines().next().unwrap().expect("expect");
        println!("{:?}", first_line);

        let result = parse_bed_file(&first_line);

        if let Some((ctg, st, en)) = result {
            println!("ctg: {}", ctg);
            println!("st: {}", st);
            println!("en: {}", en);
            assert_eq!(st, 7915738);
        } else {
            panic!("Failed to parse BED record");
        }
    }

    #[rstest]
    fn test_read_bed_vec(path_to_bed_file: &str, path_to_bed_file_gzipped: &str) {
        let result1 = read_bed_vec(path_to_bed_file);
        assert_eq!(result1.len(), 20);

        let result2 = read_bed_vec(path_to_bed_file_gzipped);
        assert_eq!(result2.len(), 20);
    }

    #[rstest]
    fn test_read_narrow_peak_vec() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let path_to_narrow_peak = format!("{}{}", path_to_crate, "/tests/data/dummy.narrowPeak");
        let result1 = read_narrow_peak_vec(&path_to_narrow_peak);
        assert_eq!(result1.len(), 1);

        let path_to_narrow_peak_gzipped =
            format!("{}{}", path_to_crate, "/tests/data/dummy.narrowPeak.gz");

        let result2 = read_narrow_peak_vec(&path_to_narrow_peak_gzipped);
        assert_eq!(result2.len(), 1);
    }

    #[rstest]
    fn test_read_narrow_peak_chrom_sizes() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let path_to_narrow_peak = format!("{}{}", path_to_crate, "/tests/data/dummy.narrowPeak");
        let _result1 = read_chromosome_sizes(path_to_narrow_peak.as_str());
    }

    #[rstest]
    fn test_read_narrow_peak_core_counts() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let path_to_narrow_peak = format!("{}{}", path_to_crate, "/tests/data/dummy.narrowPeak");
        let chrom_sizes = read_chromosome_sizes(path_to_narrow_peak.as_str()).unwrap();
        let narrow_peak_vec: Vec<Chromosome> = read_narrow_peak_vec(path_to_narrow_peak.as_str());
        let stepsize = 1;

        for chromosome in narrow_peak_vec.iter() {
            let current_chrom_size = *chrom_sizes.get(&chromosome.chrom).unwrap() as i32;
            let _result = core_counts(
                &chromosome.starts,
                &chromosome.ends,
                current_chrom_size,
                stepsize,
            );
        }
    }

    #[rstest]
    fn test_read_narrow_peak_starts_counts() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let path_to_narrow_peak = format!("{}{}", path_to_crate, "/tests/data/dummy.narrowPeak");
        let chrom_sizes = read_chromosome_sizes(path_to_narrow_peak.as_str()).unwrap();
        let narrow_peak_vec: Vec<Chromosome> = read_narrow_peak_vec(path_to_narrow_peak.as_str());
        let stepsize = 1;
        let smooth_size = 1;

        for chromosome in narrow_peak_vec.iter() {
            let current_chrom_size = *chrom_sizes.get(&chromosome.chrom).unwrap() as i32;
            let _result = start_end_counts(
                &chromosome.starts,
                current_chrom_size,
                smooth_size,
                stepsize,
            );
        }
    }

    #[rstest]
    fn test_read_bed_vec_length(path_to_sorted_small_bed_file: &str) {
        let chromosomes: Vec<Chromosome> = read_bed_vec(path_to_sorted_small_bed_file);
        let num_chromosomes = chromosomes.len();

        assert_eq!(num_chromosomes, 5);
    }

    #[rstest]
    fn test_read_bam_header(path_to_small_bam_file: &str) {
        let chromosomes: Vec<Chromosome> = read_bam_header(path_to_small_bam_file);
        let num_chromosomes = chromosomes.len();
        println!("Number of chroms: {}", num_chromosomes);
        assert_eq!(num_chromosomes, 1);
    }

    #[rstest]
    fn test_process_bam(
        path_to_small_bam_file: &str,
    ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/tests/hg38.chrom.sizes");
        let chromsizerefpath = chromsizerefpath.as_str();
        let combinedbedpath = path_to_small_bam_file;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bw";
        let filetype = "bam";
        let num_threads = 2;
        let score = false;
        let stepsize = 1;
        let zoom = 0;

        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        )
        .expect("Uniwig main failed!");

        Ok(())
    }

    #[rstest]
    fn test_process_bam_to_bed(
        path_to_small_bam_file: &str,
    ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/tests/hg38.chrom.sizes");
        let chromsizerefpath = chromsizerefpath.as_str();
        let combinedbedpath = path_to_small_bam_file;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bed";
        let filetype = "bam";
        let num_threads = 2;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        )
        .expect("Uniwig main failed!");

        Ok(())
    }

    #[rstest]
    fn test_run_uniwig_main_wig_type() -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        // This test uses the bed file to determine chromsizes for speed
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!("{}{}", path_to_crate, "/tests/data/test5.bed");
        let combinedbedpath = tempbedpath.as_str();

        let chromsizerefpath = combinedbedpath;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 5;
        let output_type = "wig";
        let filetype = "bed";
        let num_threads = 6;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        )
        .expect("Uniwig main failed!");

        Ok(())
    }

    #[rstest]
    fn test_run_uniwig_main_npy_type() -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        // This test uses the bed file to determine chromsizes for speed
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!("{}{}", path_to_crate, "/tests/data/test5.bed");
        let combinedbedpath = tempbedpath.as_str();

        let chromsizerefpath = combinedbedpath;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 2;
        let output_type = "npy";
        let filetype = "bed";
        let num_threads = 6;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        )
        .expect("Uniwig main failed!");
        Ok(())
    }

    #[rstest]
    fn test_reading_chrom_sizes() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        // Read from sizes file
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/tests/hg38.chrom.sizes");
        let chrom_sizes = read_chromosome_sizes(chromsizerefpath.as_str()).unwrap();
        let chrom_name = String::from("chr13");
        let current_chrom_size = chrom_sizes[&chrom_name.clone()] as i32;
        assert_eq!(current_chrom_size, 114364328);

        // Read from BED file
        let tempbedpath = format!("{}{}", path_to_crate, "/tests/data/test5.bed");
        let combinedbedpath = tempbedpath.as_str();
        let chrom_sizes = read_chromosome_sizes(combinedbedpath).unwrap();
        let chrom_name = String::from("chr1");
        let current_chrom_size = chrom_sizes[&chrom_name.clone()] as i32;
        assert_eq!(current_chrom_size, 32);
    }

    #[rstest]
    fn test_uniwig_mismatched_chrom_sizes(_path_to_bed_file: &str) {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        // Read from sizes file
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/tests/hg38.chrom.sizes");

        // Read from BED file that contains chromosomes not in size file
        let tempbedpath = format!("{}{}", path_to_crate, "/tests/data/test_unknown_chrom.bed");
        let combinedbedpath = tempbedpath.as_str();
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 5;
        let output_type = "npy";
        let filetype = "bed";
        let num_threads: i32 = 6;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        let result = uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            &chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        );

        assert!(result.is_ok());
    }

    #[rstest]
    fn test_uniwig_write_bw(_path_to_bed_file: &str) {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        // Read from sizes file
        let directory_bed_graphs: String = format!("{}{}", path_to_crate, "/tests/data");
        let chrom_sizes: String = format!("{}{}", path_to_crate, "/tests/data/dummy.chrom.sizes");
        let num_threads = 2;
        let zoom = 0;

        write_bw_files(
            directory_bed_graphs.as_str(),
            chrom_sizes.as_str(),
            num_threads,
            zoom,
        );
    }

    #[rstest]
    fn test_uniwig_wiggle_output(
        _path_to_dummy_bed_file: &str,
        _path_to_dummy_chromsizes: &str,
        _path_to_start_wig_output: &str,
        _path_to_core_wig_output: &str,
    ) {
        let chromsizerefpath = _path_to_dummy_chromsizes;
        let combinedbedpath = _path_to_dummy_bed_file;
        let test_output_path = _path_to_start_wig_output;
        let core_test_output_path = _path_to_core_wig_output;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let mut bwfileheader_path = path.into_os_string().into_string().unwrap();
        bwfileheader_path.push_str("/final/");

        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "wig";
        let filetype = "bed";
        let num_threads: i32 = 2;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        let result = uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            &chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        );

        assert!(result.is_ok());

        // Test _start.wig output
        let path = PathBuf::from(&tempdir.path());
        let mut final_start_file_path = path.into_os_string().into_string().unwrap();
        final_start_file_path.push_str("/final/_start.wig");
        let final_start_file_path = final_start_file_path.as_str();

        let file1 = File::open(final_start_file_path).unwrap();
        let file2 = File::open(test_output_path).unwrap();

        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            let line1 = lines1.next().transpose().unwrap();
            let line2 = lines2.next().transpose().unwrap();

            match (line1, line2) {
                (Some(line1), Some(line2)) => {
                    assert_eq!(line1, line2);
                }
                (None, None) => {
                    break; // Both files reached the end
                }
                _ => {
                    panic!("FILES ARE NOT EQUAL!!!")
                }
            }
        }

        // Test _core.wig output
        let path = PathBuf::from(&tempdir.path());
        let mut final_core_file_path = path.into_os_string().into_string().unwrap();
        final_core_file_path.push_str("/final/_core.wig");
        let final_core_file_path = final_core_file_path.as_str();

        let file1 = File::open(final_core_file_path).unwrap();
        let file2 = File::open(core_test_output_path).unwrap();

        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            let line1 = lines1.next().transpose().unwrap();
            let line2 = lines2.next().transpose().unwrap();

            match (line1, line2) {
                (Some(line1), Some(line2)) => {
                    assert_eq!(line1, line2);
                }
                (None, None) => {
                    break; // Both files reached the end
                }
                _ => {
                    panic!("FILES ARE NOT EQUAL!!!")
                }
            }
        }
    }

    #[rstest]
    fn test_uniwig_bedgraph_output(
        _path_to_dummy_bed_file: &str,
        _path_to_dummy_chromsizes: &str,
        _path_to_start_bedgraph_output: &str,
        _path_to_core_bedgraph_output: &str,
    ) {
        let chromsizerefpath = _path_to_dummy_chromsizes;
        let combinedbedpath = _path_to_dummy_bed_file;
        let test_output_path = _path_to_start_bedgraph_output;
        let core_test_output_path = _path_to_core_bedgraph_output;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let mut bwfileheader_path = path.into_os_string().into_string().unwrap();
        bwfileheader_path.push_str("/final/");

        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bedgraph";
        let filetype = "bed";
        let num_threads: i32 = 2;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        let result = uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            &chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        );

        assert!(result.is_ok());

        // Test _start.wig output
        let path = PathBuf::from(&tempdir.path());
        let mut final_start_file_path = path.into_os_string().into_string().unwrap();
        final_start_file_path.push_str("/final/_start.bedGraph");
        let final_start_file_path = final_start_file_path.as_str();

        let file1 = File::open(final_start_file_path).unwrap();
        let file2 = File::open(test_output_path).unwrap();

        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            let line1 = lines1.next().transpose().unwrap();
            let line2 = lines2.next().transpose().unwrap();

            match (line1, line2) {
                (Some(line1), Some(line2)) => {
                    assert_eq!(line1, line2);
                }
                (None, None) => {
                    break; // Both files reached the end
                }
                _ => {
                    panic!("FILES ARE NOT EQUAL!!!")
                }
            }
        }

        // Test _core.wig output
        let path = PathBuf::from(&tempdir.path());
        let mut final_core_file_path = path.into_os_string().into_string().unwrap();
        final_core_file_path.push_str("/final/_core.bedGraph");
        let final_core_file_path = final_core_file_path.as_str();

        let file1 = File::open(final_core_file_path).unwrap();
        let file2 = File::open(core_test_output_path).unwrap();

        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        loop {
            let line1 = lines1.next().transpose().unwrap();
            let line2 = lines2.next().transpose().unwrap();

            match (line1, line2) {
                (Some(line1), Some(line2)) => {
                    assert_eq!(line1, line2);
                }
                (None, None) => {
                    break; // Both files reached the end
                }
                _ => {
                    panic!("FILES ARE NOT EQUAL!!!")
                }
            }
        }
    }

    #[rstest]
    fn test_process_narrowpeak(
        path_to_dummy_narrowpeak: &str,
    ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/tests/hg38.chrom.sizes");
        let chromsizerefpath = chromsizerefpath.as_str();
        let combinedbedpath = path_to_dummy_narrowpeak;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bw";
        let filetype = "narrowpeak";
        let num_threads = 2;
        let score = true;
        let stepsize = 1;
        let zoom = 2;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        )
        .expect("Uniwig main failed!");

        Ok(())
    }

    #[rstest]
    fn test_process_bed_to_bw(
        _path_to_dummy_bed_file: &str,
    ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let chromsizerefpath: String = format!("{}{}", path_to_crate, "/tests/hg38.chrom.sizes");
        let chromsizerefpath = chromsizerefpath.as_str();
        let combinedbedpath = _path_to_dummy_bed_file;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        let mut bwfileheader_path = path.into_os_string().into_string().unwrap();
        bwfileheader_path.push_str("/final/");
        let bwfileheader = bwfileheader_path.as_str();

        let smoothsize: i32 = 1;
        let output_type = "bw";
        let filetype = "bed";
        let num_threads = 2;
        let score = true;
        let stepsize = 1;
        let zoom = 1;
        let vec_count_type = vec!["start", "end", "core"];

        uniwig_main(
            vec_count_type,
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            bwfileheader,
            output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        )
        .expect("Uniwig main failed!");

        Ok(())
    }
}
