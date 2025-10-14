#![allow(non_snake_case)]
use gtars::bbcache::client::BBClient;
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
    "tests/data/tokenizers/peaks.bed"
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
    "tests/data/tokenizers/peaks.bed.gz"
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

#[fixture]
fn path_to_bed_gz_from_bb() -> &'static str {
    "tests/data/6b2e163a1d4319d99bd465c6c78a9741.bed.gz"
}

#[fixture]
fn bbid() -> &'static str {
    "6b2e163a1d4319d99bd465c6c78a9741"
}

#[fixture]
fn bsid() -> &'static str {
    "gse127562"
}

#[fixture]
fn path_to_bedset() -> &'static str {
    "tests/data/bedset"
}

mod tests {
    use super::*;
    use gtars::igd::create::{
        create_igd_f, gdata_t, igd_add, igd_saveT, igd_save_db, igd_t, parse_bed,
    };
    use gtars::igd::search::{
        get_file_info_tsv, get_igd_info, get_tsv_path, igd_search, igd_t_from_disk,
    };

    use gtars::uniwig::{uniwig_main, Chromosome};

    use gtars::uniwig::counting::{core_counts, start_end_counts};
    use gtars::uniwig::reading::{
        create_chrom_vec_default_score, create_chrom_vec_scores, parse_bedlike_file,
        read_bam_header, read_chromosome_sizes,
    };

    use gtars::uniwig::utils::npy_to_wig;
    use gtars::uniwig::writing::write_bw_files;

    // use gtars::bbcache::client::BBClient;

    use byteorder::{LittleEndian, ReadBytesExt};
    // use flate2::read::GzDecoder;
    use std::collections::HashMap;
    use std::collections::HashSet;
    use std::fs;
    use std::fs::{read_dir, OpenOptions};
    use std::io::{Seek, SeekFrom};

    

    // UNIWIG TESTS
    #[rstest]
    fn test_uniwig_parsed_bed_file(path_to_bed_file: &str) {
        let path = Path::new(path_to_bed_file);
        let file = File::open(path).unwrap();

        let mut reader = BufReader::new(file);
        let first_line = reader.by_ref().lines().next().unwrap().expect("expect");
        println!("{:?}", first_line);

        let result = parse_bedlike_file(&first_line);

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
    fn test_create_chrom_vec_default_score(path_to_bed_file: &str, path_to_bed_file_gzipped: &str) {
        let result1 = create_chrom_vec_default_score(path_to_bed_file);
        assert_eq!(result1.len(), 20);

        let result2 = create_chrom_vec_default_score(path_to_bed_file_gzipped);
        assert_eq!(result2.len(), 20);
    }

    #[rstest]
    fn test_create_chrom_vec_scores() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let path_to_narrow_peak = format!("{}{}", path_to_crate, "/tests/data/dummy.narrowPeak");
        let result1 = create_chrom_vec_scores(&path_to_narrow_peak);
        assert_eq!(result1.len(), 1);

        let path_to_narrow_peak_gzipped =
            format!("{}{}", path_to_crate, "/tests/data/dummy.narrowPeak.gz");

        let result2 = create_chrom_vec_scores(&path_to_narrow_peak_gzipped);
        assert_eq!(result2.len(), 1);
    }

    #[rstest]
    fn test_read_scored_core_counts() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let path_to_narrow_peak = format!("{}{}", path_to_crate, "/tests/data/dummy.narrowPeak");
        let chrom_sizes = read_chromosome_sizes(path_to_narrow_peak.as_str()).unwrap();
        let narrow_peak_vec: Vec<Chromosome> =
            create_chrom_vec_scores(path_to_narrow_peak.as_str());
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
    fn test_read_scored_starts_counts() {
        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let path_to_narrow_peak = format!("{}{}", path_to_crate, "/tests/data/dummy.narrowPeak");
        let chrom_sizes = read_chromosome_sizes(path_to_narrow_peak.as_str()).unwrap();
        let narrow_peak_vec: Vec<Chromosome> =
            create_chrom_vec_scores(path_to_narrow_peak.as_str());
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
        let chromosomes: Vec<Chromosome> =
            create_chrom_vec_default_score(path_to_sorted_small_bed_file);
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
    fn test_run_uniwig_main_directory_type() -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        // This test uses the bed file to determine chromsizes for speed
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!("{}{}", path_to_crate, "/tests/data/dir_of_files/dir_beds/");
        let combinedbedpath = tempbedpath.as_str();

        //let chromsizerefpath = combinedbedpath;

        let chromsizerefpath = format!(
            "{}{}",
            path_to_crate, "/tests/data/dir_of_files/dummy.chrom.sizes"
        );
        let chromsizerefpath = chromsizerefpath.as_str();
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        //let bwfileheader = "/home/drc/Downloads/gtars_uniwig_30june2025/output/";

        let smoothsize: i32 = 2;
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
    fn test_run_uniwig_main_directory_narrowpeaks_type(
    ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        // This test uses the bed file to determine chromsizes for speed
        let path_to_crate = env!("CARGO_MANIFEST_DIR");

        let tempbedpath = format!(
            "{}{}",
            path_to_crate, "/tests/data/dir_of_files/dir_narrowpeaks/"
        );
        let combinedbedpath = tempbedpath.as_str();

        //let chromsizerefpath = combinedbedpath;

        let chromsizerefpath = format!(
            "{}{}",
            path_to_crate, "/tests/data/dir_of_files/dummy.chrom.sizes"
        );
        let chromsizerefpath = chromsizerefpath.as_str();
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let bwfileheader_path = path.into_os_string().into_string().unwrap();
        let bwfileheader = bwfileheader_path.as_str();

        //let bwfileheader = "/home/drc/Downloads/gtars_uniwig_30june2025/output/";

        let smoothsize: i32 = 2;
        let output_type = "wig";
        let filetype = "narrowpeak";
        let num_threads = 6;
        let score = true;
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

        let original_bedgraph_path = format!("{}/tests/data/test1.bedGraph", path_to_crate);
        let chrom_sizes_path = format!("{}/tests/data/dummy.chrom.sizes", path_to_crate);

        let temp_dir = tempfile::tempdir().unwrap();
        let temp_bedgraph_path = temp_dir.path().join("test1.bedGraph");

        fs::copy(&original_bedgraph_path, &temp_bedgraph_path)
            .expect("Failed to copy .bedGraph file to temporary directory");

        let num_threads = 2;
        let zoom = 0;

        write_bw_files(
            temp_bedgraph_path.to_str().expect("Invalid temp path"), // Use the path in the temp directory
            chrom_sizes_path.as_str(),
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
        let chromsizerefpath: String =
            format!("{}{}", path_to_crate, "/tests/data/dummy.chrom.sizes");
        let chromsizerefpath = chromsizerefpath.as_str();
        let combinedbedpath = path_to_dummy_narrowpeak;

        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());

        // For some reason, you cannot chain .as_string() to .unwrap() and must create a new line.
        let mut bwfileheader_path = path.into_os_string().into_string().unwrap();
        bwfileheader_path.push_str("/final/");
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
        let chromsizerefpath: String =
            format!("{}{}", path_to_crate, "/tests/data/dummy.chrom.sizes");
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

    #[rstest]
    fn test_npy_to_wig(
        _path_to_dummy_bed_file: &str,
        _path_to_dummy_chromsizes: &str,
    ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        let chromsizerefpath = _path_to_dummy_chromsizes;
        let combinedbedpath = _path_to_dummy_bed_file;
        let tempdir = tempfile::tempdir()?; // use `?` for idiomatic error handling
        let path = PathBuf::from(tempdir.path());

        let smoothsize = 1;
        let wig_output_type = "wig";
        let npy_output_type = "npy";
        let filetype = "bed";
        let num_threads = 6;
        let score = false;
        let stepsize = 1;
        let zoom = 0;
        let vec_count_type = vec!["start", "end", "core"];

        // Generate npy output
        let npyfileheader_path = format!("{}/npyfinal/", path.display());
        let npyfileheader = npyfileheader_path.as_str();

        let _ = uniwig_main(
            vec_count_type.clone(),
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            npyfileheader,
            npy_output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        );

        // Generate wig output
        let wigfileheader_path = format!("{}/wigfinal/", path.display());
        let wigfileheader = wigfileheader_path.as_str();

        let _ = uniwig_main(
            vec_count_type.clone(),
            smoothsize,
            combinedbedpath,
            chromsizerefpath,
            wigfileheader,
            wig_output_type,
            filetype,
            num_threads,
            score,
            stepsize,
            zoom,
            false,
            true,
            1.0,
        );

        // Run npy_to_wig
        let genwigfileheader_path = format!("{}/genwigfinal/", path.display());
        let genwigfileheader = genwigfileheader_path.as_str();

        let npy_header_path = Path::new(npyfileheader);
        let gen_wig_header_path = Path::new(genwigfileheader);
        let _ = npy_to_wig(npy_header_path, gen_wig_header_path);

        // Compare output directories
        let ref_wig_header_path = Path::new(wigfileheader);

        let mut files1: Vec<_> = read_dir(ref_wig_header_path)?
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();
        let mut files2: Vec<_> = read_dir(gen_wig_header_path)?
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();

        files1.sort();
        files2.sort();

        assert_eq!(files1, files2, "Directory file names differ");

        for file_name in files1 {
            let path1 = gen_wig_header_path.join(&file_name);
            let path2 = ref_wig_header_path.join(&file_name);

            let mut f1 = File::open(&path1)?;
            let mut f2 = File::open(&path2)?;

            let mut buf1 = Vec::new();
            let mut buf2 = Vec::new();

            f1.read_to_end(&mut buf1)?;
            f2.read_to_end(&mut buf2)?;

            assert_eq!(
                buf1,
                buf2,
                "File contents differ between:\n  {}\nand\n  {}",
                path1.display(),
                path2.display()
            );
        }

        Ok(())
    }

    #[rstest]
    fn test_bbcache_local(
        _path_to_bed_gz_from_bb: &str,
        _bbid: &str,
        _path_to_bedset: &str,
    ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
        fn cleaned_subfolders(subfolder: PathBuf) {
            let subdirs: Vec<_> = read_dir(&subfolder)
                .unwrap_or_else(|e| {
                    panic!("Failed to read directory {}: {}", subfolder.display(), e)
                })
                .filter_map(Result::ok)
                .filter(|entry| entry.path().is_dir())
                .collect();

            // Assert no subdirectories exist
            assert!(
                subdirs.is_empty(),
                "Subfolders found in {}: {:?}",
                subfolder.display(),
                subdirs.iter().map(|e| e.path()).collect::<Vec<_>>()
            );
        }
        let tempdir = tempfile::tempdir()?;
        let cache_folder = PathBuf::from(tempdir.path());

        let mut bbc =
            BBClient::new(Some(cache_folder.clone()), None).expect("Failed to create BBClient");

        let bed_id = bbc
            .add_local_bed_to_cache(PathBuf::from(_path_to_bed_gz_from_bb), Some(false))
            .unwrap();
        assert_eq!(&bed_id, _bbid);

        let bedset_id = bbc
            .add_local_folder_as_bedset(PathBuf::from(_path_to_bedset))
            .unwrap();
        assert!(bbc.seek(&bedset_id).is_ok());

        bbc.remove(&bedset_id)
            .expect("Failed to remove bedset file and its bed files");
        let bedset_subfolder = cache_folder.join("bedsets");
        cleaned_subfolders(bedset_subfolder);

        bbc.remove(_bbid).expect("Failed to remove cached bed file");
        let bedfile_subfolder = cache_folder.join("bedfiles");
        cleaned_subfolders(bedfile_subfolder);
        Ok(())
    }

    // This test should be mocked and not use bedbase. Commented for now.
    // #[rstest]
    // fn test_bbcache_bedbase(
    //     _path_to_bed_gz_from_bb: &str,
    //     _bbid: &str,
    //     _bsid: &str,
    // ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
    //     fn read_gzip_file(path: impl AsRef<std::path::Path>) -> Vec<u8> {
    //         let file = File::open(path).expect("Failed to open file");
    //         let mut decoder = GzDecoder::new(BufReader::new(file));
    //         let mut contents = Vec::new();
    //         decoder
    //             .read_to_end(&mut contents)
    //             .expect("Failed to read decompressed contents");
    //         contents
    //     }
    //
    //     let tempdir = tempfile::tempdir()?;
    //     let cache_folder = PathBuf::from(tempdir.path());
    //
    //     let mut bbc =
    //         BBClient::new(Some(cache_folder.clone()), None).expect("Failed to create BBClient");
    //
    //     let _rs = bbc.load_bed(_bbid).expect("Failed to load bed file");
    //
    //     assert!(bbc.seek(_bbid).is_ok());
    //
    //     let cached_bed_path = bbc.seek(_bbid).expect("Failed to seek cached bed file");
    //     let cached_content = read_gzip_file(&cached_bed_path);
    //     let comparison_content = read_gzip_file(_path_to_bed_gz_from_bb);
    //     assert_eq!(
    //         cached_content, comparison_content,
    //         "Cached content does not match the original content"
    //     );

    //     let bedset = bbc.load_bedset(_bsid).unwrap();
    //     assert!(bbc.seek(_bsid).is_ok());
    //     for rs in bedset.region_sets {
    //         let bed_id = rs.identifier();
    //         assert!(bbc.seek(&bed_id.clone()).is_ok());
    //         let bed_in_set = bbc.load_bed(&bed_id).unwrap();
    //         assert_eq!(bed_id, bed_in_set.identifier());
    //     }
    //     Ok(())
    // }
}
