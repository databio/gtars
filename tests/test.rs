#![allow(non_snake_case)]
use gtars::bbcache::client::BBClient;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use rstest::*;

// #[fixture]
// fn path_to_data() -> &'static str {
//     "tests/data"
// }

// #[fixture]
// fn path_to_bed_file() -> &'static str {
//     "tests/data/tokenizers/peaks.bed"
// }

// #[fixture]
// fn path_to_sorted_small_bed_file() -> &'static str {
//     "tests/data/test_sorted_small.bed"
// }

// #[fixture]
// fn path_to_small_bam_file() -> &'static str {
//     "tests/data/test_chr22_small.bam"
//     //"/home/drc/Downloads/bam files for rust test/test1_sort_dedup.bam"
// }

// #[fixture]
// fn path_to_chrom_sizes_file() -> &'static str {
//     "tests/hg38.chrom.sizes"
// }

// #[fixture]
// fn path_to_bed_file_gzipped() -> &'static str {
//     "tests/data/tokenizers/peaks.bed.gz"
// }

// #[fixture]
// fn path_to_dummy_bed_file() -> &'static str {
//     "tests/data/dummy.bed"
// }

// #[fixture]
// fn path_to_dummy_chromsizes() -> &'static str {
//     "tests/data/dummy.chrom.sizes"
// }

// #[fixture]
// fn path_to_dummy_narrowpeak() -> &'static str {
//     "tests/data/dummy.narrowPeak"
// }

// #[fixture]
// fn path_to_start_wig_output() -> &'static str {
//     "tests/data/out/_start.wig"
// }

// #[fixture]
// fn path_to_core_wig_output() -> &'static str {
//     "tests/data/out/_core.wig"
// }

// #[fixture]
// fn path_to_start_bedgraph_output() -> &'static str {
//     "tests/data/out/_start.bedGraph"
// }

// #[fixture]
// fn path_to_core_bedgraph_output() -> &'static str {
//     "tests/data/out/_core.bedGraph"
// }

// #[fixture]
// fn path_to_bed_gz_from_bb() -> &'static str {
//     "tests/data/6b2e163a1d4319d99bd465c6c78a9741.bed.gz"
// }

// #[fixture]
// fn bbid() -> &'static str {
//     "6b2e163a1d4319d99bd465c6c78a9741"
// }

// #[fixture]
// fn bsid() -> &'static str {
//     "gse127562"
// }

// #[fixture]
// fn path_to_bedset() -> &'static str {
//     "tests/data/bedset"
// }

// mod tests {
   

//     #[rstest]
//     fn test_npy_to_wig(
//         _path_to_dummy_bed_file: &str,
//         _path_to_dummy_chromsizes: &str,
//     ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
//         let chromsizerefpath = _path_to_dummy_chromsizes;
//         let combinedbedpath = _path_to_dummy_bed_file;
//         let tempdir = tempfile::tempdir()?; // use `?` for idiomatic error handling
//         let path = PathBuf::from(tempdir.path());

//         let smoothsize = 1;
//         let wig_output_type = "wig";
//         let npy_output_type = "npy";
//         let filetype = "bed";
//         let num_threads = 6;
//         let score = false;
//         let stepsize = 1;
//         let zoom = 0;
//         let vec_count_type = vec!["start", "end", "core"];

//         // Generate npy output
//         let npyfileheader_path = format!("{}/npyfinal/", path.display());
//         let npyfileheader = npyfileheader_path.as_str();

//         let _ = uniwig_main(
//             vec_count_type.clone(),
//             smoothsize,
//             combinedbedpath,
//             chromsizerefpath,
//             npyfileheader,
//             npy_output_type,
//             filetype,
//             num_threads,
//             score,
//             stepsize,
//             zoom,
//             false,
//             true,
//             1.0,
//         );

//         // Generate wig output
//         let wigfileheader_path = format!("{}/wigfinal/", path.display());
//         let wigfileheader = wigfileheader_path.as_str();

//         let _ = uniwig_main(
//             vec_count_type.clone(),
//             smoothsize,
//             combinedbedpath,
//             chromsizerefpath,
//             wigfileheader,
//             wig_output_type,
//             filetype,
//             num_threads,
//             score,
//             stepsize,
//             zoom,
//             false,
//             true,
//             1.0,
//         );

//         // Run npy_to_wig
//         let genwigfileheader_path = format!("{}/genwigfinal/", path.display());
//         let genwigfileheader = genwigfileheader_path.as_str();

//         let npy_header_path = Path::new(npyfileheader);
//         let gen_wig_header_path = Path::new(genwigfileheader);
//         let _ = npy_to_wig(npy_header_path, gen_wig_header_path);

//         // Compare output directories
//         let ref_wig_header_path = Path::new(wigfileheader);

//         let mut files1: Vec<_> = read_dir(ref_wig_header_path)?
//             .map(|entry| entry.unwrap().file_name().into_string().unwrap())
//             .collect();
//         let mut files2: Vec<_> = read_dir(gen_wig_header_path)?
//             .map(|entry| entry.unwrap().file_name().into_string().unwrap())
//             .collect();

//         files1.sort();
//         files2.sort();

//         assert_eq!(files1, files2, "Directory file names differ");

//         for file_name in files1 {
//             let path1 = gen_wig_header_path.join(&file_name);
//             let path2 = ref_wig_header_path.join(&file_name);

//             let mut f1 = File::open(&path1)?;
//             let mut f2 = File::open(&path2)?;

//             let mut buf1 = Vec::new();
//             let mut buf2 = Vec::new();

//             f1.read_to_end(&mut buf1)?;
//             f2.read_to_end(&mut buf2)?;

//             assert_eq!(
//                 buf1,
//                 buf2,
//                 "File contents differ between:\n  {}\nand\n  {}",
//                 path1.display(),
//                 path2.display()
//             );
//         }

//         Ok(())
//     }

//     #[rstest]
//     fn test_bbcache_local(
//         _path_to_bed_gz_from_bb: &str,
//         _bbid: &str,
//         _path_to_bedset: &str,
//     ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
//         fn cleaned_subfolders(subfolder: PathBuf) {
//             let subdirs: Vec<_> = read_dir(&subfolder)
//                 .unwrap_or_else(|e| {
//                     panic!("Failed to read directory {}: {}", subfolder.display(), e)
//                 })
//                 .filter_map(Result::ok)
//                 .filter(|entry| entry.path().is_dir())
//                 .collect();

//             // Assert no subdirectories exist
//             assert!(
//                 subdirs.is_empty(),
//                 "Subfolders found in {}: {:?}",
//                 subfolder.display(),
//                 subdirs.iter().map(|e| e.path()).collect::<Vec<_>>()
//             );
//         }
//         let tempdir = tempfile::tempdir()?;
//         let cache_folder = PathBuf::from(tempdir.path());

//         let mut bbc =
//             BBClient::new(Some(cache_folder.clone()), None).expect("Failed to create BBClient");

//         let bed_id = bbc
//             .add_local_bed_to_cache(PathBuf::from(_path_to_bed_gz_from_bb), Some(false))
//             .unwrap();
//         assert_eq!(&bed_id, _bbid);

//         let bedset_id = bbc
//             .add_local_folder_as_bedset(PathBuf::from(_path_to_bedset))
//             .unwrap();
//         assert!(bbc.seek(&bedset_id).is_ok());

//         bbc.remove(&bedset_id)
//             .expect("Failed to remove bedset file and its bed files");
//         let bedset_subfolder = cache_folder.join("bedsets");
//         cleaned_subfolders(bedset_subfolder);

//         bbc.remove(_bbid).expect("Failed to remove cached bed file");
//         let bedfile_subfolder = cache_folder.join("bedfiles");
//         cleaned_subfolders(bedfile_subfolder);
//         Ok(())
//     }

//     // This test should be mocked and not use bedbase. Commented for now.
//     // #[rstest]
//     // fn test_bbcache_bedbase(
//     //     _path_to_bed_gz_from_bb: &str,
//     //     _bbid: &str,
//     //     _bsid: &str,
//     // ) -> Result<(), Box<(dyn std::error::Error + 'static)>> {
//     //     fn read_gzip_file(path: impl AsRef<std::path::Path>) -> Vec<u8> {
//     //         let file = File::open(path).expect("Failed to open file");
//     //         let mut decoder = GzDecoder::new(BufReader::new(file));
//     //         let mut contents = Vec::new();
//     //         decoder
//     //             .read_to_end(&mut contents)
//     //             .expect("Failed to read decompressed contents");
//     //         contents
//     //     }
//     //
//     //     let tempdir = tempfile::tempdir()?;
//     //     let cache_folder = PathBuf::from(tempdir.path());
//     //
//     //     let mut bbc =
//     //         BBClient::new(Some(cache_folder.clone()), None).expect("Failed to create BBClient");
//     //
//     //     let _rs = bbc.load_bed(_bbid).expect("Failed to load bed file");
//     //
//     //     assert!(bbc.seek(_bbid).is_ok());
//     //
//     //     let cached_bed_path = bbc.seek(_bbid).expect("Failed to seek cached bed file");
//     //     let cached_content = read_gzip_file(&cached_bed_path);
//     //     let comparison_content = read_gzip_file(_path_to_bed_gz_from_bb);
//     //     assert_eq!(
//     //         cached_content, comparison_content,
//     //         "Cached content does not match the original content"
//     //     );

//     //     let bedset = bbc.load_bedset(_bsid).unwrap();
//     //     assert!(bbc.seek(_bsid).is_ok());
//     //     for rs in bedset.region_sets {
//     //         let bed_id = rs.identifier();
//     //         assert!(bbc.seek(&bed_id.clone()).is_ok());
//     //         let bed_in_set = bbc.load_bed(&bed_id).unwrap();
//     //         assert_eq!(bed_id, bed_in_set.identifier());
//     //     }
//     //     Ok(())
//     // }
// }
