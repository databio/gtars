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

// }
