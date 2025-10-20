#![allow(nonstandard_style)]

pub mod create;
pub mod search;

#[cfg(feature = "bloom")]
pub mod igdbloom; // Bloom filter work is only available as library functions NOT CLI

pub mod consts {
    pub const IGD_CMD: &str = "igd";
    pub const IGD_CREATE: &str = "create";
    pub const IGD_SEARCH: &str = "search";
}

#[cfg(test)]
mod tests {

    use rstest::rstest;

    use crate::create::{create_igd_f, gdata_t, igd_add, igd_save_db, igd_saveT, igd_t, parse_bed};
    // Import get_igd_info if it is defined in another module
    use crate::search::{
        get_file_info_tsv, get_igd_info, get_tsv_path, igd_search, igd_t_from_disk,
    };

    use std::collections::HashMap;
    use std::path::{Path, PathBuf};
    use tempfile::NamedTempFile;

    use byteorder::{LittleEndian, ReadBytesExt};
    use std::collections::HashSet;
    use std::fs::OpenOptions;
    use std::io::{BufReader, Read, Seek, SeekFrom};

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
        let testfilelists = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/igd_file_list_01/")
            .to_string_lossy()
            .to_string();

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
        db_path_unwrapped.push('/');
        let db_output_path = db_path_unwrapped.clone();

        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let testfilelists = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/igd_file_list_01/")
            .to_string_lossy()
            .to_string();

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
                for _j in 0..tmpi {
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
                        idx,
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
    fn test_igd_create_removes_temp_dir() {
        let tempdir = tempfile::tempdir().unwrap();
        let path = PathBuf::from(&tempdir.path());
        let mut db_path_unwrapped = path.into_os_string().into_string().unwrap();
        db_path_unwrapped.push('/');
        let db_output_path = db_path_unwrapped.clone();

        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let testfilelists = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests/data/igd_file_list_01/")
            .to_string_lossy()
            .to_string();

        let demo_name = String::from("demo");

        let _igd_saved = create_igd_f(&db_output_path, &testfilelists, &demo_name);

        let temp_folder = format!("{}{}", db_output_path, "data0/");
        let path = Path::new(&temp_folder);

        // Assert path does not exist
        assert!(!path.exists());
    }

    #[rstest]
    #[case(
        "/../tests/data/igd_file_list_01/",
        "/../tests/data/igd_query_files/query1.bed",
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
        db_path_unwrapped.push('/');
        let db_output_path = db_path_unwrapped.clone();

        let path_to_crate = env!("CARGO_MANIFEST_DIR");
        let testfilelists = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests")
            .join(input.trim_start_matches('/'))
            .to_string_lossy()
            .to_string();

        let demo_name = String::from("demo");

        let _igd_saved = create_igd_f(&db_output_path, &testfilelists, &demo_name);

        println!("dboutput_path {}", db_output_path);

        db_path_unwrapped.push_str("/demo.igd");

        let queryfile = PathBuf::from(path_to_crate)
            .parent()
            .unwrap()
            .join("tests")
            .join(query_file.trim_start_matches('/'))
            .to_string_lossy()
            .to_string();
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
}
