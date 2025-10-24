pub mod client;
pub mod consts;
pub mod utils;

#[cfg(test)]
mod tests {
    use super::client::BBClient;
    use rstest::{fixture, rstest};
    use std::fs::read_dir;
    use std::path::PathBuf;

    #[fixture]
    fn path_to_bed_gz_from_bb() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/6b2e163a1d4319d99bd465c6c78a9741.bed.gz")
    }

    #[fixture]
    fn bbid() -> PathBuf {
        "6b2e163a1d4319d99bd465c6c78a9741".into()
    }

    #[fixture]
    fn path_to_bedset() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data/bedset")
    }

    #[rstest]
    fn test_bbcache_local(
        path_to_bed_gz_from_bb: PathBuf,
        bbid: PathBuf,
        path_to_bedset: PathBuf,
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

        let mut bbc = BBClient::builder()
            .with_cache_folder(cache_folder.clone())
            .finish()?;

        let bed_id = bbc
            .add_local_bed_to_cache(path_to_bed_gz_from_bb, Some(false))
            .unwrap();
        assert_eq!(&bed_id, &bbid.to_string_lossy());

        let bedset_id = bbc.add_local_folder_as_bedset(path_to_bedset).unwrap();
        assert!(bbc.seek(&bedset_id).is_ok());

        bbc.remove(&bedset_id)
            .expect("Failed to remove bedset file and its bed files");
        let bedset_subfolder = cache_folder.join("bedsets");
        cleaned_subfolders(bedset_subfolder);

        bbc.remove(&bbid.to_string_lossy())
            .expect("Failed to remove cached bed file");
        let bedfile_subfolder = cache_folder.join("bedfiles");
        cleaned_subfolders(bedfile_subfolder);
        Ok(())
    }
}
