use vocab::create_count_map;
use rstest::{fixture, rstest};

#[fixture]
fn data_path() -> &'static str {
    "tests/data/files"
}

#[fixture]
fn universe_path() -> &'static str {
    "tests/data/universe.bed"
}

#[rstest]
fn total_overlaps(data_path: &'static str, universe_path: &'static str) {
    let cnt_map = create_count_map(data_path, universe_path).unwrap();
    // all of our data only intersects with three regions in the universe
    assert_eq!(cnt_map.len(), 3);
}

#[rstest]
fn overlap_count(data_path: &'static str, universe_path: &'static str) {
    let cnt_map = create_count_map(data_path, universe_path).unwrap();

    // two duplicate files, so each counted overlap should be 2
    for (_, cnt) in cnt_map {
        assert_eq!(cnt, 2);
    }
}
