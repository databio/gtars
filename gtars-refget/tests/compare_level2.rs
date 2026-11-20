//! Tests for `compare_level2_json`, the store-free comparison of two level-2
//! seqcol bodies.
//!
//! The fixtures in `tests/data/comparison/` are copied from the refget
//! (Python) repo's API compliance suite: `collections/*.fa.json` are level-2
//! bodies and `expected/compare_<a>_<b>.json` the spec comparison results.
//!
//! One deliberate edit: where only one element is shared (`name_length_pairs`
//! in the four `pair_swap` fixtures), `a_and_b_same_order` is `null`, not
//! refget Python's `true`. The seqcol spec says order is undefined (null) with
//! fewer than 2 overlapping elements; refget Python's `_compare_elements`
//! checks `A_count + B_count < 1` instead.

use gtars_refget::compare_level2_json;
use serde_json::{Value, json};
use std::fs;
use std::path::Path;

const DATA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/comparison");

fn read_json(path: &Path) -> Value {
    serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap()
}

#[test]
fn matches_refget_compliance_fixtures() {
    let mut checked = 0;
    for entry in fs::read_dir(Path::new(DATA).join("expected")).unwrap() {
        let path = entry.unwrap().path();
        let fname = path.file_name().unwrap().to_str().unwrap();
        // compare_<a>.fa_<b>.fa.json
        let stem = fname
            .strip_prefix("compare_")
            .unwrap()
            .strip_suffix(".json")
            .unwrap();
        let (a_name, b_name) = stem.split_once(".fa_").unwrap();
        let a = read_json(&Path::new(DATA).join(format!("collections/{a_name}.fa.json")));
        let b = read_json(&Path::new(DATA).join(format!("collections/{b_name}.json")));
        let expected = read_json(&path);
        let digest = |side: &str| expected["digests"][side].as_str().map(String::from);

        let result = compare_level2_json(&a, &b, digest("a"), digest("b")).unwrap();
        assert_eq!(
            serde_json::to_value(&result).unwrap(),
            expected,
            "fixture {fname}"
        );
        checked += 1;
    }
    assert_eq!(checked, 15);
}

#[test]
fn self_comparison_is_all_same_order() {
    let a = read_json(&Path::new(DATA).join("collections/base.fa.json"));
    let result = compare_level2_json(&a, &a, None, None).unwrap();
    let v = serde_json::to_value(&result).unwrap();
    assert_eq!(v["digests"], json!({"a": null, "b": null}));
    assert!(result.attributes.a_only.is_empty() && result.attributes.b_only.is_empty());
    for attr in &result.attributes.a_and_b {
        assert_eq!(
            result.array_elements.a_and_b_same_order[attr],
            Some(true),
            "{attr}"
        );
    }
}

#[test]
fn one_shared_element_has_null_order() {
    let a = json!({"names": ["chr1", "chr2"], "lengths": [10, 20]});
    let b = json!({"names": ["chr1", "chrM"], "lengths": [10, 30]});
    let v = serde_json::to_value(compare_level2_json(&a, &b, None, None).unwrap()).unwrap();
    assert_eq!(
        v["array_elements"]["a_and_b_count"],
        json!({"names": 1, "lengths": 1})
    );
    assert_eq!(
        v["array_elements"]["a_and_b_same_order"],
        json!({"names": null, "lengths": null})
    );
}

#[test]
fn unbalanced_duplicates_have_null_order() {
    // A has 4 twice, B once: overlap counts the min, order is undefined.
    let a = json!({"lengths": [8, 4, 4]});
    let b = json!({"lengths": [8, 4]});
    let r = compare_level2_json(&a, &b, None, None).unwrap();
    assert_eq!(r.array_elements.a_and_b_count["lengths"], 2);
    assert_eq!(r.array_elements.a_and_b_same_order["lengths"], None);

    // Balanced duplicates still get an order verdict.
    let c = json!({"lengths": [4, 8, 4]});
    let r = compare_level2_json(&a, &c, None, None).unwrap();
    assert_eq!(r.array_elements.a_and_b_count["lengths"], 3);
    assert_eq!(r.array_elements.a_and_b_same_order["lengths"], Some(false));
}

#[test]
fn object_elements_compare_by_content() {
    // Same pairs, different key order inside the objects, swapped order.
    let a = json!({"name_length_pairs": [
        {"name": "chr1", "length": 10},
        {"name": "chr2", "length": 20},
    ]});
    let b = json!({"name_length_pairs": [
        {"length": 20, "name": "chr2"},
        {"length": 10, "name": "chr1"},
        {"length": 5, "name": "chrM"},
    ]});
    let r = compare_level2_json(&a, &b, Some("A".into()), None).unwrap();
    assert_eq!(r.array_elements.a_and_b_count["name_length_pairs"], 2);
    assert_eq!(
        r.array_elements.a_and_b_same_order["name_length_pairs"],
        Some(false)
    );
    assert_eq!(r.array_elements.b_count["name_length_pairs"], 3);
    assert_eq!(r.digests.a.as_deref(), Some("A"));
    assert_eq!(r.digests.b, None);
}

#[test]
fn string_and_number_are_different_elements() {
    // Matches refget Python, where "8" != 8.
    let a = json!({"lengths": [8, 4]});
    let b = json!({"lengths": ["8", "4"]});
    let r = compare_level2_json(&a, &b, None, None).unwrap();
    assert_eq!(r.array_elements.a_and_b_count["lengths"], 0);
}

#[test]
fn attribute_sets_are_split() {
    let a = json!({"names": ["x"], "lengths": [1]});
    let b = json!({"names": ["x"], "sequences": ["SQ.a"]});
    let r = compare_level2_json(&a, &b, None, None).unwrap();
    assert_eq!(r.attributes.a_only, vec!["lengths"]);
    assert_eq!(r.attributes.b_only, vec!["sequences"]);
    assert_eq!(r.attributes.a_and_b, vec!["names"]);
    assert!(!r.array_elements.a_and_b_count.contains_key("lengths"));
}

#[test]
fn rejects_non_object_and_non_array() {
    assert!(compare_level2_json(&json!([1]), &json!({}), None, None).is_err());
    let err = compare_level2_json(&json!({"names": "chr1"}), &json!({}), None, None).unwrap_err();
    assert!(err.to_string().contains("names"));
}
