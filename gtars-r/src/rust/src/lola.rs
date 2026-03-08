use extendr_api::prelude::*;

use gtars_core::models::RegionSet;
use gtars_igd::igd::Igd;
use gtars_lola::database::{RegionDB, RegionSetAnno};
use gtars_lola::enrichment::run_lola;
use gtars_lola::models::{Direction, LolaConfig};
use gtars_lola::universe;

// =========================================================================
// Helper macros
// =========================================================================

macro_rules! with_regiondb {
    ($ptr:expr, $db:ident, $body:block) => {{
        let ext_ptr = <ExternalPtr<RegionDB>>::try_from($ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid RegionDB pointer".into()))?;
        let $db = &*ext_ptr;
        $body
    }};
}

macro_rules! with_regionset {
    ($ptr:expr, $rs:ident, $body:block) => {{
        let ext_ptr = <ExternalPtr<RegionSet>>::try_from($ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer".into()))?;
        let $rs = &*ext_ptr;
        $body
    }};
}

// =========================================================================
// Helper functions
// =========================================================================

/// Extract a list of RegionSet pointers from an R list.
fn extract_region_sets(user_sets: List) -> extendr_api::Result<Vec<RegionSet>> {
    let mut sets = Vec::new();
    for item in user_sets.iter() {
        let (_, val) = item;
        let ext_ptr = <ExternalPtr<RegionSet>>::try_from(val)
            .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer in list".into()))?;
        sets.push((*ext_ptr).clone());
    }
    Ok(sets)
}

/// Convert LOLA results to an R list (data.frame-like structure).
fn results_to_list(results: &[gtars_lola::models::LolaResult]) -> List {
    let n = results.len();
    let mut user_set = Vec::with_capacity(n);
    let mut db_set = Vec::with_capacity(n);
    let mut p_value_log = Vec::with_capacity(n);
    let mut odds_ratio = Vec::with_capacity(n);
    let mut support = Vec::with_capacity(n);
    let mut rnk_pv = Vec::with_capacity(n);
    let mut rnk_or = Vec::with_capacity(n);
    let mut rnk_sup = Vec::with_capacity(n);
    let mut max_rnk = Vec::with_capacity(n);
    let mut mean_rnk = Vec::with_capacity(n);
    let mut b_vec = Vec::with_capacity(n);
    let mut c_vec = Vec::with_capacity(n);
    let mut d_vec = Vec::with_capacity(n);
    let mut filename = Vec::with_capacity(n);

    let mut q_value: Vec<Option<f64>> = Vec::with_capacity(n);

    for r in results {
        user_set.push((r.user_set + 1) as i32); // 1-based for R
        db_set.push((r.db_set + 1) as i32);
        p_value_log.push(r.p_value_log);
        odds_ratio.push(r.odds_ratio);
        support.push(r.support as i32);
        rnk_pv.push(r.rnk_pv as i32);
        rnk_or.push(r.rnk_or as i32);
        rnk_sup.push(r.rnk_sup as i32);
        max_rnk.push(r.max_rnk as i32);
        mean_rnk.push(r.mean_rnk);
        b_vec.push(r.b as i32);
        c_vec.push(r.c as i32);
        d_vec.push(r.d as i32);
        filename.push(r.filename.clone());
        q_value.push(r.q_value);
    }

    // Convert Option<f64> to Rfloat (NA for None)
    let q_value_r: Vec<Rfloat> = q_value
        .iter()
        .map(|q| match q {
            Some(v) => Rfloat::from(*v),
            None => Rfloat::na(),
        })
        .collect();

    list!(
        userSet = user_set,
        dbSet = db_set,
        pValueLog = p_value_log,
        oddsRatio = odds_ratio,
        support = support,
        rnkPV = rnk_pv,
        rnkOR = rnk_or,
        rnkSup = rnk_sup,
        maxRnk = max_rnk,
        meanRnk = mean_rnk,
        b = b_vec,
        c = c_vec,
        d = d_vec,
        filename = filename,
        qValue = q_value_r
    )
}

// =========================================================================
// Database loading
// =========================================================================

/// Load a LOLA region database from a folder.
/// @export
/// @param db_location Path to the LOLA database folder
/// @param collections Optional character vector of collection names to load
/// @param limit Optional integer limit on files per collection
#[extendr(r_name = "load_region_db")]
pub fn r_load_region_db(
    db_location: &str,
    collections: Nullable<Vec<String>>,
    limit: Nullable<i32>,
) -> extendr_api::Result<Robj> {
    let path = std::path::Path::new(db_location);

    let coll_filter: Option<Vec<String>> = match collections {
        Nullable::NotNull(c) => Some(c),
        Nullable::Null => None,
    };

    let coll_refs: Option<Vec<&str>> = coll_filter
        .as_ref()
        .map(|v| v.iter().map(|s| s.as_str()).collect());

    let lim = match limit {
        Nullable::NotNull(l) => Some(l as usize),
        Nullable::Null => None,
    };

    let db = RegionDB::from_lola_folder(path, coll_refs.as_deref(), lim)
        .map_err(|e| extendr_api::Error::Other(format!("Failed to load RegionDB: {}", e)))?;

    Ok(ExternalPtr::new(db).into())
}

/// Load a region database from BED file paths and optional metadata.
/// @export
/// @param bed_files Character vector of paths to BED files
/// @param filenames Optional character vector of display names for each file
#[extendr(r_name = "load_region_db_from_beds")]
pub fn r_load_region_db_from_beds(
    bed_files: Vec<String>,
    filenames: Nullable<Vec<String>>,
) -> extendr_api::Result<Robj> {
    let names: Vec<String> = match filenames {
        Nullable::NotNull(n) => n,
        Nullable::Null => bed_files
            .iter()
            .map(|p| {
                std::path::Path::new(p)
                    .file_name()
                    .map(|f| f.to_string_lossy().into_owned())
                    .unwrap_or_else(|| p.clone())
            })
            .collect(),
    };

    let mut region_sets = Vec::new();
    let mut igd_inputs = Vec::new();
    let mut region_anno = Vec::new();

    for (i, bed_path) in bed_files.iter().enumerate() {
        let path = std::path::Path::new(bed_path);
        let rs = RegionSet::try_from(path)
            .map_err(|e| extendr_api::Error::Other(format!("Failed to read {}: {}", bed_path, e)))?;

        let regions: Vec<(String, i32, i32)> = rs
            .regions
            .iter()
            .map(|r| (r.chr.clone(), r.start as i32, r.end as i32))
            .collect();

        let name = names.get(i).cloned().unwrap_or_default();
        igd_inputs.push((name.clone(), regions));
        region_sets.push(rs);
        region_anno.push(RegionSetAnno {
            filename: name,
            ..Default::default()
        });
    }

    let igd = Igd::from_region_sets(igd_inputs);
    let db = RegionDB::from_igd_with_regions(igd, region_sets, region_anno);

    Ok(ExternalPtr::new(db).into())
}

// Note: merge_region_dbs is not exposed to R because ExternalPtr doesn't
// support taking ownership (needed for RegionDB::merge). Users should
// reload from disk or merge at the R level.

/// List region set filenames in a database.
/// @export
/// @param db ExternalPtr to RegionDB
/// @param collections Optional character vector filter
#[extendr(r_name = "list_region_sets")]
pub fn r_list_region_sets(
    db: Robj,
    collections: Nullable<Vec<String>>,
) -> extendr_api::Result<Vec<String>> {
    with_regiondb!(db, db_ref, {
        let coll_filter: Option<Vec<String>> = match collections {
            Nullable::NotNull(c) => Some(c),
            Nullable::Null => None,
        };
        let coll_refs: Option<Vec<&str>> = coll_filter
            .as_ref()
            .map(|v| v.iter().map(|s| s.as_str()).collect());
        Ok(db_ref.list_region_sets(coll_refs.as_deref()))
    })
}

// =========================================================================
// Core LOLA
// =========================================================================

/// Run LOLA enrichment analysis.
/// @export
/// @param user_sets_list R list of RegionSet external pointers
/// @param universe_ptr ExternalPtr to universe RegionSet
/// @param db_ptr ExternalPtr to RegionDB
/// @param min_overlap Minimum base-pair overlap (default 1)
/// @param direction "enrichment" or "depletion"
#[extendr(r_name = "run_lola")]
pub fn r_run_lola(
    user_sets_list: List,
    universe_ptr: Robj,
    db_ptr: Robj,
    min_overlap: i32,
    direction: &str,
) -> extendr_api::Result<Robj> {
    let user_sets = extract_region_sets(user_sets_list)?;

    with_regionset!(universe_ptr, universe, {
        with_regiondb!(db_ptr, db, {
            let dir = match direction {
                "depletion" | "less" => Direction::Depletion,
                _ => Direction::Enrichment,
            };

            let config = LolaConfig {
                min_overlap,
                direction: dir,
            };

            let results = run_lola(&db.igd, &user_sets, universe, &config)
                .map_err(|e| extendr_api::Error::Other(format!("LOLA error: {}", e)))?;

            Ok(results_to_list(&results).into())
        })
    })
}

// =========================================================================
// Universe handling
// =========================================================================

/// Check universe appropriateness for user sets.
/// @export
/// @param user_sets_list R list of RegionSet external pointers
/// @param universe_ptr ExternalPtr to universe RegionSet
#[extendr(r_name = "check_universe")]
pub fn r_check_universe(
    user_sets_list: List,
    universe_ptr: Robj,
) -> extendr_api::Result<Robj> {
    let user_sets = extract_region_sets(user_sets_list)?;

    with_regionset!(universe_ptr, universe, {
        let report = universe::check_universe_appropriateness(&user_sets, universe);

        let mut us_index = Vec::new();
        let mut total = Vec::new();
        let mut in_univ = Vec::new();
        let mut coverage = Vec::new();
        let mut m2m = Vec::new();
        let mut warnings: Vec<String> = Vec::new();

        for ur in &report.user_set_reports {
            us_index.push((ur.user_set_index + 1) as i32);
            total.push(ur.total_regions as i32);
            in_univ.push(ur.regions_in_universe as i32);
            coverage.push(ur.coverage);
            m2m.push(ur.many_to_many_count as i32);
            if !ur.warnings.is_empty() {
                warnings.extend(ur.warnings.clone());
            }
        }

        Ok(list!(
            userSet = us_index,
            totalRegions = total,
            regionsInUniverse = in_univ,
            coverage = coverage,
            manyToMany = m2m,
            warnings = warnings
        )
        .into())
    })
}

/// Redefine user sets in terms of universe regions.
/// @export
/// @param user_sets_list R list of RegionSet external pointers
/// @param universe_ptr ExternalPtr to universe RegionSet
#[extendr(r_name = "redefine_user_sets")]
pub fn r_redefine_user_sets(
    user_sets_list: List,
    universe_ptr: Robj,
) -> extendr_api::Result<Robj> {
    let user_sets = extract_region_sets(user_sets_list)?;

    with_regionset!(universe_ptr, universe, {
        let redefined = universe::redefine_user_sets(&user_sets, universe);

        let result_list: Vec<Robj> = redefined
            .into_iter()
            .map(|rs| ExternalPtr::new(rs).into())
            .collect();

        Ok(List::from_values(result_list).into())
    })
}

/// Build a restricted universe from user sets.
/// @export
/// @param user_sets_list R list of RegionSet external pointers
#[extendr(r_name = "build_restricted_universe")]
pub fn r_build_restricted_universe(user_sets_list: List) -> extendr_api::Result<Robj> {
    let user_sets = extract_region_sets(user_sets_list)?;
    let restricted = universe::build_restricted_universe(&user_sets);
    Ok(ExternalPtr::new(restricted).into())
}

// =========================================================================
// Module registration
// =========================================================================

extendr_module! {
    mod lola;
    fn r_load_region_db;
    fn r_load_region_db_from_beds;
    fn r_list_region_sets;
    fn r_run_lola;
    fn r_check_universe;
    fn r_redefine_user_sets;
    fn r_build_restricted_universe;
}
