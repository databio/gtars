use extendr_api::prelude::*;

use gtars_core::models::{RegionSet, RegionSetList};
use gtars_igd::igd::Igd;
use gtars_lola::database::{RegionDB, RegionSetAnno};
use gtars_lola::enrichment::run_lola;
use gtars_lola::models::{Direction, LolaConfig};
use gtars_lola::output::{annotate_results, apply_fdr_correction};
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

macro_rules! with_regionsetlist {
    ($ptr:expr, $rsl:ident, $body:block) => {{
        let ext_ptr = <ExternalPtr<RegionSetList>>::try_from($ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid RegionSetList pointer".into()))?;
        let $rsl = &*ext_ptr;
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
    use gtars_lola::output::results_to_columns;

    let c = results_to_columns(results);

    // R-specific conversions: 1-based indices, i32 casts, NA types
    let user_set: Vec<i32> = c.user_set.iter().map(|&v| (v + 1) as i32).collect();
    let db_set: Vec<i32> = c.db_set.iter().map(|&v| (v + 1) as i32).collect();
    let support: Vec<i32> = c.support.iter().map(|&v| v as i32).collect();
    let rnk_pv: Vec<i32> = c.rnk_pv.iter().map(|&v| v as i32).collect();
    let rnk_or: Vec<i32> = c.rnk_or.iter().map(|&v| v as i32).collect();
    let rnk_sup: Vec<i32> = c.rnk_sup.iter().map(|&v| v as i32).collect();
    let max_rnk: Vec<i32> = c.max_rnk.iter().map(|&v| v as i32).collect();
    let b_vec: Vec<i32> = c.b.iter().map(|&v| v as i32).collect();
    let c_vec: Vec<i32> = c.c.iter().map(|&v| v as i32).collect();
    let d_vec: Vec<i32> = c.d.iter().map(|&v| v as i32).collect();
    let db_set_size: Vec<i32> = c.db_set_size.iter().map(|&v| v as i32).collect();

    let q_value_r: Vec<Rfloat> = c.q_value.iter()
        .map(|q| match q { Some(v) => Rfloat::from(*v), None => Rfloat::na() })
        .collect();

    let to_rstr = |v: &[Option<String>]| -> Vec<Rstr> {
        v.iter()
            .map(|s| match s {
                Some(val) => Rstr::from(val.as_str()),
                None => Rstr::na(),
            })
            .collect()
    };

    list!(
        userSet = user_set,
        dbSet = db_set,
        collection = to_rstr(&c.collection),
        pValueLog = c.p_value_log,
        oddsRatio = c.odds_ratio,
        support = support,
        rnkPV = rnk_pv,
        rnkOR = rnk_or,
        rnkSup = rnk_sup,
        maxRnk = max_rnk,
        meanRnk = c.mean_rnk,
        b = b_vec,
        c = c_vec,
        d = d_vec,
        description = to_rstr(&c.description),
        cellType = to_rstr(&c.cell_type),
        tissue = to_rstr(&c.tissue),
        antibody = to_rstr(&c.antibody),
        treatment = to_rstr(&c.treatment),
        dataSource = to_rstr(&c.data_source),
        filename = c.filename,
        qValue = q_value_r,
        size = db_set_size
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

            let mut results = run_lola(&db.igd, &user_sets, universe, &config)
                .map_err(|e| extendr_api::Error::Other(format!("LOLA error: {}", e)))?;

            annotate_results(&mut results, db);
            apply_fdr_correction(&mut results);

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
        let universe_igd = Igd::from_single_region_set(universe);
        let report = universe::check_universe_appropriateness(&user_sets, &universe_igd);

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
        let universe_igd = Igd::from_single_region_set(universe);
        let redefined = universe::redefine_user_sets(&user_sets, universe, &universe_igd);

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
// RegionDB accessors
// =========================================================================

/// Get per-file annotations from a RegionDB as a data.frame-like list.
///
/// Returns a list with columns: filename, cellType, description, tissue,
/// dataSource, antibody, treatment, collection. This matches the structure
/// of R LOLA's regionDB$regionAnno.
/// @export
/// @param db ExternalPtr to RegionDB
#[extendr(r_name = "regiondb_anno")]
pub fn r_regiondb_anno(db: Robj) -> extendr_api::Result<Robj> {
    with_regiondb!(db, db_ref, {
        let n = db_ref.region_anno.len();
        let mut filename = Vec::with_capacity(n);
        let mut cell_type: Vec<Option<String>> = Vec::with_capacity(n);
        let mut description: Vec<Option<String>> = Vec::with_capacity(n);
        let mut tissue: Vec<Option<String>> = Vec::with_capacity(n);
        let mut data_source: Vec<Option<String>> = Vec::with_capacity(n);
        let mut antibody: Vec<Option<String>> = Vec::with_capacity(n);
        let mut treatment: Vec<Option<String>> = Vec::with_capacity(n);
        let mut collection: Vec<Option<String>> = Vec::with_capacity(n);
        let mut size = Vec::with_capacity(n);

        for (i, a) in db_ref.region_anno.iter().enumerate() {
            filename.push(a.filename.clone());
            cell_type.push(empty_to_na(&a.cell_type));
            description.push(empty_to_na(&a.description));
            tissue.push(empty_to_na(&a.tissue));
            data_source.push(empty_to_na(&a.data_source));
            antibody.push(empty_to_na(&a.antibody));
            treatment.push(empty_to_na(&a.treatment));
            collection.push(empty_to_na(&a.collection));
            size.push(db_ref.region_sets[i].len() as i32);
        }

        let to_rstr = |v: &[Option<String>]| -> Vec<Rstr> {
            v.iter()
                .map(|s| match s {
                    Some(val) => Rstr::from(val.as_str()),
                    None => Rstr::na(),
                })
                .collect()
        };

        Ok(list!(
            filename = filename,
            cellType = to_rstr(&cell_type),
            description = to_rstr(&description),
            tissue = to_rstr(&tissue),
            dataSource = to_rstr(&data_source),
            antibody = to_rstr(&antibody),
            treatment = to_rstr(&treatment),
            collection = to_rstr(&collection),
            size = size
        )
        .into())
    })
}

/// Get collection-level annotations from a RegionDB as a data.frame-like list.
/// @export
/// @param db ExternalPtr to RegionDB
#[extendr(r_name = "regiondb_collection_anno")]
pub fn r_regiondb_collection_anno(db: Robj) -> extendr_api::Result<Robj> {
    with_regiondb!(db, db_ref, {
        let n = db_ref.collection_anno.len();
        let mut collector = Vec::with_capacity(n);
        let mut date = Vec::with_capacity(n);
        let mut source = Vec::with_capacity(n);
        let mut description = Vec::with_capacity(n);
        let mut collection_name = Vec::with_capacity(n);

        for a in &db_ref.collection_anno {
            collector.push(a.collector.clone());
            date.push(a.date.clone());
            source.push(a.source.clone());
            description.push(a.description.clone());
            collection_name.push(a.collection_name.clone());
        }

        Ok(list!(
            collectionname = collection_name,
            collector = collector,
            date = date,
            source = source,
            description = description
        )
        .into())
    })
}

/// Get a single RegionSet from a RegionDB by 1-based index.
/// @export
/// @param db ExternalPtr to RegionDB
/// @param index 1-based integer index into the region sets
#[extendr(r_name = "regiondb_region_set")]
pub fn r_regiondb_region_set(db: Robj, index: i32) -> extendr_api::Result<Robj> {
    with_regiondb!(db, db_ref, {
        let i = (index - 1) as usize; // convert to 0-based
        if i >= db_ref.region_sets.len() {
            return Err(extendr_api::Error::Other(format!(
                "Index {} out of range (database has {} region sets)",
                index,
                db_ref.region_sets.len()
            )));
        }
        Ok(ExternalPtr::new(db_ref.region_sets[i].clone()).into())
    })
}

// =========================================================================
// RegionSetList
// =========================================================================

/// Create a RegionSetList from a list of RegionSet external pointers.
/// @export
/// @param sets_list R list of RegionSet external pointers
#[extendr(r_name = "regionsetlist_from_sets")]
pub fn r_regionsetlist_from_sets(sets_list: List) -> extendr_api::Result<Robj> {
    let sets = extract_region_sets(sets_list)?;
    Ok(ExternalPtr::new(RegionSetList::new(sets)).into())
}

/// Extract region sets from a RegionDB by 1-based indices as a RegionSetList.
/// @export
/// @param db ExternalPtr to RegionDB
/// @param indices Integer vector of 1-based indices
#[extendr(r_name = "regionsetlist_from_db")]
pub fn r_regionsetlist_from_db(db: Robj, indices: Vec<i32>) -> extendr_api::Result<Robj> {
    with_regiondb!(db, db_ref, {
        let idx: Vec<usize> = indices.iter().map(|&i| (i - 1) as usize).collect();
        let rsl = db_ref.get_region_set_list(&idx);
        Ok(ExternalPtr::new(rsl).into())
    })
}

/// Get the number of region sets in a RegionSetList.
/// @export
/// @param rsl ExternalPtr to RegionSetList
#[extendr(r_name = "regionsetlist_length")]
pub fn r_regionsetlist_length(rsl: Robj) -> extendr_api::Result<i32> {
    with_regionsetlist!(rsl, rsl_ref, {
        Ok(rsl_ref.len() as i32)
    })
}

/// Get a single RegionSet from a RegionSetList by 1-based index.
/// @export
/// @param rsl ExternalPtr to RegionSetList
/// @param index 1-based integer index
#[extendr(r_name = "regionsetlist_get")]
pub fn r_regionsetlist_get(rsl: Robj, index: i32) -> extendr_api::Result<Robj> {
    with_regionsetlist!(rsl, rsl_ref, {
        let i = (index - 1) as usize;
        match rsl_ref.get(i) {
            Some(rs) => Ok(ExternalPtr::new(rs.clone()).into()),
            None => Err(extendr_api::Error::Other(format!(
                "Index {} out of range (list has {} region sets)",
                index,
                rsl_ref.len()
            ))),
        }
    })
}

/// Flatten a RegionSetList into a single RegionSet (no merging).
/// @export
/// @param rsl ExternalPtr to RegionSetList
#[extendr(r_name = "regionsetlist_concat")]
pub fn r_regionsetlist_concat(rsl: Robj) -> extendr_api::Result<Robj> {
    with_regionsetlist!(rsl, rsl_ref, {
        Ok(ExternalPtr::new(rsl_ref.concat()).into())
    })
}

/// Get the names from a RegionSetList, or NULL if no names.
/// @export
/// @param rsl ExternalPtr to RegionSetList
#[extendr(r_name = "regionsetlist_names")]
pub fn r_regionsetlist_names(rsl: Robj) -> extendr_api::Result<Robj> {
    with_regionsetlist!(rsl, rsl_ref, {
        match &rsl_ref.names {
            Some(names) => Ok(names.clone().into()),
            None => Ok(().into()),
        }
    })
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
    fn r_regiondb_anno;
    fn r_regiondb_collection_anno;
    fn r_regiondb_region_set;
    fn r_regionsetlist_from_sets;
    fn r_regionsetlist_from_db;
    fn r_regionsetlist_length;
    fn r_regionsetlist_get;
    fn r_regionsetlist_concat;
    fn r_regionsetlist_names;
}
