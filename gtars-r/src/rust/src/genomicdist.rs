use std::collections::HashMap;

use extendr_api::prelude::*;

use gtars_core::models::{Region, RegionSet};
use gtars_genomicdist::models::{GenomeAssembly, TssIndex};
use gtars_overlaprs::RegionSetOverlaps;
use gtars_genomicdist::{
    calc_dinucl_freq_per_region, calc_gc_content, calc_summary_signal, chrom_karyotype_key,
    consensus, genome_partition_list, calc_expected_partitions, calc_partitions,
    CoordinateMode, GenomicDistAnnotation, GenomicIntervalSetStatistics, GeneModel, IntervalRanges,
    PartitionList, SignalMatrix, Strand, StrandedRegionSet, DINUCL_ORDER,
};

// =========================================================================
// Helper macros
// =========================================================================

macro_rules! with_regionset {
    ($ptr:expr, $rs:ident, $body:block) => {{
        let ext_ptr = <ExternalPtr<RegionSet>>::try_from($ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer".into()))?;
        let $rs = &*ext_ptr;
        $body
    }};
}

macro_rules! with_assembly {
    ($ptr:expr, $asm:ident, $body:block) => {{
        let ext_ptr = <ExternalPtr<GenomeAssembly>>::try_from($ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid GenomeAssembly pointer".into()))?;
        let $asm = &*ext_ptr;
        $body
    }};
}

macro_rules! with_partitions {
    ($ptr:expr, $pl:ident, $body:block) => {{
        let ext_ptr = <ExternalPtr<PartitionList>>::try_from($ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid PartitionList pointer".into()))?;
        let $pl = &*ext_ptr;
        $body
    }};
}

// =========================================================================
// Helper functions
// =========================================================================

/// Convert an R signed integer to u32 with a clear error on negative values.
fn checked_u32(value: i32, name: &str) -> extendr_api::Result<u32> {
    u32::try_from(value).map_err(|_| {
        extendr_api::Error::Other(format!("{} must be non-negative, got {}", name, value))
    })
}

/// Construct RegionSet from R vectors (0-based half-open coordinates)
fn regionset_from_vecs(chrs: Vec<String>, starts: Vec<i32>, ends: Vec<i32>) -> extendr_api::Result<RegionSet> {
    let regions: Vec<Region> = chrs
        .into_iter()
        .zip(starts.into_iter().zip(ends.into_iter()))
        .map(|(chr, (start, end))| Ok(Region {
            chr,
            start: checked_u32(start, "start")?,
            end: checked_u32(end, "end")?,
            rest: None,
        }))
        .collect::<extendr_api::Result<Vec<Region>>>()?;
    Ok(RegionSet::from(regions))
}

/// Build HashMap<String, u32> from parallel name/size vectors
fn chrom_sizes_from_vecs(names: Vec<String>, sizes: Vec<i32>) -> extendr_api::Result<HashMap<String, u32>> {
    names
        .into_iter()
        .zip(sizes.into_iter())
        .map(|(name, size)| Ok((name, checked_u32(size, "chrom_size")?)))
        .collect()
}

// =========================================================================
// 0. RegionSet core
// =========================================================================

/// Load a RegionSet from a BED/narrowPeak/gzip file
/// @export
/// @param bed_path Path to a BED, narrowPeak, or gzipped BED file
#[extendr(r_name = "load_regionset")]
pub fn r_load_regionset(bed_path: &str) -> extendr_api::Result<Robj> {
    let rs = RegionSet::try_from(bed_path)
        .map_err(|e| extendr_api::Error::Other(format!("Loading RegionSet: {}", e)))?;
    Ok(ExternalPtr::new(rs).into())
}

/// Create a RegionSet from R vectors (0-based half-open coordinates)
/// @export
/// @param chrs Character vector of chromosome names
/// @param starts Integer vector of start positions (0-based)
/// @param ends Integer vector of end positions (half-open)
#[extendr(r_name = "regionset_from_vectors")]
pub fn r_regionset_from_vectors(chrs: Vec<String>, starts: Vec<i32>, ends: Vec<i32>) -> extendr_api::Result<Robj> {
    let rs = regionset_from_vecs(chrs, starts, ends)?;
    Ok(ExternalPtr::new(rs).into())
}

/// Extract chr/start/end vectors from a RegionSet pointer
/// @export
/// @param rs_ptr External pointer to a RegionSet
#[extendr(r_name = "regionset_to_vectors")]
pub fn r_regionset_to_vectors(rs_ptr: Robj) -> extendr_api::Result<List> {
    with_regionset!(rs_ptr, rs, {
        let chrs: Vec<String> = rs.regions.iter().map(|r| r.chr.clone()).collect();
        let starts: Vec<f64> = rs.regions.iter().map(|r| r.start as f64).collect();
        let ends: Vec<f64> = rs.regions.iter().map(|r| r.end as f64).collect();
        Ok(list!(chr = chrs, start = starts, end = ends))
    })
}

/// Get the number of regions in a RegionSet
/// @export
/// @param rs_ptr External pointer to a RegionSet
#[extendr(r_name = "regionset_length")]
pub fn r_regionset_length(rs_ptr: Robj) -> extendr_api::Result<i32> {
    with_regionset!(rs_ptr, rs, {
        Ok(rs.regions.len() as i32)
    })
}

// =========================================================================
// 1. RegionSet Statistics
// =========================================================================

/// Calculate region widths (end - start)
/// @export
/// @param rs_ptr External pointer to a RegionSet
#[extendr(r_name = "r_calc_widths")]
pub fn r_calc_widths(rs_ptr: Robj) -> extendr_api::Result<Vec<f64>> {
    with_regionset!(rs_ptr, rs, {
        Ok(rs.calc_widths().into_iter().map(|w| w as f64).collect())
    })
}

/// Calculate distances between consecutive regions on each chromosome
/// @export
/// @param rs_ptr External pointer to a RegionSet
#[extendr(r_name = "r_calc_neighbor_distances")]
pub fn r_calc_neighbor_distances(rs_ptr: Robj) -> extendr_api::Result<Vec<f64>> {
    with_regionset!(rs_ptr, rs, {
        let dists = rs
            .calc_neighbor_distances()
            .map_err(|e| extendr_api::Error::Other(format!("{}", e)))?;
        Ok(dists.into_iter().map(|d| d as f64).collect())
    })
}

/// Calculate distance from each region to its nearest neighbor
/// @export
/// @param rs_ptr External pointer to a RegionSet
#[extendr(r_name = "r_calc_nearest_neighbors")]
pub fn r_calc_nearest_neighbors(rs_ptr: Robj) -> extendr_api::Result<Vec<f64>> {
    with_regionset!(rs_ptr, rs, {
        let dists = rs
            .calc_nearest_neighbors()
            .map_err(|e| extendr_api::Error::Other(format!("{}", e)))?;
        Ok(dists.into_iter().map(|d| d as f64).collect())
    })
}

/// Calculate per-chromosome statistics
/// @export
/// @param rs_ptr External pointer to a RegionSet
#[extendr(r_name = "r_chromosome_statistics")]
pub fn r_chromosome_statistics(rs_ptr: Robj) -> extendr_api::Result<List> {
    with_regionset!(rs_ptr, rs, {
        let stats = rs.chromosome_statistics();
        let mut entries: Vec<_> = stats.into_iter().collect();
        entries.sort_by(|(a, _), (b, _)| {
            chrom_karyotype_key(a).cmp(&chrom_karyotype_key(b))
        });

        let mut chr_names: Vec<String> = Vec::new();
        let mut n_regions: Vec<f64> = Vec::new();
        let mut start_pos: Vec<f64> = Vec::new();
        let mut end_pos: Vec<f64> = Vec::new();
        let mut min_len: Vec<f64> = Vec::new();
        let mut max_len: Vec<f64> = Vec::new();
        let mut mean_len: Vec<f64> = Vec::new();
        let mut median_len: Vec<f64> = Vec::new();

        for (chr, s) in &entries {
            chr_names.push(chr.clone());
            n_regions.push(s.number_of_regions as f64);
            start_pos.push(s.start_nucleotide_position as f64);
            end_pos.push(s.end_nucleotide_position as f64);
            min_len.push(s.minimum_region_length as f64);
            max_len.push(s.maximum_region_length as f64);
            mean_len.push(s.mean_region_length);
            median_len.push(s.median_region_length);
        }

        Ok(list!(
            chromosome = chr_names,
            number_of_regions = n_regions,
            start_nucleotide_position = start_pos,
            end_nucleotide_position = end_pos,
            minimum_region_length = min_len,
            maximum_region_length = max_len,
            mean_region_length = mean_len,
            median_region_length = median_len
        ))
    })
}

/// Calculate region distribution across chromosome bins
/// @export
/// @param rs_ptr External pointer to a RegionSet
/// @param n_bins Number of bins (default 250)
#[extendr(r_name = "r_region_distribution")]
pub fn r_region_distribution(rs_ptr: Robj, n_bins: i32, chrom_names: Robj, chrom_lengths: Robj) -> extendr_api::Result<List> {
    let n_bins_u32 = checked_u32(n_bins, "n_bins")?;
    with_regionset!(rs_ptr, rs, {
        let dist = if !chrom_names.is_null() && !chrom_lengths.is_null() {
            let names: Vec<String> = chrom_names.as_str_iter()
                .ok_or_else(|| extendr_api::Error::Other("chrom_names must be character".into()))?
                .map(|s| s.to_string())
                .collect();
            let lengths: Vec<f64> = chrom_lengths.as_real_iter()
                .ok_or_else(|| extendr_api::Error::Other("chrom_lengths must be numeric".into()))?
                .copied()
                .collect();
            let chrom_sizes: HashMap<String, u32> = names
                .into_iter()
                .zip(lengths.into_iter().map(|v| v as u32))
                .collect();
            rs.region_distribution_with_chrom_sizes(n_bins_u32, &chrom_sizes)
        } else {
            rs.region_distribution_with_bins(n_bins_u32)
        };
        let mut bins: Vec<_> = dist.into_values().collect();
        bins.sort_by(|a, b| {
            chrom_karyotype_key(&a.chr)
                .cmp(&chrom_karyotype_key(&b.chr))
                .then(a.start.cmp(&b.start))
        });

        let mut chrs: Vec<String> = Vec::new();
        let mut starts: Vec<f64> = Vec::new();
        let mut ends: Vec<f64> = Vec::new();
        let mut counts: Vec<f64> = Vec::new();
        let mut rids: Vec<f64> = Vec::new();

        for bin in &bins {
            chrs.push(bin.chr.clone());
            starts.push(bin.start as f64);
            ends.push(bin.end as f64);
            counts.push(bin.n as f64);
            rids.push(bin.rid as f64);
        }

        Ok(list!(
            chr = chrs,
            start = starts,
            end = ends,
            n = counts,
            rid = rids
        ))
    })
}

// =========================================================================
// 2. GC Content & Dinucleotide Freq
// =========================================================================

/// Load a genome assembly from a FASTA file
/// @export
/// @param fasta_path Path to a FASTA file
#[extendr(r_name = "load_genome_assembly")]
pub fn r_load_genome_assembly(fasta_path: &str) -> extendr_api::Result<Robj> {
    let assembly = GenomeAssembly::try_from(fasta_path)
        .map_err(|e| extendr_api::Error::Other(format!("Loading FASTA: {}", e)))?;
    Ok(ExternalPtr::new(assembly).into())
}

/// Calculate GC content for each region
/// @export
/// @param rs_ptr External pointer to a RegionSet
/// @param assembly_ptr External pointer to a GenomeAssembly
/// @param ignore_unk_chroms Skip regions on chromosomes not in the assembly
#[extendr(r_name = "r_calc_gc_content")]
pub fn r_calc_gc_content(
    rs_ptr: Robj,
    assembly_ptr: Robj,
    ignore_unk_chroms: bool,
) -> extendr_api::Result<Vec<f64>> {
    with_regionset!(rs_ptr, rs, {
        with_assembly!(assembly_ptr, asm, {
            calc_gc_content(rs, asm, ignore_unk_chroms)
                .map_err(|e| extendr_api::Error::Other(format!("{}", e)))
        })
    })
}

/// Calculate per-region dinucleotide frequencies (percentages)
/// @export
/// @param rs_ptr External pointer to a RegionSet
/// @param assembly_ptr External pointer to a GenomeAssembly
#[extendr(r_name = "r_calc_dinucl_freq")]
pub fn r_calc_dinucl_freq(rs_ptr: Robj, assembly_ptr: Robj) -> extendr_api::Result<List> {
    with_regionset!(rs_ptr, rs, {
        with_assembly!(assembly_ptr, asm, {
            let (labels, matrix) = calc_dinucl_freq_per_region(rs, asm)
                .map_err(|e| extendr_api::Error::Other(format!("{}", e)))?;

            // Build column names (uppercase to match GD: AA, AC, AG, ...)
            let col_names: Vec<String> = DINUCL_ORDER
                .iter()
                .map(|d| d.to_string().unwrap_or_else(|_| "??".into()).to_uppercase())
                .collect();

            // Build per-dinucleotide columns (each is a Vec<f64> of length n_regions)
            let n = labels.len();
            let mut columns: Vec<Vec<f64>> = vec![Vec::with_capacity(n); 16];
            for row in &matrix {
                for (col_idx, &val) in row.iter().enumerate() {
                    columns[col_idx].push(val);
                }
            }

            // Assemble named list: region + 16 dinucleotide columns
            let mut pairs: Vec<(&str, Robj)> = Vec::with_capacity(17);
            let region_robj = labels.into_iter().collect_robj();
            pairs.push(("region", region_robj));
            // Need to keep column Robjs alive
            let col_robjs: Vec<Robj> = columns
                .into_iter()
                .map(|v| v.into_iter().collect_robj())
                .collect();
            for (i, name) in col_names.iter().enumerate() {
                pairs.push((name.as_str(), col_robjs[i].clone()));
            }

            Ok(List::from_pairs(pairs))
        })
    })
}

// =========================================================================
// 3. Interval Ranges
// =========================================================================

/// Trim regions to chromosome boundaries
/// @export
/// @param rs_ptr External pointer to a RegionSet
/// @param chrom_names Character vector of chromosome names
/// @param chrom_sizes Integer vector of chromosome sizes
#[extendr(r_name = "r_trim")]
pub fn r_trim(
    rs_ptr: Robj,
    chrom_names: Vec<String>,
    chrom_sizes: Vec<i32>,
) -> extendr_api::Result<Robj> {
    let sizes = chrom_sizes_from_vecs(chrom_names, chrom_sizes)?;
    with_regionset!(rs_ptr, rs, {
        let trimmed = rs.trim(&sizes);
        Ok(ExternalPtr::new(trimmed).into())
    })
}

/// Generate promoter regions relative to each region's start
/// @export
/// @param rs_ptr External pointer to a RegionSet
/// @param upstream Bases upstream of start
/// @param downstream Bases downstream of start
#[extendr(r_name = "r_promoters")]
pub fn r_promoters(rs_ptr: Robj, upstream: i32, downstream: i32) -> extendr_api::Result<Robj> {
    let upstream_u32 = checked_u32(upstream, "upstream")?;
    let downstream_u32 = checked_u32(downstream, "downstream")?;
    with_regionset!(rs_ptr, rs, {
        let result = rs.promoters(upstream_u32, downstream_u32);
        Ok(ExternalPtr::new(result).into())
    })
}

/// Merge overlapping and adjacent intervals per chromosome
/// @export
/// @param rs_ptr External pointer to a RegionSet
#[extendr(r_name = "r_reduce")]
pub fn r_reduce(rs_ptr: Robj) -> extendr_api::Result<Robj> {
    with_regionset!(rs_ptr, rs, {
        let result = rs.reduce();
        Ok(ExternalPtr::new(result).into())
    })
}

/// Set difference: remove portions of A that overlap with B
/// @export
/// @param rs_ptr_a External pointer to RegionSet A
/// @param rs_ptr_b External pointer to RegionSet B
#[extendr(r_name = "r_setdiff")]
pub fn r_setdiff(rs_ptr_a: Robj, rs_ptr_b: Robj) -> extendr_api::Result<Robj> {
    let ext_a = <ExternalPtr<RegionSet>>::try_from(rs_ptr_a)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (a)".into()))?;
    let ext_b = <ExternalPtr<RegionSet>>::try_from(rs_ptr_b)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (b)".into()))?;
    let result = ext_a.setdiff(&*ext_b);
    Ok(ExternalPtr::new(result).into())
}

/// Pairwise intersection by index position
/// @export
/// @param rs_ptr_a External pointer to RegionSet A
/// @param rs_ptr_b External pointer to RegionSet B
#[extendr(r_name = "r_pintersect")]
pub fn r_pintersect(rs_ptr_a: Robj, rs_ptr_b: Robj) -> extendr_api::Result<Robj> {
    let ext_a = <ExternalPtr<RegionSet>>::try_from(rs_ptr_a)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (a)".into()))?;
    let ext_b = <ExternalPtr<RegionSet>>::try_from(rs_ptr_b)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (b)".into()))?;
    let result = ext_a.pintersect(&*ext_b);
    Ok(ExternalPtr::new(result).into())
}

/// Combine two region sets without merging
/// @export
/// @param rs_ptr_a External pointer to RegionSet A
/// @param rs_ptr_b External pointer to RegionSet B
#[extendr(r_name = "r_concat")]
pub fn r_concat(rs_ptr_a: Robj, rs_ptr_b: Robj) -> extendr_api::Result<Robj> {
    let ext_a = <ExternalPtr<RegionSet>>::try_from(rs_ptr_a)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (a)".into()))?;
    let ext_b = <ExternalPtr<RegionSet>>::try_from(rs_ptr_b)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (b)".into()))?;
    let result = ext_a.concat(&*ext_b);
    Ok(ExternalPtr::new(result).into())
}

/// Merge two region sets into a minimal non-overlapping set
/// @export
/// @param rs_ptr_a External pointer to RegionSet A
/// @param rs_ptr_b External pointer to RegionSet B
#[extendr(r_name = "r_union")]
pub fn r_union(rs_ptr_a: Robj, rs_ptr_b: Robj) -> extendr_api::Result<Robj> {
    let ext_a = <ExternalPtr<RegionSet>>::try_from(rs_ptr_a)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (a)".into()))?;
    let ext_b = <ExternalPtr<RegionSet>>::try_from(rs_ptr_b)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (b)".into()))?;
    let result = ext_a.union(&*ext_b);
    Ok(ExternalPtr::new(result).into())
}

/// Nucleotide-level Jaccard similarity between two region sets
/// @export
/// @param rs_ptr_a External pointer to RegionSet A
/// @param rs_ptr_b External pointer to RegionSet B
#[extendr(r_name = "r_jaccard")]
pub fn r_jaccard(rs_ptr_a: Robj, rs_ptr_b: Robj) -> extendr_api::Result<f64> {
    let ext_a = <ExternalPtr<RegionSet>>::try_from(rs_ptr_a)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (a)".into()))?;
    let ext_b = <ExternalPtr<RegionSet>>::try_from(rs_ptr_b)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (b)".into()))?;
    Ok(ext_a.jaccard(&*ext_b))
}

/// Shift all regions by a fixed offset
/// @export
/// @param rs_ptr External pointer to a RegionSet
/// @param offset Integer offset (positive = downstream, negative = upstream)
#[extendr(r_name = "r_shift")]
pub fn r_shift(rs_ptr: Robj, offset: i32) -> extendr_api::Result<Robj> {
    with_regionset!(rs_ptr, rs, {
        let result = rs.shift(offset as i64);
        Ok(ExternalPtr::new(result).into())
    })
}

/// Generate flanking regions
/// @export
/// @param rs_ptr External pointer to a RegionSet
/// @param width Flank width in base pairs
/// @param use_start If TRUE, flank upstream of start; if FALSE, downstream of end
/// @param both If TRUE, flank on both sides of the anchor
#[extendr(r_name = "r_flank")]
pub fn r_flank(rs_ptr: Robj, width: i32, use_start: bool, both: bool) -> extendr_api::Result<Robj> {
    let width_u32 = checked_u32(width, "width")?;
    with_regionset!(rs_ptr, rs, {
        let result = rs.flank(width_u32, use_start, both);
        Ok(ExternalPtr::new(result).into())
    })
}

/// Resize regions to a fixed width
/// @export
/// @param rs_ptr External pointer to a RegionSet
/// @param width New width in base pairs
/// @param fix Anchor point: "start", "end", or "center"
#[extendr(r_name = "r_resize")]
pub fn r_resize(rs_ptr: Robj, width: i32, fix: &str) -> extendr_api::Result<Robj> {
    let width_u32 = checked_u32(width, "width")?;
    with_regionset!(rs_ptr, rs, {
        let result = rs.resize(width_u32, fix);
        Ok(ExternalPtr::new(result).into())
    })
}

/// Narrow regions by specifying a relative sub-range
/// @export
/// @param rs_ptr External pointer to a RegionSet
/// @param start 1-based relative start (NA to omit)
/// @param end 1-based relative end (NA to omit)
/// @param width Width of the sub-range (NA to omit)
#[extendr(r_name = "r_narrow")]
pub fn r_narrow(rs_ptr: Robj, start: Robj, end: Robj, width: Robj) -> extendr_api::Result<Robj> {
    with_regionset!(rs_ptr, rs, {
        let s = if start.is_na() { None } else { Some(checked_u32(i32::try_from(start).unwrap_or(1), "start")?) };
        let e = if end.is_na() { None } else { Some(checked_u32(i32::try_from(end).unwrap_or(1), "end")?) };
        let w = if width.is_na() { None } else { Some(checked_u32(i32::try_from(width).unwrap_or(1), "width")?) };
        let result = rs.narrow(s, e, w);
        Ok(ExternalPtr::new(result).into())
    })
}

/// Break regions into non-overlapping disjoint pieces
/// @export
/// @param rs_ptr External pointer to a RegionSet
#[extendr(r_name = "r_disjoin")]
pub fn r_disjoin(rs_ptr: Robj) -> extendr_api::Result<Robj> {
    with_regionset!(rs_ptr, rs, {
        let result = rs.disjoin();
        Ok(ExternalPtr::new(result).into())
    })
}

/// Return gaps between regions per chromosome
/// @export
/// @param rs_ptr External pointer to a RegionSet
#[extendr(r_name = "r_gaps")]
pub fn r_gaps(rs_ptr: Robj) -> extendr_api::Result<Robj> {
    with_regionset!(rs_ptr, rs, {
        let result = rs.gaps();
        Ok(ExternalPtr::new(result).into())
    })
}

/// Range-level intersection of two region sets
/// @export
/// @param rs_ptr_a External pointer to RegionSet A
/// @param rs_ptr_b External pointer to RegionSet B
#[extendr(r_name = "r_intersect")]
pub fn r_intersect(rs_ptr_a: Robj, rs_ptr_b: Robj) -> extendr_api::Result<Robj> {
    let ext_a = <ExternalPtr<RegionSet>>::try_from(rs_ptr_a)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (a)".into()))?;
    let ext_b = <ExternalPtr<RegionSet>>::try_from(rs_ptr_b)
        .map_err(|_| extendr_api::Error::Other("Invalid RegionSet pointer (b)".into()))?;
    let result = ext_a.intersect(&*ext_b);
    Ok(ExternalPtr::new(result).into())
}

/// Compute consensus regions from a list of RegionSet pointers
/// @export
/// @param rs_list An R list of external pointers to RegionSets
#[extendr(r_name = "r_consensus")]
pub fn r_consensus(rs_list: List) -> extendr_api::Result<List> {
    let mut sets: Vec<RegionSet> = Vec::with_capacity(rs_list.len());
    for (i, item) in rs_list.values().enumerate() {
        let ext = <ExternalPtr<RegionSet>>::try_from(item).map_err(|_| {
            extendr_api::Error::Other(format!("Invalid RegionSet pointer at index {}", i + 1))
        })?;
        sets.push((*ext).clone());
    }
    let result = consensus(&sets);
    let chrs: Vec<String> = result.iter().map(|r| r.chr.clone()).collect();
    let starts: Vec<f64> = result.iter().map(|r| r.start as f64).collect();
    let ends: Vec<f64> = result.iter().map(|r| r.end as f64).collect();
    let counts: Vec<f64> = result.iter().map(|r| r.count as f64).collect();
    Ok(list!(chr = chrs, start = starts, end = ends, count = counts))
}

/// Find all overlapping (queryHits, subjectHits) pairs between two RegionSets
/// @export
/// @param query_ptr External pointer to query RegionSet
/// @param subject_ptr External pointer to subject RegionSet
/// @param minoverlap Minimum overlap in base pairs (default 1)
#[extendr(r_name = "r_find_overlaps")]
pub fn r_find_overlaps(
    query_ptr: Robj,
    subject_ptr: Robj,
    minoverlap: i32,
) -> extendr_api::Result<List> {
    let ext_q = <ExternalPtr<RegionSet>>::try_from(query_ptr)
        .map_err(|_| extendr_api::Error::Other("Invalid query RegionSet pointer".into()))?;
    let ext_s = <ExternalPtr<RegionSet>>::try_from(subject_ptr)
        .map_err(|_| extendr_api::Error::Other("Invalid subject RegionSet pointer".into()))?;
    let indices = RegionSetOverlaps::find_overlaps(&*ext_q, &*ext_s, Some(minoverlap));
    // Flatten Vec<Vec<usize>> to (queryHits, subjectHits) pair vecs, 1-based for R
    let mut query_hits: Vec<i32> = Vec::new();
    let mut subject_hits: Vec<i32> = Vec::new();
    for (q_idx, hits) in indices.iter().enumerate() {
        for &s_idx in hits {
            query_hits.push(q_idx as i32 + 1);
            subject_hits.push(s_idx as i32 + 1);
        }
    }
    Ok(list!(queryHits = query_hits, subjectHits = subject_hits))
}

/// Count the number of subject regions overlapping each query region
/// @export
/// @param query_ptr External pointer to query RegionSet
/// @param subject_ptr External pointer to subject RegionSet
/// @param minoverlap Minimum overlap in base pairs (default 1)
#[extendr(r_name = "r_count_overlaps")]
pub fn r_count_overlaps(
    query_ptr: Robj,
    subject_ptr: Robj,
    minoverlap: i32,
) -> extendr_api::Result<Vec<i32>> {
    let ext_q = <ExternalPtr<RegionSet>>::try_from(query_ptr)
        .map_err(|_| extendr_api::Error::Other("Invalid query RegionSet pointer".into()))?;
    let ext_s = <ExternalPtr<RegionSet>>::try_from(subject_ptr)
        .map_err(|_| extendr_api::Error::Other("Invalid subject RegionSet pointer".into()))?;
    let counts = RegionSetOverlaps::count_overlaps(&*ext_q, &*ext_s, Some(minoverlap));
    Ok(counts.into_iter().map(|c| c as i32).collect())
}

// =========================================================================
// 4. Partitions
// =========================================================================

/// Build a PartitionList from RegionSet pointers (genes, exons, optional UTRs)
/// @export
/// @param genes_ptr External pointer to genes RegionSet
/// @param exons_ptr External pointer to exons RegionSet
/// @param three_utr_ptr External pointer to 3'UTR RegionSet (or NULL)
/// @param five_utr_ptr External pointer to 5'UTR RegionSet (or NULL)
/// @param core_prom Core promoter size in bp
/// @param prox_prom Proximal promoter size in bp
/// @param chrom_names Chromosome names for trim (empty to skip)
/// @param chrom_sizes_vec Chromosome sizes for trim (parallel to chrom_names)
#[extendr(r_name = "r_partition_list_from_regions")]
pub fn r_partition_list_from_regions(
    genes_ptr: Robj,
    exons_ptr: Robj,
    three_utr_ptr: Robj,
    five_utr_ptr: Robj,
    core_prom: i32,
    prox_prom: i32,
    chrom_names: Vec<String>,
    chrom_sizes_vec: Vec<i32>,
) -> extendr_api::Result<Robj> {
    let genes_ext = <ExternalPtr<RegionSet>>::try_from(genes_ptr)
        .map_err(|_| extendr_api::Error::Other("Invalid genes RegionSet pointer".into()))?;
    let exons_ext = <ExternalPtr<RegionSet>>::try_from(exons_ptr)
        .map_err(|_| extendr_api::Error::Other("Invalid exons RegionSet pointer".into()))?;

    let three_utr: Option<RegionSet> = if three_utr_ptr.is_null() {
        None
    } else {
        let ext = <ExternalPtr<RegionSet>>::try_from(three_utr_ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid 3'UTR RegionSet pointer".into()))?;
        Some(ext.regions.clone().into())
    };

    let five_utr: Option<RegionSet> = if five_utr_ptr.is_null() {
        None
    } else {
        let ext = <ExternalPtr<RegionSet>>::try_from(five_utr_ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid 5'UTR RegionSet pointer".into()))?;
        Some(ext.regions.clone().into())
    };

    let model = GeneModel {
        genes: StrandedRegionSet::unstranded(genes_ext.reduce()),
        exons: StrandedRegionSet::unstranded(exons_ext.reduce()),
        three_utr: three_utr
            .map(|rs| StrandedRegionSet::unstranded(rs.reduce()))
            .filter(|srs| !srs.is_empty()),
        five_utr: five_utr
            .map(|rs| StrandedRegionSet::unstranded(rs.reduce()))
            .filter(|srs| !srs.is_empty()),
    };

    let sizes = if chrom_names.is_empty() {
        None
    } else {
        Some(chrom_sizes_from_vecs(chrom_names, chrom_sizes_vec)?)
    };
    let core_prom_u32 = checked_u32(core_prom, "core_prom")?;
    let prox_prom_u32 = checked_u32(prox_prom, "prox_prom")?;
    let pl = genome_partition_list(&model, core_prom_u32, prox_prom_u32, sizes.as_ref());
    Ok(ExternalPtr::new(pl).into())
}

/// Build a PartitionList from RegionSet pointers with strand information
/// @export
/// @param genes_chrs Gene chromosome names
/// @param genes_starts Gene start positions (0-based)
/// @param genes_ends Gene end positions (half-open)
/// @param genes_strands Gene strand strings ("+", "-", "*")
/// @param exons_chrs Exon chromosome names
/// @param exons_starts Exon start positions (0-based)
/// @param exons_ends Exon end positions (half-open)
/// @param exons_strands Exon strand strings
/// @param three_utr_chrs 3'UTR chrs (empty vector if absent)
/// @param three_utr_starts 3'UTR starts
/// @param three_utr_ends 3'UTR ends
/// @param three_utr_strands 3'UTR strands
/// @param five_utr_chrs 5'UTR chrs (empty vector if absent)
/// @param five_utr_starts 5'UTR starts
/// @param five_utr_ends 5'UTR ends
/// @param five_utr_strands 5'UTR strands
/// @param core_prom Core promoter size in bp
/// @param prox_prom Proximal promoter size in bp
/// @param chrom_names Chromosome names for trim (empty to skip)
/// @param chrom_sizes_vec Chromosome sizes for trim (parallel to chrom_names)
#[extendr(r_name = "r_partition_list_from_regions_stranded")]
pub fn r_partition_list_from_regions_stranded(
    genes_chrs: Vec<String>,
    genes_starts: Vec<i32>,
    genes_ends: Vec<i32>,
    genes_strands: Vec<String>,
    exons_chrs: Vec<String>,
    exons_starts: Vec<i32>,
    exons_ends: Vec<i32>,
    exons_strands: Vec<String>,
    three_utr_chrs: Vec<String>,
    three_utr_starts: Vec<i32>,
    three_utr_ends: Vec<i32>,
    three_utr_strands: Vec<String>,
    five_utr_chrs: Vec<String>,
    five_utr_starts: Vec<i32>,
    five_utr_ends: Vec<i32>,
    five_utr_strands: Vec<String>,
    core_prom: i32,
    prox_prom: i32,
    chrom_names: Vec<String>,
    chrom_sizes_vec: Vec<i32>,
) -> extendr_api::Result<Robj> {
    let parse_strands = |sv: Vec<String>| -> Vec<Strand> {
        sv.iter()
            .map(|s| Strand::from_char(s.chars().next().unwrap_or('.')))
            .collect()
    };

    // Do NOT reduce here — genome_partition_list handles all reductions.
    // Pre-reducing genes would collapse overlapping genes before promoter
    // construction, losing individual gene promoters (R computes promoters
    // from raw genes, then reduces the promoters).
    let genes_rs = regionset_from_vecs(genes_chrs, genes_starts, genes_ends)?;
    let genes_srs = StrandedRegionSet::new(genes_rs, parse_strands(genes_strands));

    let exons_rs = regionset_from_vecs(exons_chrs, exons_starts, exons_ends)?;
    let exons_srs = StrandedRegionSet::new(exons_rs, parse_strands(exons_strands));

    let three_utr = if three_utr_chrs.is_empty() {
        None
    } else {
        let rs = regionset_from_vecs(three_utr_chrs, three_utr_starts, three_utr_ends)?;
        let srs = StrandedRegionSet::new(rs, parse_strands(three_utr_strands));
        if srs.is_empty() { None } else { Some(srs) }
    };

    let five_utr = if five_utr_chrs.is_empty() {
        None
    } else {
        let rs = regionset_from_vecs(five_utr_chrs, five_utr_starts, five_utr_ends)?;
        let srs = StrandedRegionSet::new(rs, parse_strands(five_utr_strands));
        if srs.is_empty() { None } else { Some(srs) }
    };

    let model = GeneModel {
        genes: genes_srs,
        exons: exons_srs,
        three_utr,
        five_utr,
    };

    let sizes = if chrom_names.is_empty() {
        None
    } else {
        Some(chrom_sizes_from_vecs(chrom_names, chrom_sizes_vec)?)
    };
    let core_prom_u32 = checked_u32(core_prom, "core_prom")?;
    let prox_prom_u32 = checked_u32(prox_prom, "prox_prom")?;
    let pl = genome_partition_list(&model, core_prom_u32, prox_prom_u32, sizes.as_ref());
    Ok(ExternalPtr::new(pl).into())
}

/// Build a PartitionList from a GTF file
/// @export
/// @param gtf_path Path to a GTF or GTF.gz file
/// @param filter_protein_coding Keep only protein-coding genes
/// @param convert_ensembl_ucsc Prepend "chr" to bare chromosome names
/// @param core_prom Core promoter size in bp
/// @param prox_prom Proximal promoter size in bp
/// @param chrom_names Chromosome names for trim (empty to skip)
/// @param chrom_sizes_vec Chromosome sizes for trim (parallel to chrom_names)
#[extendr(r_name = "r_partition_list_from_gtf")]
pub fn r_partition_list_from_gtf(
    gtf_path: &str,
    filter_protein_coding: bool,
    convert_ensembl_ucsc: bool,
    core_prom: i32,
    prox_prom: i32,
    chrom_names: Vec<String>,
    chrom_sizes_vec: Vec<i32>,
) -> extendr_api::Result<Robj> {
    let model = GeneModel::from_gtf(gtf_path, filter_protein_coding, convert_ensembl_ucsc)
        .map_err(|e| extendr_api::Error::Other(format!("Loading GTF: {}", e)))?;
    let sizes = if chrom_names.is_empty() {
        None
    } else {
        Some(chrom_sizes_from_vecs(chrom_names, chrom_sizes_vec)?)
    };
    let core_prom_u32 = checked_u32(core_prom, "core_prom")?;
    let prox_prom_u32 = checked_u32(prox_prom, "prox_prom")?;
    let pl = genome_partition_list(&model, core_prom_u32, prox_prom_u32, sizes.as_ref());
    Ok(ExternalPtr::new(pl).into())
}

/// Assign regions to genomic partitions
/// @export
/// @param rs_ptr External pointer to query RegionSet
/// @param partition_ptr External pointer to PartitionList
/// @param bp_proportion If TRUE, count base pairs; if FALSE, count regions
#[extendr(r_name = "r_calc_partitions")]
pub fn r_calc_partitions(
    rs_ptr: Robj,
    partition_ptr: Robj,
    bp_proportion: bool,
) -> extendr_api::Result<List> {
    with_regionset!(rs_ptr, rs, {
        with_partitions!(partition_ptr, pl, {
            let result = calc_partitions(rs, pl, bp_proportion);
            let names: Vec<String> = result.counts.iter().map(|(n, _)| n.clone()).collect();
            let counts: Vec<f64> = result.counts.iter().map(|(_, c)| *c as f64).collect();
            Ok(list!(
                partition = names,
                count = counts,
                total = result.total as f64
            ))
        })
    })
}

/// Calculate observed vs expected partition enrichment
/// @export
/// @param rs_ptr External pointer to query RegionSet
/// @param partition_ptr External pointer to PartitionList
/// @param chrom_names Character vector of chromosome names
/// @param chrom_sizes Integer vector of chromosome sizes
/// @param bp_proportion If TRUE, count base pairs; if FALSE, count regions
#[extendr(r_name = "r_calc_expected_partitions")]
pub fn r_calc_expected_partitions(
    rs_ptr: Robj,
    partition_ptr: Robj,
    chrom_names: Vec<String>,
    chrom_sizes: Vec<i32>,
    bp_proportion: bool,
) -> extendr_api::Result<List> {
    let sizes = chrom_sizes_from_vecs(chrom_names, chrom_sizes)?;
    with_regionset!(rs_ptr, rs, {
        with_partitions!(partition_ptr, pl, {
            let result = calc_expected_partitions(rs, pl, &sizes, bp_proportion);
            let partitions: Vec<String> = result.rows.iter().map(|r| r.partition.clone()).collect();
            let observed: Vec<f64> = result.rows.iter().map(|r| r.observed).collect();
            let expected: Vec<f64> = result.rows.iter().map(|r| r.expected).collect();
            let log10_oe: Vec<f64> = result.rows.iter().map(|r| r.log10_oe).collect();
            let pvals: Vec<f64> = result.rows.iter().map(|r| r.chi_sq_pval).collect();
            Ok(list!(
                partition = partitions,
                observed = observed,
                expected = expected,
                log10OE = log10_oe,
                pvalue = pvals
            ))
        })
    })
}

// =========================================================================
// 5. Signal
// =========================================================================

/// Calculate summary signal (overlap query with signal matrix)
/// @export
/// @param rs_ptr External pointer to query RegionSet
/// @param signal_region_ids Character vector of region IDs (chr_start_end format)
/// @param condition_names Character vector of condition/column names
/// @param values_flat Numeric vector of signal values (row-major, length = n_regions * n_conditions)
/// @param n_regions Number of signal regions
/// @param n_conditions Number of conditions
#[extendr(r_name = "r_calc_summary_signal")]
pub fn r_calc_summary_signal(
    rs_ptr: Robj,
    signal_region_ids: Vec<String>,
    condition_names: Vec<String>,
    values_flat: Vec<f64>,
    _n_regions: i32,
    n_conditions: i32,
) -> extendr_api::Result<List> {
    with_regionset!(rs_ptr, rs, {
        let n_cond = n_conditions as usize;

        // Parse region IDs into RegionSet
        let regions: Vec<Region> = signal_region_ids
            .iter()
            .filter_map(|id| {
                let parts: Vec<&str> = id.splitn(3, '_').collect();
                if parts.len() == 3 {
                    let chr = parts[0].to_string();
                    let start = parts[1].parse::<u32>().ok()?;
                    let end = parts[2].parse::<u32>().ok()?;
                    Some(Region {
                        chr,
                        start,
                        end,
                        rest: None,
                    })
                } else {
                    None
                }
            })
            .collect();

        let signal_matrix = SignalMatrix {
            regions: RegionSet::from(regions),
            condition_names: condition_names.clone(),
            n_conditions: n_cond,
            values: values_flat,
        };

        let result = calc_summary_signal(rs, &signal_matrix, CoordinateMode::GRanges)
            .map_err(|e| extendr_api::Error::Other(format!("{}", e)))?;

        // Build signal summary matrix as list of vectors
        let region_labels: Vec<String> = result.signal_matrix.iter().map(|(l, _)| l.clone()).collect();
        let mut signal_cols: Vec<Vec<f64>> = vec![Vec::new(); n_cond];
        for (_label, signals) in &result.signal_matrix {
            for (j, &val) in signals.iter().enumerate() {
                if j < n_cond {
                    signal_cols[j].push(val);
                }
            }
        }

        // Build matrixStats
        let stat_conditions: Vec<String> = result.matrix_stats.iter().map(|s| s.condition.clone()).collect();
        let lower_whiskers: Vec<f64> = result.matrix_stats.iter().map(|s| s.lower_whisker).collect();
        let lower_hinges: Vec<f64> = result.matrix_stats.iter().map(|s| s.lower_hinge).collect();
        let medians: Vec<f64> = result.matrix_stats.iter().map(|s| s.median).collect();
        let upper_hinges: Vec<f64> = result.matrix_stats.iter().map(|s| s.upper_hinge).collect();
        let upper_whiskers: Vec<f64> = result.matrix_stats.iter().map(|s| s.upper_whisker).collect();

        // Build the signal matrix as a named list
        let mut signal_list = List::from_values(&signal_cols);
        let _ = signal_list.set_names(condition_names.clone());

        Ok(list!(
            condition_names = condition_names,
            region_labels = region_labels,
            signal_matrix = signal_list,
            matrixStats = list!(
                condition = stat_conditions,
                lower_whisker = lower_whiskers,
                lower_hinge = lower_hinges,
                median = medians,
                upper_hinge = upper_hinges,
                upper_whisker = upper_whiskers
            )
        ))
    })
}

// =========================================================================
// 6. TSS / Feature Distances
// =========================================================================

/// Calculate absolute distance from each query region to nearest feature midpoint.
/// Returns NA for regions on chromosomes without features.
/// @export
/// @param query_ptr External pointer to query RegionSet
/// @param features_ptr External pointer to features RegionSet
#[extendr(r_name = "r_calc_tss_distances")]
pub fn r_calc_tss_distances(query_ptr: Robj, features_ptr: Robj) -> extendr_api::Result<Doubles> {
    with_regionset!(query_ptr, query, {
        let features_ext = <ExternalPtr<RegionSet>>::try_from(features_ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid features RegionSet pointer".into()))?;
        let index = TssIndex::from_region_set((*features_ext).clone(), CoordinateMode::GRanges)
            .map_err(|e| extendr_api::Error::Other(format!("Building TssIndex: {}", e)))?;
        let dists = index
            .calc_tss_distances(query, CoordinateMode::GRanges)
            .map_err(|e| extendr_api::Error::Other(format!("{}", e)))?;
        let result: Doubles = dists
            .into_iter()
            .map(|d| {
                if d == u32::MAX {
                    Rfloat::na()
                } else {
                    Rfloat::from(d as f64)
                }
            })
            .collect();
        Ok(result)
    })
}

/// Calculate signed distance from each query region to nearest feature.
/// Returns NA for regions on chromosomes without features.
/// @export
/// @param query_ptr External pointer to query RegionSet
/// @param features_ptr External pointer to features RegionSet
#[extendr(r_name = "r_calc_feature_distances")]
pub fn r_calc_feature_distances(
    query_ptr: Robj,
    features_ptr: Robj,
) -> extendr_api::Result<Doubles> {
    with_regionset!(query_ptr, query, {
        let features_ext = <ExternalPtr<RegionSet>>::try_from(features_ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid features RegionSet pointer".into()))?;
        let index = TssIndex::from_region_set((*features_ext).clone(), CoordinateMode::GRanges)
            .map_err(|e| extendr_api::Error::Other(format!("Building TssIndex: {}", e)))?;
        let dists = index
            .calc_feature_distances(query, CoordinateMode::GRanges)
            .map_err(|e| extendr_api::Error::Other(format!("{}", e)))?;
        let result: Doubles = dists
            .into_iter()
            .map(|d| {
                if d == i64::MAX {
                    Rfloat::na()
                } else {
                    Rfloat::from(d as f64)
                }
            })
            .collect();
        Ok(result)
    })
}

// =========================================================================
// 7. Binary Asset Loading
// =========================================================================

/// Load a GDA binary file and return a GenomicDistAnnotation pointer.
///
/// The returned pointer gives access to the gene model, partition list,
/// and TSS index contained in the pre-compiled binary.
/// @export
/// @param path Path to a .bin GDA file (produced by `gtars prep`)
#[extendr(r_name = "load_gda_bin")]
pub fn r_load_gda_bin(path: &str) -> extendr_api::Result<Robj> {
    let ann = GenomicDistAnnotation::load_bin(path)
        .map_err(|e| extendr_api::Error::Other(format!("Loading GDA binary: {}", e)))?;
    Ok(ExternalPtr::new(ann).into())
}

/// Extract the gene model from a GDA annotation pointer.
/// @export
/// @param gda_ptr External pointer to a GenomicDistAnnotation
#[extendr(r_name = "gda_gene_model")]
pub fn r_gda_gene_model(gda_ptr: Robj) -> extendr_api::Result<Robj> {
    let ext_ptr = <ExternalPtr<GenomicDistAnnotation>>::try_from(gda_ptr)
        .map_err(|_| extendr_api::Error::Other("Invalid GenomicDistAnnotation pointer".into()))?;
    let model = ext_ptr.gene_model.clone();
    Ok(ExternalPtr::new(model).into())
}

/// Build a PartitionList from a GDA annotation pointer.
/// @export
/// @param gda_ptr External pointer to a GenomicDistAnnotation
/// @param core_prom Core promoter size in bp
/// @param prox_prom Proximal promoter size in bp
/// @param chrom_names Chromosome names for trim (empty to skip)
/// @param chrom_sizes_vec Chromosome sizes for trim (parallel to chrom_names)
#[extendr(r_name = "gda_partition_list")]
pub fn r_gda_partition_list(
    gda_ptr: Robj,
    core_prom: i32,
    prox_prom: i32,
    chrom_names: Vec<String>,
    chrom_sizes_vec: Vec<i32>,
) -> extendr_api::Result<Robj> {
    let ext_ptr = <ExternalPtr<GenomicDistAnnotation>>::try_from(gda_ptr)
        .map_err(|_| extendr_api::Error::Other("Invalid GenomicDistAnnotation pointer".into()))?;
    let sizes = if chrom_names.is_empty() {
        None
    } else {
        Some(chrom_sizes_from_vecs(chrom_names, chrom_sizes_vec)?)
    };
    let core_prom_u32 = checked_u32(core_prom, "core_prom")?;
    let prox_prom_u32 = checked_u32(prox_prom, "prox_prom")?;
    let pl = genome_partition_list(&ext_ptr.gene_model, core_prom_u32, prox_prom_u32, sizes.as_ref());
    Ok(ExternalPtr::new(pl).into())
}

/// Derive a TssIndex from a GDA annotation (using strand-aware TSS positions).
/// @export
/// @param gda_ptr External pointer to a GenomicDistAnnotation
#[extendr(r_name = "gda_tss_index")]
pub fn r_gda_tss_index(gda_ptr: Robj) -> extendr_api::Result<Robj> {
    let ext_ptr = <ExternalPtr<GenomicDistAnnotation>>::try_from(gda_ptr)
        .map_err(|_| extendr_api::Error::Other("Invalid GenomicDistAnnotation pointer".into()))?;
    let model = &ext_ptr.gene_model;
    let tss_regions: Vec<Region> = model.genes.inner.regions.iter()
        .zip(model.genes.strands.iter())
        .map(|(r, strand)| {
            let tss_pos = match strand {
                Strand::Minus => r.end.saturating_sub(1),
                _ => r.start,
            };
            Region { chr: r.chr.clone(), start: tss_pos, end: tss_pos + 1, rest: None }
        })
        .collect();
    let tss_rs = RegionSet { regions: tss_regions, header: None, path: None };
    let index = TssIndex::try_from(tss_rs)
        .map_err(|e| extendr_api::Error::Other(format!("Building TssIndex from GDA: {}", e)))?;
    Ok(ExternalPtr::new(index).into())
}

/// Load a SignalMatrix from a packed binary file.
/// @export
/// @param path Path to a .bin signal matrix file (produced by `gtars prep`)
#[extendr(r_name = "load_signal_matrix_bin")]
pub fn r_load_signal_matrix_bin(path: &str) -> extendr_api::Result<Robj> {
    let sm = SignalMatrix::load_bin(path)
        .map_err(|e| extendr_api::Error::Other(format!("Loading signal matrix binary: {}", e)))?;
    Ok(ExternalPtr::new(sm).into())
}

/// Load a SignalMatrix from a TSV file.
/// @export
/// @param path Path to a TSV signal matrix file
#[extendr(r_name = "load_signal_matrix_tsv")]
pub fn r_load_signal_matrix_tsv(path: &str) -> extendr_api::Result<Robj> {
    let sm = SignalMatrix::from_tsv(path)
        .map_err(|e| extendr_api::Error::Other(format!("Loading signal matrix TSV: {}", e)))?;
    Ok(ExternalPtr::new(sm).into())
}

/// Calculate summary signal from a SignalMatrix pointer.
/// @export
/// @param rs_ptr External pointer to query RegionSet
/// @param sm_ptr External pointer to a SignalMatrix
#[extendr(r_name = "r_calc_summary_signal_from_matrix")]
pub fn r_calc_summary_signal_from_matrix(
    rs_ptr: Robj,
    sm_ptr: Robj,
) -> extendr_api::Result<List> {
    with_regionset!(rs_ptr, rs, {
        let sm_ext = <ExternalPtr<SignalMatrix>>::try_from(sm_ptr)
            .map_err(|_| extendr_api::Error::Other("Invalid SignalMatrix pointer".into()))?;

        let result = calc_summary_signal(rs, &*sm_ext, CoordinateMode::GRanges)
            .map_err(|e| extendr_api::Error::Other(format!("{}", e)))?;

        let n_cond = result.condition_names.len();
        let region_labels: Vec<String> = result.signal_matrix.iter().map(|(l, _)| l.clone()).collect();
        let mut signal_cols: Vec<Vec<f64>> = vec![Vec::new(); n_cond];
        for (_label, signals) in &result.signal_matrix {
            for (j, &val) in signals.iter().enumerate() {
                if j < n_cond {
                    signal_cols[j].push(val);
                }
            }
        }

        let stat_conditions: Vec<String> = result.matrix_stats.iter().map(|s| s.condition.clone()).collect();
        let lower_whiskers: Vec<f64> = result.matrix_stats.iter().map(|s| s.lower_whisker).collect();
        let lower_hinges: Vec<f64> = result.matrix_stats.iter().map(|s| s.lower_hinge).collect();
        let medians: Vec<f64> = result.matrix_stats.iter().map(|s| s.median).collect();
        let upper_hinges: Vec<f64> = result.matrix_stats.iter().map(|s| s.upper_hinge).collect();
        let upper_whiskers: Vec<f64> = result.matrix_stats.iter().map(|s| s.upper_whisker).collect();

        let mut signal_list = List::from_values(&signal_cols);
        let _ = signal_list.set_names(result.condition_names.clone());

        Ok(list!(
            condition_names = result.condition_names,
            region_labels = region_labels,
            signal_matrix = signal_list,
            matrixStats = list!(
                condition = stat_conditions,
                lower_whisker = lower_whiskers,
                lower_hinge = lower_hinges,
                median = medians,
                upper_hinge = upper_hinges,
                upper_whisker = upper_whiskers
            )
        ))
    })
}

// =========================================================================
// Module registration
// =========================================================================

extendr_module! {
    mod genomicdist;
    // RegionSet core
    fn r_load_regionset;
    fn r_regionset_from_vectors;
    fn r_regionset_to_vectors;
    fn r_regionset_length;
    // Statistics
    fn r_calc_widths;
    fn r_calc_neighbor_distances;
    fn r_calc_nearest_neighbors;
    fn r_chromosome_statistics;
    fn r_region_distribution;
    // GC / Dinucleotide
    fn r_load_genome_assembly;
    fn r_calc_gc_content;
    fn r_calc_dinucl_freq;
    // Interval Ranges
    fn r_trim;
    fn r_promoters;
    fn r_reduce;
    fn r_setdiff;
    fn r_pintersect;
    fn r_concat;
    fn r_union;
    fn r_jaccard;
    fn r_shift;
    fn r_flank;
    fn r_resize;
    fn r_narrow;
    fn r_disjoin;
    fn r_gaps;
    fn r_intersect;
    fn r_find_overlaps;
    fn r_count_overlaps;
    fn r_consensus;
    // Partitions
    fn r_partition_list_from_regions;
    fn r_partition_list_from_regions_stranded;
    fn r_partition_list_from_gtf;
    fn r_calc_partitions;
    fn r_calc_expected_partitions;
    // Signal
    fn r_calc_summary_signal;
    // TSS / Feature distances
    fn r_calc_tss_distances;
    fn r_calc_feature_distances;
    // Binary asset loading
    fn r_load_gda_bin;
    fn r_gda_gene_model;
    fn r_gda_partition_list;
    fn r_gda_tss_index;
    fn r_load_signal_matrix_bin;
    fn r_load_signal_matrix_tsv;
    fn r_calc_summary_signal_from_matrix;
}

// =========================================================================
// Unit tests (pure Rust, no R runtime needed)
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // --- checked_u32 helper tests ---

    #[test]
    fn test_checked_u32_zero() {
        assert_eq!(checked_u32(0, "x").unwrap(), 0u32);
    }

    #[test]
    fn test_checked_u32_positive() {
        assert_eq!(checked_u32(42, "x").unwrap(), 42u32);
    }

    #[test]
    fn test_checked_u32_max_i32() {
        assert_eq!(checked_u32(i32::MAX, "x").unwrap(), i32::MAX as u32);
    }

    #[test]
    fn test_checked_u32_negative_one() {
        let err = checked_u32(-1, "start").unwrap_err();
        let msg = format!("{}", err);
        assert!(msg.contains("start"), "error should name the parameter");
        assert!(msg.contains("-1"), "error should include the value");
    }

    #[test]
    fn test_checked_u32_min_i32() {
        let err = checked_u32(i32::MIN, "end").unwrap_err();
        let msg = format!("{}", err);
        assert!(msg.contains("end"));
        assert!(msg.contains(&i32::MIN.to_string()));
    }

    // --- regionset_from_vecs tests ---

    #[test]
    fn test_regionset_from_vecs_valid() {
        let rs = regionset_from_vecs(
            vec!["chr1".into(), "chr2".into()],
            vec![0, 100],
            vec![50, 200],
        )
        .unwrap();
        assert_eq!(rs.regions.len(), 2);
        assert_eq!(rs.regions[0].start, 0);
        assert_eq!(rs.regions[0].end, 50);
        assert_eq!(rs.regions[1].start, 100);
        assert_eq!(rs.regions[1].end, 200);
    }

    #[test]
    fn test_regionset_from_vecs_negative_start() {
        let result = regionset_from_vecs(
            vec!["chr1".into()],
            vec![-5],
            vec![100],
        );
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("start"));
        assert!(msg.contains("-5"));
    }

    #[test]
    fn test_regionset_from_vecs_negative_end() {
        let result = regionset_from_vecs(
            vec!["chr1".into()],
            vec![0],
            vec![-1],
        );
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("end"));
    }

    // --- chrom_sizes_from_vecs tests ---

    #[test]
    fn test_chrom_sizes_from_vecs_valid() {
        let sizes = chrom_sizes_from_vecs(
            vec!["chr1".into(), "chr2".into()],
            vec![248956422, 242193529],
        )
        .unwrap();
        assert_eq!(sizes["chr1"], 248956422);
        assert_eq!(sizes["chr2"], 242193529);
    }

    #[test]
    fn test_chrom_sizes_from_vecs_negative() {
        let result = chrom_sizes_from_vecs(
            vec!["chr1".into()],
            vec![-100],
        );
        assert!(result.is_err());
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("chrom_size"));
    }
}
