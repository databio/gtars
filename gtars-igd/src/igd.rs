//! Unified IGD (Integrated Genome Database) implementation.
//!
//! Provides a single `Igd` struct that supports both construction and querying
//! in-memory, with optional disk persistence. This replaces the split between
//! `igd_t` (creation-only) and `igd_t_from_disk` (query-only) in the legacy API.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write as IoWrite};
use std::path::{Path, PathBuf};

use byteorder::{LittleEndian, ReadBytesExt};
use gtars_core::consts::{BED_FILE_EXTENSION, GZ_FILE_EXTENSION};
use gtars_core::models::RegionSet;
use gtars_core::utils::get_dynamic_reader;

use crate::create::MAX_CHROM_NAME_LEN;

/// A single genomic interval record within the IGD index.
#[derive(Default, Clone, Debug)]
pub struct Record {
    /// Index of the source file (DB set) this interval came from.
    pub file_idx: u32,
    /// Region start position (0-based).
    pub start: i32,
    /// Region end position (exclusive).
    pub end: i32,
    /// BED score value (preserved for backward compat; LOLA ignores this).
    pub value: i32,
}

/// A tile (bin) within a contig, holding records that start in this tile's range.
#[derive(Default, Clone, Debug)]
pub struct Tile {
    /// Genomic interval records in this tile, sorted by start position after finalization.
    pub records: Vec<Record>,
}

/// A contig (chromosome) within the IGD index.
#[derive(Default, Clone, Debug)]
pub struct Contig {
    /// Chromosome name (e.g., "chr1").
    pub name: String,
    /// Tiles for this contig. Tile i covers positions [i*nbp, (i+1)*nbp).
    pub tiles: Vec<Tile>,
}

/// Metadata about a source file indexed into the IGD.
#[derive(Default, Clone, Debug)]
pub struct FileInfo {
    /// Filename of the source BED file.
    pub filename: String,
    /// Number of regions from this file.
    pub num_regions: u32,
    /// Average region width in base pairs.
    pub avg_region_width: f64,
}

/// Unified IGD index: can be built in memory or loaded from disk, and queried directly.
#[derive(Debug)]
pub struct Igd {
    /// Tile/bin size in base pairs (default 16384 = 2^14).
    pub nbp: i32,
    /// Per-chromosome data.
    pub contigs: Vec<Contig>,
    /// Per-file metadata (filename, region count, avg width).
    pub file_info: Vec<FileInfo>,
    /// Chromosome name → index into `contigs`.
    chrom_index: HashMap<String, usize>,
    /// Whether tiles have been sorted (required before querying).
    finalized: bool,
}

impl Igd {
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Create a new empty IGD with default tile size (16384bp).
    pub fn new() -> Self {
        Self::with_tile_size(16384)
    }

    /// Create a new empty IGD with a custom tile size.
    pub fn with_tile_size(nbp: i32) -> Self {
        Igd {
            nbp,
            contigs: Vec::new(),
            file_info: Vec::new(),
            chrom_index: HashMap::new(),
            finalized: false,
        }
    }

    /// Add a single interval from a given file index. The IGD must not yet be finalized.
    ///
    /// Intervals with `start >= end` or negative coordinates are silently skipped.
    ///
    /// # Panics
    ///
    /// Panics if called after [`finalize`](Self::finalize).
    pub fn add(&mut self, chrom: &str, start: i32, end: i32, value: i32, file_idx: u32) {
        assert!(
            !self.finalized,
            "Cannot add intervals after finalization"
        );
        if start < 0 || end < 0 || start >= end {
            return;
        }

        let n1 = start / self.nbp;
        let n2 = (end - 1) / self.nbp;
        let needed_tiles = (n2 + 1) as usize;

        // Get or create contig
        let ctg_idx = if let Some(&idx) = self.chrom_index.get(chrom) {
            idx
        } else {
            let idx = self.contigs.len();
            self.contigs.push(Contig {
                name: chrom.to_string(),
                tiles: Vec::new(),
            });
            self.chrom_index.insert(chrom.to_string(), idx);
            idx
        };

        let contig = &mut self.contigs[ctg_idx];

        // Expand tiles if needed
        if contig.tiles.len() < needed_tiles {
            contig.tiles.resize_with(needed_tiles, Tile::default);
        }

        // Add record to each spanned tile
        let record = Record {
            file_idx,
            start,
            end,
            value,
        };

        for i in n1..=n2 {
            contig.tiles[i as usize].records.push(record.clone());
        }
    }

    /// Finalize the IGD: sort all tile records by start position.
    /// Must be called after all intervals are added and before any queries.
    pub fn finalize(&mut self) {
        if self.finalized {
            return;
        }
        for contig in &mut self.contigs {
            for tile in &mut contig.tiles {
                tile.records.sort_by_key(|r| r.start);
            }
        }
        self.finalized = true;
    }

    /// Build an IGD from a directory of BED files.
    pub fn from_bed_dir(path: &Path) -> anyhow::Result<Self> {
        let mut bed_files: Vec<PathBuf> = Vec::new();

        // Collect BED/gz files (validation happens during the parse pass)
        for entry in fs::read_dir(path)? {
            let p = entry?.path();
            if let Some(ext) = p.extension().and_then(|e| e.to_str()) {
                if (ext == BED_FILE_EXTENSION.trim_start_matches('.')
                    || ext == GZ_FILE_EXTENSION.trim_start_matches('.'))
                    && p.is_file()
                {
                    bed_files.push(p);
                }
            }
        }

        bed_files.sort(); // deterministic order

        let mut igd = Igd::new();
        let mut file_infos: Vec<FileInfo> = Vec::new();

        for bed_path in &bed_files {
            let reader = match get_dynamic_reader(bed_path) {
                Ok(r) => r,
                Err(_) => continue,
            };
            let mut count: u32 = 0;
            let mut total_width: u64 = 0;
            let mut has_valid_line = false;
            let file_idx = file_infos.len();

            for line in reader.lines() {
                let line = match line {
                    Ok(l) => l,
                    Err(_) => continue,
                };
                if let Some((chrom, start, end, score)) = Self::parse_bed_line(&line) {
                    has_valid_line = true;
                    if start >= 0 {
                        igd.add(&chrom, start, end, score, file_idx as u32);
                        count += 1;
                        total_width += (end - start) as u64;
                    }
                }
            }

            // Skip files with no parseable BED lines
            if !has_valid_line {
                continue;
            }

            let filename = bed_path
                .file_name()
                .map(|f| f.to_string_lossy().into_owned())
                .unwrap_or_default();

            file_infos.push(FileInfo {
                filename,
                num_regions: count,
                avg_region_width: if count > 0 {
                    total_width as f64 / count as f64
                } else {
                    0.0
                },
            });
        }

        igd.file_info = file_infos;
        igd.finalize();
        Ok(igd)
    }

    /// Build an IGD from an iterator of (filename, regions) pairs.
    /// Each region is (chrom, start, end). For BEDbase integration where
    /// regions come from an API rather than disk.
    pub fn from_region_sets<I>(sets: I) -> Self
    where
        I: IntoIterator<Item = (String, Vec<(String, i32, i32)>)>,
    {
        let mut igd = Igd::new();
        let mut file_infos: Vec<FileInfo> = Vec::new();

        for (file_idx, (filename, regions)) in sets.into_iter().enumerate() {
            let mut count: u32 = 0;
            let mut total_width: u64 = 0;

            for (chrom, start, end) in &regions {
                if *start < *end {
                    igd.add(chrom, *start, *end, 0, file_idx as u32);
                    count += 1;
                    total_width += (*end - *start) as u64;
                }
            }

            file_infos.push(FileInfo {
                filename,
                num_regions: count,
                avg_region_width: if count > 0 {
                    total_width as f64 / count as f64
                } else {
                    0.0
                },
            });
        }

        igd.file_info = file_infos;
        igd.finalize();
        igd
    }

    /// Build an IGD from a vec of gtars-core RegionSets with associated filenames.
    pub fn from_named_region_sets(sets: &[(String, &RegionSet)]) -> Self {
        let mut igd = Igd::new();
        let mut file_infos: Vec<FileInfo> = Vec::with_capacity(sets.len());

        for (file_idx, (filename, region_set)) in sets.iter().enumerate() {
            let mut count: u32 = 0;
            let mut total_width: u64 = 0;

            for region in &region_set.regions {
                if region.start < region.end {
                    let start = region.start as i32;
                    let end = region.end as i32;
                    igd.add(&region.chr, start, end, 0, file_idx as u32);
                    count += 1;
                    total_width += (end - start) as u64;
                }
            }

            file_infos.push(FileInfo {
                filename: filename.clone(),
                num_regions: count,
                avg_region_width: if count > 0 {
                    total_width as f64 / count as f64
                } else {
                    0.0
                },
            });
        }

        igd.file_info = file_infos;
        igd.finalize();
        igd
    }

    /// Load an IGD from a pre-built `.igd` file on disk.
    pub fn from_igd_file(path: &Path) -> anyhow::Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Read header
        let mut buf4 = [0u8; 4];
        reader.read_exact(&mut buf4)?;
        let nbp = i32::from_le_bytes(buf4);
        reader.read_exact(&mut buf4)?;
        let g_type = i32::from_le_bytes(buf4);
        reader.read_exact(&mut buf4)?;
        let n_ctg = i32::from_le_bytes(buf4);

        // Read tiles-per-contig
        let mut n_tiles: Vec<i32> = Vec::with_capacity(n_ctg as usize);
        for _ in 0..n_ctg {
            reader.read_exact(&mut buf4)?;
            n_tiles.push(i32::from_le_bytes(buf4));
        }

        // Read region counts per tile
        let mut n_cnts: Vec<Vec<i32>> = Vec::with_capacity(n_ctg as usize);
        for i in 0..n_ctg as usize {
            let k = n_tiles[i];
            let mut counts = Vec::with_capacity(k as usize);
            for _ in 0..k {
                counts.push(reader.read_i32::<LittleEndian>()?);
            }
            n_cnts.push(counts);
        }

        // Read chromosome names
        let mut chrom_names: Vec<String> = Vec::with_capacity(n_ctg as usize);
        for _ in 0..n_ctg {
            let mut buf = [0u8; MAX_CHROM_NAME_LEN];
            reader.read_exact(&mut buf)?;
            let name = String::from_utf8_lossy(&buf)
                .trim_matches('\0')
                .to_string();
            chrom_names.push(name);
        }

        // Determine record size
        let rec_size: usize = if g_type == 0 { 12 } else { 16 };

        // Read region data into tiles
        let mut contigs: Vec<Contig> = Vec::with_capacity(n_ctg as usize);
        let mut chrom_index: HashMap<String, usize> = HashMap::new();

        for i in 0..n_ctg as usize {
            let mut tiles: Vec<Tile> = Vec::with_capacity(n_tiles[i] as usize);

            for j in 0..n_tiles[i] as usize {
                let count = n_cnts[i][j];
                let mut records = Vec::with_capacity(count as usize);

                for _ in 0..count {
                    let idx = reader.read_i32::<LittleEndian>()?;
                    let start = reader.read_i32::<LittleEndian>()?;
                    let end = reader.read_i32::<LittleEndian>()?;
                    let value = if rec_size == 16 {
                        reader.read_i32::<LittleEndian>()?
                    } else {
                        0
                    };
                    records.push(Record {
                        file_idx: idx as u32,
                        start,
                        end,
                        value,
                    });
                }

                tiles.push(Tile { records });
            }

            chrom_index.insert(chrom_names[i].clone(), i);
            contigs.push(Contig {
                name: chrom_names[i].clone(),
                tiles,
            });
        }

        // Load file info from the companion .tsv file
        let tsv_path = path.with_extension("tsv");
        let file_info = if tsv_path.exists() {
            Self::load_file_info_tsv(&tsv_path)?
        } else {
            Vec::new()
        };

        Ok(Igd {
            nbp,
            contigs,
            file_info,
            chrom_index,
            finalized: true, // disk data is already sorted
        })
    }

    /// Save the IGD to a `.igd` binary file and companion `.tsv` metadata.
    ///
    /// # Panics
    ///
    /// Panics if called before [`finalize`](Self::finalize).
    pub fn save(&self, path: &Path) -> anyhow::Result<()> {
        assert!(self.finalized, "Must finalize before saving");

        // Write main .igd file
        let mut buffer = Vec::new();
        let g_type: i32 = 1; // gType=1 means 16-byte records with value field

        // Header
        buffer.write_all(&self.nbp.to_le_bytes())?;
        buffer.write_all(&g_type.to_le_bytes())?;
        buffer.write_all(&(self.contigs.len() as i32).to_le_bytes())?;

        // Tiles per contig
        for contig in &self.contigs {
            buffer.write_all(&(contig.tiles.len() as i32).to_le_bytes())?;
        }

        // Region counts per tile
        for contig in &self.contigs {
            for tile in &contig.tiles {
                buffer.write_all(&(tile.records.len() as i32).to_le_bytes())?;
            }
        }

        // Chromosome names (40 bytes each, null-padded)
        for contig in &self.contigs {
            let mut name_bytes = contig.name.as_bytes().to_vec();
            name_bytes.resize(MAX_CHROM_NAME_LEN, 0);
            buffer.write_all(&name_bytes)?;
        }

        // Region data (sorted within each tile)
        for contig in &self.contigs {
            for tile in &contig.tiles {
                for rec in &tile.records {
                    buffer.write_all(&(rec.file_idx as i32).to_le_bytes())?;
                    buffer.write_all(&rec.start.to_le_bytes())?;
                    buffer.write_all(&rec.end.to_le_bytes())?;
                    buffer.write_all(&rec.value.to_le_bytes())?;
                }
            }
        }

        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(path, &buffer)?;

        // Write companion .tsv
        let tsv_path = path.with_extension("tsv");
        let mut tsv = String::new();
        tsv.push_str("Index\tFile\tNumber of Regions\tAvg size\n");
        for (i, fi) in self.file_info.iter().enumerate() {
            tsv.push_str(&format!(
                "{}\t{}\t{}\t{:.2}\n",
                i, fi.filename, fi.num_regions, fi.avg_region_width
            ));
        }
        fs::write(tsv_path, tsv)?;

        Ok(())
    }

    // -----------------------------------------------------------------------
    // Querying
    // -----------------------------------------------------------------------

    /// Query a single interval, incrementing hit counts per file.
    /// Returns the number of overlaps found.
    ///
    /// `hits` must have length >= `self.file_info.len()`.
    /// `min_overlap` is the minimum number of overlapping base pairs required (default 1).
    ///
    /// Intervals with `start >= end` or `end <= 0` return 0 immediately.
    /// Negative `start` is clamped to 0.
    ///
    /// # Panics
    ///
    /// Panics if called before [`finalize`](Self::finalize).
    pub fn count_overlaps(
        &self,
        chrom: &str,
        start: i32,
        end: i32,
        min_overlap: i32,
        hits: &mut [u64],
    ) -> u32 {
        assert!(self.finalized, "Must finalize before querying");

        if start >= end || end <= 0 {
            return 0;
        }
        let start = start.max(0);

        let ctg_idx = match self.chrom_index.get(chrom) {
            Some(&idx) => idx,
            None => return 0,
        };

        let contig = &self.contigs[ctg_idx];
        let n_tiles = contig.tiles.len() as i32;

        let n1 = start / self.nbp;
        let mut n2 = (end - 1) / self.nbp;

        if n1 >= n_tiles {
            return 0;
        }
        n2 = n2.min(n_tiles - 1);

        let mut total_overlaps: u32 = 0;

        // First tile (n1): binary search + backward scan
        let tile = &contig.tiles[n1 as usize];
        if !tile.records.is_empty() && end > tile.records[0].start {
            // Binary search: find rightmost record with start < end
            let mut tl: i32 = 0;
            let mut tr: i32 = tile.records.len() as i32 - 1;

            while tl < tr - 1 {
                let tm = (tl + tr) / 2;
                if tile.records[tm as usize].start < end {
                    tl = tm;
                } else {
                    tr = tm;
                }
            }
            if tile.records[tr as usize].start < end {
                tl = tr;
            }

            // Scan backward checking overlap
            for i in (0..=tl).rev() {
                let rec = &tile.records[i as usize];
                let overlap_bp = rec.end.min(end) - rec.start.max(start);
                if overlap_bp >= min_overlap {
                    hits[rec.file_idx as usize] += 1;
                    total_overlaps += 1;
                }
            }
        }

        // Subsequent tiles (n1+1 through n2)
        if n2 > n1 {
            let mut bd = self.nbp * (n1 + 1);

            for j in (n1 + 1)..=n2 {
                let tile = &contig.tiles[j as usize];
                if tile.records.is_empty() {
                    bd += self.nbp;
                    continue;
                }

                if end > tile.records[0].start {
                    // Skip records that start before this tile's boundary
                    // (they were already counted in a previous tile)
                    let mut ts: i32 = 0;
                    while ts < tile.records.len() as i32
                        && tile.records[ts as usize].start < bd
                    {
                        ts += 1;
                    }

                    // Binary search for rightmost record with start < end
                    let mut tl: i32 = 0;
                    let mut tr: i32 = tile.records.len() as i32 - 1;

                    while tl < tr - 1 {
                        let tm = (tl + tr) / 2;
                        if tile.records[tm as usize].start < end {
                            tl = tm;
                        } else {
                            tr = tm;
                        }
                    }
                    if tile.records[tr as usize].start < end {
                        tl = tr;
                    }

                    // Scan from ts to tl checking overlap
                    for i in (ts..=tl).rev() {
                        let rec = &tile.records[i as usize];
                        let overlap_bp = rec.end.min(end) - rec.start.max(start);
                        if overlap_bp >= min_overlap {
                            hits[rec.file_idx as usize] += 1;
                            total_overlaps += 1;
                        }
                    }
                }

                bd += self.nbp;
            }
        }

        total_overlaps
    }

    /// Query all regions in a set, returning total pairwise hit counts per file.
    /// Each query region can contribute multiple hits to the same file.
    pub fn count_set_overlaps(&self, regions: &RegionSet, min_overlap: i32) -> Vec<u64> {
        let mut hits = vec![0u64; self.file_info.len()];
        for region in &regions.regions {
            self.count_overlaps(
                &region.chr,
                region.start as i32,
                region.end as i32,
                min_overlap,
                &mut hits,
            );
        }
        hits
    }

    /// Count the number of query regions that overlap each file (binary per query).
    ///
    /// Unlike `count_set_overlaps` which counts total pairwise overlaps,
    /// this counts at most 1 per query region per file. This matches
    /// R LOLA's `countOverlaps()` semantics where `support = sum(countOverlaps > 0)`.
    pub fn count_region_hits(&self, regions: &RegionSet, min_overlap: i32) -> Vec<u64> {
        let n_files = self.file_info.len();
        let mut totals = vec![0u64; n_files];
        let mut per_region = vec![0u64; n_files];

        for region in &regions.regions {
            // Zero out per-region hits
            for h in per_region.iter_mut() {
                *h = 0;
            }

            self.count_overlaps(
                &region.chr,
                region.start as i32,
                region.end as i32,
                min_overlap,
                &mut per_region,
            );

            // Binary: if this region hit file i at all, count it once
            for i in 0..n_files {
                if per_region[i] > 0 {
                    totals[i] += 1;
                }
            }
        }
        totals
    }

    /// Query all regions given as (chrom, start, end) tuples.
    pub fn count_regions_overlaps(
        &self,
        regions: &[(String, i32, i32)],
        min_overlap: i32,
    ) -> Vec<u64> {
        let mut hits = vec![0u64; self.file_info.len()];
        for (chrom, start, end) in regions {
            self.count_overlaps(chrom, *start, *end, min_overlap, &mut hits);
        }
        hits
    }

    /// Build an IGD from a single RegionSet (for two-set overlap queries).
    ///
    /// Stores the original region index in the `value` field of each record,
    /// enabling `find_overlaps_regionset` to return subject indices.
    pub fn from_single_region_set(rs: &RegionSet) -> Self {
        let mut igd = Igd::new();
        igd.file_info = vec![FileInfo {
            filename: String::new(),
            num_regions: rs.regions.len() as u32,
            avg_region_width: if rs.regions.is_empty() {
                0.0
            } else {
                rs.regions.iter().map(|r| (r.end - r.start) as f64).sum::<f64>()
                    / rs.regions.len() as f64
            },
        }];

        for (i, region) in rs.regions.iter().enumerate() {
            igd.add(
                &region.chr,
                region.start as i32,
                region.end as i32,
                i as i32, // store original index in value field
                0,
            );
        }

        igd.finalize();
        igd
    }

    /// Find all overlapping (query_idx, subject_idx) pairs between a query
    /// RegionSet and the subject indexed in this IGD.
    ///
    /// The IGD must have been built with `from_single_region_set` so that
    /// `record.value` stores the original subject region index.
    ///
    /// # Panics
    ///
    /// Panics if called before [`finalize`](Self::finalize).
    pub fn find_overlaps_regionset(
        &self,
        query: &RegionSet,
        min_overlap: i32,
    ) -> Vec<(u32, u32)> {
        assert!(self.finalized, "Must finalize before querying");

        let mut pairs: Vec<(u32, u32)> = Vec::new();

        for (q_idx, region) in query.regions.iter().enumerate() {
            let ctg_idx = match self.chrom_index.get(&region.chr) {
                Some(&idx) => idx,
                None => continue,
            };

            let start = region.start as i32;
            let end = region.end as i32;
            let contig = &self.contigs[ctg_idx];
            let n_tiles = contig.tiles.len() as i32;

            let n1 = start / self.nbp;
            let mut n2 = (end - 1) / self.nbp;
            if n1 >= n_tiles {
                continue;
            }
            n2 = n2.min(n_tiles - 1);

            // Use a set to deduplicate subject hits (records span multiple tiles)
            let mut seen_subjects = std::collections::HashSet::new();

            // First tile
            let tile = &contig.tiles[n1 as usize];
            if !tile.records.is_empty() && end > tile.records[0].start {
                let mut tl: i32 = 0;
                let mut tr: i32 = tile.records.len() as i32 - 1;
                while tl < tr - 1 {
                    let tm = (tl + tr) / 2;
                    if tile.records[tm as usize].start < end {
                        tl = tm;
                    } else {
                        tr = tm;
                    }
                }
                if tile.records[tr as usize].start < end {
                    tl = tr;
                }
                for i in (0..=tl).rev() {
                    let rec = &tile.records[i as usize];
                    let overlap_bp = rec.end.min(end) - rec.start.max(start);
                    if overlap_bp >= min_overlap && seen_subjects.insert(rec.value as u32) {
                        pairs.push((q_idx as u32, rec.value as u32));
                    }
                }
            }

            // Subsequent tiles
            if n2 > n1 {
                let mut bd = self.nbp * (n1 + 1);
                for j in (n1 + 1)..=n2 {
                    let tile = &contig.tiles[j as usize];
                    if tile.records.is_empty() {
                        bd += self.nbp;
                        continue;
                    }
                    if end > tile.records[0].start {
                        let mut ts: i32 = 0;
                        while ts < tile.records.len() as i32
                            && tile.records[ts as usize].start < bd
                        {
                            ts += 1;
                        }
                        let mut tl: i32 = 0;
                        let mut tr: i32 = tile.records.len() as i32 - 1;
                        while tl < tr - 1 {
                            let tm = (tl + tr) / 2;
                            if tile.records[tm as usize].start < end {
                                tl = tm;
                            } else {
                                tr = tm;
                            }
                        }
                        if tile.records[tr as usize].start < end {
                            tl = tr;
                        }
                        for i in (ts..=tl).rev() {
                            let rec = &tile.records[i as usize];
                            let overlap_bp = rec.end.min(end) - rec.start.max(start);
                            if overlap_bp >= min_overlap
                                && seen_subjects.insert(rec.value as u32)
                            {
                                pairs.push((q_idx as u32, rec.value as u32));
                            }
                        }
                    }
                    bd += self.nbp;
                }
            }
        }

        pairs
    }

    /// Count the number of subject regions overlapping each query region.
    ///
    /// Returns a `Vec<u32>` of length `query.regions.len()` where each entry
    /// is the number of distinct subject regions overlapping that query region.
    ///
    /// The IGD must have been built with `from_single_region_set`.
    ///
    /// # Panics
    ///
    /// Panics if called before [`finalize`](Self::finalize).
    pub fn count_overlaps_per_query(
        &self,
        query: &RegionSet,
        min_overlap: i32,
    ) -> Vec<u32> {
        assert!(self.finalized, "Must finalize before querying");

        let mut counts = vec![0u32; query.regions.len()];

        for (q_idx, region) in query.regions.iter().enumerate() {
            let ctg_idx = match self.chrom_index.get(&region.chr) {
                Some(&idx) => idx,
                None => continue,
            };

            let start = region.start as i32;
            let end = region.end as i32;
            let contig = &self.contigs[ctg_idx];
            let n_tiles = contig.tiles.len() as i32;

            let n1 = start / self.nbp;
            let mut n2 = (end - 1) / self.nbp;
            if n1 >= n_tiles {
                continue;
            }
            n2 = n2.min(n_tiles - 1);

            let mut seen = std::collections::HashSet::new();

            // First tile
            let tile = &contig.tiles[n1 as usize];
            if !tile.records.is_empty() && end > tile.records[0].start {
                let mut tl: i32 = 0;
                let mut tr: i32 = tile.records.len() as i32 - 1;
                while tl < tr - 1 {
                    let tm = (tl + tr) / 2;
                    if tile.records[tm as usize].start < end {
                        tl = tm;
                    } else {
                        tr = tm;
                    }
                }
                if tile.records[tr as usize].start < end {
                    tl = tr;
                }
                for i in (0..=tl).rev() {
                    let rec = &tile.records[i as usize];
                    let overlap_bp = rec.end.min(end) - rec.start.max(start);
                    if overlap_bp >= min_overlap && seen.insert(rec.value) {
                        counts[q_idx] += 1;
                    }
                }
            }

            // Subsequent tiles
            if n2 > n1 {
                let mut bd = self.nbp * (n1 + 1);
                for j in (n1 + 1)..=n2 {
                    let tile = &contig.tiles[j as usize];
                    if tile.records.is_empty() {
                        bd += self.nbp;
                        continue;
                    }
                    if end > tile.records[0].start {
                        let mut ts: i32 = 0;
                        while ts < tile.records.len() as i32
                            && tile.records[ts as usize].start < bd
                        {
                            ts += 1;
                        }
                        let mut tl: i32 = 0;
                        let mut tr: i32 = tile.records.len() as i32 - 1;
                        while tl < tr - 1 {
                            let tm = (tl + tr) / 2;
                            if tile.records[tm as usize].start < end {
                                tl = tm;
                            } else {
                                tr = tm;
                            }
                        }
                        if tile.records[tr as usize].start < end {
                            tl = tr;
                        }
                        for i in (ts..=tl).rev() {
                            let rec = &tile.records[i as usize];
                            let overlap_bp = rec.end.min(end) - rec.start.max(start);
                            if overlap_bp >= min_overlap && seen.insert(rec.value) {
                                counts[q_idx] += 1;
                            }
                        }
                    }
                    bd += self.nbp;
                }
            }
        }

        counts
    }

    /// Number of source files indexed.
    pub fn num_files(&self) -> usize {
        self.file_info.len()
    }

    /// Number of contigs (chromosomes).
    pub fn num_contigs(&self) -> usize {
        self.contigs.len()
    }

    /// Total number of records across all tiles.
    /// Note: intervals spanning multiple tiles are counted once per tile.
    pub fn total_records(&self) -> usize {
        self.contigs
            .iter()
            .flat_map(|c| &c.tiles)
            .map(|t| t.records.len())
            .sum()
    }

    // -----------------------------------------------------------------------
    // Internal helpers
    // -----------------------------------------------------------------------

    /// Parse a BED line into (chrom, start, end, score).
    fn parse_bed_line(line: &str) -> Option<(String, i32, i32, i32)> {
        let mut fields = line.split('\t');
        let chrom = fields.next()?;
        let start: i32 = fields.next()?.parse().ok()?;
        let end: i32 = fields.next()?.parse().ok()?;

        if chrom.len() >= 40 || end <= 0 {
            return None;
        }

        let _ = fields.next(); // skip col4 (name)
        let score: i32 = fields
            .next()
            .and_then(|s| s.parse().ok())
            .unwrap_or(-1);

        Some((chrom.to_string(), start, end, score))
    }

    /// Load file info from a companion .tsv file.
    fn load_file_info_tsv(tsv_path: &Path) -> anyhow::Result<Vec<FileInfo>> {
        let file = File::open(tsv_path)?;
        let reader = BufReader::new(file);
        let mut infos = Vec::new();

        for (i, line) in reader.lines().enumerate() {
            if i == 0 {
                continue; // skip header
            }
            let line = line?;
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 4 {
                continue;
            }

            infos.push(FileInfo {
                filename: fields[1].trim().to_string(),
                num_regions: fields[2].trim().parse().unwrap_or(0),
                avg_region_width: fields[3].trim().parse().unwrap_or(0.0),
            });
        }

        Ok(infos)
    }
}

impl Default for Igd {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data")
    }

    #[test]
    fn test_igd_build_and_query_basic() {
        // Build IGD manually with known intervals
        let mut igd = Igd::new();
        igd.file_info = vec![
            FileInfo {
                filename: "file0.bed".into(),
                num_regions: 2,
                avg_region_width: 100.0,
            },
            FileInfo {
                filename: "file1.bed".into(),
                num_regions: 1,
                avg_region_width: 50.0,
            },
        ];

        // file0: two regions on chr1
        igd.add("chr1", 100, 200, 0, 0);
        igd.add("chr1", 300, 400, 0, 0);
        // file1: one region on chr1 overlapping file0's first region
        igd.add("chr1", 150, 250, 0, 1);

        igd.finalize();

        // Query region that overlaps file0[0] and file1[0]
        let mut hits = vec![0u64; 2];
        let n = igd.count_overlaps("chr1", 120, 180, 1, &mut hits);
        assert_eq!(n, 2);
        assert_eq!(hits[0], 1); // file0
        assert_eq!(hits[1], 1); // file1

        // Query region that overlaps only file0[1]
        let mut hits = vec![0u64; 2];
        let n = igd.count_overlaps("chr1", 350, 380, 1, &mut hits);
        assert_eq!(n, 1);
        assert_eq!(hits[0], 1);
        assert_eq!(hits[1], 0);

        // Query region that overlaps nothing
        let mut hits = vec![0u64; 2];
        let n = igd.count_overlaps("chr1", 500, 600, 1, &mut hits);
        assert_eq!(n, 0);
        assert_eq!(hits[0], 0);
        assert_eq!(hits[1], 0);
    }

    #[test]
    fn test_igd_min_overlap() {
        let mut igd = Igd::new();
        igd.file_info = vec![FileInfo {
            filename: "file0.bed".into(),
            num_regions: 1,
            avg_region_width: 100.0,
        }];

        // file0: region [100, 200)
        igd.add("chr1", 100, 200, 0, 0);
        igd.finalize();

        // Query [190, 250) — only 10bp overlap with [100,200)
        let mut hits = vec![0u64; 1];
        igd.count_overlaps("chr1", 190, 250, 1, &mut hits);
        assert_eq!(hits[0], 1); // 10bp >= 1bp

        let mut hits = vec![0u64; 1];
        igd.count_overlaps("chr1", 190, 250, 10, &mut hits);
        assert_eq!(hits[0], 1); // 10bp >= 10bp

        let mut hits = vec![0u64; 1];
        igd.count_overlaps("chr1", 190, 250, 11, &mut hits);
        assert_eq!(hits[0], 0); // 10bp < 11bp
    }

    #[test]
    fn test_igd_multi_tile_spanning() {
        // Test an interval that spans multiple tiles
        let mut igd = Igd::new();
        igd.file_info = vec![FileInfo {
            filename: "file0.bed".into(),
            num_regions: 1,
            avg_region_width: 20000.0,
        }];

        // Interval spanning tiles 0 and 1 (tile size = 16384)
        igd.add("chr1", 10000, 20000, 0, 0);
        igd.finalize();

        // Query in tile 0
        let mut hits = vec![0u64; 1];
        igd.count_overlaps("chr1", 11000, 12000, 1, &mut hits);
        assert_eq!(hits[0], 1);

        // Query in tile 1
        let mut hits = vec![0u64; 1];
        igd.count_overlaps("chr1", 17000, 18000, 1, &mut hits);
        assert_eq!(hits[0], 1);

        // Query spanning both tiles
        let mut hits = vec![0u64; 1];
        igd.count_overlaps("chr1", 15000, 19000, 1, &mut hits);
        assert_eq!(hits[0], 1); // should count only once per file
    }

    #[test]
    fn test_igd_unknown_chrom() {
        let mut igd = Igd::new();
        igd.file_info = vec![FileInfo {
            filename: "file0.bed".into(),
            num_regions: 1,
            avg_region_width: 100.0,
        }];
        igd.add("chr1", 100, 200, 0, 0);
        igd.finalize();

        let mut hits = vec![0u64; 1];
        let n = igd.count_overlaps("chrZ", 100, 200, 1, &mut hits);
        assert_eq!(n, 0);
    }

    #[test]
    fn test_igd_from_bed_dir() {
        let bed_dir = test_data_dir().join("igd_file_list_01");
        let igd = Igd::from_bed_dir(&bed_dir).unwrap();

        assert_eq!(igd.num_files(), 1);
        assert_eq!(igd.num_contigs(), 3); // chr1, chr2, chr3

        // Query all 8 regions from the file against itself — each should hit
        let mut hits = vec![0u64; 1];
        // chr1: [1,100), [200,300), [32768,32868), [49152,49352)
        igd.count_overlaps("chr1", 1, 100, 1, &mut hits);
        igd.count_overlaps("chr1", 200, 300, 1, &mut hits);
        igd.count_overlaps("chr1", 32768, 32868, 1, &mut hits);
        igd.count_overlaps("chr1", 49152, 49352, 1, &mut hits);
        // chr2: [1,100), [200,300)
        igd.count_overlaps("chr2", 1, 100, 1, &mut hits);
        igd.count_overlaps("chr2", 200, 300, 1, &mut hits);
        // chr3: [32768,32868), [49152,49352)
        igd.count_overlaps("chr3", 32768, 32868, 1, &mut hits);
        igd.count_overlaps("chr3", 49152, 49352, 1, &mut hits);

        assert_eq!(hits[0], 8); // all 8 regions overlap
    }

    #[test]
    fn test_igd_from_region_sets() {
        let sets = vec![
            (
                "set1.bed".to_string(),
                vec![
                    ("chr1".to_string(), 100, 200),
                    ("chr1".to_string(), 300, 400),
                ],
            ),
            (
                "set2.bed".to_string(),
                vec![("chr1".to_string(), 150, 350)],
            ),
        ];

        let igd = Igd::from_region_sets(sets);
        assert_eq!(igd.num_files(), 2);

        let mut hits = vec![0u64; 2];
        igd.count_overlaps("chr1", 160, 170, 1, &mut hits);
        assert_eq!(hits[0], 1); // set1: [100,200) overlaps
        assert_eq!(hits[1], 1); // set2: [150,350) overlaps
    }

    #[test]
    fn test_igd_save_and_reload() {
        let mut igd = Igd::new();
        igd.file_info = vec![
            FileInfo {
                filename: "file0.bed".into(),
                num_regions: 2,
                avg_region_width: 100.0,
            },
            FileInfo {
                filename: "file1.bed".into(),
                num_regions: 1,
                avg_region_width: 100.0,
            },
        ];

        igd.add("chr1", 100, 200, 5, 0);
        igd.add("chr1", 300, 400, 10, 0);
        igd.add("chr2", 50, 150, 0, 1);
        igd.finalize();

        // Save
        let tmpdir = tempfile::tempdir().unwrap();
        let igd_path = tmpdir.path().join("test.igd");
        igd.save(&igd_path).unwrap();

        // Reload
        let igd2 = Igd::from_igd_file(&igd_path).unwrap();
        assert_eq!(igd2.num_contigs(), 2);
        assert_eq!(igd2.num_files(), 2);

        // Query reloaded IGD should match original
        let mut hits_orig = vec![0u64; 2];
        let mut hits_reload = vec![0u64; 2];

        igd.count_overlaps("chr1", 150, 350, 1, &mut hits_orig);
        igd2.count_overlaps("chr1", 150, 350, 1, &mut hits_reload);

        assert_eq!(hits_orig, hits_reload);
    }

    #[test]
    fn test_igd_disk_roundtrip_matches_legacy() {
        // Build IGD from test data using new API, save, reload, and verify
        // query results match what the legacy create+search pipeline produces
        let bed_dir = test_data_dir().join("igd_file_list_01");
        let igd = Igd::from_bed_dir(&bed_dir).unwrap();

        let tmpdir = tempfile::tempdir().unwrap();
        let igd_path = tmpdir.path().join("demo.igd");
        igd.save(&igd_path).unwrap();

        let igd2 = Igd::from_igd_file(&igd_path).unwrap();

        // Query from the query file and compare
        let mut hits1 = vec![0u64; igd.num_files()];
        let mut hits2 = vec![0u64; igd2.num_files()];

        let query_regions = vec![
            ("chr1", 1, 100),
            ("chr1", 200, 300),
            ("chr1", 32768, 32868),
            ("chr1", 49152, 49352),
            ("chr2", 1, 100),
            ("chr2", 200, 300),
            ("chr3", 32768, 32868),
            ("chr3", 49152, 49352),
        ];

        for (c, s, e) in &query_regions {
            igd.count_overlaps(c, *s, *e, 1, &mut hits1);
            igd2.count_overlaps(c, *s, *e, 1, &mut hits2);
        }

        assert_eq!(hits1, hits2);
    }

    #[test]
    fn test_igd_count_set_overlaps() {
        let sets = vec![
            (
                "db1.bed".to_string(),
                vec![
                    ("chr1".to_string(), 100, 200),
                    ("chr1".to_string(), 500, 600),
                ],
            ),
            (
                "db2.bed".to_string(),
                vec![("chr1".to_string(), 150, 250)],
            ),
        ];

        let igd = Igd::from_region_sets(sets);

        // Create a query RegionSet
        let query = RegionSet::from(vec![
            gtars_core::models::Region {
                chr: "chr1".to_string(),
                start: 120,
                end: 180,
                rest: None,
            },
            gtars_core::models::Region {
                chr: "chr1".to_string(),
                start: 520,
                end: 560,
                rest: None,
            },
        ]);

        let hits = igd.count_set_overlaps(&query, 1);
        assert_eq!(hits[0], 2); // db1: both [100,200) and [500,600) hit
        assert_eq!(hits[1], 1); // db2: only [150,250) hits first query region
    }

    #[test]
    fn test_igd_pairwise_overlap_counting() {
        // Verify that IGD counts pairwise overlaps (not binary per query region).
        // If one query region overlaps 3 DB regions in the same file, that's 3 hits.
        let mut igd = Igd::new();
        igd.file_info = vec![FileInfo {
            filename: "file0.bed".into(),
            num_regions: 3,
            avg_region_width: 50.0,
        }];

        // Three overlapping DB regions in same file
        igd.add("chr1", 100, 200, 0, 0);
        igd.add("chr1", 120, 220, 0, 0);
        igd.add("chr1", 140, 240, 0, 0);
        igd.finalize();

        // One query region overlapping all three
        let mut hits = vec![0u64; 1];
        igd.count_overlaps("chr1", 150, 190, 1, &mut hits);
        assert_eq!(hits[0], 3); // pairwise: 3 overlaps, not 1
    }

    #[test]
    fn test_igd_empty() {
        let igd = Igd::new();
        // Can't query unfinalized IGD — but an empty one should work after finalize
        let mut igd = igd;
        igd.finalize();

        assert_eq!(igd.num_files(), 0);
        assert_eq!(igd.num_contigs(), 0);
    }

    fn make_region(chr: &str, start: u32, end: u32) -> gtars_core::models::Region {
        gtars_core::models::Region {
            chr: chr.to_string(),
            start,
            end,
            rest: None,
        }
    }

    #[test]
    fn test_from_single_region_set() {
        let subject = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
            make_region("chr2", 50, 150),
        ]);

        let igd = Igd::from_single_region_set(&subject);
        assert_eq!(igd.num_files(), 1);
        assert_eq!(igd.file_info[0].num_regions, 3);
    }

    #[test]
    fn test_find_overlaps_regionset_basic() {
        let subject = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 300, 400),
            make_region("chr1", 500, 600),
        ]);

        let query = RegionSet::from(vec![
            make_region("chr1", 150, 350), // overlaps subject 0 and 1
            make_region("chr1", 550, 650), // overlaps subject 2
            make_region("chr1", 700, 800), // no overlap
        ]);

        let igd = Igd::from_single_region_set(&subject);
        let mut pairs = igd.find_overlaps_regionset(&query, 1);
        pairs.sort();

        assert_eq!(pairs, vec![(0, 0), (0, 1), (1, 2)]);
    }

    #[test]
    fn test_find_overlaps_regionset_no_overlap() {
        let subject = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let query = RegionSet::from(vec![make_region("chr1", 300, 400)]);

        let igd = Igd::from_single_region_set(&subject);
        let pairs = igd.find_overlaps_regionset(&query, 1);
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_find_overlaps_regionset_min_overlap() {
        let subject = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let query = RegionSet::from(vec![make_region("chr1", 190, 300)]); // 10bp overlap

        let igd = Igd::from_single_region_set(&subject);

        let pairs1 = igd.find_overlaps_regionset(&query, 1);
        assert_eq!(pairs1.len(), 1);

        let pairs50 = igd.find_overlaps_regionset(&query, 50);
        assert!(pairs50.is_empty());
    }

    #[test]
    fn test_find_overlaps_regionset_multi_chrom() {
        let subject = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr2", 100, 200),
        ]);

        let query = RegionSet::from(vec![
            make_region("chr1", 150, 180),
            make_region("chr2", 150, 180),
            make_region("chr3", 150, 180), // no subject on chr3
        ]);

        let igd = Igd::from_single_region_set(&subject);
        let mut pairs = igd.find_overlaps_regionset(&query, 1);
        pairs.sort();

        assert_eq!(pairs, vec![(0, 0), (1, 1)]);
    }

    #[test]
    fn test_count_overlaps_per_query_basic() {
        let subject = RegionSet::from(vec![
            make_region("chr1", 100, 200),
            make_region("chr1", 150, 250),
            make_region("chr1", 500, 600),
        ]);

        let query = RegionSet::from(vec![
            make_region("chr1", 160, 180), // overlaps subject 0 and 1
            make_region("chr1", 550, 580), // overlaps subject 2
            make_region("chr1", 700, 800), // no overlap
        ]);

        let igd = Igd::from_single_region_set(&subject);
        let counts = igd.count_overlaps_per_query(&query, 1);

        assert_eq!(counts, vec![2, 1, 0]);
    }

    #[test]
    fn test_count_overlaps_per_query_empty() {
        let subject = RegionSet::from(vec![make_region("chr1", 100, 200)]);
        let query = RegionSet::from(vec![]);

        let igd = Igd::from_single_region_set(&subject);
        let counts = igd.count_overlaps_per_query(&query, 1);
        assert!(counts.is_empty());
    }

    #[test]
    fn test_find_overlaps_multi_tile_dedup() {
        // Subject region spans multiple tiles — should only appear once per query
        let subject = RegionSet::from(vec![
            make_region("chr1", 10000, 40000), // spans tiles 0, 1, 2
        ]);

        let query = RegionSet::from(vec![
            make_region("chr1", 15000, 35000), // also spans tiles
        ]);

        let igd = Igd::from_single_region_set(&subject);
        let pairs = igd.find_overlaps_regionset(&query, 1);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0], (0, 0));

        let counts = igd.count_overlaps_per_query(&query, 1);
        assert_eq!(counts, vec![1]);
    }

    #[test]
    fn test_from_region_sets_skips_invalid_intervals() {
        // Mix of valid and invalid intervals — invalid should be skipped in count/width
        let sets = vec![(
            "test.bed".to_string(),
            vec![
                ("chr1".to_string(), 100, 200), // valid: 100bp
                ("chr1".to_string(), 300, 300), // invalid: start == end
                ("chr1".to_string(), 500, 400), // invalid: start > end
                ("chr1".to_string(), 600, 700), // valid: 100bp
            ],
        )];
        let igd = Igd::from_region_sets(sets);

        assert_eq!(igd.file_info.len(), 1);
        assert_eq!(igd.file_info[0].num_regions, 2, "Only 2 valid intervals");
        assert!(
            (igd.file_info[0].avg_region_width - 100.0).abs() < 1e-6,
            "avg width should be 100.0, got {}",
            igd.file_info[0].avg_region_width
        );
    }

    #[test]
    fn test_add_negative_coordinates_no_panic() {
        let mut igd = Igd::new();
        igd.file_info = vec![FileInfo {
            filename: "test.bed".into(),
            num_regions: 0,
            avg_region_width: 0.0,
        }];

        // These should all be silently skipped, no panic
        igd.add("chr1", -100, 200, 0, 0);
        igd.add("chr1", 100, -200, 0, 0);
        igd.add("chr1", -100, -50, 0, 0);

        // Add one valid interval to verify the IGD still works
        igd.add("chr1", 100, 200, 0, 0);
        igd.finalize();

        // Query should find the one valid interval
        let mut hits = vec![0u64; 1];
        let count = igd.count_overlaps("chr1", 150, 160, 1, &mut hits);
        assert_eq!(count, 1);
    }

    #[test]
    fn test_query_negative_coordinates() {
        let sets = vec![(
            "test.bed".to_string(),
            vec![("chr1".to_string(), 100, 200)],
        )];
        let igd = Igd::from_region_sets(sets);

        let mut hits = vec![0u64; 1];
        // Negative start, positive end — should clamp start to 0 and still work
        let count = igd.count_overlaps("chr1", -50, 150, 1, &mut hits);
        assert_eq!(count, 1, "Should find overlap after clamping negative start");

        // Both negative — should return 0
        hits[0] = 0;
        let count = igd.count_overlaps("chr1", -100, -50, 1, &mut hits);
        assert_eq!(count, 0, "Both negative should return 0");
    }

    #[test]
    fn test_parse_bed_line_no_chr_prefix() {
        // Non-UCSC chromosome names should now be accepted
        let result = Igd::parse_bed_line("1\t100\t200\tname\t500");
        assert!(result.is_some(), "Should parse non-chr-prefixed chromosomes");
        let (chrom, start, end, score) = result.unwrap();
        assert_eq!(chrom, "1");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
        assert_eq!(score, 500);
    }

    #[test]
    fn test_from_bed_dir_large_coordinates() {
        // Intervals with coordinates > 321M should now be accepted
        let sets = vec![(
            "test.bed".to_string(),
            vec![("chr1".to_string(), 400_000_000, 400_001_000)],
        )];
        let igd = Igd::from_region_sets(sets);

        let mut hits = vec![0u64; 1];
        let count = igd.count_overlaps("chr1", 400_000_500, 400_000_600, 1, &mut hits);
        assert_eq!(count, 1, "Should find overlap at coordinates > 321M");
    }
}
