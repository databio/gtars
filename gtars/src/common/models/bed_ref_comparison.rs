use crate::common::models::{RegionSet};
use std::path::{Path, PathBuf};
use std::collections::{HashMap, HashSet};

use crate::common::utils::get_dynamic_reader;
use std::io::BufRead;
use std::fs;


#[derive(Clone, Debug)]
pub struct ReferenceGenomeMetadata {
    pub genome: String,
    pub digest: String,
    pub description: String,
    pub collection: HashMap<String, u32>
}

impl ReferenceGenomeMetadata {
    pub fn new(genome: String, digest: String, description: String, collection: HashMap<String, u32>) -> Self {
        ReferenceGenomeMetadata { genome, digest, description, collection }
    }
    pub fn try_from(value: &Path, name: &str) -> anyhow::Result<Self> {
        // chrom sizes file
        let path = value;

        let file_reader = get_dynamic_reader(path).expect("!Can't read file");

        let mut chrom_sizes: HashMap<String, u32> = HashMap::new();

        for line in file_reader.lines() {
            let string_line = line?;

            let parts: Vec<String> = string_line.split('\t').map(|s| s.to_string()).collect();

            chrom_sizes.insert(parts[0].to_string(), parts[1].parse::<u32>()?);

        }
        let ref_genome = ReferenceGenomeMetadata {
            genome: name.to_string(),
            digest: name.to_string(),
            description: name.to_string(),
            collection: chrom_sizes
        };

        Ok(ref_genome)

    }
}



#[derive(Debug, Default)]
pub struct ChromNameStats {
    pub xs: f64, // Extra Sequences (equivalent to sensitivity)
    pub q_and_m: f64,
    pub q_and_not_m: f64,
    pub not_q_and_m: f64,
    pub jaccard_index: f64,
    pub jaccard_index_binary: f64,
    pub passed_chrom_names: bool,
}

#[derive(Debug, Default)]
pub struct ChromLengthStats {
    pub oobr: f64,
    pub beyond_range: bool,
    pub num_of_chrom_beyond: i32,
    pub percentage_bed_chrom_beyond: f64,
    pub percentage_genome_chrom_beyond: f64,
}

#[derive(Debug, Default)]
pub struct SequenceFitStats{
    sequence_fit: Option<f64>
}

#[derive(Debug, Default)]
pub struct RatingModel{
    assigned_points: i32,
    tier_ranking: i32,
}

#[derive(Debug, Default)]
pub struct CompatibilityStats {
    chrom_name_stats: ChromNameStats,
    chrom_length_stats: Option<ChromLengthStats>,
    chrom_sequence_fit_stats: SequenceFitStats,
    compatibility: RatingModel,
}

#[derive(Debug, Default)]
pub struct CompatibilityConcise {
    pub xs: f64,
    pub oobr: Option<f64>,
    pub sequence_fit: Option<f64>,
    pub assigned_points: i32,
    pub tier_ranking: i32,
}

#[derive(Debug, Default, Clone)]
pub struct ReferenceValidator {
    pub reference_genomes: Vec<ReferenceGenomeMetadata>,
}


fn _list_chrom_files(folder_path: &Path) -> Vec<PathBuf> {
        let mut files: Vec<PathBuf> = Vec::new();

        if let Ok(entries) = fs::read_dir(folder_path) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_file() {
                    files.push(path);
                }
            }
        }
        files
}

impl ReferenceValidator {
    pub fn new(reference_genomes: Vec<ReferenceGenomeMetadata>) -> Self {
        Self { reference_genomes }
    }

    pub fn determine_compatibility(&self, rs: RegionSet) -> HashMap<String, CompatibilityConcise> {
        if self.reference_genomes.is_empty(){
            return HashMap::new();
        };

        let mut compatibility_map: HashMap<String, CompatibilityConcise> = HashMap::new();
        let bed_info: HashMap<String, u32> = rs.get_max_end_per_chr();
        println!("Determining compatibility");

        for genome in self.reference_genomes.iter(){
            compatibility_map.insert(genome.digest.clone(), get_concise_stats(genome.clone(), &bed_info));
            // println!("Determining chrom length stats for {}", genome.digest);
        };

        compatibility_map
    }

    pub fn try_from(folder_path: &Path) -> Self {
        let chrom_file_list = _list_chrom_files(folder_path);

        let mut genomes_vector: Vec<ReferenceGenomeMetadata> = Vec::new();

        for file in chrom_file_list{

            genomes_vector.push(ReferenceGenomeMetadata::try_from(&file, &file.file_name().unwrap().to_str().unwrap()).unwrap());
        };

        ReferenceValidator::new(genomes_vector)
    }
}



fn calculate_rating(xs: f64, oobr: Option<f64>, sequence_fit: Option<f64>) -> RatingModel {
    let mut points_rating: i32 = 0;
    let mut tier_ranking: i32 = 0;

    // 1. Check extra sequences sensitivity.
    // sensitivity = 1 is considered great and no points should be assigned
    if xs < 0.3 {
        points_rating += 6  // 3 + 1 + 1 + 1;
    } else if xs < 0.5 {
        points_rating += 5  // 3 + 1 + 1
    } else if xs < 0.7 {
        points_rating += 4 // 3 + 1
    } else if xs < 1.0 {
        points_rating += 3
    };
    //  2. Check OOBR and assign points based on sensitivity
    //  only assessed if no extra chroms in query bed file
    if oobr.is_some() {
        let oobr_value = oobr.unwrap();

        if oobr_value < 0.3 {
            points_rating += 6  // 3 + 1 + 1 + 1
        } else if oobr_value < 0.5 {
            points_rating += 5  // 3 + 1 + 1
        } else if oobr_value < 0.7 {
            points_rating += 4 // 3 + 1
        } else if oobr_value < 1.0 {
            points_rating += 3
        }
        // Do nothing here, points have already been added when Assessing XS if it is not == 1

    };

    if sequence_fit.is_some() {
        let sequence_fit_value = sequence_fit.unwrap();

        // 3. Check Sequence Fit - comparing lengths in queries vs lengths of queries in ref genome vs not in ref genome
        if sequence_fit_value < 0.90 {
            points_rating += 1
        };
        if sequence_fit_value < 0.60 {
            points_rating += 1
        };
        if sequence_fit_value < 0.60 {
            points_rating += 1 // TODO: is it intentional
        };
    }

    if points_rating == 0 {
        tier_ranking = 1;
    } else if points_rating >= 1 && points_rating <= 3 {
        tier_ranking = 2;
    } else if points_rating >= 4 && points_rating <= 6 {
        tier_ranking = 3;
    } else if 7 <= points_rating {
        tier_ranking = 4;
    };

    RatingModel{
        assigned_points: points_rating,
        tier_ranking,
    }

}


// pub fn caclulate_chrom_stats(genome_info: ReferenceGenomeMetadata, rs: RegionSet) -> CompatibilityStats {
pub fn caclulate_chrom_stats(genome_info: ReferenceGenomeMetadata, bed_info: &HashMap<String, u32>) -> CompatibilityStats {

    // let bed_info: HashMap<String, u32> = rs.get_max_end_per_chr();

    let mut passed_chrom_names: bool = true;

    let genome_chrom_set: HashSet<_> = genome_info.collection.keys().collect();
    let bed_chrom_set: HashSet<_> = bed_info.keys().collect();

    let chrom_intersection: HashSet<_> = genome_chrom_set.intersection(&bed_chrom_set).clone().collect();
    let chrom_union: HashSet<_> = genome_chrom_set.union(&bed_chrom_set).clone().collect();

    let q_and_m = chrom_intersection.len() as f64;
    let q_and_not_m = bed_info.len() as f64 - chrom_intersection.len() as f64;
    // Symmetric difference (keys that are in one but not both)
    let not_q_and_m = genome_chrom_set.symmetric_difference(&bed_chrom_set).count() as f64;

    let chrom_jaccard_index = (chrom_intersection.len() as f64)/chrom_union.len() as f64;
    let jaccard_binary: f64 = q_and_not_m / (q_and_m + not_q_and_m + q_and_not_m);
    let sensitivity: f64 = q_and_m / (q_and_m + q_and_not_m);

    if q_and_not_m > 0.0 {
        passed_chrom_names = false;
    }

    let chrom_stats = ChromNameStats{
        xs: sensitivity,
        q_and_m: q_and_m,
        q_and_not_m: q_and_m,
        not_q_and_m: not_q_and_m,
        jaccard_index: chrom_jaccard_index,
        jaccard_index_binary: jaccard_binary,
        passed_chrom_names: passed_chrom_names,
    };


    // Layer 2:  Check Lengths, but only if layer 1 is passing [all chroms are in ref genome]
    let mut chrom_length_stats: Option<ChromLengthStats> = None;
    if passed_chrom_names {
        let mut chroms_beyond_range: bool = false;
        let mut num_of_chrom_beyond: i32 = 0;
        let mut num_chrom_within_bounds: i32 = 0;

        for key in  bed_chrom_set.iter() {
            let key_string = key.to_owned();
            if genome_chrom_set.contains(key_string) {

                if bed_info.get(key_string).unwrap() > genome_info.collection.get(key_string).unwrap() {

                    num_of_chrom_beyond += 1;
                    chroms_beyond_range = true;
                } else {
                    num_chrom_within_bounds += 1;
                }
            }
        }
        let oobr: f64 = ( num_chrom_within_bounds / ( num_chrom_within_bounds + num_of_chrom_beyond)) as f64;

        chrom_length_stats = Some(ChromLengthStats{
                oobr,
                beyond_range: chroms_beyond_range,
                num_of_chrom_beyond: num_of_chrom_beyond,
                percentage_bed_chrom_beyond: (
                    100.0 * num_of_chrom_beyond as f64 / bed_chrom_set.len() as f64
                ),
                percentage_genome_chrom_beyond: (
                    100.0 * num_of_chrom_beyond as f64 / genome_chrom_set.len() as f64),
        });

    };

    // Layer 3 Calculate Sequence Fit if any query chrom names were present
    let mut sequence_fit: Option<f64> = None;

    if chrom_intersection.len() > 0{
        let mut bed_sum: f64 = 0.0;
        let mut ref_genome_sum:f64 = 0.0;

        for key in chrom_intersection {
            bed_sum += genome_info.collection.get(key.to_owned()).unwrap().clone() as f64;
        }
        for key in genome_info.collection.keys() {
            ref_genome_sum += genome_info.collection.get(key).unwrap().clone() as f64;
        }

         sequence_fit = Some(bed_sum / ref_genome_sum);
    }

    let oobr_option: Option<f64> = match &chrom_length_stats {
        Some(chrom_length_stats) => {
            Some(chrom_length_stats.oobr)
        }
        _ => None,
    };

    let tier_obj = calculate_rating(
        chrom_stats.xs,
        oobr_option,
        sequence_fit,
    );


    CompatibilityStats {
        chrom_name_stats: chrom_stats,
        chrom_length_stats: chrom_length_stats,
        chrom_sequence_fit_stats: SequenceFitStats{sequence_fit: sequence_fit},
        compatibility: tier_obj,
}

}

// pub fn get_concise_stats(genome_info: ReferenceGenomeMetadata, rs: RegionSet) -> CompatibilityConcise{
pub fn get_concise_stats(genome_info: ReferenceGenomeMetadata, bed_info: &HashMap<String, u32>) -> CompatibilityConcise{
    let chrom_stats:CompatibilityStats = caclulate_chrom_stats(genome_info, bed_info);

    let oobr_option: Option<f64> = match &chrom_stats.chrom_length_stats {
        Some(chrom_length_stats) => {
            Some(chrom_length_stats.oobr)
        }
        _ => None,
    };

    CompatibilityConcise{
        xs: chrom_stats.chrom_name_stats.xs,
        oobr: oobr_option,
        sequence_fit: chrom_stats.chrom_sequence_fit_stats.sequence_fit,
        assigned_points: chrom_stats.compatibility.tier_ranking,
        tier_ranking: chrom_stats.compatibility.tier_ranking,
    }
}




#[cfg(test)]
mod tests {
    use std::io::Error;
    use std::path::PathBuf;
    use super::*;
    use rstest::*;

    fn get_test_path(file_name: &str) -> anyhow::Result<PathBuf, Error> {
        let file_path: PathBuf = std::env::current_dir()
            .unwrap()
            .join("tests/data/regionset")
            .join(file_name);
        Ok(file_path)
    }

    #[rstest]
    fn open_bed_file(){

        let file_path = get_test_path("dummy.narrowPeak.bed.gz").unwrap();
        let bed_file: RegionSet = RegionSet::try_from(file_path).unwrap();

        println!("{}", bed_file);
    }

    #[rstest]
    fn open_chrom_sizes_file(){
        let file_path = Path::new("/home/bnt4me/virginia/repos/bedboss/bedboss/refgenome_validator/chrom_sizes/ensembl_hg38.chrom.sizes");
        let name: String = String::from("ensembl_hg38.chrom.sizes");

        let ref_obj: ReferenceGenomeMetadata = ReferenceGenomeMetadata::try_from(file_path, &name).unwrap();

        println!("{:?}", ref_obj);

    }

    #[rstest]
    fn run_calculations(){
        let bed_path2 = Path::new("/home/bnt4me/virginia/repos/bedboss/test/data/bed/hg38/GSM6732293_Con_liver-IP2.bed");
        let bed_path = Path::new("/home/bnt4me/Downloads/003c91fed233b4def93aa1fcb743a317.bed.gz");
        let ref_path = Path::new("/home/bnt4me/virginia/repos/bedboss/bedboss/refgenome_validator/chrom_sizes/ucsc_hg38.chrom.sizes");
        let name: String = String::from("ensembl_hg38.chrom.sizes");

        let ref_obj: ReferenceGenomeMetadata = ReferenceGenomeMetadata::try_from(ref_path, &name).unwrap();
        let bed_file: RegionSet = RegionSet::try_from(bed_path).unwrap();
        let bed_file2: RegionSet = RegionSet::try_from(bed_path2).unwrap();

        // let kj = caclulate_chrom_stats(ref_obj, bed_file);
        // let kj = get_concise_stats(ref_obj, bed_file);


        let folder_path = Path::new("/home/bnt4me/virginia/repos/bedboss/bedboss/refgenome_validator/chrom_sizes");

        let refVal = ReferenceValidator::try_from(folder_path);

        let kj = refVal.determine_compatibility( bed_file);

        println!("{:?}", kj);

        let kj2 = refVal.determine_compatibility( bed_file2);
        println!("{:?}", kj2);

    }
}