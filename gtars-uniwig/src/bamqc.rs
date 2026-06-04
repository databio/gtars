use noodles::bam;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Default)]
pub struct BamQcResult {
    pub total_reads: u64,
    pub distinct: u64,
    pub m1: u64,
    pub m2: u64,
    pub dups: u64,
    pub mito_reads: u64,
    pub nrf: f64,
    pub pbc1: f64,
    pub pbc2: f64,
}

impl BamQcResult {
    pub fn mito_rate(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            self.mito_reads as f64 / self.total_reads as f64
        }
    }

    pub fn dup_rate(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            self.dups as f64 / self.total_reads as f64
        }
    }
}

#[derive(Debug, Clone, Default)]
struct ChromResult {
    chrom: String,
    position_counts: HashMap<(i64, i64), u64>,
    total_reads: u64,
    dup_count: u64,
    mito_reads: u64,
    is_paired: bool,
}

fn is_mitochondrial(chrom: &str) -> bool {
    let lower = chrom.to_lowercase();
    lower == "chrm" || lower == "mt" || lower == "chrmt" || lower.contains("rcrsd")
}

fn process_chromosome(bam_path: &Path, chrom: &str) -> Result<ChromResult, Box<dyn Error + Send + Sync>> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let region = chrom.parse()?;
    let records = reader.query(&header, &region)?;

    let mut result = ChromResult::default();
    result.chrom = chrom.to_string();

    if is_mitochondrial(chrom) {
        for record_result in records {
            let record = record_result?;
            if !record.flags().is_unmapped() {
                result.total_reads += 1;
                result.mito_reads += 1;
                if record.flags().is_duplicate() {
                    result.dup_count += 1;
                }
            }
        }
        return Ok(result);
    }

    for record_result in records {
        let record = record_result?;

        if record.flags().is_unmapped() {
            continue;
        }

        result.total_reads += 1;

        if record.flags().is_duplicate() {
            result.dup_count += 1;
        }

        let pos = match record.alignment_start() {
            Some(p) => match p {
                Ok(position) => position.get() as i64,
                Err(_) => continue,
            },
            None => continue,
        };

        let key = if record.flags().is_segmented() {
            result.is_paired = true;
            if record.flags().is_first_segment() {
                let tlen = record.template_length() as i64;
                (pos, tlen)
            } else {
                continue;
            }
        } else {
            let qlen = record.sequence().len() as i64;
            (pos, qlen)
        };

        *result.position_counts.entry(key).or_insert(0) += 1;
    }

    Ok(result)
}

pub fn compute_bam_qc_parallel(bam_path: &Path, num_threads: usize) -> Result<BamQcResult, Box<dyn Error>> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()?;

    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let chromosomes: Vec<String> = header
        .reference_sequences()
        .iter()
        .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).to_string())
        .collect();

    let results: Vec<ChromResult> = pool.install(|| {
        chromosomes
            .par_iter()
            .filter_map(|chrom| {
                match process_chromosome(bam_path, chrom) {
                    Ok(result) => Some(result),
                    Err(_) => None,
                }
            })
            .collect()
    });

    let mut total_reads: u64 = 0;
    let mut dup_count: u64 = 0;
    let mut mito_count: u64 = 0;
    let mut is_paired_data = false;
    let mut global_counts: HashMap<(String, i64, i64), u64> = HashMap::new();

    for chrom_result in results {
        total_reads += chrom_result.total_reads;
        dup_count += chrom_result.dup_count;
        mito_count += chrom_result.mito_reads;
        if chrom_result.is_paired {
            is_paired_data = true;
        }
        for ((pos, tlen), count) in chrom_result.position_counts {
            let key = (chrom_result.chrom.clone(), pos, tlen);
            *global_counts.entry(key).or_insert(0) += count;
        }
    }

    let m_distinct = global_counts.len() as u64;
    let mut m1: u64 = 0;
    let mut m2: u64 = 0;

    for &count in global_counts.values() {
        if count == 1 {
            m1 += 1;
        } else if count == 2 {
            m2 += 1;
        }
    }

    let effective_total = if is_paired_data {
        total_reads / 2
    } else {
        total_reads
    };

    let total_f = effective_total.max(1) as f64;
    let m_distinct_f = m_distinct.max(1) as f64;
    let m2_f = m2.max(1) as f64;

    let nrf = m1 as f64 / total_f;
    let pbc1 = m1 as f64 / m_distinct_f;
    let pbc2 = m1 as f64 / m2_f;

    Ok(BamQcResult {
        total_reads: effective_total,
        distinct: m_distinct,
        m1,
        m2: m2.max(1),
        dups: dup_count,
        mito_reads: mito_count,
        nrf,
        pbc1,
        pbc2,
    })
}

pub fn compute_bam_qc(bam_path: &Path) -> Result<BamQcResult, Box<dyn Error>> {
    let mut reader = bam::io::reader::Builder.build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let mut position_counts: HashMap<(String, i64, i64), u64> = HashMap::new();
    let mut total_reads: u64 = 0;
    let mut dup_count: u64 = 0;
    let mut mito_count: u64 = 0;
    let mut is_paired_data = false;

    let reference_sequences = header.reference_sequences();

    for result in reader.records() {
        let record = result?;

        if record.flags().is_unmapped() {
            continue;
        }

        total_reads += 1;

        if record.flags().is_duplicate() {
            dup_count += 1;
        }

        let chrom_idx = match record.reference_sequence_id() {
            Some(id) => match id {
                Ok(idx) => idx,
                Err(_) => continue,
            },
            None => continue,
        };

        let chrom_name = reference_sequences
            .get_index(chrom_idx)
            .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).to_string())
            .unwrap_or_default();

        if is_mitochondrial(&chrom_name) {
            mito_count += 1;
            continue;
        }

        let pos = match record.alignment_start() {
            Some(p) => match p {
                Ok(position) => position.get() as i64,
                Err(_) => continue,
            },
            None => continue,
        };

        let key = if record.flags().is_segmented() {
            is_paired_data = true;
            if record.flags().is_first_segment() {
                let tlen = record.template_length() as i64;
                (chrom_name, pos, tlen)
            } else {
                continue;
            }
        } else {
            let qlen = record.sequence().len() as i64;
            (chrom_name, pos, qlen)
        };

        *position_counts.entry(key).or_insert(0) += 1;
    }

    let m_distinct = position_counts.len() as u64;

    let mut m1: u64 = 0;
    let mut m2: u64 = 0;

    for &count in position_counts.values() {
        if count == 1 {
            m1 += 1;
        } else if count == 2 {
            m2 += 1;
        }
    }

    let effective_total = if is_paired_data {
        total_reads / 2
    } else {
        total_reads
    };

    let total_f = effective_total.max(1) as f64;
    let m_distinct_f = m_distinct.max(1) as f64;
    let m2_f = m2.max(1) as f64;

    let nrf = m1 as f64 / total_f;
    let pbc1 = m1 as f64 / m_distinct_f;
    let pbc2 = m1 as f64 / m2_f;

    Ok(BamQcResult {
        total_reads: effective_total,
        distinct: m_distinct,
        m1,
        m2: m2.max(1),
        dups: dup_count,
        mito_reads: mito_count,
        nrf,
        pbc1,
        pbc2,
    })
}

pub fn write_bam_qc_tsv<W: Write>(result: &BamQcResult, writer: &mut W) -> Result<(), Box<dyn Error>> {
    writeln!(
        writer,
        "Total_read_pairs\tDistinct_read_pairs\tOne_read_pair\tTwo_read_pairs\tDuplicate_rate\tMitochondria_reads\tMitochondria_rate\tNRF\tPBC1\tPBC2"
    )?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        result.total_reads,
        result.distinct,
        result.m1,
        result.m2,
        result.dup_rate(),
        result.mito_reads,
        result.mito_rate(),
        result.nrf,
        result.pbc1,
        result.pbc2
    )?;
    Ok(())
}

pub fn run_bam_qc(bam_path: &Path, output_path: &Path) -> Result<BamQcResult, Box<dyn Error>> {
    let result = compute_bam_qc(bam_path)?;
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    write_bam_qc_tsv(&result, &mut writer)?;
    Ok(result)
}

pub fn run_bam_qc_parallel(bam_path: &Path, output_path: &Path, num_threads: usize) -> Result<BamQcResult, Box<dyn Error>> {
    let result = compute_bam_qc_parallel(bam_path, num_threads)?;
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    write_bam_qc_tsv(&result, &mut writer)?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn test_data_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/data")
    }

    #[test]
    fn test_is_mitochondrial() {
        assert!(is_mitochondrial("chrM"));
        assert!(is_mitochondrial("chrm"));
        assert!(is_mitochondrial("MT"));
        assert!(is_mitochondrial("chrMT"));
        assert!(is_mitochondrial("rCRSd"));
        assert!(!is_mitochondrial("chr1"));
        assert!(!is_mitochondrial("chr22"));
    }

    #[test]
    fn test_compute_bam_qc_small() {
        let bam_path = test_data_path().join("test_chr22_small.bam");
        if !bam_path.exists() {
            eprintln!("Test BAM file not found, skipping: {:?}", bam_path);
            return;
        }
        let result = compute_bam_qc(&bam_path).expect("Failed to compute BAM QC");
        assert!(result.total_reads > 0, "Expected some reads");
        assert!(result.nrf >= 0.0 && result.nrf <= 1.0, "NRF should be between 0 and 1");
        assert!(result.pbc1 >= 0.0 && result.pbc1 <= 1.0, "PBC1 should be between 0 and 1");
    }

    #[test]
    fn test_compute_bam_qc_dummy() {
        let bam_path = test_data_path().join("dummy.bam");
        if !bam_path.exists() {
            eprintln!("Test BAM file not found, skipping: {:?}", bam_path);
            return;
        }
        let result = compute_bam_qc(&bam_path).expect("Failed to compute BAM QC");
        println!("Dummy BAM results: {:?}", result);
    }
}
