//! BED file format classifier
//!
//! This module provides functionality to classify BED files according to their format
//! and compliance with UCSC and ENCODE specifications.

use crate::errors::BedClassifierError;
use gtars_core::models::RegionSet;
use regex::Regex;
use std::fmt::{self, Display};

#[cfg(feature = "bedclassifier")]
use polars::prelude::*;

///
/// Data format enumeration for BED file classification
///
#[derive(Clone, Debug, PartialEq)]
pub enum DataFormat {
    Unknown,
    UcscBed,
    UcscBedRs,
    BedLike,
    BedLikeRs,
    EncodeNarrowPeak,
    EncodeNarrowPeakRs,
    EncodeBroadPeak,
    EncodeBroadPeakRs,
    EncodeGappedPeak,
    EncodeGappedPeakRs,
    EncodeRnaElements,
    EncodeRnaElementsRs,
}

impl Display for DataFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            DataFormat::Unknown => "unknown_data_format",
            DataFormat::UcscBed => "ucsc_bed",
            DataFormat::UcscBedRs => "ucsc_bed_rs",
            DataFormat::BedLike => "bed_like",
            DataFormat::BedLikeRs => "bed_like_rs",
            DataFormat::EncodeNarrowPeak => "encode_narrowpeak",
            DataFormat::EncodeNarrowPeakRs => "encode_narrowpeak_rs",
            DataFormat::EncodeBroadPeak => "encode_broadpeak",
            DataFormat::EncodeBroadPeakRs => "encode_broadpeak_rs",
            DataFormat::EncodeGappedPeak => "encode_gappedpeak",
            DataFormat::EncodeGappedPeakRs => "encode_gappedpeak_rs",
            DataFormat::EncodeRnaElements => "encode_rna_elements",
            DataFormat::EncodeRnaElementsRs => "encode_rna_elements_rs",
        };
        write!(f, "{}", s)
    }
}

///
/// BED file classification output
///
#[derive(Clone, Debug)]
pub struct BedClassificationOutput {
    pub bed_compliance: String,
    pub data_format: DataFormat,
    pub compliant_columns: usize,
    pub non_compliant_columns: usize,
}

impl Display for BedClassificationOutput {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "BedClassificationOutput {{ bed_compliance: {}, data_format: {}, compliant_columns: {}, non_compliant_columns: {} }}",
            self.bed_compliance, self.data_format, self.compliant_columns, self.non_compliant_columns
        )
    }
}

///
/// Classify BED file format
///
/// This function analyzes the BED file structure and returns classification details.
/// It uses polars DataFrame to perform column validation similar to bedboss bedclassifier.
///
/// # Arguments
/// - region_set: Reference to a RegionSet to classify
///
/// # Returns
/// BedClassificationOutput with format classification details
///
/// # Example
/// ```no_run
/// use gtars_genomicdist::bed_classifier::classify_bed;
/// use gtars_core::models::RegionSet;
///
/// let region_set = RegionSet::try_from("file.bed").unwrap();
/// let classification = classify_bed(&region_set).unwrap();
/// println!("Format: {}", classification.data_format);
/// ```
///
#[cfg(feature = "bedclassifier")]
pub fn classify_bed(region_set: &RegionSet) -> Result<BedClassificationOutput, BedClassifierError> {
    // Get polars DataFrame
    let df = match region_set.to_polars() {
        Ok(df) => df,
        Err(_) => {
            return Ok(BedClassificationOutput {
                bed_compliance: "unknown_bed_compliance".to_string(),
                data_format: DataFormat::Unknown,
                compliant_columns: 0,
                non_compliant_columns: 0,
            });
        }
    };

    let num_cols = df.width();
    let mut compliant_columns = 0;
    let mut bed_format_named = DataFormat::UcscBed;
    let mut relaxed = false;

    // Helper function to check if a column matches a regex pattern
    let check_string_column = |df: &DataFrame, col_idx: usize, pattern: &str| -> bool {
        if let Ok(col) = df.column(&format!("column_{}", col_idx + 1)) {
            if let Ok(series) = col.cast(&DataType::String) {
                let regex = Regex::new(pattern).unwrap();
                if let Ok(ca) = series.str() {
                    return ca
                        .into_iter()
                        .all(|opt| opt.map(|s| regex.is_match(s)).unwrap_or(false));
                }
            }
        }
        false
    };

    // Helper function to check if a column is integer and meets conditions
    let check_int_column =
        |df: &DataFrame, col_idx: usize, min_val: Option<i64>, max_val: Option<i64>| -> bool {
            if let Ok(col) = df.column(&format!("column_{}", col_idx + 1)) {
                if let Ok(ca) = col.i64() {
                    return ca.into_iter().all(|opt| {
                        opt.map(|val| {
                            let mut valid = true;
                            if let Some(min) = min_val {
                                valid = valid && val >= min;
                            }
                            if let Some(max) = max_val {
                                valid = valid && val <= max;
                            }
                            valid
                        })
                        .unwrap_or(false)
                    });
                }
            }
            false
        };

    // Helper function to check if column is float or contains -1
    let check_float_or_minus_one = |df: &DataFrame, col_idx: usize| -> bool {
        if let Ok(col) = df.column(&format!("column_{}", col_idx + 1)) {
            // Try as float first
            if col.dtype().is_float() {
                return true;
            }
            // Try as int and check if all values are -1
            if let Ok(ca) = col.i64() {
                return ca
                    .into_iter()
                    .all(|opt| opt.map(|v| v == -1).unwrap_or(false));
            }
        }
        false
    };

    // RGB color pattern
    let regex_colors =
        r"^(?:\d|[1-9]\d|1\d{2}|2[0-4]\d|25[0-5])(?:,(?:\d|[1-9]\d|1\d{2}|2[0-4]\d|25[0-5])){0,2}$";

    // Validate columns
    for col_idx in 0..num_cols {
        let is_valid = match col_idx {
            0 => check_string_column(&df, col_idx, r"[A-Za-z0-9_]{1,255}"),
            1 => check_int_column(&df, col_idx, Some(0), None),
            2 => check_int_column(&df, col_idx, Some(0), None),
            3 => check_string_column(&df, col_idx, r"[\x20-\x7e]{1,255}"),
            4 => {
                // Standard check: integer 0-1000
                if check_int_column(&df, col_idx, Some(0), Some(1000)) {
                    true
                } else {
                    // Relaxed check: any non-negative integer
                    if check_int_column(&df, col_idx, Some(0), None) {
                        relaxed = true;
                        true
                    } else {
                        false
                    }
                }
            }
            5 => {
                // Check for strand: +, -, or .
                if let Ok(col) = df.column(&format!("column_{}", col_idx + 1)) {
                    if let Ok(series) = col.cast(&DataType::String) {
                        if let Ok(ca) = series.str() {
                            ca.into_iter().all(|opt| {
                                opt.map(|s| s == "+" || s == "-" || s == ".")
                                    .unwrap_or(false)
                            })
                        } else {
                            false
                        }
                    } else {
                        false
                    }
                } else {
                    false
                }
            }
            6 => check_int_column(&df, col_idx, Some(0), None),
            7 => check_int_column(&df, col_idx, Some(0), None),
            8 => check_string_column(&df, col_idx, regex_colors),
            9 => check_int_column(&df, col_idx, None, None),
            10 => check_string_column(&df, col_idx, r"^(0(,\d+)*|\d+(,\d+)*)?,?$"),
            11 => check_string_column(&df, col_idx, r"^(0(,\d+)*|\d+(,\d+)*)?,?$"),
            12 => check_float_or_minus_one(&df, col_idx),
            13 => {
                if let Ok(col) = df.column(&format!("column_{}", col_idx + 1)) {
                    if let Ok(ca) = col.i64() {
                        if let Some(first) = ca.get(0) {
                            first != -1
                        } else {
                            false
                        }
                    } else {
                        false
                    }
                } else {
                    false
                }
            }
            _ => false,
        };

        if is_valid && col_idx < 12 {
            compliant_columns += 1;
        } else {
            let nccols = num_cols - compliant_columns;

            // Check for special ENCODE formats
            if col_idx >= 6 {
                // NarrowPeak: 10 columns
                if num_cols == 10
                    && col_idx == 6
                    && check_float_or_minus_one(&df, 6)
                    && check_float_or_minus_one(&df, 7)
                    && check_float_or_minus_one(&df, 8)
                    && check_int_column(&df, 9, None, None)
                {
                    bed_format_named = if relaxed {
                        DataFormat::EncodeNarrowPeakRs
                    } else {
                        DataFormat::EncodeNarrowPeak
                    };
                    return Ok(BedClassificationOutput {
                        bed_compliance: format!("bed{}+{}", compliant_columns, nccols),
                        data_format: bed_format_named,
                        compliant_columns,
                        non_compliant_columns: nccols,
                    });
                }

                // BroadPeak or RNA elements: 9 columns
                if num_cols == 9 && col_idx == 6 {
                    if check_float_or_minus_one(&df, 6)
                        && check_float_or_minus_one(&df, 7)
                        && check_float_or_minus_one(&df, 8)
                    {
                        bed_format_named = if relaxed {
                            DataFormat::EncodeBroadPeakRs
                        } else {
                            DataFormat::EncodeBroadPeak
                        };
                        return Ok(BedClassificationOutput {
                            bed_compliance: format!("bed{}+{}", compliant_columns, nccols),
                            data_format: bed_format_named,
                            compliant_columns,
                            non_compliant_columns: nccols,
                        });
                    } else if check_float_or_minus_one(&df, 6) && check_float_or_minus_one(&df, 7) {
                        // Check if column 8 first value is not -1 (RNA elements)
                        if let Ok(col) = df.column(&format!("column_{}", 9)) {
                            if let Ok(ca) = col.i64() {
                                if let Some(first) = ca.get(0) {
                                    if first != -1 {
                                        bed_format_named = if relaxed {
                                            DataFormat::EncodeRnaElementsRs
                                        } else {
                                            DataFormat::EncodeRnaElements
                                        };
                                        return Ok(BedClassificationOutput {
                                            bed_compliance: format!(
                                                "bed{}+{}",
                                                compliant_columns, nccols
                                            ),
                                            data_format: bed_format_named,
                                            compliant_columns,
                                            non_compliant_columns: nccols,
                                        });
                                    }
                                }
                            }
                        }
                    }
                }

                // GappedPeak: 15 columns
                if num_cols == 15
                    && col_idx == 12
                    && check_float_or_minus_one(&df, 12)
                    && check_float_or_minus_one(&df, 13)
                    && check_float_or_minus_one(&df, 14)
                {
                    bed_format_named = if relaxed {
                        DataFormat::EncodeGappedPeakRs
                    } else {
                        DataFormat::EncodeGappedPeak
                    };
                    return Ok(BedClassificationOutput {
                        bed_compliance: format!("bed{}+{}", compliant_columns, nccols),
                        data_format: bed_format_named,
                        compliant_columns,
                        non_compliant_columns: nccols,
                    });
                }
            }

            // Default to BED_LIKE if validation fails
            bed_format_named = if relaxed {
                if nccols == 0 {
                    DataFormat::UcscBedRs
                } else {
                    DataFormat::BedLikeRs
                }
            } else {
                DataFormat::BedLike
            };

            return Ok(BedClassificationOutput {
                bed_compliance: format!("bed{}+{}", compliant_columns, nccols),
                data_format: bed_format_named,
                compliant_columns,
                non_compliant_columns: nccols,
            });
        }
    }

    // All columns validated successfully
    bed_format_named = if relaxed {
        DataFormat::UcscBedRs
    } else {
        DataFormat::UcscBed
    };

    Ok(BedClassificationOutput {
        bed_compliance: format!("bed{}+0", compliant_columns),
        data_format: bed_format_named,
        compliant_columns,
        non_compliant_columns: 0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn get_test_path(file_name: &str) -> std::path::PathBuf {
        std::env::current_dir()
            .unwrap()
            .join("../tests/data/regionset")
            .join(file_name)
    }

    #[cfg(feature = "bedclassifier")]
    #[test]
    fn test_classify_bed_narrowpeak() {
        let file_path = get_test_path("dummy.narrowPeak");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let classification = classify_bed(&region_set).unwrap();
        println!("Classification: {}", classification);
        println!("Bed compliance: {}", classification.bed_compliance);
        println!("Data format: {}", classification.data_format);
        println!("Compliant columns: {}", classification.compliant_columns);
        println!(
            "Non-compliant columns: {}",
            classification.non_compliant_columns
        );

        // NarrowPeak files should have 10 columns (bed6+4)
        assert!(classification.bed_compliance.starts_with("bed"));
        assert_eq!(classification.data_format, DataFormat::EncodeNarrowPeak);
    }

    #[cfg(feature = "bedclassifier")]
    #[test]
    fn test_classify_bed_basic() {
        let file_path = get_test_path("dummy_headers.bed");
        let region_set = RegionSet::try_from(file_path.to_str().unwrap()).unwrap();

        let classification = classify_bed(&region_set).unwrap();
        println!("Classification: {}", classification);
        println!("Bed compliance: {}", classification.bed_compliance);
        println!("Data format: {}", classification.data_format);

        // Basic BED files should be classified
        assert!(classification.bed_compliance.starts_with("bed"));
        assert!(classification.compliant_columns >= 3); // At least chr, start, end
        assert_eq!(classification.data_format, DataFormat::UcscBed);
    }
}
