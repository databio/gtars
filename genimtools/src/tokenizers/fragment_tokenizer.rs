use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::common::models::Region;
use crate::tokenizers::TreeTokenizer;
use anyhow::{Context, Result};

use super::Tokenizer;

pub struct FragmentTokenizer {
    tokenizer: TreeTokenizer,
}

impl TryFrom<&Path> for FragmentTokenizer {
    type Error = anyhow::Error;
    ///
    /// # Arguments
    /// - `value` - the path to the bed file
    ///
    /// # Returns
    /// A new FragmentTokenizer
    fn try_from(value: &Path) -> Result<Self> {
        let tokenizer = TreeTokenizer::try_from(value)?;
        Ok(Self { tokenizer })
    }
}

impl FragmentTokenizer {
    ///
    /// Tokenize a fragments file into a known universe.
    ///
    /// A `fragments.tsv.gz` file is a tab-separated file with the following columns:
    /// | Column Number | Name        | Description                                                                                                      |
    /// |---------------|-------------|------------------------------------------------------------------------------------------------------------------|
    /// | 1             | chrom       | Reference genome chromosome of fragment                                                                           |
    /// | 2             | chromStart  | Adjusted start position of fragment on chromosome.                                                                |
    /// | 3             | chromEnd    | Adjusted end position of fragment on chromosome. The end position is exclusive, so represents the position immediately following the fragment interval. |
    /// | 4             | barcode     | The 10x Barcode of this fragment. This corresponds to the CB tag attached to the corresponding BAM file records for this fragment.                   |
    /// | 5             | readSupport | The total number of read pairs associated with this fragment. This includes the read pair marked unique and all duplicate read pairs.                 |
    ///
    /// For performance reasons, this will only tokenize cells onto disk as a `.gtok` file, since storing all cells in memory is not feasible
    ///
    ///
    /// # Arguments
    /// - `fragments_file_path` - the path to the fragments file
    pub fn tokenize_fragments(&self, fragments_file_path: &Path) -> Result<()> {
        let is_gzipped = fragments_file_path.extension() == Some(OsStr::new("gz"));
        let file = File::open(fragments_file_path).with_context(|| "Failed to open bed file.")?;

        let file: Box<dyn Read> = match is_gzipped {
            true => Box::new(GzDecoder::new(file)),
            false => Box::new(file),
        };

        let reader = BufReader::new(file);

        let file_map: HashMap<String, File> = HashMap::new();

        for (line_num, line) in reader.lines().enumerate() {
            let line = line
                .with_context(|| format!("Failed parsing line {} in fragments file", line_num))?;

            let fields: Vec<&str> = line.split('\t').collect();

            if fields.len() != 5 {
                anyhow::bail!(
                    "Detected improper number of fields in fragments file line {}: {}",
                    line_num,
                    line
                );
            }

            let chr = fields[0];
            let start = fields[1].parse::<u32>().with_context(|| {
                format!(
                    "Failed to parse start position in fragments file line {}: {}",
                    line_num, line
                )
            })?;
            let end = fields[2].parse::<u32>().with_context(|| {
                format!(
                    "Failed to parse end position in fragments file line {}: {}",
                    line_num, line
                )
            })?;
            let barcode = fields[3];
            let read_support = fields[4].parse::<u32>().with_context(|| {
                format!(
                    "Failed to parse read support in fragments file line {}: {}",
                    line_num, line
                )
            })?;

            let r = Region {
                chr: chr.to_string(),
                start,
                end,
            };

            // get actual tokens
            let tokens = self.tokenizer.tokenize_region(&r);

            // get file path
            let file_handle = file_map.get(barcode);
            match file_handle {
                Some(fh) => {
                    // append to file
                }
                None => {}
            }
        }

        Ok(())
    }
}
