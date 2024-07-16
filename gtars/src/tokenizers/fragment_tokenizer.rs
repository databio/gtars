use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

use crate::common::models::{Region, TokenizedRegionSet};
use crate::common::utils::get_dynamic_reader;
use crate::io::consts::{GTOK_HEADER, GTOK_U32_FLAG};

use anyhow::{Context, Result};

use super::Tokenizer;

pub struct FragmentTokenizer<T>
where
    T: Tokenizer,
{
    pub tokenizer: T,
}

impl<T> FragmentTokenizer<T>
where
    T: Tokenizer,
{
    pub fn new(tokenizer: T) -> Self {
        Self { tokenizer }
    }
}

impl<T> FragmentTokenizer<T>
where
    T: Tokenizer,
{
    fn parse_fragment_file_line(line: String) -> Result<(String, u32, u32, String, u32)> {
        let fields: Vec<&str> = line.split_whitespace().collect();

        if fields.len() < 4 {
            anyhow::bail!("Detected improper number of fields detected");
        }

        let chr = fields[0];
        let start = fields[1].parse::<u32>()?;
        let end = fields[2].parse::<u32>()?;
        let barcode = fields[3];
        let _read_support = fields[4].parse::<u32>()?;

        Ok((
            chr.to_string(),
            start,
            end,
            barcode.to_string(),
            _read_support,
        ))
    }

    fn init_gtok_file(filename: &str) -> Result<()> {
        // make sure the path exists
        let path = std::path::Path::new(filename);

        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        } else {
            anyhow::bail!("Failed to create parent directories for gtok file!")
        }

        let file = File::create(filename).with_context(|| "Failed to create gtok file!")?;
        let mut writer = BufWriter::new(file);

        writer
            .write_all(GTOK_HEADER)
            .with_context(|| "Failed to write GTOK header to file!")?;

        // assume large and write u32 flag
        writer
            .write_all(&GTOK_U32_FLAG.to_le_bytes())
            .with_context(|| "Failed to write GTOK size flag to file!")?;

        Ok(())
    }

    fn append_tokens_to_gtok_file(writer: &mut BufWriter<File>, tokens: &[u32]) -> Result<()> {
        for token in tokens.iter() {
            writer
                .write_all(&token.to_le_bytes())
                .with_context(|| "Failed to write token to gtok file")?;
        }

        Ok(())
    }

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
    pub fn tokenize_fragments_to_gtoks(
        &self,
        fragments_file_path: &Path,
        out_path: &Path,
    ) -> Result<()> {
        let reader = get_dynamic_reader(fragments_file_path)?;

        let mut file_map: HashMap<String, BufWriter<File>> = HashMap::new();

        for (line_num, line) in reader.lines().enumerate() {
            let line = line
                .with_context(|| format!("Failed parsing line {} in fragments file", line_num))?;

            let (chr, start, end, barcode, _read_support) = Self::parse_fragment_file_line(line)
                .with_context(|| format!("Failed parsing line {} in fragments file", line_num))?;

            let r = Region {
                chr: chr.to_string(),
                start,
                end,
            };

            // get actual tokens
            let tokens = self.tokenizer.tokenize_region(&r);

            // get file path
            let file = file_map.get_mut(&barcode);

            // determine if we need to create a new file or append to an existing one
            match file {
                Some(fh) => {
                    // append to file
                    Self::append_tokens_to_gtok_file(fh, &tokens.ids)?;
                }
                None => {
                    // create new file
                    let new_file_name = out_path.join(format!("{}.gtok", barcode));
                    let new_file_name = new_file_name.to_str().unwrap();
                    Self::init_gtok_file(new_file_name)?;

                    // append to file -- need to open file again in append mode
                    let file = OpenOptions::new()
                        .append(true)
                        .open(new_file_name)
                        .with_context(|| "Failed to open gtok file for appending")?;

                    let mut writer = BufWriter::new(file);
                    Self::append_tokens_to_gtok_file(&mut writer, &tokens.ids)?;

                    // insert into map
                    file_map.insert(barcode.to_string(), writer);
                }
            }
        }

        Ok(())
    }

    ///
    /// Tokenize a fragments file into a known universe. This version will filter out any cell barcodes that are not in the provided filter.
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
    pub fn tokenize_fragments_to_gtoks_with_filter(
        &self,
        fragments_file_path: &Path,
        out_path: &Path,
        filter: Vec<String>,
    ) -> Result<()> {
        let reader = get_dynamic_reader(fragments_file_path)?;

        let mut file_map: HashMap<String, BufWriter<File>> = HashMap::new();
        let filter: HashSet<String> = HashSet::from_iter(filter);

        for (line_num, line) in reader.lines().enumerate() {
            if line_num % 10_000 == 0 {
                print!("Processed {} lines\r", line_num);
                std::io::stdout().flush().unwrap();
            }

            let line = line
                .with_context(|| format!("Failed parsing line {} in fragments file", line_num))?;

            let (chr, start, end, barcode, _read_support) = Self::parse_fragment_file_line(line)
                .with_context(|| format!("Failed parsing line {} in fragments file", line_num))?;

            if !filter.contains(&barcode) {
                // skip! -- barcode not in filter
                continue;
            }

            let r = Region {
                chr: chr.to_string(),
                start,
                end,
            };

            // get actual tokens
            let tokens = self.tokenizer.tokenize_region(&r);

            // get file path
            let file = file_map.get_mut(&barcode);
            match file {
                Some(fh) => {
                    // append to file
                    Self::append_tokens_to_gtok_file(fh, &tokens.ids)?;
                }
                None => {
                    // create new file
                    let new_file_name = out_path.join(format!("{}.gtok", barcode));
                    let new_file_name = new_file_name.to_str().unwrap();
                    Self::init_gtok_file(new_file_name)?;

                    // append to file -- need to open file again in append mode
                    let file = OpenOptions::new()
                        .append(true)
                        .open(new_file_name)
                        .with_context(|| "Failed to open gtok file for appending")?;

                    let mut writer = BufWriter::new(file);
                    Self::append_tokens_to_gtok_file(&mut writer, &tokens.ids)?;

                    // insert into map
                    file_map.insert(barcode.to_string(), writer);
                }
            }
        }

        Ok(())
    }

    ///
    /// Tokenize a fragments file into a vector of TokenizedRegionSets. A fragments file represents a collection of cells,
    /// therefore we need to demultiplex the reads and tokenize each cell separately.
    ///
    /// This function will consume more memory than the `tokenize_fragments_to_gtoks` function, because it
    /// will store all cells and their tokens in memory.
    ///
    /// # Arguments
    /// - `fragments_file_path` - the path to the fragments file
    pub fn tokenize_fragments(
        &self,
        fragments_file_path: &Path,
    ) -> Result<Vec<TokenizedRegionSet>> {
        let reader = get_dynamic_reader(fragments_file_path)?;

        let mut barcode_ids_map: HashMap<String, Vec<u32>> = HashMap::new();

        for (line_num, line) in reader.lines().enumerate() {
            if line_num % 10_000 == 0 {
                print!("Processed {} lines\r", line_num);
                std::io::stdout().flush().unwrap();
            }

            let line = line
                .with_context(|| format!("Failed parsing line {} in fragments file", line_num))?;

            if line_num % 10_000 == 0 {
                print!("Processed {} lines\r", line_num);
                std::io::stdout().flush().unwrap();
            }

            let (chr, start, end, barcode, _read_support) = Self::parse_fragment_file_line(line)
                .with_context(|| format!("Failed parsing line {} in fragments file", line_num))?;

            let r = Region {
                chr: chr.to_string(),
                start,
                end,
            };

            // get actual tokens
            let tokens = self.tokenizer.tokenize_region(&r);

            let barcode_tokens = barcode_ids_map.entry(barcode).or_default();

            barcode_tokens.extend(tokens.ids);
        }

        Ok(barcode_ids_map
            .into_values()
            .map(|ids| TokenizedRegionSet::new(ids, self.tokenizer.get_universe()))
            .collect())
    }

    pub fn tokenize_fragments_with_filter(
        &self,
        fragments_file_path: &Path,
        filter: Vec<String>,
    ) -> Result<Vec<TokenizedRegionSet>> {
        let reader = get_dynamic_reader(fragments_file_path)?;

        let mut barcode_ids_map: HashMap<String, Vec<u32>> = HashMap::new();
        let filter: HashSet<String> = HashSet::from_iter(filter);

        for (line_num, line) in reader.lines().enumerate() {
            if line_num % 10_000 == 0 {
                print!("Processed {} lines\r", line_num);
                std::io::stdout().flush().unwrap();
            }

            let line = line
                .with_context(|| format!("Failed parsing line {} in fragments file", line_num))?;

            let (chr, start, end, barcode, _read_support) = Self::parse_fragment_file_line(line)
                .with_context(|| format!("Failed parsing line {} in fragments file", line_num))?;

            if !filter.contains(&barcode) {
                // skip! -- barcode not in filter
                continue;
            }

            let r = Region {
                chr: chr.to_string(),
                start,
                end,
            };

            // get actual tokens
            let tokens = self.tokenizer.tokenize_region(&r);

            let barcode_tokens = barcode_ids_map.entry(barcode).or_default();

            barcode_tokens.extend(tokens.ids);
        }

        Ok(barcode_ids_map
            .into_values()
            .map(|ids| TokenizedRegionSet::new(ids, self.tokenizer.get_universe()))
            .collect())
    }
}
