// Universe stuff
pub const CHR_COL_NAME: &str = "chr";
pub const START_COL_NAME: &str = "start";
pub const END_COL_NAME: &str = "end";
pub const DELIMITER: char = '\t';

pub const BED_FILE_EXTENSION: &str = "bed";
pub const BAM_FILE_EXTENSION: &str = "bam";
pub const GZ_FILE_EXTENSION: &str = "gz";
pub const IGD_FILE_EXTENSION: &str = "igd";

// Special tokens
pub mod special_tokens {
    pub const PAD_CHR: &str = "chrPAD";
    pub const PAD_START: u8 = 0;
    pub const PAD_END: u8 = 0;
    pub const MASK_CHR: &str = "chrMASK";
    pub const MASK_START: u8 = 0;
    pub const MASK_END: u8 = 0;
    pub const UNKNOWN_CHR: &str = "chrUNK";
    pub const UNKNOWN_START: u8 = 0;
    pub const UNKNOWN_END: u8 = 0;
    pub const EOS_CHR: &str = "chrEOS";
    pub const EOS_START: u8 = 0;
    pub const EOS_END: u8 = 0;
    pub const BOS_CHR: &str = "chrBOS";
    pub const BOS_START: u8 = 0;
    pub const BOS_END: u8 = 0;
    pub const CLS_CHR: &str = "chrCLS";
    pub const CLS_START: u8 = 0;
    pub const CLS_END: u8 = 0;
    pub const SEP_CHR: &str = "chrSEP";
    pub const SEP_START: u8 = 0;
    pub const SEP_END: u8 = 0;
}

// GTOK stuff
pub const GTOK_EXT: &str = "gtok";
