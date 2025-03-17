use crate::common::models::Region;

#[derive(Clone, Debug)]
pub struct SpecialTokens {
    pub unk: Region,
    pub pad: Region,
    pub mask: Region,
    pub cls: Region,
    pub eos: Region,
    pub bos: Region,
    pub sep: Region,
}

impl Default for SpecialTokens {
    fn default() -> Self {
        SpecialTokens {
            unk: Region {
                chr: "chrUNK".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            pad: Region {
                chr: "chrPAD".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            mask: Region {
                chr: "chrMASK".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            cls: Region {
                chr: "chrCLS".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            eos: Region {
                chr: "chrEOS".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            bos: Region {
                chr: "chrBOS".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
            sep: Region {
                chr: "chrSEP".to_string(),
                start: 0,
                end: 0,
                rest: None,
            },
        }
    }
}

impl From<SpecialTokens> for Vec<Region> {
    fn from(val: SpecialTokens) -> Self {
        vec![
            val.unk, val.pad, val.mask, val.cls, val.eos, val.bos, val.sep,
        ]
    }
}