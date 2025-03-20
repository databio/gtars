use crate::{
    common::models::Region,
    tokenizers::config::{SpecialToken, SpecialTokenAssignment},
};

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

impl From<Vec<SpecialTokenAssignment>> for SpecialTokens {
    fn from(value: Vec<SpecialTokenAssignment>) -> Self {
        let mut special_tokens = SpecialTokens::default();

        for token in value {
            let region_str = token.token;
            let token_assignment = token.name;

            let parts = region_str.split(":").collect::<Vec<&str>>();

            // todo: Make this better
            if parts.len() != 2 {
                println!("Invalid region string: {}... Skipping", region_str);
                continue;
            }

            let chr = parts[0].to_string();

            let start_end = parts[1].split("-").collect::<Vec<&str>>();
            if start_end.len() != 2 {
                println!("Invalid start-end string: {}... Skipping", region_str);
                continue;
            }

            let start = start_end[0].parse::<u32>().unwrap();
            let end = start_end[1].parse::<u32>().unwrap();

            let region = Region {
                chr,
                start,
                end,
                rest: None,
            };
            match token_assignment {
                SpecialToken::Unk => special_tokens.unk = region,
                SpecialToken::Pad => special_tokens.pad = region,
                SpecialToken::Mask => special_tokens.mask = region,
                SpecialToken::Cls => special_tokens.cls = region,
                SpecialToken::Bos => special_tokens.bos = region,
                SpecialToken::Eos => special_tokens.eos = region,
                SpecialToken::Sep => special_tokens.sep = region,
            }
        }

        special_tokens
    }
}

impl From<SpecialTokens> for Vec<Region> {
    fn from(val: SpecialTokens) -> Self {
        vec![
            val.unk, val.pad, val.mask, val.cls, val.eos, val.bos, val.sep,
        ]
    }
}

impl From<&SpecialTokens> for Vec<Region> {
    fn from(val: &SpecialTokens) -> Self {
        vec![
            val.unk.clone(),
            val.pad.clone(),
            val.mask.clone(),
            val.cls.clone(),
            val.eos.clone(),
            val.bos.clone(),
            val.sep.clone(),
        ]
    }
}
