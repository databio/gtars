use crate::config::{SpecialToken, SpecialTokenAssignment};

#[derive(Clone, Debug)]
pub struct SpecialTokens {
    pub unk: String,
    pub pad: String,
    pub mask: String,
    pub cls: String,
    pub eos: String,
    pub bos: String,
    pub sep: String,
}

impl Default for SpecialTokens {
    fn default() -> Self {
        SpecialTokens {
            unk: "<unk>".to_string(),
            pad: "<pad>".to_string(),
            mask: "<mask>".to_string(),
            cls: "<cls>".to_string(),
            eos: "<eos>".to_string(),
            bos: "<bos>".to_string(),
            sep: "<sep>".to_string(),
        }
    }
}

impl From<Vec<SpecialTokenAssignment>> for SpecialTokens {
    fn from(value: Vec<SpecialTokenAssignment>) -> Self {
        let mut special_tokens = SpecialTokens::default();

        for token in value {
            let token_assignment = token.name;
            match token_assignment {
                SpecialToken::Unk => special_tokens.unk = token.token,
                SpecialToken::Pad => special_tokens.pad = token.token,
                SpecialToken::Mask => special_tokens.mask = token.token,
                SpecialToken::Cls => special_tokens.cls = token.token,
                SpecialToken::Bos => special_tokens.bos = token.token,
                SpecialToken::Eos => special_tokens.eos = token.token,
                SpecialToken::Sep => special_tokens.sep = token.token,
            }
        }

        special_tokens
    }
}

impl From<SpecialTokens> for Vec<String> {
    fn from(val: SpecialTokens) -> Self {
        vec![
            val.unk, val.pad, val.mask, val.cls, val.eos, val.bos, val.sep,
        ]
    }
}

impl From<&SpecialTokens> for Vec<String> {
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
