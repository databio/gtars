use std::fmt::Display;

#[derive(Eq, PartialEq, Hash)]
pub enum SpecialToken {
    Unk,
    Pad,
    Mask,
    Cls,
    Bos,
    Eos,
    Sep,
}

impl Display for SpecialToken {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let display_str = match self {
            SpecialToken::Unk => "UNK",
            SpecialToken::Pad => "PAD",
            SpecialToken::Mask => "MASK",
            SpecialToken::Cls => "CLS",
            SpecialToken::Bos => "BOS",
            SpecialToken::Eos => "EOS",
            SpecialToken::Sep => "SEP",
        };
        write!(f, "{}", display_str)
    }
}
