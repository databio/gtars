pub struct Encoding {
    pub input_ids: Vec<u32>,
    pub attention_mask: Vec<u8>,
    pub special_tokens_mask: Vec<u8>,
}

pub struct BatchEncoding {
    pub encodings: Vec<Encoding>,
}
