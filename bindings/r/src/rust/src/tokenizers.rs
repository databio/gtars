use extendr_api::prelude::*;
use gtars::tokenizers::{BitsTree, Universe};

#[extendr]
#[derive(Debug)]
struct Tokenizer {
    pub bits_tree: BitsTree
}

#[extendr]
impl Tokenizer {
    fn new(chrs: &[String], starts: &[usize], ends: &[usize]) -> Self {
        let universe = Universe::from_vectors(chrs, starts, ends);
        let bits_tree = BitsTree::from(universe);
        Tokenizer { bits_tree }
    }
}

extendr_module! {
    mod tokenizers;
    impl Tokenizer;
}