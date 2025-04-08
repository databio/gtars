use extendr_api::prelude::*;
use gtars::tokenizers::{BitsTree, GTokenize, Universe};

pub fn extract_chrs_starts_ends_from_robjs(
    chrs: Robj,
    starts: Robj,
    ends: Robj,
) -> Result<(Vec<String>, Vec<usize>, Vec<usize>)> {
    // extract the Vec<String>
    let chrs: Vec<String> = if let Some(chrs) = chrs.as_string_vector() {
        chrs
    } else {
        return Err("Could not extract chrs as a slice of Strings".into());
    };

    // extract the starts as Vec<usize>
    let starts: Vec<usize> = if let Some(starts) = starts.as_integer_vector() {
        starts.into_iter().map(|i| i as usize).collect()
    } else {
        return Err("Could not extract starts as a slice of usize.".into());
    };

    let ends: Vec<usize> = if let Some(ends) = ends.as_integer_vector() {
        ends.into_iter().map(|i| i as usize).collect()
    } else {
        return Err("Could not extract ends as a slice of usize.".into());
    };

    Ok((chrs, starts, ends))
}

#[extendr]
#[derive(Debug)]
struct Tokenizer {
    pub bits_tree: BitsTree,
}

#[extendr]
impl Tokenizer {
    pub fn new(chrs: Robj, starts: Robj, ends: Robj) -> Result<Self> {
        let (chrs, starts, ends) = extract_chrs_starts_ends_from_robjs(chrs, starts, ends)?;

        let universe = Universe::from_vectors(&chrs, &starts, &ends);
        let bits_tree = BitsTree::from(universe);
        Ok(Tokenizer { bits_tree })
    }

    ///
    /// Tokenize a list of regions (passed as chrs, starts, and ends) into the universe.
    ///
    /// This will return two vectors:
    ///  1. The universe hits, and
    ///  2. the query hits
    ///
    pub fn tokenize(&self, gr: S4) {

        let seqnames = gr.get_slot("seqnames").unwrap();
        let seqnames_s4 = S4::try_from(seqnames).unwrap();
        // TODO: get seqnames from Rle S4 object somehow...

        // get "ranges", slot -- coerce back to S4
        let ranges = gr.get_slot("ranges").unwrap();
        let ranges_s4 = S4::try_from(ranges).unwrap();

        // get "starts", slot -- coerce to integer slice
        let starts = ranges_s4.get_slot("start").unwrap().as_integer_slice().unwrap();

        // get "widths", slot -- coerce to integer slice
        // then use widths to build ends (GRanges doesnt store ends, actually!)
        let widths = ranges_s4.get_slot("width").unwrap().as_integer_slice().unwrap();
        let ends = starts.iter().zip(widths)
            .map(|(s, w)| {
                s + w
            })
            .collect::<Vec<i32>>();
        
        // let starts = gr.get_slot("ranges").unwrap()

        // let mut query_hits = Vec::new();
        // let mut subject_hits = Vec::new();

        // let timing_start = std::time::Instant::now();

        // chrs.iter()
        //     .enumerate()
        //     .zip(starts.iter())
        //     .zip(ends.iter())
        //     .for_each(|(((query_indx, chr), &start), &end)| {
        //     if let Some(tree) = self.bits_tree.get_tree_for_chr(chr) {
        //         tree.find(start as u32, end as u32)
        //         .for_each(|iv| {
        //             query_hits.push(query_indx + 1);
        //             // val is that actual id of the region in the universe
        //             subject_hits.push((iv.val + 1) as usize);
        //         });
        //     }
        //     });

        // let timing_duration = timing_start.elapsed();
        // println!("Time taken for tokenization: {:?}", timing_duration);

        // let result = vec![
        //     query_hits.into_robj(), 
        //     subject_hits.into_robj()
        // ];

        // Ok(result)

    }

    pub fn universe_length(&self) -> usize {
        self.bits_tree.get_vocab_size()
    }
}

extendr_module! {
    mod tokenizers;
    impl Tokenizer;
}