use extendr_api::prelude::*;
use gtars::tokenizers::{BitsTree, Universe};

#[extendr]
#[derive(Debug)]
struct Tokenizer {
    pub bits_tree: BitsTree
}

#[extendr]
impl Tokenizer {
    fn new(chrs: Robj, starts: Robj, ends: Robj) -> Result<Self> {

        // extract the Vec<String>
        let chrs: Vec<String> = if let Some(chrs) = chrs.as_string_vector() {
            chrs
        } else {
            return Err("Could not extract chrs as a slice of Strings".into())
        };

        // extract the starts as Vec<usize>
        let starts: Vec<usize> = if let Some(starts) = starts.as_integer_slice() {
            starts.iter().map(|i| *i as usize).collect()
        } else {
            return Err("Could not extract starts as a slice of usize.".into())
        };

        let ends: Vec<usize> = if let Some(ends) = ends.as_integer_slice() {
            ends.iter().map(|i| *i as usize).collect()
        } else {
            return Err("Could not extract ends as a slice of usize.".into())
        };

        let universe = Universe::from_vectors(&chrs, &starts, &ends);
        let bits_tree = BitsTree::from(universe);
        Ok(Tokenizer { bits_tree })
    }
}

extendr_module! {
    mod tokenizers;
    impl Tokenizer;
}

// the trait bound `&[std::string::String]: extendr_api::TryFrom<extendr_api::Robj>` is not satisfied
// the following other types implement trait `std::convert::From<T>`:
//   `&[u8]` implements `std::convert::From<&bstr::bstr::BStr>`
//   `&[u8]` implements `std::convert::From<&winnow::stream::bstr::BStr>`
//   `&[u8]` implements `std::convert::From<&winnow::stream::bytes::Bytes>`
//   `&[u8]` implements `std::convert::From<webpki::subject_name::ip_address::IpAddrRef<'_>>`
//   `[T; 1000]` implements `std::convert::From<generic_array::GenericArray<T, typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UTerm, typenum::bit::B1>, typenum::bit::B1>, typenum::bit::B1>, typenum::bit::B1>, typenum::bit::B1>, typenum::bit::B0>, typenum::bit::B1>, typenum::bit::B0>, typenum::bit::B0>, typenum::bit::B0>>>`
//   `[T; 100]` implements `std::convert::From<generic_array::GenericArray<T, typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UTerm, typenum::bit::B1>, typenum::bit::B1>, typenum::bit::B0>, typenum::bit::B0>, typenum::bit::B1>, typenum::bit::B0>, typenum::bit::B0>>>`
//   `[T; 1024]` implements `std::convert::From<generic_array::GenericArray<T, typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UInt<typenum::uint::UTerm, typenum::bit::B1>, typenum::bit::B0>, typenum::bit::B0>, typenum::bit::B0>, typenum::bit::B0>, typenum::bit::B0>, typenum::bit::B0>, typenum::bit::B0>, typenum::bit::B0>, typenum::bit::B0>, typenum::bit::B0>>>`
//   `[T; 10]` implements `std::convert::From<(T, T, T, T, T, T, T, T, T, T)>`
// and 93 others
// required for `extendr_api::Robj` to implement `std::convert::Into<&[std::string::String]>`
// required for `&[std::string::String]` to implement `extendr_api::TryFrom<extendr_api::Robj>`
// required for `extendr_api::Robj` to implement `extendr_api::TryInto<&[std::string::String]>`