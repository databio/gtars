use std::path::Path;

use special_tokens::SpecialTokens;

use super::{tokenizer_impl::TokenizerError, universe::Universe};

pub mod special_tokens;

///
/// Prepare the universe and special tokens. This function will build
/// the universe struct and prepare the special tokens if they are provided.
///
/// Doing these together is necessary, because the special tokens contribute
/// to the universe/vocab.
///
/// # Arguments:
/// - config: the tokenizer config
///
pub fn prepare_universe_and_special_tokens<P: AsRef<Path>>(
    universe_file: P,
    special_tokens: SpecialTokens,
) -> Result<(Universe, SpecialTokens), TokenizerError> {
    let mut universe = Universe::try_from(universe_file.as_ref())?;
    universe.add_special_tokens(&special_tokens);
    Ok((universe, special_tokens))
}
