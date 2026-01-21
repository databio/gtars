from typing import Any, Dict, List, Optional, Union

class Tokenizer:
    def __init__(self, path: str) -> None: ...
    @classmethod
    def from_config(cls, cfg: str) -> "Tokenizer":
        """
        Loads a tokenizer from a configuration file.
        Args:
            cfg (str): The path to the configuration file.
        Returns:
            Tokenizer: An instance of the Tokenizer class.
        """

    @classmethod
    def from_bed(cls, path: str) -> "Tokenizer":
        """
        Loads a tokenizer from a BED file.
        Args:
            path (str): The path to the BED file.
        Returns:
            Tokenizer: An instance of the Tokenizer class.
        """

    @classmethod
    def from_pretrained(cls, model_name: str) -> "Tokenizer":
        """
        Loads a pretrained tokenizer from model on Hugging Face.
        Args:
            model_name (str): The name of the pretrained model.
        Returns:
            Tokenizer: An instance of the Tokenizer class.
        """

    def tokenize(self, regions: Any) -> List[str]:
        """
        Tokenizes the input regions into a list of tokens.
        Args:
            regions (Any): The input regions to tokenize.
        Returns:
            List[str]: A list of tokens.
        """

    def encode(self, tokens: Any) -> Union[List[int], int]:
        """
        Encodes the input tokens into a list of token IDs.
        Args:
            tokens (Any): The input tokens to encode.
        Returns:
            Union[List[int], int]: A list of token IDs or a single token ID.
        """

    def decode(self, ids: Any) -> Union[List[str], str]:
        """
        Decodes the input token IDs into a list of tokens.
        Args:
            ids (Any): The input token IDs to decode.
        Returns:
            Union[List[str], str]: A list of tokens or a single token.
        """

    def convert_ids_to_tokens(self, id: Any) -> Union[List[str], str]:
        """
        Converts the input token IDs into a list of tokens.
        Args:
            id (Any): The input token IDs to convert.
        Returns:
            Union[List[str], str]: A list of tokens or a single token.
        """

    def convert_tokens_to_ids(self, region: Any) -> Union[List[int], int]:
        """
        Converts the input tokens into a list of token IDs.
        Args:
            region (Any): The input tokens to convert.
        Returns:
            Union[List[int], int]: A list of token IDs or a single token ID.
        """

    @property
    def unk_token(self) -> str: ...
    @property
    def pad_token(self) -> str: ...
    @property
    def mask_token(self) -> str: ...
    @property
    def cls_token(self) -> str: ...
    @property
    def bos_token(self) -> str: ...
    @property
    def eos_token(self) -> str: ...
    @property
    def sep_token(self) -> str: ...
    @property
    def pad_token_id(self) -> int: ...
    @property
    def mask_token_id(self) -> int: ...
    @property
    def cls_token_id(self) -> int: ...
    @property
    def bos_token_id(self) -> int: ...
    @property
    def eos_token_id(self) -> int: ...
    @property
    def sep_token_id(self) -> int: ...
    @property
    def unk_token_id(self) -> int: ...
    @property
    def vocab_size(self) -> int: ...
    @property
    def special_tokens_map(self) -> Dict[str, Optional[str]]:
        """
        Returns a dictionary mapping special tokens to their string representations.
        Returns:
            Dict[str, Optional[str]]: A dictionary mapping special tokens to their string representations.
        """

    def get_vocab(self) -> Dict[str, int]:
        """
        Returns the vocabulary of the tokenizer.
        Returns:
            Dict[str, int]: A dictionary mapping tokens to their IDs.
        """

    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...
    def __call__(self, regions: Any) -> Any: ...

class Universe:
    """
    A Python wrapper for the Universe class in Rust.
    """

    def __init__(self, universe: Any) -> None:
        """
        Initializes the Universe object.
        Args:
            universe (Any): The underlying Rust Universe object.
        """

    def len(self) -> int:
        """
        Returns the number of regions in the Universe.
        Returns:
            int: The number of regions.
        """

    def is_empty(self) -> bool:
        """
        Checks if the Universe is empty.
        Returns:
            bool: True if the Universe is empty, False otherwise.
        """

    def __len__(self) -> int:
        """
        Returns the number of regions in the Universe.
        Returns:
            int: The number of regions.
        """

    def __repr__(self) -> str:
        """
        Returns a string representation of the Universe.
        Returns:
            str: A string describing the Universe.
        """

def tokenize_fragment_file(file: str, tokenizer: Tokenizer) -> Dict[str, List[int]]:
    """
    Tokenizes a fragment file using the specified tokenizer.
    Args:
        file (str): The path to the fragment file.
        tokenizer (Tokenizer): The tokenizer to use for tokenization.
    Returns:
        Dict[str, List[int]]: A dictionary mapping cell barcodes to lists of token IDs.
    """
