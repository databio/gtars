from typing import Any, Dict, List, Optional, Union

class Tokenizer:
    def __init__(self, path: str) -> None: ...
    
    @classmethod
    def from_config(cls, cfg: str) -> 'Tokenizer':
        """
        Loads a tokenizer from a configuration file.
        Args:
            cfg (str): The path to the configuration file.
        Returns:
            Tokenizer: An instance of the Tokenizer class.
        """
    
    @classmethod
    def from_bed(cls, path: str) -> 'Tokenizer':
        """
        Loads a tokenizer from a BED file.
        Args:
            path (str): The path to the BED file.
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
    def get_unk_token(self) -> str: ...

    @property
    def get_pad_token(self) -> str: ...
    
    @property
    def get_mask_token(self) -> str: ...
    
    @property
    def get_cls_token(self) -> str: ...
    
    @property
    def get_bos_token(self) -> str: ...
    
    @property
    def get_eos_token(self) -> str: ...
    
    @property
    def get_sep_token(self) -> str: ...
    
    @property
    def get_pad_token_id(self) -> int: ...
    
    @property
    def get_mask_token_id(self) -> int: ...
    
    @property
    def get_cls_token_id(self) -> int: ...
    
    @property
    def get_bos_token_id(self) -> int: ...
    
    @property
    def get_eos_token_id(self) -> int: ...
    
    @property
    def get_sep_token_id(self) -> int: ...
    
    @property
    def get_unk_token_id(self) -> int: ...
    
    @property
    def get_vocab_size(self) -> int: ...
    
    @property
    def get_special_tokens_map(self) -> Dict[str, Optional[str]]:
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

def create_instances(
        sequences: Union[List[int], List[List[int]]],
        window_size: int,
        algorithm: str,
    ) -> List[Dict[str, Union[int, List[int]]]]:
        """
        Creates training instances for a given sequence or list of sequences.

        Args:
            sequences (Union[List[int], List[List[int]]]): A sequence or list of sequences of token IDs.
            window_size (int): The size of the context window.
            algorithm (str): The algorithm to use ('cbow' or 'sg').

        Returns:
            List[Dict[str, Union[int, List[int]]]]: A list of dictionaries representing the training instances.
        """