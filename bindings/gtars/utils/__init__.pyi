from typing import List
from ..models import TokenizedRegionPointer

def write_tokens_to_gtok(filename: str, tokens: List[str]):
    """
    Write tokens to a GTOK file.

    :param filename: The filename of the GTOK file.
    :param tokens: The tokens to write.
    """
    pass

def read_tokens_from_gtok(filename: str) -> List[int]:
    """
    Read tokens from a GTOK file.

    :param filename: The filename of the GTOK file.
    """
    pass

def write_tokens_to_gtokp(filename: str, tokens: List[TokenizedRegionPointer]):
    """
    Write tokens with positional information to a GTOKP file.

    :param filename: The filename of the GTOK file.
    :param tokens: The tokens to write.
    """
    pass

def read_tokens_from_gtokp(filename: str, unzip: bool) -> List[TokenizedRegionPointer]:
    """
    Read tokens from a GTOKP file with included positional information.

    :param filename: The filename of the GTOK file.
    :param unzip: Whether to return the positional information separately.
    """
    pass