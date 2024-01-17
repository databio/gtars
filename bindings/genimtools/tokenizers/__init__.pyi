from typing import List

class Region:
    chr: str
    start: int
    end: int

class TokenizedRegion:
    chr: str
    start: int
    end: int
    id: int

class TokenizedRegionSet:
    def regions(self) -> List[TokenizedRegion]:
        """
        Get the list of regions in this set.
        """
        
    def ids(self) -> List[int]:
        """
        Get the list of IDs in this set.
        """
    
    def ids_as_strs(self) -> List[str]:
        """
        Get the list of IDs in this set as strings. This
        is specifically meant for geniml workflows which requires
        strings
        """

class TreeTokenizer:
    def __init__(self, universe: str) -> None:
        """
        Instatiate a new TreeTokenizer. This tokenizer takes advantage of
        interval trees to compute genomic interval overlaps efficiently.
        
        :param universe: The universe of characters to use for tokenization.
        """
        pass
    
    def tokenize(self, regions: List[Region]) -> TokenizedRegionSet:
        """
        Tokenize a list of regions.
        """
    
    def tokenize_bed_file(self, bed_file: str) -> TokenizedRegionSet:
        """
        Tokenize a bed file directly.

        This was added to create a more performant tokenization strategy
        that could tokenize directly from disk in rust instead of using
        pandas.
        """