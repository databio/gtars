from typing import List, Tuple, Iterator

class Universe:
    """
    A Universe object represents a set of regions.
    
    All tokenizers are created from a Universe object.
    """

    @property
    def regions(self) -> List[Region]:
        """
        The regions in the universe.
        """
    
    def insert_token(self, region: Region) -> None:
        """
        Insert a token into the universe.

        :param region: The region to insert.
        """
        new_id = len(self.regions) + 1
        self.region_to_id[region] = new_id
        self.id_to_region[new_id] = region

    def convert_region_to_id(self, region: Region) -> int:
        """
        Convert a region to its corresponding ID.

        :param region: The region to convert.

        :return: The ID of the region, or None if the region is not found.
        """
        return self.region_to_id.get(region)

    def convert_chr_start_end_to_id(self, chrom: str, start: int, end: int) -> int:
        """
        Convert chromosome, start, and end positions to the corresponding ID.

        :param chrom: The chromosome name.
        :param start: The start position.
        :param end: The end position.

        :return: The ID of the region, or None if the region is not found.
        """
        region = Region(chrom, start, end)
        return self.convert_region_to_id(region)

    def convert_id_to_region(self, id: int) -> Region:
        """
        Convert an ID to its corresponding region.

        :param id: The ID to convert.

        :return: The region corresponding to the ID, or None if the ID is not found.
        """
        return self.id_to_region.get(id)

    def __len__(self) -> int:
        """
        Get the number of regions in the universe.

        :return: The number of regions.
        """
        return len(self.regions)

    def is_empty(self) -> bool:
        """
        Check if the universe is empty.

        :return: True if the universe is empty, False otherwise.
        """
        return len(self.regions) == 0

    def __repr__(self) -> str:
        """
        Get a string representation of the universe.

        :return: The string representation.
        """
        return f"Universe with {len(self)} regions"

class Region:
    def __new__(cls, chrom: str, start: int, end: int) -> Region:
        """
        Construct a new Region object.

        :param chrom: The chromosome name.
        :param start: The start position.
        :param end: The end position.
        """
    
    @property
    def chr(self) -> str:
        """
        The chromosome name for this region.
        """

    @property
    def start(self) -> int:
        """
        The start position for this region.
        """

    @property
    def end(self) -> int:
        """
        The end position for this region.
        """

    def __repr__(self) -> str: ...

class TokenizedRegion:
    """
    A TokenizedRegion object represents a tokenized region.
    """

    @property
    def chr(self) -> str:
        """
        The chromosome name for this region.
        """
    
    @property
    def start(self) -> int:
        """
        The start position for this region.
        """
    
    @property
    def end(self) -> int:
        """
        The end position for this region.
        """

    @property
    def id(self) -> int:
        """
        The integer representation of the tokenized region.
        """
    
    @property
    def universe(self) -> Universe:
        """
        The universe object.
        """
    
    def to_region(self) -> Region:
        """
        Convert the tokenized region back to the original region.
        """
    
    def __repr__(self) -> str: ...

class RegionSet:
    def __new__(cls, path: str) -> RegionSet:
        """
        Construct a new RegionSet object.

        :param path: The path to the BED file.
        """
    
    def __repr__(self) -> str: ...

    def __len__(self) -> int: ...

    def __iter__(self) -> Iterator[Region]: ...

    def __next__(self) -> Region: ...

    def __getitem__(self, indx: int) -> Region: ...

class TokenizedRegionSet:
    def __new__(cls, regions: List[Region], tokens: List[int]) -> TokenizedRegionSet:
        """
        Construct a new TokenizedRegionSet object.

        :param regions: The original regions.
        :param tokens: The tokenized regions.
        """
    
    @property
    def ids(self) -> List[int]:
        """
        Integer representation of the tokenized regions.
        """
    
    @property
    def universe(self) -> Universe:
        """
        The universe object.
        """

    def to_bit_vector(self) -> List[int]:
        """
        Convert the tokenized regions to a bit vector.
        """
    
    def to_regions(self) -> List[Region]:
        """
        Convert the tokenized regions back to the original regions.
        """
    
    def to_ids(self) -> List[int]:
        """
        Get the integer representations of the tokenized regions.
        """
    
    def ids_as_strs(self) -> List[str]:
        """
        Get the integer representations of the tokenized regions as strings. This
        is useful for applications that require string representations of the
        tokenized regions.
        """
    
    def __len__(self) -> int: ...

    def __repr__(self) -> str: ...

class TreeTokenizer:
    def __new__(cls, path: str) -> TreeTokenizer:
        """
        Construct a new TreeTokenize from a universe file.

        :param path: The path to the universe file. This should be a BED file.
        """
    
    def unknown_token(self) -> Region:
        """
        Get the unknown token.
        """

    def padding_token(self) -> Region:
        """
        Get the padding token.
        """

    def mask_token(self) -> Region:
        """
        Get the mask token.
        """

    def cls_token(self) -> Region:
        """
        Get the CLS token.
        """

    def bos_token(self) -> Region:
        """
        Get the BOS token.
        """

    def eos_token(self) -> Region:
        """
        Get the EOS token.
        """

    def sep_token(self) -> Region:
        """
        Get the SEP token.
        """
    
    def unknown_token_id(self) -> int:
        """
        Get the ID of the unknown token.
        """

    def padding_token_id(self) -> int:
        """
        Get the ID of the padding token.
        """

    def mask_token_id(self) -> int:
        """
        Get the ID of the mask token.
        """

    def cls_token_id(self) -> int:
        """
        Get the ID of the CLS token.
        """

    def bos_token_id(self) -> int:
        """
        Get the ID of the BOS token.
        """

    def eos_token_id(self) -> int:
        """
        Get the ID of the EOS token.
        """

    def sep_token_id(self) -> int:
        """
        Get the ID of the SEP token.
        """

    def vocab_size(self) -> int:
        """
        Get the vocabulary size.
        """

    def tokenize(self, regions: List[Region]) -> List[Region]:
        """
        Tokenize a list of regions. This will only return the tokenized regions.

        :param regions: The regions to tokenize.

        :return: The tokenized regions as a list.
        """

    def tokenize_bed_file(self, path: str) -> List[Region]:
        """
        Tokenize a BED file directly.

        :param path: The path to the BED file.

        :return: The tokenized regions as a list.
        """

    def encode(self, regions: List[Region]) -> List[int]:
        """
        Encode a list of regions. This will return the integer representation of the tokenized regions.

        :param regions: The regions to encode.

        :return: The integer representation of the tokenized regions.
        """

    def decode(self, ids: List[int]) -> List[Region]:
        """
        Decode a list of integer representations of the tokenized regions.

        :param ids: The integer representations of the tokenized regions.

        :return: The decoded regions.
        """

    def vocab(self) -> List[Tuple[Region, int]]:
        """
        Get the vocabulary.

        :return: The vocabulary as a list of tuples.
        """
    
    @property
    def universe(self) -> Universe:
        """
        The universe object.
        """

    def __call__(self, regions: List[Region]) -> TokenizedRegionSet:
        """
        Tokenize a list of regions.

        :param regions: The regions to tokenize.

        :return: A TokenizedRegionSet object.
        """

    def __len__(self) -> int:
        """
        Get the vocabulary size.
        """

    def __repr__(self) -> str:
        """
        Get a string representation of the tokenizer.
        """

class FragmentTokenizer:
    def __new__(cls, path: str) -> FragmentTokenizer:
        """
        Construct a new FragmentTokenizer from a universe file.

        :param path: The path to the universe file. This should be a BED file.
        """
    
    def tokenize_fragments_to_gtoks(self, file_path: str, out_path: str = None, filter: List[str] = None) -> None:
        """
        Tokenize a file containing fragments.

        :param file_path: The path to the file containing fragments.
        :param out_path: The path to the output file. If None, the output is written to the standard output.
        :param filter: A list of chromosomes to filter. If None, all chromosomes are included.
        """