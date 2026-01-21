from typing import List, Dict

class Region:
    chr: str
    start: int
    end: int

    rest: str

    def __repr__(self): ...
    def __str__(self): ...
    def __len__(self) -> int:
        """
        Size of the region
        """
        ...

class ChromosomeStatistics:
    chromosome: str
    number_of_regions: int
    minimum_region_length: int
    maximum_region_length: int
    mean_region_length: float
    median_region_length: float
    start_nucleotide_position: int
    end_nucleotide_position: int

class RegionSet:
    regions: List[Region]
    header: str
    path: str

    def __init__(self, path: str) -> "RegionSet":
        """
        :param path: path to the bed file
        """
        ...

    @classmethod
    def from_regions(cls, regions: List[Region]) -> "RegionSet":
        """
        :param regions: list of regions
        """
        ...

    @property
    def identifier(self) -> str:
        """
        Get RegionSet identifier

        :return: bed digest/identifier
        """
        ...

    @property
    def path(self) -> str:
        """
        Get RegionSet identifier

        :return: bed file path
        """
        ...

    @property
    def header(self) -> str:
        """
        Header of the bed file
        """
        ...

    @property
    def file_digest(self) -> str:
        """
        Digest of whole bed file (with all columns)
        """
        ...

    def to_bed(self, path: str) -> None:
        """
        Save RegionSet as bed file

        :param path: path to the bed file that has to be saved
        """
        ...

    def to_bed_gz(self, path: str) -> None:
        """
        Save RegionSet as bed.gz file

        :param path: path to the bed file that has to be saved
        """
        ...

    def to_bigbed(self, out_path: str, chrom_size: str) -> None:
        """
        Save RegionSet as bigBed file

        :param out_path: path bigbed file path
        :param chrom_size: path to the chrom sizes file
        """
        ...

    def sort(self) -> None:
        """
        Sort the regions
        """
        ...

    def file_digest(self) -> str:
        """
        Get full file digest (all columns). It differs from identifier which is calculated only from first three columns.
        """
        ...

    def region_widths(self) -> List[int]:
        """
        Get list of region widths
        """
        ...

    def mean_region_width(self) -> int:
        """
        Mean width of the regions
        """
        ...

    def get_max_end_per_chr(self) -> Dict[str, int]:
        """
        Get Max end coordinate of nucleotide for each chromosome
        """
        ...

    def get_nucleotide_length(self) -> int:
        """
        Get total number of nucleotides in RegionSet
        """
        ...

    def chromosome_statistics(self) -> Dict[str, ChromosomeStatistics]:
        """
        Get a dictionary of ChromosomeStatistics for each chromosome in the RegionSet
        """
        ...

    def __len__(self) -> int:
        """
        Size of the regionset
        """
        ...

    def __iter__(self): ...
    def __getitem__(self, indx: int): ...
    def __repr__(self): ...
    def __str__(self): ...

class GenomeAssembly:
    def __init__(self, path: str) -> "GenomeAssembly":
        """
        :param path: path to the fasta file
        """
        ...

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...


class TssIndex:
    def __init__(self, path: str) -> "TssIndex":
        """
        :param path: path to the bed file with tss regions
        """
        ...

    def calc_tss_distances(self, rs: RegionSet) -> List[int]:
        """
        :param rs: RegionSet

        :return : list of region distances to the nearest TSS
        """
        ...

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...