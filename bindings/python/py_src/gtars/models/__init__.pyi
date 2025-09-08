from typing import List, Dict, Union

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
    def __len__(self) -> int:
        """
        Size of the regionset
        """
        ...
    def __iter__(self): ...
    def __getitem__(self, indx: int): ...
    def __repr__(self): ...
    def __str__(self): ...

class ReferenceGenomeMetadata:
    def __init__(self, genome: str, digest: str, description: str, collection: Dict[str, int]) -> "ReferenceGenomeMetadata":
        """
        Initialize ReferenceGenomeMetadata

        :param genome: Genome name of the reference
        :param digest: Digest of the reference
        :param description: Description of the reference
        :collection: dict of chr name: length
        """
        ...

    @classmethod
    def from_path(path: str, name: str) -> "ReferenceGenomeMetadata":
        """
        Initialize ReferenceGenomeMetadata from chrom_sizes file

        :param path: Path to the fasta file
        :param name: Genome/digest. It will be saved to digest/name/description
        """
        ...

    @property
    def name() -> str:
        """
        Get name of the refernce genome
        """
        ...


    @property
    def digets() -> str:
        """
        Get digest of the refernce genome
        """
        ...

    @property
    def description() -> str:
        """
        Get description of the refernce genome
        """
        ...

    @property
    def collection() -> Dict[str, int]:
        """
        Get chrom with leght of the refernce genome
        """
        ...

class CompatibilityConcise:
    @property
    def xs() -> float:
        """
        xs
        """
        ...

    @property
    def oobr() -> Union[float, None]:
        """
        oobr
        """
        ...
    @property
    def sequence_fit() -> Union[float, None]:
        """
        sequence fit
        """
        ...

    @property
    def description() -> str:
        """
        Get description of the refernce genome
        """
        ...

    @property
    def assigned_points() -> int:
        """
        assigned_points
        """
        ...

        @property
    def tier_ranking() -> int:
        """
        tier_ranking
        """
        ...


class ReferenceValidator:
    def __init__(self, ref_genome: List[ReferenceGenomeMetadata]) -> "ReferenceValidator":
        """
        Initialize Reference Validator by providing list of Reference genomes metadata objects

        :param ref_genome: list of refGenomeMetadata objects
        """
        ...

    @classmethod
    def from_path() -> "ReferenceValidator":
        ...

    def determine_compatibility(self, rs: RegionSet) -> Dict[str,CompatibilityConcise]:
        """
        Determine compatibility of the bed file

        :param rs: RegionSet object
        """
        ...
