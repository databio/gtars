from typing import List, Dict, Optional

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
        Get list of region widths (alias for widths())
        """
        ...

    def widths(self) -> List[int]:
        """
        Get list of region widths (end - start for each region)
        """
        ...

    def neighbor_distances(self) -> List[Optional[float]]:
        """
        Distances between consecutive regions on each chromosome.
        Returns None for missing values.
        """
        ...

    def nearest_neighbors(self) -> List[Optional[float]]:
        """
        Distance from each region to its nearest neighbor.
        Returns None for regions with no neighbor on the same chromosome.
        """
        ...

    def distribution(self, n_bins: int = 250) -> List[Dict[str, object]]:
        """
        Region distribution across genomic bins.

        :param n_bins: number of bins (default 250)
        :return: list of dicts with keys: chr, start, end, n, rid
        """
        ...

    def trim(self, chrom_sizes: Dict[str, int]) -> "RegionSet":
        """
        Clamp regions to chromosome boundaries.
        Regions on unknown chromosomes are dropped.
        """
        ...

    def promoters(self, upstream: int, downstream: int) -> "RegionSet":
        """
        Generate promoter regions relative to each region's start.
        """
        ...

    def reduce(self) -> "RegionSet":
        """
        Merge overlapping and adjacent intervals per chromosome.
        """
        ...

    def setdiff(self, other: "RegionSet") -> "RegionSet":
        """
        Set difference: remove portions overlapping with other.
        """
        ...

    def pintersect(self, other: "RegionSet") -> "RegionSet":
        """
        Pairwise intersection by index position.
        """
        ...

    def concat(self, other: "RegionSet") -> "RegionSet":
        """
        Combine two region sets without merging.
        """
        ...

    def union(self, other: "RegionSet") -> "RegionSet":
        """
        Merge two region sets into minimal non-overlapping set.
        """
        ...

    def jaccard(self, other: "RegionSet") -> float:
        """
        Nucleotide-level Jaccard similarity (|intersection| / |union|).
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

    def subset_by_overlaps(self, other: "RegionSet") -> "RegionSet":
        """
        Return a new RegionSet containing only regions from self that
        overlap at least one region in other.

        Builds an AIList index from other and queries each region in self.
        """
        ...

    def count_overlaps(self, other: "RegionSet") -> List[int]:
        """
        Count how many regions in other overlap each region in self.

        Returns a list of integers with one entry per region.
        """
        ...

    def any_overlaps(self, other: "RegionSet") -> List[bool]:
        """
        Check whether each region in self overlaps any region in other.

        Returns a list of booleans with one entry per region.
        """
        ...

    def find_overlaps(self, other: "RegionSet") -> List[List[int]]:
        """
        Find indices into other that overlap each region in self.

        Returns a list of lists, where each inner list contains the
        0-based indices of regions in other that overlap the
        corresponding region in self.
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

    @staticmethod
    def from_regionset(rs: RegionSet) -> "TssIndex":
        """
        Create a TssIndex from a RegionSet.
        """
        ...

    def calc_tss_distances(self, rs: RegionSet) -> List[int]:
        """
        Unsigned distance from each query region to nearest TSS.
        """
        ...

    def feature_distances(self, rs: RegionSet) -> List[Optional[float]]:
        """
        Signed distance from each query region to nearest feature.
        Returns None for regions on chromosomes with no features.
        """
        ...

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...


class GeneModel:
    @staticmethod
    def from_gtf(
        path: str,
        filter_protein_coding: bool = True,
        convert_ensembl_ucsc: bool = True,
    ) -> "GeneModel":
        """
        Load a gene model from a GTF file.
        """
        ...

    @property
    def n_genes(self) -> int: ...

    @property
    def n_exons(self) -> int: ...

    def __repr__(self) -> str: ...


class PartitionList:
    @staticmethod
    def from_gene_model(
        gene_model: GeneModel,
        core_prom: int,
        prox_prom: int,
        chrom_sizes: Optional[Dict[str, int]] = None,
    ) -> "PartitionList":
        """
        Build partitions from a gene model.
        """
        ...

    @staticmethod
    def from_gtf(
        path: str,
        core_prom: int,
        prox_prom: int,
        filter_protein_coding: bool = True,
        convert_ensembl_ucsc: bool = True,
        chrom_sizes: Optional[Dict[str, int]] = None,
    ) -> "PartitionList":
        """
        Build partitions directly from a GTF file.
        """
        ...

    def partition_names(self) -> List[str]:
        """
        List of partition category names.
        """
        ...

    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...


class SignalMatrix:
    @staticmethod
    def load_bin(path: str) -> "SignalMatrix":
        """
        Load a SignalMatrix from a packed binary file.
        """
        ...

    @staticmethod
    def from_tsv(path: str) -> "SignalMatrix":
        """
        Load a SignalMatrix from a TSV file.
        """
        ...

    @property
    def condition_names(self) -> List[str]: ...

    @property
    def n_conditions(self) -> int: ...

    @property
    def n_regions(self) -> int: ...

    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...


class GenomicDistAnnotation:
    @staticmethod
    def load_bin(path: str) -> "GenomicDistAnnotation":
        """
        Load a GenomicDistAnnotation from a GDA binary file.
        """
        ...

    @staticmethod
    def from_gtf(
        path: str,
        filter_protein_coding: bool = True,
        convert_ensembl_ucsc: bool = True,
    ) -> "GenomicDistAnnotation":
        """
        Create a GenomicDistAnnotation from a GTF file.
        """
        ...

    def gene_model(self) -> GeneModel:
        """
        Extract the gene model.
        """
        ...

    def partition_list(
        self,
        core_prom: int,
        prox_prom: int,
        chrom_sizes: Optional[Dict[str, int]] = None,
    ) -> PartitionList:
        """
        Build a PartitionList from the gene model.
        """
        ...

    def tss_index(self) -> TssIndex:
        """
        Derive a strand-aware TssIndex from the gene model.
        """
        ...

    def __repr__(self) -> str: ...