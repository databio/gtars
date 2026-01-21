from typing import Optional, List, Dict
from gtars.models import RegionSet, GenomeAssembly

def calc_gc_content(
    rs: RegionSet, genome: GenomeAssembly, ignore_unk_chroms: Optional[bool]
) -> List[float]:
    """
    Calculate GC content for RegionSet (bed file)

    :param rs: RegionSet object
    :param genome: GenomeAssembly object (reference genome)
    :param ignore_unk_chroms: ignore unknown chromosome for provided genome assembly

    :return: List[float]
    """
    ...

def calc_dinucleotide_frequency(
    rs: RegionSet, genome: GenomeAssembly
) -> Dict[str, int]:
    """
    Calculate Dinucleotide frequencies for RegionSet (bed file)

    :param rs: RegionSet object
    :param genome: GenomeAssembly object (reference genome)
    :return: Dict[str, int]
    """
    ...
