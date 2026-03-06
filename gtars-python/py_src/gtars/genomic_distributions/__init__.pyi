from typing import Optional, List, Dict
from gtars.models import RegionSet, GenomeAssembly, PartitionList, SignalMatrix

def calc_gc_content(
    rs: RegionSet, genome: GenomeAssembly, ignore_unk_chroms: Optional[bool] = False
) -> List[float]:
    """
    Calculate GC content for each region.

    :param rs: RegionSet object
    :param genome: GenomeAssembly object (reference genome)
    :param ignore_unk_chroms: ignore unknown chromosomes
    :return: GC content per region (0 to 1)
    """
    ...

def calc_dinucleotide_frequency(
    rs: RegionSet, genome: GenomeAssembly
) -> Dict[str, int]:
    """
    Calculate dinucleotide frequencies across all regions.

    :param rs: RegionSet object
    :param genome: GenomeAssembly object (reference genome)
    :return: dict mapping dinucleotide names to counts
    """
    ...

def calc_partitions(
    rs: RegionSet, partition_list: PartitionList, bp_proportion: bool = False
) -> Dict[str, object]:
    """
    Classify regions into genomic partitions.

    :param rs: query RegionSet
    :param partition_list: PartitionList to classify against
    :param bp_proportion: if True, count base pairs; if False, count regions
    :return: dict with keys: partition (list), count (list), total (int)
    """
    ...

def calc_expected_partitions(
    rs: RegionSet,
    partition_list: PartitionList,
    chrom_sizes: Dict[str, int],
    bp_proportion: bool = False,
) -> Dict[str, object]:
    """
    Observed vs expected partition enrichment with chi-square test.

    :param rs: query RegionSet
    :param partition_list: PartitionList to classify against
    :param chrom_sizes: chromosome sizes for expected calculation
    :param bp_proportion: if True, count base pairs; if False, count regions
    :return: dict with keys: partition, observed, expected, log10OE, pvalue
    """
    ...

def calc_summary_signal(
    rs: RegionSet, signal_matrix: SignalMatrix
) -> Dict[str, object]:
    """
    Overlap query regions with a signal matrix and compute summary statistics.

    :param rs: query RegionSet
    :param signal_matrix: SignalMatrix to overlap with
    :return: dict with keys: condition_names, region_labels, signal_matrix, matrix_stats
    """
    ...

def median_abs_distance(distances: List[float]) -> Optional[float]:
    """
    Median of absolute values, filtering sentinel values.

    :param distances: list of distances (NaN/inf treated as missing)
    :return: median absolute distance, or None if all missing
    """
    ...

def consensus(region_sets: List[RegionSet]) -> List[Dict[str, object]]:
    """
    Compute consensus regions from multiple region sets.

    :param region_sets: list of RegionSet objects
    :return: list of dicts with keys: chr, start, end, count
    """
    ...
