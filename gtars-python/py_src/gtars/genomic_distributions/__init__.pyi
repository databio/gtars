from typing import Optional, Union, List, Dict
from gtars.models import RegionSet, GenomeAssembly, BinaryGenomeAssembly, PartitionList, SignalMatrix

def calc_gc_content(
    rs: RegionSet,
    genome: Union[GenomeAssembly, BinaryGenomeAssembly],
    ignore_unk_chroms: Optional[bool] = False,
) -> List[float]:
    """
    Calculate GC content for each region.

    :param rs: RegionSet object
    :param genome: GenomeAssembly or BinaryGenomeAssembly (reference genome)
    :param ignore_unk_chroms: ignore unknown chromosomes
    :return: GC content per region (0 to 1)
    """
    ...

def calc_dinucl_freq(
    rs: RegionSet,
    genome: Union[GenomeAssembly, BinaryGenomeAssembly],
    raw_counts: bool = False,
    ignore_unk_chroms: bool = False,
) -> Dict[str, object]:
    """
    Per-region dinucleotide frequencies, matching R GenomicDistributions ``calcDinuclFreq``.

    :param rs: RegionSet object
    :param genome: GenomeAssembly or BinaryGenomeAssembly (reference genome)
    :param raw_counts: if False (default), return percentages (0–100) per row;
        if True, return raw integer-valued counts
    :param ignore_unk_chroms: skip regions on chromosomes not in the assembly
    :return: dict with keys:
        - ``region_labels``: list of ``chr_start_end`` strings (one per region)
        - ``dinucleotides``: list of 16 dinucleotide names in canonical order
        - ``frequencies``: list of lists — one row per region, 16 values per row
          matching ``dinucleotides`` order

    For pooled global counts across all regions, sum each column of the
    raw-counts matrix.
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
