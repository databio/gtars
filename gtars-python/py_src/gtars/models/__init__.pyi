from typing import List, Dict, Optional, Tuple

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

class SpacingStats:
    """
    Summary statistics over the distribution of inter-region spacings.

    Returned by ``RegionSet.inter_peak_spacing()``. All float fields are
    NaN when ``n_gaps`` is 0 (empty, singleton-per-chromosome, or
    all-overlapping inputs).
    """
    n_gaps: int
    mean: float
    median: float
    std: float
    iqr: float
    log_mean: float
    log_std: float

    def __repr__(self) -> str: ...

class ClusterStats:
    """
    Cluster-level summary statistics at a given stitching radius.

    Returned by ``RegionSet.peak_clusters(radius_bp, min_cluster_size)``.

    ``min_cluster_size`` applies uniformly to every size-dependent field
    except ``max_cluster_size`` (which is always the biggest cluster in
    the input regardless of filter). The identity
    ``n_clusters * mean_cluster_size == n_clustered_peaks`` holds at any
    threshold.

    **Default** ``min_cluster_size = 2``: fields describe "clusters of at
    least 2 peaks". ``n_clusters`` excludes singletons,
    ``n_clustered_peaks`` is peaks with at least one neighbor within
    ``radius_bp``, ``mean_cluster_size`` is their average size, and
    ``fraction_clustered`` is ``n_clustered_peaks / total_peaks``. This
    matches the typical "how clustered are my peaks?" use case.

    **Pass** ``min_cluster_size=1`` to get the simple-average view —
    ``mean_cluster_size`` becomes ``total_peaks / n_clusters``, but
    ``n_clustered_peaks`` degenerates to ``total_peaks`` and
    ``fraction_clustered`` to ``1.0`` (tautologies under this threshold).

    ``mean_cluster_size`` and ``fraction_clustered`` are NaN for empty
    inputs. ``mean_cluster_size`` is also NaN when no clusters meet the
    threshold (e.g. all singletons with ``min_cluster_size=2``).
    """
    radius_bp: int
    n_clusters: int
    n_clustered_peaks: int
    mean_cluster_size: float
    max_cluster_size: int
    fraction_clustered: float

    def __repr__(self) -> str: ...

class DensityVector:
    """
    Dense per-window peak count vector and its chromosome binning.

    Returned by ``RegionSet.density_vector(chrom_sizes, n_bins)``. Unlike
    ``distribution()`` (which returns only non-empty bins), this carries
    the full zero-padded vector, ordered by karyotypic chromosome order
    and then bin index.

    Notes
    -----
    ``n_bins`` is the **target bin count for the longest chromosome in
    chrom_sizes**, not the total length of ``counts``. Bin width is
    derived as ``max(chrom_sizes.values()) // n_bins`` (minimum 1 bp),
    and every chromosome is tiled with windows of that width. The length
    of ``counts`` is ``sum(ceil(chrom_size / bin_width))`` across
    chromosomes in ``chrom_sizes``, which can substantially exceed
    ``n_bins`` when many chromosomes are present. To target a specific
    bin width in bp, pass ``n_bins = max_chrom_len // desired_width``.

    The last bin on each chromosome is narrower than ``bin_width``
    whenever ``chrom_size`` is not an exact multiple of ``bin_width``.
    Chromosomes shorter than ``bin_width`` (common with UCSC alt /
    random / unplaced contigs like ``chrUn_*``, ``*_random``, ``*_alt``)
    reduce to a single bin whose effective width equals the chromosome
    size rather than ``bin_width``. Entries in ``counts`` are therefore
    counts per bin, not counts per ``bin_width`` bp — bins of different
    effective widths are not directly comparable as densities when
    ``chrom_sizes`` contains contigs significantly shorter than
    ``bin_width``.
    """
    n_bins: int
    bin_width: int
    counts: List[int]
    chrom_offsets: List[Tuple[str, int]]

    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...

class DensityHomogeneity:
    """
    Summary of how evenly peaks are distributed across genome windows.

    Returned by ``RegionSet.density_homogeneity(chrom_sizes, n_bins)``.
    Poisson-distributed peaks give cv ≈ 1; clustered sets give cv >> 1;
    evenly-spread sets give cv << 1.

    Notes
    -----
    See ``DensityVector`` for the definition of ``n_bins`` (target bin
    count for the longest chromosome, not the total window count) and
    for the treatment of chromosomes shorter than the derived bin width.
    Short contigs each contribute a narrow single-bin entry which
    dilutes ``mean_count``, inflates ``n_windows``, and raises ``gini``.

    The Gini coefficient is biased high for very sparse count
    distributions (many zero-count windows). Check ``n_nonzero_windows``
    before interpreting Gini on sparse peak sets.
    """
    bin_width: int
    n_windows: int
    n_nonzero_windows: int
    mean_count: float
    variance: float
    cv: float
    gini: float

    def __repr__(self) -> str: ...

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

    def neighbor_distances(self) -> List[int]:
        """
        Signed gaps between consecutive regions on each chromosome.

        Output length may be shorter than input region count — chromosomes with
        fewer than 2 regions are skipped (no neighbors to measure against).
        Output is NOT aligned 1:1 with input regions.
        """
        ...

    def nearest_neighbors(self) -> List[int]:
        """
        Distance from each region to its nearest neighbor on the same chromosome.

        Output length may be shorter than input region count — chromosomes with
        only one region are skipped. Output is NOT aligned 1:1 with input regions.
        """
        ...

    def distribution(
        self,
        n_bins: int = 250,
        chrom_sizes: Optional[Dict[str, int]] = None,
    ) -> List[Dict[str, object]]:
        """
        Region distribution across genomic bins.

        :param n_bins: number of bins for the longest chromosome (default 250)
        :param chrom_sizes: optional mapping of chromosome name to length. When provided,
            per-chromosome bin sizes are derived from the reference genome
            (bin_size = chrom_size / n_bins per chrom). This produces outputs that are
            comparable across BED files and aligned with reference genome positions.
            When absent, bin size is derived from the BED file's observed max end
            coordinate — outputs will NOT be comparable across files.

            When ``chrom_sizes`` is provided, regions on chromosomes not listed in
            ``chrom_sizes`` are skipped, as are regions whose midpoint falls beyond
            the stated chromosome size (common with assembly mismatches). The summed
            bin counts may therefore be lower than the input region count.
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

    def coverage(self, other: "RegionSet") -> float:
        """
        Fraction of self's base pairs covered by other (after merging overlaps).
        Returns a value in [0.0, 1.0].
        """
        ...

    def overlap_coefficient(self, other: "RegionSet") -> float:
        """
        Overlap coefficient: |intersection_bp| / min(|self_bp|, |other_bp|).
        Returns a value in [0.0, 1.0].
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

    def intersect_all(self, other: "RegionSet") -> "RegionSet":
        """
        All-vs-all genomic intersection.

        For each pair of overlapping regions between self and other,
        compute intersection coordinates [max(a.start, b.start), min(a.end, b.end)].
        Returns a RegionSet of all intersection fragments.

        Unlike pintersect (which pairs by index position), this finds ALL
        overlapping pairs across the two sets.
        """
        ...

    def subtract(self, other: "RegionSet") -> "RegionSet":
        """
        Genomic subtraction (alias for setdiff).

        Remove portions of self that overlap with other.
        Both inputs are reduced before subtraction.
        """
        ...

    def closest(self, other: "RegionSet") -> List[tuple]:
        """
        Find nearest region in other for each region in self.

        Returns list of (self_index, other_index, distance) tuples.
        Distance is 0 for overlapping regions.
        Regions on chromosomes absent in other are omitted.
        """
        ...

    def cluster(self, max_gap: int = 0) -> List[int]:
        """
        Cluster nearby regions.

        Assign a cluster ID to each region. Regions within max_gap
        distance on the same chromosome are assigned the same cluster.
        Returns cluster IDs in original region order.
        """
        ...

    def chromosome_statistics(self) -> Dict[str, ChromosomeStatistics]:
        """
        Get a dictionary of ChromosomeStatistics for each chromosome in the RegionSet
        """
        ...

    def gaps(self, chrom_sizes: Dict[str, int]) -> "RegionSet":
        """
        Return gaps between regions per chromosome, bounded by chrom sizes.

        Emits leading gaps (from 0), inter-region gaps, trailing gaps
        (to ``chrom_size``), and full-chromosome gaps for any chromosome
        in ``chrom_sizes`` that has no regions. Regions on chromosomes
        not listed in ``chrom_sizes`` are skipped.
        """
        ...

    def inter_peak_spacing(self) -> SpacingStats:
        """
        Summary statistics over the distribution of inter-region spacings.

        Wraps ``neighbor_distances()`` and reduces to scalar mean / median /
        std / IQR / log-mean / log-std. Useful as a descriptive scalar
        summary of how regularly peaks are spaced.
        """
        ...

    def peak_clusters(
        self,
        radius_bp: int,
        min_cluster_size: int = 2,
    ) -> ClusterStats:
        """
        Cluster-level summary statistics at a given stitching radius.

        Wraps ``cluster(radius_bp)``: groups regions into single-linkage
        connected components where two regions link if the bp gap between
        ``prev.end`` and ``next.start`` is at most ``radius_bp``.
        Chromosome-scoped (no cross-chromosome linking).

        :param radius_bp: stitching radius in bp.
        :param min_cluster_size: minimum cluster size that qualifies
            clusters for inclusion in the size-dependent fields of the
            returned ``ClusterStats`` (``n_clusters``,
            ``n_clustered_peaks``, ``mean_cluster_size``, and
            ``fraction_clustered``). ``max_cluster_size`` is always
            unfiltered.

            **Default 2**: every field describes "clusters with at least
            2 peaks". Matches typical enhancer-clustering or
            super-enhancer-stitching analyses and gives
            ``fraction_clustered`` a useful meaning (fraction of peaks
            with at least one neighbor within ``radius_bp``).

            **Pass 1** to include singletons and get the simple-average
            view — ``mean_cluster_size`` becomes ``total_peaks /
            n_clusters`` but ``n_clustered_peaks`` degenerates to
            ``total_peaks`` and ``fraction_clustered`` to ``1.0``.

            The identity ``n_clusters * mean_cluster_size ==
            n_clustered_peaks`` holds at any threshold.
        """
        ...

    def density_vector(
        self,
        chrom_sizes: Dict[str, int],
        n_bins: int,
    ) -> DensityVector:
        """
        Dense zero-padded per-window peak count vector.

        Unlike ``distribution()`` (which returns only non-empty bins), this
        returns the full zero-padded count vector with one entry per window
        on every chromosome in ``chrom_sizes``, ordered by karyotypic
        chromosome order and bin index. Suitable for ML feature extraction.

        :param chrom_sizes: mapping of chromosome name to length in bp.
        :param n_bins: **target bin count for the longest chromosome in
            chrom_sizes** — not the length of the returned ``counts``.
            Bin width is derived as ``max(chrom_sizes.values()) //
            n_bins`` (minimum 1 bp); shorter chromosomes get
            proportionally fewer bins. Total window count is
            ``sum(ceil(size / bin_width))`` across ``chrom_sizes`` and
            can substantially exceed ``n_bins``.
        :return: a ``DensityVector`` — see its docstring for per-chromosome
            bin-width behavior and the treatment of short contigs (the
            last bin on each chromosome and all bins on chromosomes
            shorter than ``bin_width`` are narrower than ``bin_width``,
            so counts are per bin, not per ``bin_width`` bp).
        """
        ...

    def density_homogeneity(
        self,
        chrom_sizes: Dict[str, int],
        n_bins: int,
    ) -> DensityHomogeneity:
        """
        Summary statistics over the dense per-window count vector.

        Returns mean, population variance, coefficient of variation, Gini
        coefficient, and the count of nonzero windows. Useful as a scalar
        measure of peak spatial evenness.

        :param chrom_sizes: mapping of chromosome name to length in bp.
        :param n_bins: **target bin count for the longest chromosome**,
            not the total window count — see
            :py:meth:`RegionSet.density_vector` for the full semantic.
            Short contigs in ``chrom_sizes`` each contribute a narrow
            single-bin entry which dilutes ``mean_count``, inflates
            ``n_windows``, and raises ``gini``.
        :return: a ``DensityHomogeneity``. Note: Gini is biased high for
            very sparse count distributions; consult
            ``n_nonzero_windows`` before interpreting Gini on sparse
            peak sets.
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

class BinaryGenomeAssembly:
    def __init__(self, path: str) -> "BinaryGenomeAssembly":
        """
        Open a .fab binary FASTA file for zero-copy mmap access.

        Create .fab files with ``gtars prep --fasta <file>``.

        :param path: path to the .fab file
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
