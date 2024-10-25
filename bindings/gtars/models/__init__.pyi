from typing import List

class Region:
    chr: str
    start: int
    end: int

class RegionSet:
    regions: List[Region]

class TokenizedRegionPointer:
    id: int
    chrom_id: int
    source_start: int
    source_end: int