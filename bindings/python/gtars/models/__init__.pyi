from typing import List

class Region:
    chr: str
    start: int
    end: int

    rest: str

class RegionSet:
    regions: List[Region]
    header: str
    path: str

    def __init__(path: str):
        """
        :param path: path to the bed file
        """
        ...
