from typing import Optional

def prune_universe(data: str, universe: str, min_count: Optional[int], output: Optional[str]) -> None:
    """
    Given a universe and a path to its corresponding data, "prune" the universe to remove regions that
    are not present in the data. This is useful for reducing the size of the universe to remove unused
    regions. It also can be used to remove regions that fall below a certain frequency threshold.

    :param str data: The path to the data file.
    :param str universe: The path to the universe file.
    :param int min_count: The minimum number of times a region must appear in the data to be kept.
    :param str output: The path to the output file.
    """