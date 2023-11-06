from typing import Optional

def prune_universe(universe: str, data: str, min_count: Optional[int], output: Optional[str]) -> None:
    """
    Prune a universe file.

    Pruning a universe file removes all items that appear less than min_count. This is useful for
    removing universe regions that are never seen during training or ultra-rare regions.

    :param str universe: Path to the universe file.
    :param str data: Path to the data (folder of bed files).
    :param int min_count: Minimum number of times a region must appear to be kept.
    :param str output: Path to the output file.
    """