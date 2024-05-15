def prune_universe(data: str, universe: str, min_count: int = None, output: int = None):
    """
    Prune the universe based on the data. This is useful for reducing the total number of tokens in the universe.

    :param data: The data to use for pruning. It is a path to a directory containing the data.
    :param universe: The universe to prune. It should be a path to a BED file.
    :param min_count: The minimum number of times a token should appear in the data to be included in the pruned universe.
    :param output: The output file. If not provided, the pruned universe will be saved in the same directory as the input universe.
    """