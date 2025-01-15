from typing import List

class PyCoverage:
    def __init__(self, path: str):
        """
        Initialize a PyCoverage object.

        :param path: The path to the BigWig file.
        """
        self.bw_file = path

    def stats(self, chr: str, start: int, end: int, stat_type: str) -> float:
        """
        Get statistics from the BigWig file.

        :param chr: The chromosome name.
        :param start: The start position.
        :param end: The end position.
        :param stat_type: The type of statistic to compute.
        :return: The computed statistic.
        """
        pass

    def __repr__(self) -> str:
        """
        Return a string representation of the PyCoverage object.
        """
        return f"Coverage({self.bw_file})"

    def __str__(self) -> str:
        """
        Return a string representation of the PyCoverage object.
        """
        return f"Coverage({self.bw_file})"

def write_tokens_to_gtok(filename: str, tokens: List[str]):
    """
    Write tokens to a GTOK file.

    :param filename: The filename of the GTOK file.
    :param tokens: The tokens to write.
    """
    pass

def read_tokens_from_gtok(filename: str) -> List[str]:
    """
    Read tokens from a GTOK file.

    :param filename: The filename of the GTOK file.
    """
    pass

def read_tokens_from_gtok_as_strings(filename: str) -> List[str]:
    """
    Read tokens from a GTOK file and return them as strings.

    :param filename: The filename of the GTOK file.
    """
    pass

def read_big_wig_to_coverage(filename: str) -> "PyCoverage":
    """
    Read a BigWig file and return a PyCoverage object.

    :param filename: The filename of the BigWig file.
    """
    pass