from typing import List

class Interval:
    """
    Represents a range of values.
    """

    def __new__(cls, start: int, end: int) -> Interval:
        """
        Create a new Interval object.

        :param start: The start of the interval.
        :param end: The end of the interval.
        """
    
    @property
    def start(self) -> int:
        """
        The start of the interval.
        """
    
    @property
    def end(self) -> int:
        """
        The end of the interval.
        """
    
    def __repr__(self) -> str: ...

class AIList:
    """
    The augmented interval list (AILIST) object.

    This object will compute region overlaps very efficiently.
    """

    def __new__(cls, intervals: List[Interval], minimum_coverage_length: int = None) -> AIList:
        """
        Create a new AIList object.

        :param intervals: A list of intervals.
        :param minimum_coverage_length: The minimum length of the coverage.
        """
    
    def query(self, interval: Interval) -> List[Interval]:
        """
        Query the AIList object for overlapping intervals.

        :param interval: The interval to query.
        """