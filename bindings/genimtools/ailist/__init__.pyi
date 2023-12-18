from typing import List

class Interval:
    start: int
    end: int
    def __init__(self, start: int, end: int):
        """
        Create a new Interval.
        """

class AIList:
    def __init__(self, intervals: List[Interval], minimum_coverage_length: int) -> AIList:
        """
        Create a new AIList.

        :param intervals: The list of intervals to use.
        :param minimum_coverage_length: The minimum number of intervals that are covered before being decomposed into another list.
        """

    def query(self, interval: Interval) -> List[Interval]:
        """
        Query the AIList for the intervals that overlap with the given interval.

        :param interval: The interval to query.
        """
