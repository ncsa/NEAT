"""
A class describing one type of variant for NEAT.
"""
import numpy as np

from .base_variant import BaseVariant


class Duplication(BaseVariant):
    """
    A duplication is when a segment of the chromosome repeats itself for N number of bases. For this variant type,
    the insertion point for the duplication is at the end of the original segment.
        position1 = The position of the first base of the segment being duplicated
        length = The length of the segment being duplicated
        position2 = position1 + length, i.e., the position where the duplication starts.

    :param position1: The position of the first base of the segment being duplicated
    :param length: The length of the segment being duplicated
    :param position2: The position where the duplication starts. position1 + length or
        position1 - length, for reverse strand duplication.
    :param is_input: True if from an input vcf, in which case this variant will get priority.
    :param kwargs: can be used to store data from input vars or unused variables from the base class.
    """

    def __init__(self,
                 position1: int,
                 length: int,
                 position2: int,
                 genotype: np.ndarray = None,
                 qual_score: int = None,
                 is_input: bool = False,
                 **kwargs):

        self.position1 = position1
        self.length = length
        self.position2 = position2
        # Note that 'stream' is essentially orientation, but that is reserved for official designations.
        self.stream = 1
        if self.position2 == position1 - length:
            self.stream = -1
        elif not self.position2 == position1 + length:
            raise ValueError("Second position of a dup must = first position - length or first position + length")
        self.genotype = genotype
        self.qual_score = qual_score
        self.is_input = is_input
        self.metadata = kwargs

    def contains(self, position: int) -> bool:
        """
        Checks whether the given position is within the variant.

        :param position: The given position, relative to the reference
        :return: True or False
        """
        return self.position1 <= position < self.position1 + self.length

    def __lt__(self, position: int) -> bool:
        if self.stream == 1:
            return position < self.position2
        else:
            return position < self.position2 - self.length

    def __gt__(self, position: int) -> bool:
        if self.stream == 1:
            return position > self.position2 + self.length
        else:
            return position > self.position2

    def __le__(self, position: int) -> bool:
        if self.stream == 1:
            return position <= self.position2
        else:
            return position <= self.position2 - self.length

    def __ge__(self, position: int) -> bool:
        if self.stream == 1:
            return position >= self.position2 + self.length
        else:
            return position >= self.position2

    def __eq__(self, other) -> bool:
        if other.type == self.__class__:
            return other.position1 == self.position1 and other.length == self.length
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.position2}, {self.length})'
