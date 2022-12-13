"""
A class describing one type of variant for NEAT.
"""

import numpy as np

from Bio.Seq import Seq
from .base_variant import BaseVariant


class Deletion(BaseVariant):
    """
        A deletion is when any number of consecutive bases are absent in the reads, compared to the reference. This can
        happen for several biological reasons. The Mutation model sets the parameters and boundaries for these variants.
        We will use the convention that the length is a negative value for a deletion.

        :param position1: Position of the variant. This includes the base before the first base actually deleted,
            similar to VCF notation.
        :param length: The length of the deletion.
        :param genotype: the list of ploids with a 1 indicating where this variant is found
        :param qual_score: the quality score for this variant, as read or generated
        :param is_input: True if from an input vcf, in which case this variant will get priority.
        :param kwargs: can be used to store data from input vars or unused variables from the base class.
        """

    def __init__(self,
                 position1: int,
                 length: int,
                 genotype: np.ndarray = None,
                 qual_score: int = None,
                 is_input: bool = False,
                 **kwargs):
        self.position1 = position1
        self.length = length
        self.genotype = genotype
        self.qual_score = qual_score
        self.is_input = is_input
        self.metadata = kwargs

    def contains(self, other) -> bool:
        """
        Checks whether the given position1 is within the variant.

        :param other: The given position, relative to the reference
        :return: True or False
        """
        if np.array_equal(other.genotype, self.genotype):
            return self.position1 <= other.position1 < self.position1 + self.length
        else:
            return False

    def __lt__(self, position: int) -> bool:
        return position < self.position1

    def __gt__(self, position: int) -> bool:
        return position > self.position1 + self.length

    def __le__(self, position: int) -> bool:
        return position <= self.position1

    def __ge__(self, position: int) -> bool:
        return position >= self.position1 + self.length

    def __eq__(self, other) -> bool:
        if other.type == self.__class__:
            return other.position == self.position1 and other.length == self.length
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.position1}, {self.length})'
