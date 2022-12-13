"""
A class describing one type of variant for NEAT.
"""
import numpy as np

from Bio.Seq import Seq
from .base_variant import BaseVariant


class Inversion(BaseVariant):
    """
    An inversion is a type of variant where a segment of the chromosome is reversed from its orientation
    in the reference.

    :param position1: start position of the variant. Is also the start position of the original segment.
    :param length: Length of the inversion
    :param orientation: For inversion, should = -1, be definition.
    :param is_input: True if from an input vcf, in which case this variant will get priority.
    :param kwargs: can be used to store data from input vars or unused variables from the base class.
    """
    def __init__(self,
                 position1: int,
                 length: int,
                 orientation: int,
                 genotype: np.ndarray = None,
                 qual_score: int = None,
                 is_input: bool = False,
                 **kwargs):

        self.position1 = position1
        self.length = length
        self.orientation = orientation
        if not self.orientation == -1:
            raise ValueError("For an inversion, the orientation must be -1")
        self.genotype = genotype
        self.qual_score = qual_score
        self.is_input = is_input
        self.metadata = kwargs

    def contains(self, position: int) -> bool:
        """
        Checks whether the given position1 is within the variant.

        :param position: The given position1, relative to the reference
        :return: True or False
        """
        return self.position1 <= position < self.position1 + self.length

    def __lt__(self, position):
        return position < self.position1

    def __gt__(self, position):
        return position > self.position1 + self.length

    def __le__(self, position):
        return position <= self.position1

    def __ge__(self, position):
        return position >= self.position1 + self.length

    def __eq__(self, other):
        if other.type == self.__class__:
            return other.position1 == self.position1 and other.length == self.length
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.position1}, {self.length})'
