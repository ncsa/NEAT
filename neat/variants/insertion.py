"""
A class describing one type of variant for NEAT.
"""
import numpy as np

from Bio.Seq import Seq
from .base_variant import BaseVariant


class Insertion(BaseVariant):
    """
    An insertion of N number of bases into a nucleotide. This is frequently due to a polymerase slip, and lobeled
    a duplication, but here we assume N number of bases (based on the model) can be inserted randomly, without
    needing to duplicate any existing bases.

    :param position1: The location, on the reference, where the variant begins. Note that this follows the VCF
        convention of the Ref being the same as the first base of the alternate.
    :param length: The number of bases to be inserted.
    :param alt: The ref + inserted sequence
    :param genotype: the list of ploids with a 1 indicating where this variant is found
    :param qual_score: the quality score for this variant, as read or generated
    :param is_input: True if from an input vcf, in which case this variant will get priority.
    :param kwargs: can be used to store data from input vars or unused variables from the base class.
    """

    def __init__(self,
                 position1: int,
                 length: int,
                 alt: str or Seq = None,
                 genotype: np.ndarray = None,
                 qual_score: int | None = None,
                 is_input: bool = False,
                 **kwargs):

        self.position1 = position1
        self.length = length
        self.alt = alt
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
            return other.position == self.position1 and other.alt == self.alt and other.length == self.length
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.position1}, {self.alt})'
