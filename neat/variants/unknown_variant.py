"""
A class describing a placeholder type of variant for NEAT, used for input variants of unsupported type.
"""

from Bio.Seq import Seq
import numpy as np
from .base_variant import BaseVariant


__all__ = [
    "UnknownVariant"
]


class UnknownVariant(BaseVariant):
    """
    An unknown variant type from an input vcf. This will just hold the data from the vcf for later output.

    :param position1: The location of the variant and the base on the reference that is substituted.
    :param is_input: True if from an input vcf, in which case this variant will get priority.
    :param kwargs: can be used to store data from input vars or unused variables from the base class.
    """
    def __init__(self,
                 position1: int,
                 genotype: np.ndarray,
                 qual_score: int,
                 is_input: bool = False,
                 **kwargs):

        self.position1 = position1
        self.genotype = genotype
        self.qual_score = qual_score
        self.is_input = is_input
        self.metadata = kwargs

    def __lt__(self, position):
        return position < self.position1

    def __gt__(self, position):
        return position > self.position1

    def __le__(self, position):
        return position <= self.position1

    def __ge__(self, position):
        return position >= self.position1

    def __eq__(self, other):
        if other.type == self.__class__:
            return other.position1 == self.position1 and other.alt == self.alt and self.genotype == other.genotype
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.position1})'
