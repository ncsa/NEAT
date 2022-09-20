"""
A class describing one type of variant for NEAT.
"""
import numpy as np

from Bio.Seq import Seq
from .base_variant import BaseVariant


__all__ = [
    "SingleNucleotideVariant"
]


class SingleNucleotideVariant(BaseVariant):
    """
    A Single Nucleotide Variant (SNV, also known as a Single Nucleotide Polymorphism, or SNP) is when one base
    is replaced by another base. This is the most common type of variant.

    :param position1: The location of the variant and the base on the reference that is substituted.
    :param alt: The new base that replaced the original.
    :param genotype: the list of ploids with a 1 indicating where this variant is found
    :param qual_score: the quality score for this variant, as read or generated
    :param is_input: True if from an input vcf, in which case this variant will get priority.
    :param kwargs: can be used to store data from input vars or unused variables from the base class.
    """
    def __init__(self,
                 position1: int,
                 alt: str or Seq = None,
                 genotype: np.ndarray = None,
                 qual_score: int | None = None,
                 is_input: bool = False,
                 **kwargs):

        self.position1 = position1
        self.alt = alt

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
        if type(other) == self.__class__:
            return other.position1 == self.position1 and other.alt == self.alt
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.position1}, {self.alt})'
