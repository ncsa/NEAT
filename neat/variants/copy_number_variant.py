"""
A class describing one type of variant for NEAT.
"""
import numpy as np

from Bio.Seq import Seq
from .base_variant import BaseVariant


class CopyNumberVariant(BaseVariant):
    """
    A copy number variation occurs when a segment of repeated DNA either shrinks or grows. This can have numerous
    biological consequences, and is an indicator for various types of cancer.

    :param position1: location of the variant.
    """

    def __init__(self,
                 position1: int,
                 length: int = None,
                 position2: int = None,
                 orientation: int = None,
                 alt: str or Seq = None,
                 genotype: np.ndarray = None,
                 qual_score: int | None = None,
                 is_input:  bool = False,
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
            return other.position == self.position1
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.position1})'
