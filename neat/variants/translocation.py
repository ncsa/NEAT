"""
A class describing one type of variant for NEAT.
"""
import numpy as np

from Bio.Seq import Seq
from .base_variant import BaseVariant


class Translocation(BaseVariant):
    """
    A translocation occurs when one segment of a chromosome switches places with a segment from another chromosome (or
    another copy of the same chromosome). Translocations are complex variants, difficult for software to properly
    identify, and are often associated with complex problems, such as cancer or genetic defects.

    :param kwargs: can be used to store data from input vars.
    """

    def __init__(self,
                 position1: int,
                 length: int,
                 position2: int,
                 orientation: int,
                 alt: str or Seq,
                 genotype: np.ndarray = None,
                 qual_score: int | None = None,
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
            return other.position == self.position1
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.position1})'

