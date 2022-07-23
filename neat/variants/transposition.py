"""
A class describing one type of variant for NEAT.
"""
import numpy as np

from Bio.Seq import Seq
from .base_variant import BaseVariant


class Transposition(BaseVariant):
    """
    A transposon is a type of variant that can move from one location to another, creating a transposition variant.
    A transposon can even move between chromosomes. Which chromosome the transposon moves to is handled outside of this
    class. Instead, the following parameters define the transposition:

    :param position1: The location of the first base of the variant transposable element.
    :param length: the length of the transposon
    :param position2: The location of the first base as it appears in the reference.
    :param orientation: 1 for positive orientation, -1 for negative.
    :param kwargs: can be used to store data from input vars or unused variables from the base class.
    """

    def __init__(self,
                 position1: int,
                 length: int,
                 position2: int,
                 orientation: int,
                 genotype: np.ndarray = None,
                 qual_score: int = None,
                 is_input: bool = False,
                 **kwargs):

        self.position1 = position1
        self.length = length
        self.position2 = position2
        self.orientation = orientation
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
            return other.position1 == self.position1 and other.alt == self.alt
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.position1})'
