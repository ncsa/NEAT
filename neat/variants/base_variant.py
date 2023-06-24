"""
Module with definitions of variant classes.
"""

__all__ = [
    "BaseVariant"
]

import abc


class BaseVariant(abc.ABC):
    """
    A template for a variant to include in NEAT
    """

    @abc.abstractmethod
    def __init__(self, position1, length, position2, orientation, alt, genotype, qual_score, is_input, **kwargs):
        """
        Initialize the object.

        These parameters mean different things, depending on the variant type
        The standard metadata for a variant is given below, to help fill out the vcf:
            {"REF": str - reference string
             "ID": str - id string
             "FILTER": str - filter string for variant
             "INFO": str - info string
             "FORMAT": str - format string
             "NEAT_sample": str - sample string
             "NEAT_cancer_sample": str (optional)}

        """
        self.position1 = position1
        self.length = length
        self.position2 = position2
        self.orientation = orientation
        self.alt = alt
        self.genotype = genotype
        self.qual_score = qual_score
        self.is_input = is_input
        self.metadata = kwargs

    @abc.abstractmethod
    def __lt__(self, position):
        """
        Check if position is to the 5' side of this variant
        :param position: The position to check relative to this variant
        :return: True or False
        """

    @abc.abstractmethod
    def __gt__(self, position):
        """
        Check if position is to the 3' side of this variant.

        :param position: The position to check relative to this variant
        :return: True or False
        """

    @abc.abstractmethod
    def __le__(self, position):
        """
        Check if position is at or to the 5' side of this variant

        :param position: The position to check relative to this variant
        :return: True or False
        """

    @abc.abstractmethod
    def __ge__(self, position):
        """
        Check if position is at or to the 3' side of the end of this variant.

        :param position: The position to check relative to this variant
        :return: True or False
        """

    @abc.abstractmethod
    def __eq__(self, other):
        """
        Check if another variant is equal to this variant.

        :param other: Another variant object
        :return: True or False
        """

    @abc.abstractmethod
    def __repr__(self):
        """
        A string representation of this variant

        :return: A string representation of this variant
        """
        prefix = f'{self.__class__.__name__}'

    def get_qual_score(self):
        if self.qual_score:
            return self.qual_score
        else:
            return self.metadata['QUAL']

    def get_alt(self):
        if self.alt:
            return self.alt
        else:
            return self.metadata['ALT']

    def get_0_location(self):
        return int(self.position1)

    def get_1_location(self):
        return self.position1 + 1
