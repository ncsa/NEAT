"""
Classes for the mutation models used for each variant type to be included,
along with helper functions, used in this simulation. Every Variant type in variants > variant_types
must have a corresponding model in order to be fully implemented.
"""

import logging

from numpy.random import Generator

from .default_sequencing_error_model import *

__all__ = [
    "FragmentLengthModel"
]

_LOG = logging.getLogger(__name__)


class FragmentLengthModel:
    """
    A model of the fragment length based on mean and standard deviation of the dataset. Used both
    to generate random fragment lengths and random read lengths. Since a read is essentially a fragment as well,
    and the stastistical models used in NEAT are similar, we'll use fragment to mean read here.

    :param fragment_mean: the mean of the collection of fragment lengths derived from data
    :param fragment_std: the standard deviation of the collection of fragment lengths derived from data
    :param rng: the random number generator for the run
    """

    def __init__(
            self,
            fragment_mean: float,
            fragment_std: float
    ):
        self.fragment_mean = fragment_mean
        self.fragment_st_dev = fragment_std

    def generate_fragments(
            self,
            number_of_fragments: int,
            rng: Generator
    ) -> list:
        """
        Generates a number of fragments based on the total length needed, and the mean and standard deviation of the set

        :param number_of_fragments: The number of fragments needed.
        :param rng: the random number generator to use
        :return: A list of fragment random fragment lengths sampled from the model.
        """
        # Due to natural variation in genome lengths, it's difficult to harden this code against all the possible
        # inputs. In order to harden this code against infinite loops caused by fragment means that the wrong size for
        # the genome, we introduce a small number of standard fragments, to ensure enough variability that our code can
        # complete. a few small fragments should increase the variability of the set. Most of these are too small
        # to create a read, so they become spacers instead.
        extra_fragments = [10, 11, 12, 13, 14, 28, 31]
        # generates a distribution, assuming normality, then rounds the result and converts to ints
        dist = np.round(rng.normal(self.fragment_mean, self.fragment_st_dev, size=number_of_fragments)).astype(int)
        dist = [abs(x) for x in dist]
        # We'll append enough fragments to pad out distribution and add variability. Don't know if the cost of doing
        # this is worth it though.
        dist.extend(extra_fragments)
        rng.shuffle(dist)  # this shuffle mixes extra fragments in.
        return dist[:number_of_fragments]
