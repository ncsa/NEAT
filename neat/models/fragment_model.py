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

        Every returned length is a draw from the model itself. Earlier versions spliced a fixed
        set of tiny lengths ([10, 11, 12, 13, 14, 28, 31]) into each batch as anti-infinite-loop
        padding for fragment means that did not suit the genome. Those values were never real
        samples, and while the ordinary read-length fragment floor hid them, the lower floor used
        when short inserts are kept (generate_reads.MIN_SHORT_INSERT) admitted the 28 and 31 bp
        entries — putting a spike of artificial, adapter-heavy inserts into the output. The
        samplers now bound their own retries instead (see generate_reads._sample_fragments).

        :param number_of_fragments: The number of fragments needed.
        :param rng: the random number generator to use
        :return: A list of fragment random fragment lengths sampled from the model.
        """
        # generates a distribution, assuming normality, then rounds the result and converts to ints
        dist = np.round(rng.normal(self.fragment_mean, self.fragment_st_dev, size=number_of_fragments)).astype(int)
        return np.abs(dist).tolist()
