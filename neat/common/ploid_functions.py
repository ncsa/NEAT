"""
Functions associated with ploidy, to help model polyploid organisms.
"""

import logging
import numpy as np

from numpy.random import Generator, default_rng

__all__ = [
    'pick_ploids',
    'get_genotype_string'
]

_LOG = logging.getLogger(__name__)


def pick_ploids(ploidy: int,
                homozygous_frequency: float,
                number_alts: int = 1,
                rng: Generator = default_rng()) -> np.ndarray:
    """
    Applies alts to random ploids. Picks at least one, maybe more.
    :param ploidy: how many copies of each chromosome this organism has
    :param homozygous_frequency: the freqency of homozygous mutations.
    :param number_alts: If there is more than one alternate, this will assign a random alternate to the ploids
    :param rng: the random number generator for this run
    :return: a list of strings representing the genotype of each ploid.
    """
    # number of ploids to make this mutation on (always at least 1)
    how_many = 1
    if rng.random() < homozygous_frequency:
        # We'll consider this one to be homozygous
        how_many = ploidy
    else:
        if ploidy == 1:
            # special case where heteroyzgous makes no sense
            how_many = 1
        else:
            # if it's polyploid, we'll consider it to be on roughly half the ploids
            # TODO may need to improve the modeling for polyploid, maybe
            how_many = ploidy//2

    # wp is just the temporary genotype list, a hat tip to the old version of NEAT
    wp = np.zeros(ploidy)
    while how_many > 0:
        x = rng.choice(range(ploidy))
        # pick a random alternate. in VCF terminology, 0 = REF, 1 = ALT1, 2 = ALT2, etc
        wp[x] = rng.choice(range(1, number_alts + 1))
        how_many -= 1

    return np.array(wp)


def get_genotype_string(genotype_list: np.ndarray):
    """
    This returns the string form of the genotype array, for printing in a vcf

    :param genotype_list: The list to convert to string
    :return: the string form
    """
    return "/".join([str(int(x)) for x in genotype_list])
