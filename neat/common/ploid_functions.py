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
    :param number_alts: If there is more than one alt, this will assign a random alt to the ploids
    :param rng: the random number generator for this run
    :return: a list of strings representing the genotype of each ploid.
    """
    # number of ploids to make this mutation on (always at least 1)
    if rng.random() < homozygous_frequency:
        if ploidy <= 2:
            # If it's homozygous with a ploidy of 2 it's on both.
            how_many = ploidy
        else:
            # if it's polyploid, it's on more than one ploid
            # TODO may need to improve the modeling for polyploid
            how_many = 1
            for i in range(ploidy):
                # Not totally sure how to model this, so I'm counting each
                # ploid as a separate homozygous event. That doesn't exactly make
                # sense though, so we'll improve this later.
                if rng.random() < homozygous_frequency:
                    how_many += 1
                else:
                    break
    else:
        how_many = 1

    # wp is just the temporary genotype list, a hat tip to the old version of NEAT
    wp = np.zeros(ploidy)
    while how_many > 0:
        x = rng.choice(range(ploidy))
        # pick a random alt. in VCF terminology, 0 = REF, 1 = ALT1, 2 = ALT2, etc
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
