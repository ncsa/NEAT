"""
Classes for the variant models included in NEAT.
Every Variant type in variants > variant_types must have a corresponding model in order to be fully implemented.
"""

import re
import logging
import abc

from numpy.random import Generator
from Bio.Seq import Seq


from ..common import TRINUC_IND
from .default_mutation_model import *
from .default_sequencing_error_model import *

__all__ = [
    "InsertionModel",
    "DeletionModel",
    "SnvModel",
]

_LOG = logging.getLogger(__name__)


class VariantModel(abc.ABC):
    _type = ...
    _description = ...

    @abc.abstractmethod
    def __init__(self):
        ...


class InsertionModel(VariantModel):
    """
    An insertion type mutation, which is usually when the DNA replication process slips and
    repeats a region of the chromosome, but can be random or happen by other means.

    Also includes a summary of the dataset, and a method to use an rng to fetch
    an insertion length or list of insertion lengths

    :param insert_len_model: keys are possible lengths of inserts, values are weights of each
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
                we'll need the rng to perform certain methods.
    """
    _type = Insertion
    _description = "An insertion of N nucleotides into a chromosome."

    def __init__(self,
                 insert_len_model: dict[int: float, ...],
                 rng: Generator = None):
        # Creating probabilities from the weights
        tot = sum(insert_len_model.values())
        self.insertion_len_model = {key: val / tot for key, val in insert_len_model.items()}
        self.rng = rng

    def get_insertion_length(self, size: int = None) -> int | list[int, ...]:
        """
        Get size number of inserts lengths. Size == 1 results in an int return, else a list of ints.

        :param size: Number of insert lengths to generate. Default is one, which returns an int value.
                     Greater than 1 returns a list of ints.
        :return: int or list of ints.
        """
        return self.rng.choice(a=list(self.insertion_len_model),
                               p=[*self.insertion_len_model.values()],
                               size=size, shuffle=False)


class DeletionModel(VariantModel):
    """
    This type is a deletion of some length. This is when a nucleotide or series of
    nucleotides is simply omitted during DNA replication.

    :param deletion_len_model: keys are possible lengths of deletion, values are probabilities of those values
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
            we'll need the rng to perform certain methods.
    """
    _type = Deletion
    _description = "A deletion of a random number of bases"

    def __init__(self,
                 deletion_len_model: dict[int: float, ...],
                 rng: Generator = None):
        # Creating probabilities from the weights
        tot = sum(deletion_len_model.values())
        self.deletion_len_model = {key: val/tot for key, val in deletion_len_model.items()}
        self.rng = rng

    def get_deletion_length(self, size: int = None) -> int | list[int, ...]:
        """
        Get size number of inserts lengths. Size == 1 results in an int return, else a list of ints.

        :param size: Number of insert lengths to generate. Default is one, which returns an int value.
                     Greater than 1 returns a list of ints.
        :return: int or list of ints.
        """
        return self.rng.choice(a=[*self.deletion_len_model],
                               p=[*self.deletion_len_model.values()],
                               size=size, shuffle=False)


class SnvModel(VariantModel):
    """
    This type is a substitution of a single base in a DNA strand, also called a single nucleotide variant (SNV)
    or single nucleotide polymorphism (SNP). This is when a nucleotide or series of
    nucleotides is simply omitted during DNA replication. This is the most common type of
    mutation.

    The original NEAT researchers found a correlation between trinucleotide context of the SNP
    and the transition rate. The summary input parameters describe that relationship.
    The default model assumes no trinculeotide bias in SNPs and chose a base to transition to at random.

    :param trinuc_trans_matrices: 3D array of matrices representing the transition probabilities of each trinucleotide combination,
                                  where each matrix represents the probability of a particular
    :param trinuc_mutation_bias: The bias of each of the 64 possible trinucleotide combinations,
                              as derived from input data. Our base assumption will be no bias
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
            we'll need the rng to perform certain methods.
    """
    _type = SingleNucleotideVariant
    _description = "Substitution"

    def __init__(self,
                 trinuc_trans_matrices: np.ndarray = None,
                 trinuc_mutation_bias: np.ndarray = None,
                 rng: Generator = None):

        self.trinuc_trans_matrices = trinuc_trans_matrices
        self.trinuc_mutation_bias = trinuc_mutation_bias
        self.no_bias = False
        if not np.any(self.trinuc_trans_matrices) and not trinuc_mutation_bias:
            self.no_bias = True
        self.trinuc_bias_map = None
        self.rng = rng

        # Some local variables for modeling
        self.local_trinuc_bias: np.array = None
        self.local_sequence: Seq or None = None

    def map_local_trinuc_bias(self,
                              sequence: Seq,
                              ngaps: np.ndarray):
        """
        Create a map of a given input sequence, showing the most likely places within the sequence
        for a substitution to occur. N regions are set to 0, so no SNV will happen in those locations.

        :param sequence: A sequence of bases to create a bias model for.
        :param ngaps: A list of dictionaries, each describing an ngap.
        :return: A list of the bias factors by position.
        """
        # If the sequence is unchanged, we don't want to process again.
        if not self.local_sequence == sequence:

            # start by assuming no bias. Each of the 64 valid trinucleotide combinations mutate with equal frequency.
            self.local_trinuc_bias = np.ones(len(sequence), dtype=float)

            # If there are ngaps we'll set the bias rate to zero in those areas
            self.local_trinuc_bias[np.where(ngaps == 0)] = 0.0

            # If the model was set up with no bias, then we skip the biasing part
            if not self.no_bias:
                # Update the map bias at the central position for that trinuc
                for trinuc in ALL_TRINUCS:
                    for match in re.finditer(trinuc, str(sequence)):
                        self.local_trinuc_bias[match.start() + 1] = self.trinuc_mutation_bias[TRINUC_IND[trinuc]]

            # Now we normalize the bias
            self.local_trinuc_bias = self.local_trinuc_bias / sum(self.local_trinuc_bias)
            self.local_sequence = sequence

    def sample_trinucs(self) -> int:
        """
        Thus functon takes a trinuc map (as generated by map_trinuc_bias) or part of a trinuc
        map and determines a random location within that map, weighted by the bias (if any)

        :return: the index of the chosen position
        """
        return int(self.rng.choice(a=np.arange(len(self.local_trinuc_bias)), p=self.local_trinuc_bias))
