"""
Class for the mutation model, which simulates a mutation in the genome of one of the possible variant types.
"""

import logging
import sys

from numpy.random import Generator
from Bio.Seq import Seq

from neat import variants

from ..common import ALLOWED_NUCL, NUC_IND, DINUC_IND
from .variant_models import InsertionModel, DeletionModel, SnvModel
from .default_mutation_model import *
from .default_sequencing_error_model import *

__all__ = [
    "MutationModel"
]

_LOG = logging.getLogger(__name__)


class MutationModel(SnvModel, InsertionModel, DeletionModel):
    """
    A mutation model. Stores various characteristics of the mutation module used in NEAT.
    Because this class in instantiating specific mutation types, kwargs will need to be employed to
    give the mutations all the parameters they need.

    :param avg_mut_rate: average mutation rate for the modeled data set
    :param homozygous_freq: How frequently reads are homozygous (hard coded as 0.010)
                            (note, we need to investigate the equivalent with polyploid organisms)
    :param variant_probs: A list of probabilities for the possible variant types. Note that SNV chance must
        always equal 1 - sum(all other variant probs), or a RunTime error will result. If the length of this list
        doesn't match the length of the list of possible variant types, a RunTime error will result.
        The current list is: (Probability of mutation being an insertion,
                              Probability of a mutation being a deletion,
                              Probability of the mutation being a single nucleotide variant)
    :param is_cancer: Whether the model is for cancer
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
            we'll need the rng to perform certain methods. Must be set for runs.
    :param trinuc_trans_matrices: The transition matrices for the trinuc
        patterns.
    :param trinuc_mut_bias: The bias for each possible trinucleotide, as measured in the
        input dataset.
    :param insert_len_model: The model for the insertion length
    :param deletion_len_model: The model for teh deletion length
    """

    def __init__(self,
                 avg_mut_rate: float = default_avg_mut_rate,
                 homozygous_freq: float = default_homozygous_freq,
                 variant_probs: dict[variants: float, ...] = default_variant_probs,
                 transition_matrix: np.ndarray = default_mutation_sub_matrix,
                 is_cancer: bool = False,
                 rng: Generator = None,
                 # Any new parameters needed for new models should go below
                 trinuc_trans_matrices: np.ndarray = default_trinuc_trans_matrices,
                 trinuc_mut_bias: np.ndarray = default_trinuc_mut_bias,
                 insert_len_model: dict[int: float] = default_insertion_len_model,
                 deletion_len_model: dict[int: float] = default_deletion_len_model):

        # Any new mutation types will need to be instantiated in the mutation model here
        SnvModel.__init__(self,
                          trinuc_trans_matrices=trinuc_trans_matrices,
                          trinuc_mutation_bias=trinuc_mut_bias)
        InsertionModel.__init__(self, insert_len_model=insert_len_model)
        DeletionModel.__init__(self, deletion_len_model=deletion_len_model)

        self.avg_mut_rate = avg_mut_rate
        self.homozygous_freq = homozygous_freq

        if not np.isclose(sum(variant_probs.values()), 1):
            _LOG.error("Probabilities do not add up to 1.")
            sys.exit(1)

        self.variant_probs = variant_probs
        self.transition_matrix = transition_matrix
        self.is_cancer = is_cancer
        self.rng = rng
        self.all_dels = []
        self.all_ins = []

    def get_mutation_type(self) -> variants:
        """
        Picks one of the mutation types at random using a weighted list from the model.
        Note that the order of mutation types is Insertion, Deletion, SNV. To update the model selection if any
        new variant types are added, you'll need to import the type and add it to the return of this method

        :return: One of the defined variant type classes.
        """
        return self.rng.choice(a=[*self.variant_probs],
                               p=[*self.variant_probs.values()])

    def is_homozygous(self) -> bool:
        """
        Randomly samples from the homozygous frequency to get either a true or false.

        :return: True or False
        """
        return True if self.rng.random() <= self.homozygous_freq else False

    """
    Each new variant will need a generation method here).
    """
    def generate_snv(self, trinucleotide: Seq, reference_location: int) -> SingleNucleotideVariant:
        """
        This takes a location on the sequence and a location within the reference and returns a new SNV

        :param trinucleotide: The trinuc of interest for this variant
        :param reference_location: The same position, relative to the reference,
            used to retrieve the current reference base.
        :return: A randomly generated variant
        """
        # First determine which matrix to use
        transition_matrix = self.trinuc_trans_matrices[DINUC_IND[trinucleotide[0] + "_" + trinucleotide[2]]]
        # then determine the trans probs based on the middle nucleotide
        transition_probs = transition_matrix[NUC_IND[trinucleotide[1]]]
        # Creating probabilities from the weights
        transition_sum = sum(transition_probs)
        transition_probs = [x/transition_sum for x in transition_probs]
        # Now pick a random alternate, weighted by the probabilities
        alt = self.rng.choice(ALLOWED_NUCL, p=transition_probs)
        temp_snv = SingleNucleotideVariant(reference_location, alt=alt)
        self.all_ins.append(temp_snv)
        return temp_snv

    def generate_insertion(self, location: int, ref: Seq) -> Insertion:
        """
        This method generates an insertion object, based on the insertion model

        :param location: The location of the variant, relative to the reference
        :param ref: The reference for which to generate the variant
        :return:
        """
        # Note that insertion length model is based on the number of bases inserted. We add 1 to the length
        # to get the length of the VCF version of the variant.
        length = self.get_insertion_length() + 1
        insertion_string = ''.join(self.rng.choice(ALLOWED_NUCL, size=length))
        alt = ref + insertion_string
        return Insertion(location, length, alt)

    def generate_deletion(self, location) -> Deletion:
        """
        Takes a location and returns a deletion object

        :param location:
        :return:
        """
        # Note that the deletion length model is based on the number of bases deleted,
        # so we add 1 to account for the common base between ref and alternate.
        length = self.get_deletion_length() + 1
        # Plus one so we make sure to grab the first base too.
        # Note: if we happen to go past the end of the sequence, it will just be shorter.
        temp_del = Deletion(location, length)
        self.all_dels.append(temp_del)
        return temp_del
