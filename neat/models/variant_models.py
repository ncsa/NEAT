"""
Classes for the variant models included in NEAT.
Every Variant type in variants > variant_types must have a corresponding model in order to be fully implemented.
"""
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
    """
    _type = Insertion
    _description = "An insertion of N nucleotides into a chromosome."

    def __init__(self, insert_len_model: dict[int, float]):
        # Creating probabilities from the weights
        tot = sum(insert_len_model.values())
        self.insertion_len_model = {key: val / tot for key, val in insert_len_model.items()}

    def get_insertion_length(self, rng: Generator, size: int = None) -> int | list[int, ...]:
        """
        Get size number of inserts lengths. Size == 1 results in an int return, else a list of ints.

        :param rng: The random number generator for the run
        :param size: Number of insert lengths to generate. Default is one, which returns an int value.
                     Greater than 1 returns a list of ints.
        :return: int or list of ints.
        """
        return rng.choice(
            a=list(self.insertion_len_model),
            p=[*self.insertion_len_model.values()],
            size=size,
            shuffle=False
        )


class DeletionModel(VariantModel):
    """
    This type is a deletion of some length. This is when a nucleotide or series of
    nucleotides is simply omitted during DNA replication.

    :param deletion_len_model: keys are possible lengths of deletion, values are probabilities of those values
    """
    _type = Deletion
    _description = "A deletion of a random number of bases"

    def __init__(self, deletion_len_model: dict[int, float]):
        # Creating probabilities from the weights
        tot = sum(deletion_len_model.values())
        self.deletion_len_model = {key: val/tot for key, val in deletion_len_model.items()}

    def get_deletion_length(self, rng: Generator, size: int = None) -> int | list[int]:
        """
        Get size number of inserts lengths. Size == 1 results in an int return, else a list of ints.

        :param rng: The random number generator for the run
        :param size: Number of insert lengths to generate. Default is one, which returns an int value.
                     Greater than 1 returns a list of ints.
        :return: int or list of ints.
        """
        return rng.choice(
            a=[*self.deletion_len_model],
            p=[*self.deletion_len_model.values()],
            size=size,
            shuffle=False
        )


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
    """
    _type = SingleNucleotideVariant
    _description = "Substitution"

    def __init__(
            self,
            trinuc_trans_matrices: np.ndarray = None,
            trinuc_mutation_bias: np.ndarray = None,
    ):

        self.trinuc_trans_matrices = trinuc_trans_matrices
        self.trinuc_mutation_bias = trinuc_mutation_bias
        self.no_bias = False
        if not np.any(self.trinuc_trans_matrices) and not trinuc_mutation_bias:
            self.no_bias = True
        self.trinuc_bias_map = None

        # Some local variables for modeling
        self.local_trinuc_bias: np.ndarray | None = None
        self.local_sequence: Seq | None = None

    def map_local_trinuc_bias(
            self,
            sequence: Seq,
            ngaps: np.ndarray
    ):
        """
        Create a map of a given input sequence, showing the most likely places within the sequence
        for a substitution to occur. N regions are set to 0, so no SNV will happen in those locations.

        :param sequence: A sequence of bases to create a bias model for.
        :param ngaps: A list of dictionaries, each describing an ngap.
        :return: A list of the bias factors by position.
        """
        # Cache-hit short-circuit. The old equality check `self.local_sequence == sequence`
        # did a full bytewise compare on every call (O(N) per call); identity check is
        # O(1) and covers the case where the caller reuses the same Seq object. A new
        # slice each call will miss either way, so we don't pay extra for the miss.
        if self.local_sequence is sequence:
            return

        n = len(sequence)
        # Start by assuming no bias. Each of the 64 valid trinucleotide combinations
        # mutate with equal frequency. N-gap positions become 0.
        self.local_trinuc_bias = np.ones(n, dtype=float)
        self.local_trinuc_bias[np.where(ngaps == 0)] = 0.0

        # If the model was set up with no bias, skip the trinuc weighting entirely.
        if not self.no_bias and n >= 3:
            # Vectorized trinucleotide scan: replaces 64 regex passes (one per trinuc)
            # over the sequence with a single O(N) numpy expression. Each center
            # position i in [1, n-1) gets its bias from the trinuc spanning
            # [i-1, i, i+1], encoded as base_left*16 + base_mid*4 + base_right
            # (matching the TRINUC_IND layout). Centers whose trinuc contains any
            # non-ACGT base keep the initial value (1.0 normally, 0.0 if in an N-gap).
            seq_bytes = np.frombuffer(str(sequence).upper().encode(), dtype=np.uint8)
            nuc_idx = np.full(n, -1, dtype=np.int8)
            nuc_idx[seq_bytes == ord('A')] = 0
            nuc_idx[seq_bytes == ord('C')] = 1
            nuc_idx[seq_bytes == ord('G')] = 2
            nuc_idx[seq_bytes == ord('T')] = 3
            tri_left = nuc_idx[:-2]
            tri_mid = nuc_idx[1:-1]
            tri_right = nuc_idx[2:]
            valid = (tri_left >= 0) & (tri_mid >= 0) & (tri_right >= 0)
            # Use int32 for the multiplied terms so we don't overflow int8.
            trinuc_idx = (
                tri_left.astype(np.int32) * 16
                + tri_mid.astype(np.int32) * 4
                + tri_right.astype(np.int32)
            )
            bias_table = np.asarray(self.trinuc_mutation_bias, dtype=float)
            # Where the trinuc is valid, look up the bias; otherwise preserve current value.
            center_view = self.local_trinuc_bias[1:-1]
            looked_up = bias_table[np.where(valid, trinuc_idx, 0)]
            self.local_trinuc_bias[1:-1] = np.where(valid, looked_up, center_view)

        # Now we normalize the bias
        total = self.local_trinuc_bias.sum()
        if total > 0:
            self.local_trinuc_bias = self.local_trinuc_bias / total
        self.local_sequence = sequence

    def sample_trinucs(self, rng: Generator) -> int:
        """
        Thus functon takes a trinuc map (as generated by map_trinuc_bias) or part of a trinuc
        map and determines a random location within that map, weighted by the bias (if any)

        :param rng: The random number generator for the run
        :return: the index of the chosen position
        """
        return int(rng.choice(a=np.arange(len(self.local_trinuc_bias)), p=self.local_trinuc_bias))
