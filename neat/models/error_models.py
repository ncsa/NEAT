"""
Classes for the error models used in the simulation. This will generate errors of the type contained in
variant_models.py, each of which is based on a type from the variants submodule.
Also contained is a simple container for storing errors, and a subclass for the simulation type to be performed.
Currently, only the traditional type is supported, but we are working on adding a Markov chain based model in the near
future.
"""

import logging

from numpy.random import Generator
from Bio.Seq import Seq
from Bio import SeqRecord

from neat import variants

from ..common import ALLOWED_NUCL, NUC_IND
from .default_mutation_model import *
from .default_sequencing_error_model import *
from .variant_models import InsertionModel, DeletionModel, SnvModel

__all__ = [
    "ErrorContainer"
]

_LOG = logging.getLogger(__name__)


class TraditionalModel:
    """
    This is the traditional NEAT model for generating quality scores. This will be replaced by
    a Markov model or an optional markov model

    :param transition_matrix: 2x2 matrix that gives the probability of each base transitioning to another.
    :param quality_scores: numpy array of ints of the PHRED quality scores possible from the sequencing machine
    :param qual_score_probs: At each position along the read_length, this gives the mean and standard deviation of
        quality scores read from the dataset used to construct the model.
    :param is_uniform: Some machines use uniform quality scores. This makes simulation a little easier.
    """

    def __init__(self,
                 transition_matrix: np.ndarray = default_error_transition_matrix,
                 quality_scores: np.ndarray = default_quality_scores,
                 qual_score_probs: np.ndarray = default_qual_score_probs,
                 is_uniform: bool = False,
                 ):

        self.transition_matrix = transition_matrix
        self.quality_scores = quality_scores
        self.quality_score_probabilities = qual_score_probs
        self.is_uniform = is_uniform

        # pre-compute the error rate for each quality score. This is the inverse of the phred score equation
        self.quality_score_error_rate: dict[int, float] = {x: 10. ** (-x / 10) for x in self.quality_scores}

        self.uniform_quality_score = None
        if self.is_uniform:
            # Set score to the lowest of the max of the quality scores and the input avg error.
            self.uniform_quality_score = min([max(self.quality_scores), int(-10. * np.log10(self.average_error) + 0.5)])

    def get_quality_scores(self,
                           input_read_length: int,
                           orig_read_len: int,
                           rng
                           ) -> list:
        """
        Takes a read_length and rng and returns an array of quality scores

        :param input_read_length: The desired length of the quality score array
        :param orig_read_len: the original read length for the model
        :param rng: random number generator.
        :return: An array of quality scores.
        """
        if self.uniform_quality_score:
            return [self.uniform_quality_score] * input_read_length
        else:
            quality_index_map = self.quality_index_remap(input_read_length, orig_read_len)
            temp_qual_array = []
            for i in quality_index_map:
                score = rng.normal(
                    self.quality_score_probabilities[i][0],
                    scale=self.quality_score_probabilities[i][1]
                )
                # make sure score is in range and an int
                score = round(score)
                if score > 42:
                    score = 42
                if score < 1:
                    score = 1

                temp_qual_array.append(score)

            return temp_qual_array

    def quality_index_remap(self, input_read_length, original_read_length):
        """
        Adjusts the quality map to the suitable read length.

        :param input_read_length: The desired length for the current read.
        :param original_read_length: The set read-length for the model.
        :return: An index map from the default read length to the new one.
        """
        if input_read_length == original_read_length:
            return np.arange(original_read_length)
        else:
            # This is basically a way to evenly spread the distribution across the number of bases in the read
            return np.array([max([0, original_read_length * n // input_read_length]) for n in range(input_read_length)])


class SequencingErrorModel(SnvModel, DeletionModel, InsertionModel, TraditionalModel):
    """
    This is a SequencingErrorModel class, based on the old SequencingError. This covers both errors and quality
    scores, since they are intimately related. There are three types of
    errors possible: substitutions, insertions, and deletions, similar to mutations, but
    simpler. Note that the three input probabilities must add up to 1.0 and the length the list of indel lengths
    must be equal to the length of its corresponding probabilities.

    :param read_length: The read length derived from real data.
    :param rescale_qualities: If set to true, NEAT will attempt to rescale the qualities based on the input error
        model, rather than using the qualities derived from the real data.
    :param variant_probs: Probability dict for each valid variant type
    :param indel_len_model: Similar to mutation model, but simpler because errors tend not to be complicated. The
        three possible variant types for errors are Insertion, Deletion, and SNV
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
        we'll need the rng to perform certain methods.
    :param avg_seq_error: A float giving the average rate of sequencing errors,
        either defined by data or user input.
    """

    def __init__(self,
                 read_length: int = default_read_length,
                 variant_probs: dict[variants: float] = default_error_variant_probs,
                 indel_len_model: dict[int: float] = default_indel_len_model,
                 insertion_model: np.ndarray = default_insertion_model,
                 rescale_qualities: bool = False,
                 avg_seq_error: float = default_avg_seq_error,
                 rng: Generator = None):

        SnvModel.__init__(self)
        InsertionModel.__init__(self, indel_len_model)
        DeletionModel.__init__(self, indel_len_model)

        self.variant_probs = variant_probs

        self.read_length = read_length
        self.read_length = read_length
        self.rescale_qualities = rescale_qualities
        self.insertion_model = insertion_model
        self.average_error = avg_seq_error
        self.rng = rng

    def get_sequencing_errors(self,
                              read_length: int,
                              padding: int,
                              reference_segment: SeqRecord,
                              quality_scores: np.ndarray,
                              rng
                              ):
        """
        Inserts errors of type substitution, insertion, or deletion into read_data, and assigns a quality score
        based on the container model.
        :param read_length: The length of the read to generate errors for.
        :param padding: this is the amount of space we have in the read for deletions.
        :param reference_segment: The section of the reference from which the read is drawn
        :param quality_scores: Array of quality scores for the read
        :return: Modified sequence and associated quality scores
        :param rng: random number generator.
        """

        error_indexes = []
        introduced_errors = []

        # The use case here would be someone running a simulation where they want no sequencing errors.
        # No need to run any loops in this case.
        if self.average_error == 0:
            return introduced_errors
        else:
            for i in range(read_length):
                if rng.random() < self.quality_score_error_rate[quality_scores[i]]:
                    error_indexes.append(i)

        total_indel_length = 0
        # To prevent deletion collisions
        del_blacklist = []

        for index in error_indexes[::-1]:
            # determine error type. Most will be SNVs
            error_type = SingleNucleotideVariant

            # Not too sure about how realistic it is to model errors as indels, but I'm leaving the code in for now.

            # This is to prevent deletion error collisions and to keep there from being too many indel errors.
            if 0 < index < self.read_length - max(
                    self.deletion_len_model) and total_indel_length > self.read_length // 4:
                error_type = self.rng.choice(a=list(self.variant_probs), p=list(self.variant_probs.values()))

            # Deletion error
            if error_type == Deletion:
                deletion_length = self.get_deletion_length()
                if padding - deletion_length < 0:
                    # No space in this read to add this deletion
                    continue
                deletion_reference = reference_segment.seq[index: index + deletion_length + 1]
                deletion_alternate = deletion_reference[0]
                introduced_errors.append(
                    ErrorContainer(Deletion, index, deletion_length, deletion_reference, deletion_alternate)
                )
                total_indel_length += deletion_length

                del_blacklist.extend(list(range(index, index + deletion_length)))
                padding -= deletion_length

            elif error_type == Insertion:
                insertion_length = self.get_insertion_length()
                insertion_reference = reference_segment[index]
                insert_string = ''.join(self.rng.choice(ALLOWED_NUCL, size=insertion_length))
                insertion_alternate = insertion_reference + insert_string
                introduced_errors.append(
                    ErrorContainer(Insertion, index, insertion_length, insertion_reference, insertion_alternate)
                )
                total_indel_length += insertion_length

            # Insert substitution error
            # Programmer note: if you add new error types, they can be added as elifs above, leaving the final
            # else dedicated to SNVs.
            else:
                snv_reference = reference_segment[index]
                nuc_index = NUC_IND[snv_reference]
                # take the zero index because this returns a list of length 1.
                snv_alt = rng.choice(ALLOWED_NUCL, p=self.transition_matrix[nuc_index])
                introduced_errors.append(
                    ErrorContainer(SingleNucleotideVariant, index, 1, snv_reference, snv_alt)
                )

        # Remove blacklisted errors
        for i in range(len(introduced_errors) - 1, -1, -1):
            if introduced_errors[i].location in del_blacklist:
                del introduced_errors[i]

        return introduced_errors, max(padding, 0)

    def generate_quality_scores(
            self,
            inp_read_len: int,
    ):
        if self.uniform_quality_score:
            return np.array([self.uniform_quality_score] * inp_read_len)
        else:

            temp_quality_array = self.get_quality_scores(
                inp_read_len,
                self.read_length,
                self.rng
            )

            if self.rescale_qualities:
                # Note that rescaling won't have much effect on binned quality scores.
                # First we rescale the quality score literals then convert back to phred score (with the 0.5
                # causing borderline cases to take the next highest number).
                rescaled_quals = [max([0, int(-10. * np.log10(self.average_error *
                                                              self.quality_score_error_rate[n]) + 0.5)])
                                  for n in temp_quality_array]
                # Now rebin the quality scores.
                temp_qual_array = np.array(rescaled_quals)
            else:
                temp_qual_array = np.array(temp_quality_array)

            return temp_qual_array[:inp_read_len]


class ErrorContainer:
    """
    Holds data for a single error

    :param error_type - the type of error this is
    :param location - the index of the start position of the variant in 0-based coordinates
    :param length - the length of the error
    :param ref - the reference sequence of the error includes base before insertion or deletion, as applicable,
        which is the same notation used in a VCF file.
    :param alt - the alternate sequence of the error (i.e., the error itself)
    """
    def __init__(self,
                 error_type: VariantTypes,
                 location: int,
                 length: int,
                 ref: str or Seq,
                 alt: str or Seq):
        self.error_type = error_type
        self.location = location
        self.length = length
        self.ref = ref
        self.alt = alt
