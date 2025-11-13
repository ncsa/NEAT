"""
Classes for the error models used in the simulation. This will generate errors of the type contained in
variant_models.py, each of which is based on a type from the variants submodule.
Also contained is a simple container for storing errors, and a subclass for the simulation type to be performed.
Currently, only the traditional type is supported, but we are working on adding a Markov chain based model in the near
future.
"""

import logging

from Bio.Seq import Seq
from Bio import SeqRecord
from numpy import median

from neat import variants

from ..common import ALLOWED_NUCL, NUC_IND
# Default values
from .default_mutation_model import *
from .default_sequencing_error_model import *
from .default_quality_score_model import *

from .variant_models import InsertionModel, DeletionModel, SnvModel

__all__ = [
    "ErrorContainer",
    "SequencingErrorModel",
    "TraditionalQualityModel"
]

_LOG = logging.getLogger(__name__)


class TraditionalQualityModel:
    """
    This is the traditional NEAT model for generating quality scores. This will be replaced by
    a Markov model or an optional markov model

    :param transition_matrix: 2x2 matrix that gives the probability of each base transitioning to another.
    :param quality_scores: numpy array of ints of the PHRED quality scores possible from the sequencing machine
    :param qual_score_probs: At each position along the read_length, this gives the mean and standard deviation of
        quality scores read from the dataset used to construct the model.
    :param is_uniform: Some machines use uniform quality scores. This makes simulation a little easier.
    """

    def __init__(
            self,
            average_error: float = default_avg_seq_error,
            transition_matrix: np.ndarray = default_error_transition_matrix,
            quality_scores: np.ndarray = default_quality_scores,
            qual_score_probs: np.ndarray = default_qual_score_probs,
            is_uniform: bool = False
    ):

        self.transition_matrix = transition_matrix
        self.quality_scores = quality_scores
        self.quality_score_probabilities = qual_score_probs
        self.is_uniform = is_uniform
        self.average_error = average_error

        # pre-compute the error rate for each quality score. This is the inverse of the phred score equation
        self.quality_score_error_rate: dict[int, float] = {x: 10. ** (-x / 10) for x in self.quality_scores}

        self.uniform_quality_score = None
        if self.is_uniform:
            # Set score to the lowest of the max of the quality scores and the input avg error.
            self.uniform_quality_score = min([max(self.quality_scores), int(-10. * np.log10(self.average_error) + 0.5)])

    def get_quality_scores(
            self,
            model_read_length: int,
            length: int,
            rng
    ) -> np.ndarray:
        """
        Takes a length and rng and returns an array of quality scores

        :param model_read_length: the original read length for the model
        :param length: The desired length of the quality score array
        :param rng: random number generator.
        :return: An array of quality scores.
        """
        if self.uniform_quality_score:
            return np.array([self.uniform_quality_score] * length)
        else:
            if length == model_read_length:
                quality_index_map = np.arange(model_read_length)
            else:
                # This is basically a way to evenly spread the distribution across the number of bases in the read
                quality_index_map = np.array(
                    [max([0, model_read_length * n // length]) for n in range(length)]
                )

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

            return np.array(temp_qual_array)


class MarkovQualityModel:
    def __init__(self):
        # TODO
        pass

    def get_quality_scores(self):
        # TODO
        pass


class SequencingErrorModel(SnvModel, DeletionModel, InsertionModel):
    """
    This is a SequencingErrorModel class, based on the old SequencingError. Instead of handling everything, we've
    refactored this to only handle sequencing errors. There are three types of errors possible: substitutions,
    insertions, and deletions, similar to mutations, but simpler and with different underlying models. Note that the
    three input probabilities must add up to 1.0 and the length the list of indel lengths must be equal to the length
    of its corresponding probabilities.

    :param read_length: The read length derived from real data.
    :param rescale_qualities: If set to true, NEAT will attempt to rescale the qualities based on the input error
        model, rather than using the qualities derived from the real data.
    :param variant_probs: Probability dict for each valid variant type.
    :param indel_len_model: Similar to mutation model, but simpler because errors tend not to be complicated. The
        three possible variant types for errors are Insertion, Deletion, and SNV.
    :param insertion_model: The specific parameters for the machine's insertion error rates.
    :param transition_matrix: The matrix for the transition for SNVs from one base to another.
    :param rescale_qualities: If qualities should be shifted because of a manually entered error rate override.
    :param avg_seq_error: A float giving the average rate of sequencing errors,
        either defined by data or user input.
    """

    def __init__(
            self,
            read_length: int = default_read_length,
            variant_probs: dict[variants: float] = default_error_variant_probs,
            indel_len_model: dict[int: float] = default_indel_len_model,
            insertion_model: np.ndarray = default_insertion_model,
            transition_matrix: np.ndarray = default_error_transition_matrix,
            rescale_qualities: bool = False,
            avg_seq_error: float = default_avg_seq_error,
    ):

        SnvModel.__init__(self)
        InsertionModel.__init__(self, indel_len_model)
        DeletionModel.__init__(self, indel_len_model)

        self.variant_probs = variant_probs

        self.read_length = read_length
        self.rescale_qualities = rescale_qualities
        self.insertion_model = insertion_model
        self.transition_matrix = transition_matrix
        self.average_error = avg_seq_error

    def get_sequencing_errors(
            self,
            padding: int,
            reference_segment: Seq,
            quality_scores: np.ndarray,
            num_errors,
            rng
    ):
        """
        Inserts errors of type substitution, insertion, or deletion into read_data, and assigns a quality score
        based on the container model.
        :param padding: this is the amount of space we have in the read for deletions.
        :param reference_segment: The section of the reference from which the read is drawn
        :param quality_scores: Array of quality scores for the read
        :param num_errors: The estimated number of errors to add.
        :param rng: random number generator.
        :return: Modified sequence and associated quality scores
        """

        error_indexes = []
        introduced_errors = []
        # pre-compute the error rate for each quality score. This is the inverse of the phred score equation
        quality_score_error_rate: dict[int, float] = {x: 10. ** (-x / 10) for x in quality_scores}

        # The use case here would be someone running a simulation where they want no sequencing errors.
        # No need to run any loops in this case.
        if self.average_error == 0:
            return introduced_errors
        else:
            i = len(quality_scores)
            while len(error_indexes) <= num_errors and i > 0:
                index = rng.choice(list(range(len(quality_scores))))
                if rng.random() < quality_score_error_rate[quality_scores[index]]:
                    error_indexes.append(index)
                i -= 1
            # This should fill in any errors to make sure we aren't coming up short
            median_score = median(quality_scores)
            while len(error_indexes) < num_errors:
                index = rng.choice(quality_scores)
                score = quality_scores[index]
                if score < median_score:
                    error_indexes.append(index)

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
                error_type = rng.choice(a=list(self.variant_probs), p=list(self.variant_probs.values()))

            # Deletion error
            if error_type == Deletion:
                deletion_length = self.get_deletion_length(rng)
                if padding - deletion_length < 0:
                    # No space in this read to add this deletion
                    continue
                deletion_reference = reference_segment[index: index + deletion_length + 1]
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
                insert_string = ''.join(rng.choice(ALLOWED_NUCL, size=insertion_length))
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
                 ref: str | Seq,
                 alt: str | Seq):
        self.error_type = error_type
        self.location = location
        self.length = length
        self.ref = ref
        self.alt = alt
