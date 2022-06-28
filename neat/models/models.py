"""
Classes for the models, along with helper functions, used in this simulation.
"""

import numpy as np
import re
import logging

from bisect import bisect_left
from numpy.random import Generator
from numpy.random import default_rng
from Bio.Seq import Seq

from ..common import ALL_TRI

__all__ = [
    "MutationModel",
    "SequencingErrorModel",
    "GcModel",
    "FragmentLengthModel",
    "Insertion",
    "Deletion",
    "Substitution"
]

_LOG = logging.getLogger(__name__)


def take_closest(my_list: np.ndarray[int, ...], my_number: float):
    """
    Helper function that assumes my_list is sorted. Returns closest value to my_number.
    This is necessary because most quality scores are binned these days.

    If two numbers are equally close, return the smaller number.
    """
    pos = bisect_left(my_list, my_number)
    if pos == 0:
        return my_list[0]
    if pos == len(my_list):
        return my_list[-1]
    before = my_list[pos - 1]
    after = my_list[pos]
    if after - my_number < my_number - before:
        return after
    else:
        return before


class VariantTypes:
    """
    Class representing the possible mutation types modeled in NEAT. Currently, there are 3, with plans to add more.
    Each subclass will have associated methods to sample the statistics of interest for that type of mutation.

    :param name: name of the mutation type
    :param name: brief description of mutation type
    """
    def __init__(self, name: str, description: str):
        self.name = name,
        self.description = description


class Insertion(VariantTypes):
    """
    An insertion type mutation, which is usually when the DNA replication process slips and
    repeats a region of the chromosome, but can be random or happen by other means.

    Also includes a summary of the dataset, and a method to use an rng to fetch
    an insertion length or list of insertion lengths

    :param insert_lengths: Potential lengths of inserts
    :param ins_len_weights: Weights of the list above
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
                we'll need the rng to perform certain methods.
    """
    def __init__(self,
                 insert_lengths: np.ndarray[int, ...],
                 ins_len_weights: np.ndarray[float, ...],
                 rng: Generator = None):
        super().__init__("Insertion", "An insertion of N nucleotides into a chromosome.")
        self.insert_lengths = insert_lengths
        self.ins_len_weights = ins_len_weights
        self.rng = rng

    def get_insertion_length(self, size: int = 1) -> int | list[int, ...]:
        """
        Get size number of inserts lengths. Size == 1 results in an int return, else a list of ints.

        :param size: Number of insert lengths to generate. Default is one, which returns an int value.
                     Greater than 1 returns a list of ints.
        :return: int or list of ints.
        """
        return self.rng.choice(a=self.insert_lengths, p=self.ins_len_weights, size=size, shuffle=False)


class Deletion(VariantTypes):
    """
    This type is a deletion of some length. This is when a nucleotide or series of
    nucleotides is simply omitted during DNA replication.

    :param deletion_lengths: Potential lengths of deletions
    :param del_len_weights: Weights of the list above
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
            we'll need the rng to perform certain methods.
    """
    def __init__(self,
                 deletion_lengths: np.ndarray[int, ...],
                 del_len_weights: np.ndarray[float, ...],
                 rng: Generator = None):
        super().__init__("Deletion", "A deletion of N bases")
        self.deletion_lengths = deletion_lengths
        self.del_len_weights = del_len_weights
        self.rng = rng

    def get_deletion_length(self, size: int = 1) -> int | list[int, ...]:
        """
        Get size number of inserts lengths. Size == 1 results in an int return, else a list of ints.

        :param size: Number of insert lengths to generate. Default is one, which returns an int value.
                     Greater than 1 returns a list of ints.
        :return: int or list of ints.
        """
        return self.rng.choice(a=self.deletion_lengths, p=self.del_len_weights, size=size, shuffle=False)


class Substitution(VariantTypes):
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
    :param trinuc_trans_bias: The bias of each of the 64 possible trinucleotide combinations,
                              as derived from input data. Our base assumption will be no bias
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
            we'll need the rng to perform certain methods.
    """
    def __init__(self,
                 trinuc_trans_matrices: np.ndarray[float, ...] = None,
                 trinuc_trans_bias: dict = None,
                 rng: Generator = None):

        super().__init__("Substitution", "A substitution of a single base")
        self.trinuc_trans_matrices = trinuc_trans_matrices
        self.trinuc_trans_bias = trinuc_trans_bias
        self.no_bias = False
        if not self.trinuc_trans_matrices and not trinuc_trans_bias:
            self.no_bias = True
        self.trinuc_bias_map = None
        self.rng = rng

    def map_trinuc_bias(self,
                        name: str,
                        sequence: Seq,
                        ngaps: list[dict]) -> np.array:
        """
        Create a map of a given input sequence, showing the most likely places within the sequence
        for a substitution to occur. This model assumes the input consists of only one of the four
        allowed nucleotides.

        :param name: name of the sequence
        :param sequence: A sequence of bases to create a bias model for.
        :param ngaps: A list of dictionaries, each describing an ngap.
        :return: A list of the bias factors by position.
        """
        # start by assuming no bias. Each of the 64 valid trinucleotide combinations mutate with equal frequency.
        default_bias = 1/64
        # Create a map initially assuming uniform mutation.
        return_map = np.full(shape=len(sequence), fill_value=default_bias)
        # If there are ngaps we'll set the bias rate to zero
        for gap in ngaps:
            if gap['chrom'] == name:
                return_map[gap['start']: gap['end']] = 0.0

        # If the model was setup with no bias, return just the default map.
        if self.no_bias:
            return return_map

        # Update the map bias at the central position for that trinuc
        for trinuc in ALL_TRI:
            for match in re.finditer(trinuc, str(sequence)):
                return_map[match.start() + 1] = self.trinuc_trans_bias[trinuc]

        return return_map

    def sample_trinucs(self,
                       trinuc_map: list) -> int:
        """
        Thus functon takes a trinuc map (as generated by map_trinuc_bias) or part of a trinuc
        map and determines a random location within that map, weighted by the bias (if any)

        :param trinuc_map: A map of the trinucleotide bias for the region of interest
        :return: the index of the chosen position
        """
        if self.no_bias:
            # choice by default assumes a uniform distribution
            return int(self.rng.choice(a=np.arange(len(trinuc_map))))
        else:
            return int(self.rng.choice(a=np.arange(len(trinuc_map)), p=trinuc_map))


class MutationModel(Substitution, Deletion, Insertion):
    """
    A mutation model. Stores various characteristics of the mutation module used in NEAT.
    Because this class in instantiating specific mutation types, kwargs will need to be employed to
    give the mutations all the parameters they need.

    :param avg_mut_rate: average mutation rate for the modeled data set
    :param homozygous_freq: How frequently reads are homozygous (hard coded as 0.010)
                            (note, we need to investigate the equivalent with polyploid organisms)
    :param insertion_chance: Probability of mutaton being an insertion.
    :param deletion_chance: Probability of a mutation being a deletion.
    :param substitution_chance: Probability of a mutation being a substitution
    :param is_cancer: Whether the model is for cancer
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
            we'll need the rng to perform certain methods.
    """
    def __init__(self,
                 avg_mut_rate: float,
                 homozygous_freq: float,
                 insertion_chance: float,
                 deletion_chance: float,
                 substitution_chance: float,
                 trinuc_trans_matrices: np.ndarray[float, ...],
                 trinuc_trans_bias: dict,
                 insertion_lengths: np.ndarray[int, ...],
                 insertion_weights: np.ndarray[float, ...],
                 deletion_lengths: np.ndarray[int, ...],
                 deletion_weights: np.ndarray[float, ...],
                 is_cancer: bool = False,
                 rng: Generator = None):

        # Any new mutation types will need to be instantiated in the mutation model here

        super(Substitution).__init__(trinuc_trans_matrices, trinuc_trans_bias)
        super(Insertion).__init__(insertion_lengths, insertion_weights)
        super(Deletion).__init__(deletion_lengths, deletion_weights)

        # We'll insert more mutation types as we go
        self.avg_mut_rate = avg_mut_rate
        self.homozygous_freq = homozygous_freq

        # Check to make sure the input for frequencies is valid
        if not np.isclose(insertion_chance + deletion_chance + substitution_chance, 1.0):
            raise ValueError("Chances for each type of mutation must add up to 1.0")

        self.insertion_chance = insertion_chance
        self.deletion_chance = deletion_chance
        self.substitution_chance = substitution_chance

        self.is_cancer = is_cancer
        self.rng = rng

    def get_mutation_type(self) -> VariantTypes:
        """
        Picks one of the mutation types at random using a weighted list from the model.
        Note that the order of mutation types is Insertion, Deletion, Substitution.

        :return: One of the defined mutation classes.
        """
        return self.rng.choice(a=repr(MutationModel.__mro__),
                               p=(self.insertion_chance, self.deletion_chance, self.substitution_chance))

    def is_homozygous(self) -> bool:
        """
        Randomly samples from the homozygous frequency to get either a true or false.

        :return: True or False
        """
        return True if self.rng.random() <= self.homozygous_freq else False


class SequencingErrorModel(Substitution, Deletion, Insertion):
    """
    This is a SequencingErrorModel class, based on the old SequencingError. This covers both errors and quality
    scores, since they are intimately related. There are three types of
    errors possible: substitutions, insertions, and deletions, similar to mutations, but
    simpler. Note that the three input probabilities must add up to 1.0 and the length the list of indel lengths
    must be equal to the length of its corresponding probabilities.

    :param avg_seq_error: A float giving the average rate of sequencing errors, either defined by data or user input.
    :param read_length: The read length derived from real data.
    :param transition_matrix: 2x2 matrix that gives the probability of each base transitioning to another.
    :param quality_scores: numpy array of ints of the PHRED quality scores possible from the sequencing machine
    :param qual_score_probs: A numpy array of arrays of floats. Each inner array has a length equal to the length
                             of the quality scores list. The number of arrays is equal to the number of positions
                             of reads in the dataset. The default model uses an array of length 101. Each row
                             tells you the probability of getting each quality score at that position along the read.
    :param rescale_qualities: If set to true, NEAT will attempt to rescale the qualities based on the input error
                              model, rather than using the qualities derived from the real data.
    :param insertion_probability: The chance an error is an insertion.
    :param deletion_probability: The chance an error is a deletion.
    :param substitution_probability: The chance an error is a substitution.
    :param indel_lengths: Similar to mutation model, but simpler because errors tend not to be complicated. One
                          list covers both insertions and deletions.
    :param indel_probs: Similar to mutation model, but simpler because errors tend not to be complicated. One
                        list covers both insertions and deletions.
    :param is_uniform: Some machines use uniform quality scores. This makes simulation a little easier.
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
            we'll need the rng to perform certain methods.
    """
    def __init__(self,
                 avg_seq_error: float,
                 read_length: int,
                 transition_matrix: np.ndarray[float, ...],
                 quality_scores: np.ndarray[int, ...],
                 qual_score_probs: np.ndarray[np.ndarray[float, ...], ...],
                 rescale_qualities: bool,
                 insertion_probability: float = 0.99,
                 deletion_probability: float = 0.004,
                 substitution_probability: float = 0.006,
                 indel_lengths: np.ndarray[int, ...] = np.array([1, 2]),
                 indel_probs: np.ndarray[float, ...] = np.array([0.999, 0.001]),
                 is_uniform: bool = False,
                 rng: Generator = None):

        if not np.isclose(insertion_probability + deletion_probability + substitution_probability, 1.0):
            raise ValueError("The values of insertion/deletion/substitution probabilities for the sequencing error "
                             "must add up to 1")

        super(Substitution).__init__()
        super(Insertion).__init__(indel_lengths, indel_probs)
        super(Deletion).__init__(indel_lengths, indel_probs)

        self.mutation_types = (Insertion, Deletion, Substitution)
        self.average_error = avg_seq_error
        self.read_length = read_length
        self.transition_matrix = transition_matrix

        # Probability that if there is an error, it is an insertion
        self.insertion_probability = insertion_probability
        # Probability that if there is an error, it is a deletion
        self.deletion_probability = deletion_probability
        # probability that if there is an error, it is a substitution
        self.substitution_probability = substitution_probability

        self.quality_scores = quality_scores
        # pre-compute the error rate for each quality score. This is the inverse of the phred score equation
        self.quality_score_error_rate: dict[int, float] = {x: 10. ** (-x / 10) for x in self.quality_scores}
        self.read_length = read_length
        self.quality_score_probability_matrix = qual_score_probs
        self.rescale_qualities = rescale_qualities
        self.is_uniform = is_uniform
        self.uniform_quality_score = None
        if self.is_uniform:
            converted_avg_err = take_closest(self.quality_scores,
                                             int(-10. * np.log10(self.average_error)))
            # Set score to the lower of the max of the quality scores and the bin closest to the input avg error.
            self.uniform_quality_score = min([max(self.quality_scores), converted_avg_err])
        self.rng = rng

    def get_sequencing_errors(self,
                              quality_scores: np.ndarray[int, ...]):
        """
        Inserts errors of type substitution, insertion, or deletion into read_data, and assigns a quality score
        based on the container model.
        :param quality_scores: Array of quality scores for the read
        :return: modified sequence and associated quality scores

        """

        sequencing_errors = []

        # The use case here would be someone running a simulation where they want no sequencing errors.
        # No need to run any loops in this case.
        if self.average_error == 0:
            return sequencing_errors
        else:
            for i in range(self.read_length):
                if self.rng.random() < self.quality_score_error_rate[quality_scores[i]]:
                    sequencing_errors.append(i)


        # number_deletions_so_far = 0
        # # don't allow indel errors to occur on subsequent positions
        # previous_indel = -2
        # # don't allow other sequencing errors to occur on bases removed by deletion errors
        # del_blacklist = []
        #
        # for index in sequencing_errors[::-1]:
        #     # determine error type
        #     is_sub = True
        #     # This check checks that we are not trying to insert indels at the end of a read,
        #     # or overlapping with a different indel.
        #     if index != 0 and \
        #             index != self.read_length - 1 - max(self.indel_lengths) and \
        #             abs(index - previous_indel) > 1:
        #         if rng.random() < self.indel_probability:
        #             is_sub = False
        #
        #     # Insert substitution error
        #     if is_sub:
        #         current_nucleotide = read_data[index]
        #         nuc_index = ALLOWED_NUCL.index(current_nucleotide)
        #         # take the zero index because this returns a list of length 1.
        #         new_nucleotide = self.substitution_model[nuc_index].sample()
        #         introduced_errors.append(('S', 1, index, current_nucleotide, new_nucleotide))
        #
        #     # insert indel error:
        #     else:
        #         # Need to take the first element because this returns a list with 1 element.
        #         indel_len = self.indel_length_model.sample()
        #
        #         # insertion error:
        #         if options.rng.random() < self.insertion_probability:
        #             current_nucleotide = read_data[index]
        #             insert = options.rng.choice(ALLOWED_NUCL, p=self.nucleotide_insertion_model, size=indel_len)
        #             new_nucleotide = current_nucleotide + "".join(insert)
        #             introduced_errors.append(('I', len(new_nucleotide) - 1, index, current_nucleotide, new_nucleotide))
        #
        #         elif index < self.read_length - 2 - number_deletions_so_far:
        #             current_nucleotide = read_data[index: index + indel_len + 1]
        #             new_nucleotide = read_data[index]
        #             number_deletions_so_far += len(current_nucleotide) - 1
        #             introduced_errors.append(('D', len(current_nucleotide) - 1, index,
        #                                       current_nucleotide, new_nucleotide))
        #             del_blacklist.extend(list(range(index + 1, index + indel_len + 1)))
        #
        #         previous_indel = index
        #
        # # Remove blacklisted errors
        # for i in range(len(introduced_errors) - 1, -1, -1):
        #     if introduced_errors[i][2] in del_blacklist:
        #         del introduced_errors[i]
        #
        # return out_qualities, introduced_errors
        return sequencing_errors

    def quality_index_remap(self, input_read_length):
        if input_read_length == self.read_length:
            return list(range(self.read_length))
        else:
            # This is basically a way to evenly spread the distribution across the number of bases in the read
            return np.array([max([0, input_read_length * n // self.read_length]) for n in range(self.read_length)])

    def get_quality_scores(self,
                           input_read_length: int) -> np.ndarray[int, ...]:
        """
        Takes a read_length and rng and returns an array of quality scores

        :param input_read_length: The desired length of the quality score array
        :return: An array of quality scores.
        """
        if self.uniform_quality_score:
            return np.array([self.uniform_quality_score] * input_read_length)
        else:
            quality_index_map = self.quality_index_remap(input_read_length)
            temp_qual_array = []
            for i in range(input_read_length):
                score = self.rng.choice(a=self.quality_scores,
                                   p=self.quality_score_probability_matrix[quality_index_map[i]])
                temp_qual_array.append(score)

        if self.rescale_qualities:
            _LOG.warning(f'Quality scores no longer exactly representative of input error probability. '
                         f'Error model scaled to match desired rate...\n'
                         f'Note that this may not have much effect on binned rates.')
            # Note that rescaling won't have much effect on binned quality scores.
            # First we rescale the quality score literals then convert back to phred score (with the 0.5
            # causing borderline cases to take the next highest number).
            rescaled_quals = [max([0, int(-10. * np.log10(self.average_error *
                                                          self.quality_score_error_rate[n]) + 0.5)])
                              for n in temp_qual_array]
            # Now rebin the quality scores.
            temp_qual_array = np.array([take_closest(self.quality_scores, n) for n in rescaled_quals])

        return temp_qual_array


class GcModel:
    """
    This model correlates GC concentration and coverage, within a given window size.
    For example, given a window size of 50, the model shows the bias at each level of
    GC concentration, from 0-50, using the index of the list as the count. For a 50-base
    window with 0 GCs, the coverage bias may be 0.25, at 20 GCs, 1.25. This is based
    on the number of reads per base within that window, on average, across the dataset.

    :param window_size: the size of the sliding window used to measure GC content.
    :param gc_bias: The coverage bias at each GC-count for the window size
    """
    def __init__(self,
                 window_size: int,
                 gc_bias: np.ndarray[float, ...]):
        # assign the model attributes.
        self.window_size = window_size
        self.gc_bias = gc_bias

    def create_coverage_bias_vector(self, sequence: Seq) -> np.ndarray[float, ...]:
        """
        Generates a vector of coverage bias per position, based on the GC concentration of the input sequence.
        :param sequence: A sequence to check coverage for and generate the vector.
        :return: A numpy array of the coverage bias, per position.
        """
        coverage_bias = []
        for i in range(0, len(sequence), self.window_size):
            subsequence = sequence[i: i + self.window_size]
            gc_count = subsequence.count('G') + subsequence.count('C')
            coverage_bias.extend([self.gc_bias[gc_count]] * self.window_size)
        return np.array(coverage_bias)


class FragmentLengthModel:
    """
    A model of the fragment length based on mean and standard deviation of the dataset.

    :param fragment_mean: the mean of the collection of fragment lengths derived from data
    :param fragment_std: the standard deviation of the collection of fragment lengths derived from data
    :param fragment_max: the largest fragment observed in the data
    :param fragment_min: the smallest fragment observed in data
    """
    def __init__(self,
                 fragment_mean: float,
                 fragment_std: float,
                 fragment_max: int,
                 fragment_min: int):
        self.fragment_mean = fragment_mean
        self.fragment_st_dev = fragment_std
        self.fragment_max = fragment_max
        self.fragment_min = fragment_min

    def generate_fragments(self,
                           count: int) -> np.ndarray[int, ...]:
        """
        Generates count number of fragments based on the mean and standard deviation.

        :param count: How many fragments to generate. Some may get filtered out, so the count is not
                      guaranteed.
        :return: A numpy list of fragment random fragment lengths sampled from the model.
        """
        # generates a distribution, assuming normality, then rounds the result and converts to ints
        dist = np.round(self.rng.normal(self.fragment_mean, self.fragment_st_dev, size=count)).astype(int)
        # filter the list to throw out outliers.
        dist = np.array([x for x in dist if self.fragment_min <= x <= self.fragment_max])
        if len(dist) < count/2:
            # This is an edge case, but just in case we had a very weird sample, we'll try to add a few more.
            add_on = self.generate_fragments(count//2)
            np.concatenate((dist, add_on))

        return dist



