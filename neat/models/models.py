"""
Classes for the mutation models used for each variant type to be included,
along with helper functions, used in this simulation. Every Variant type in variants > variant_types
must have a corresponding model in order to be fully implemented.
"""

import re
import logging
import abc

from bisect import bisect_left

import numpy as np
from numpy.random import Generator
from Bio.Seq import Seq
from Bio import SeqRecord

from neat import variants

from ..common import TRINUC_IND, ALLOWED_NUCL, NUC_IND, DINUC_IND
from .default_mutation_model import *
from .default_sequencing_error_model import *
from .default_gc_bias_model import *
from .default_fraglen_model import *
from .utils import bin_scores

__all__ = [
    "MutationModel",
    "SequencingErrorModel",
    "GcModel",
    "FragmentLengthModel",
    "InsertionModel",
    "DeletionModel",
    "SnvModel",
    "ErrorContainer"
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
        self.insert_len_model = insert_len_model
        self.rng = rng

    def get_insertion_length(self, size: int = None) -> int | list[int, ...]:
        """
        Get size number of inserts lengths. Size == 1 results in an int return, else a list of ints.

        :param size: Number of insert lengths to generate. Default is one, which returns an int value.
                     Greater than 1 returns a list of ints.
        :return: int or list of ints.
        """
        return self.rng.choice(a=list(self.insert_len_model),
                               p=[*self.insert_len_model.values()],
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
    _description = "A deletion of N bases"

    def __init__(self,
                 deletion_len_model: dict[int: float, ...],
                 rng: Generator = None):
        self.deletion_len_model = deletion_len_model
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
            raise ValueError("Probabilities do not add up to 1.")

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


class SequencingErrorModel(SnvModel, DeletionModel, InsertionModel):
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
    :param variant_probs: Probability dict for each valid variant type
    :param indel_len_model: Similar to mutation model, but simpler because errors tend not to be complicated. The
        three possible variant types for errors are Insertion, Deletion, and SNV
    :param is_uniform: Some machines use uniform quality scores. This makes simulation a little easier.
    :param rng: optional random number generator. For generating this model, no RNG is needed. But for a run,
        we'll need the rng to perform certain methods.
    """

    def __init__(self,
                 avg_seq_error: float = default_avg_seq_error,
                 read_length: int = default_read_length,
                 transition_matrix: np.ndarray = default_error_transition_matrix,
                 quality_scores: np.ndarray = default_quality_scores,
                 qual_score_probs: np.ndarray = default_qual_score_probs,
                 variant_probs: dict[variants: float] = default_error_variant_probs,
                 indel_len_model: dict[int: float] = default_indel_len_model,
                 insertion_model: np.ndarray = default_insertion_model,
                 rescale_qualities: bool = False,
                 is_uniform: bool = False,
                 rng: Generator = None):

        SnvModel.__init__(self)
        InsertionModel.__init__(self, indel_len_model)
        DeletionModel.__init__(self, indel_len_model)

        self.variant_probs = variant_probs

        self.average_error = avg_seq_error
        self.read_length = read_length
        self.transition_matrix = transition_matrix

        self.quality_scores = quality_scores
        # pre-compute the error rate for each quality score. This is the inverse of the phred score equation
        self.quality_score_error_rate: dict[int, float] = {x: 10. ** (-x / 10) for x in self.quality_scores}
        self.read_length = read_length
        self.quality_score_probabilities = qual_score_probs
        self.rescale_qualities = rescale_qualities
        self.is_uniform = is_uniform
        self.insertion_model = insertion_model
        self.uniform_quality_score = None
        if self.is_uniform:
            # bin scores returns a list, so we need the first (only) element of the list
            converted_avg_err = bin_scores(self.quality_scores,
                                           [int(-10. * np.log10(self.average_error))])[0]
            # Set score to the lowest of the max of the quality scores and the bin closest to the input avg error.
            self.uniform_quality_score = min([max(self.quality_scores), converted_avg_err])
        self.rng = rng

    def get_sequencing_errors(self,
                              read_length: int,
                              reference_segment: SeqRecord,
                              quality_scores: np.ndarray):
        """
        Inserts errors of type substitution, insertion, or deletion into read_data, and assigns a quality score
        based on the container model.
        :param read_length: The length of the read to generate errors for.
        :param reference_segment: The section of the reference from which the read is drawn
        :param quality_scores: Array of quality scores for the read
        :return: modified sequence and associated quality scores

        """

        error_indexes = []
        introduced_errors = []

        # The use case here would be someone running a simulation where they want no sequencing errors.
        # No need to run any loops in this case.
        if self.average_error == 0:
            return introduced_errors
        else:
            for i in range(read_length):
                if self.rng.random() < self.quality_score_error_rate[quality_scores[i]]:
                    error_indexes.append(i)

        num_indels_so_far = 0
        # To prevent deletion collisions
        del_blacklist = []

        for index in error_indexes[::-1]:
            # determine error type. Most will be SNVs
            error_type = SingleNucleotideVariant

            # Not too sure about how realistic it is to model errors as indels, but I'm leaving the code in for now.

            # This is to prevent deletion error collisions and to keep there from being too many indel errors.
            if 0 < index < self.read_length - max(self.deletion_len_model) and num_indels_so_far > self.read_length//2:
                error_type = self.rng.choice(a=list(self.variant_probs), p=list(self.variant_probs.values()))

            # Deletion error
            if error_type == Deletion:
                deletion_length = self.get_deletion_length()
                deletion_reference = reference_segment.seq[index: index + deletion_length + 1]
                deletion_alternate = deletion_reference[0]
                introduced_errors.append(
                    ErrorContainer(Deletion, index, deletion_length, deletion_reference, deletion_alternate)
                )
                num_indels_so_far += deletion_length
                del_blacklist.extend(list(range(index, index + deletion_length)))

            elif error_type == Insertion:
                insertion_length = self.get_insertion_length()
                insertion_reference = reference_segment[index]
                insert_string = ''.join(self.rng.choice(ALLOWED_NUCL, size=insertion_length))
                insertion_alternate = insertion_reference + insert_string
                introduced_errors.append(
                    ErrorContainer(Insertion, index, insertion_length, insertion_reference, insertion_alternate)
                )
                num_indels_so_far += insertion_length

            # Insert substitution error
            # Programmer note: if you add new error types, they can be added as elifs above, leaving the final
            # else dedicated to SNVs.
            else:
                snv_reference = reference_segment[index]
                nuc_index = NUC_IND[snv_reference]
                # take the zero index because this returns a list of length 1.
                snv_alt = self.rng.choice(ALLOWED_NUCL, p=self.transition_matrix[nuc_index])
                introduced_errors.append(
                    ErrorContainer(SingleNucleotideVariant, index, 1, snv_reference, snv_alt)
                )

        # Remove blacklisted errors
        for i in range(len(introduced_errors) - 1, -1, -1):
            if introduced_errors[i].location in del_blacklist:
                del introduced_errors[i]

        return introduced_errors

    def quality_index_remap(self, input_read_length):
        """
        Adjusts the quality map to the suitable read length.

        :param input_read_length: The desired length for the current read.
        :return: An index map from the default read length to the new one.
        """
        if input_read_length == self.read_length:
            return np.arange(self.read_length)
        else:
            # This is basically a way to evenly spread the distribution across the number of bases in the read
            return np.array([max([0, self.read_length * n // input_read_length]) for n in range(input_read_length)])

    def get_quality_scores(self,
                           input_read_length: int) -> np.ndarray:
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
                score = bin_scores(self.quality_scores,
                                   self.rng.normal(self.quality_score_probabilities[i][0],
                                                   scale=self.quality_score_probabilities[i][1]))
                temp_qual_array.append(score)

        if self.rescale_qualities:
            # Note that rescaling won't have much effect on binned quality scores.
            # First we rescale the quality score literals then convert back to phred score (with the 0.5
            # causing borderline cases to take the next highest number).
            rescaled_quals = [max([0, int(-10. * np.log10(self.average_error *
                                                          self.quality_score_error_rate[n]) + 0.5)])
                              for n in temp_qual_array]
            # Now rebin the quality scores.
            temp_qual_array = np.array(bin_scores(self.quality_scores, [n for n in rescaled_quals]))

        return temp_qual_array[:input_read_length]


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


class GcModel:
    """
    This model correlates GC concentration and coverage, within a given window size.
    For example, given a window size of 50, the model shows the bias at each level of
    GC concentration, from 0-50, using the index of the list as the count. For a 50-base
    window with 0 GCs, the coverage bias may be 0.25, at 20 GCs, 1.25. This is based
    on the number of reads per base within that window, on average, across the dataset.

    :param window_size: the size of the sliding window used to measure GC content.
    :param gc_bias: The coverage bias at each GC-count for the window size
    :param coverage: The coverage target for this particular run
    """

    def __init__(self,
                 window_size: int = None,
                 gc_bias: np.ndarray = None,
                 coverage: int = None
                 # May need rng
                 ):
        # assign the model attributes.
        self.window_size = window_size if window_size else default_window_size
        self.gc_bias = gc_bias if gc_bias else default_gc_bias
        self.coverage = coverage

    @property
    def bias_values(self):
        if self.coverage:
            return self.coverage * self.gc_bias
        else:
            return np.ones(len(self.gc_bias))

        # TODO check if these are needed
        # self.mean = self.bias_vector.mean()
        # self.deviation = self.bias_vector.std()

    def create_coverage_bias_vector(self, sequence: Seq) -> np.ndarray:
        """
        Generates a vector of coverage bias per position, based on the GC concentration of the input sequence.
        :param sequence: A sequence to check coverage for and generate the vector.
        :return: A numpy array of the coverage bias, per position.
        """
        bias_values = self.bias_values
        target_vector = [0.0] * len(sequence)
        for i in range(0, len(sequence), self.window_size):
            subsequence = sequence[i: i + self.window_size]
            gc_count = subsequence.count('G') + subsequence.count('C')
            scaling_factor = bias_values[gc_count]
            target_vector[i: i+self.window_size] = [scaling_factor] * self.window_size
        return np.array(target_vector[:len(sequence)])


class FragmentLengthModel:
    """
    A model of the fragment length based on mean and standard deviation of the dataset. Used both
    to generate random fragment lengths and random read lengths. Since a read is essentially a fragment as well,
    and the stastistical models used in NEAT are similar, we'll use fragment to mean read here.

    :param fragment_mean: the mean of the collection of fragment lengths derived from data
    :param fragment_std: the standard deviation of the collection of fragment lengths derived from data
    :param fragment_max: the largest fragment observed in the data
    :param fragment_min: the smallest fragment observed in data
    :param rng: the random number generator for the run
    """

    def __init__(self,
                 fragment_mean: float = None,
                 fragment_std: float = None,
                 fragment_max: int = None,
                 fragment_min: int = None,
                 rng: Generator = None):
        self.fragment_mean = fragment_mean if fragment_mean else default_fragment_mean
        self.fragment_st_dev = fragment_std if fragment_std else default_fragment_std
        self.fragment_max = fragment_max if fragment_max else default_fragment_max
        self.fragment_min = fragment_min if fragment_min else default_fragment_min
        self.rng = rng

    def generate_fragments(self,
                           total_length: int,
                           read_length: int,
                           coverage: int) -> list:
        """
        Generates a number of fragments based on the total length needed, and the mean and standard deviation of the set

        :param total_length: Length of the reference segment we are covering.
        :param read_length: average length of the reads
        :param coverage: the target coverage number
        :return: A list of fragment random fragment lengths sampled from the model.
        """
        # Estimate the number of fragments needed (with a 2x padding)
        number_of_fragments = int(round(total_length / read_length) * (coverage * 2))
        # generates a distribution, assuming normality, then rounds the result and converts to ints
        dist = np.round(self.rng.normal(self.fragment_mean, self.fragment_st_dev, size=number_of_fragments)).astype(int)
        # filter the list to throw out outliers.
        dist = [x for x in dist if self.fragment_min <= x <= self.fragment_max]
        # Just a sanity check to make sure our data isn't too thin:
        while number_of_fragments - len(dist) > 0:
            additional_read = self.rng.normal(loc=self.fragment_mean, scale=self.fragment_st_dev)
            if additional_read > read_length:
                distance_away = read_length - additional_read
                additional_read = read_length - distance_away
            dist.append(round(additional_read))

        return dist
