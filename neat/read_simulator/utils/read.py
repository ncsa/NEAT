"""
Class for reads. Each read is described by a name, start position, end position, quality array, mutations, errors,
and whether it is on the reverse strand.

In addition, we attach the reference sequence for later retrieval.

Methods allow comparisons between reads, based on chromosome, start and end. Also, there are methods to retrieve
both the reference sequence and the read and the actual read sequence.
"""
import logging
import numpy as np

from bisect import bisect_right
from typing import TextIO
from Bio.Seq import Seq, MutableSeq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from numpy.random import Generator

from ...common import ALLOWED_NUCL
from ...models import SequencingErrorModel, ErrorContainer, MutationModel
from ...variants import SingleNucleotideVariant, Insertion, Deletion, UnknownVariant

_LOG = logging.getLogger(__name__)


class Read:
    """
    A class representing a particular read

    :param name: name for the read
    :param raw_read: gives the read as a tuple of coordinates, e.g., (1, 110, 111, 125)
        Which gives information about the left and paired right read, or gives 0,0 for missing component
    :param reference_segment: The reference segment this read is drawn from
    :param reference_id: the ID for the reference where the segment is drawn from
    :param position: First position of the read, relative to the reference
    :param end_point: End point of the read, relative to the reference
    :param is_reverse: Whether the read is reversed
    :param is_paired: Whether this read has a proper pair.
    """
    def __init__(self,
                 name: str,
                 raw_read: tuple,
                 reference_segment: Seq,
                 reference_id: str,
                 position: int,
                 end_point: int,
                 quality_offset: int,
                 is_reverse: bool = False,
                 is_paired: bool = False):

        self.name = name
        self.raw_read = raw_read
        self.reference_segment = reference_segment
        self.reference_id = reference_id
        self.position = position
        self.end_point = end_point
        self.quality_offset = quality_offset
        self.length = end_point - position
        self.is_reverse = is_reverse
        self.is_paired = is_paired

        # These parameters won't be initialized on the initial read creation, but will
        # be generated by functions below.
        self.read_sequence: Seq = Seq("")  # initialize to empty sequence
        self.errors: list[ErrorContainer] = []  # initialize
        self.mutations: dict[int: list] = {}  # initialize
        self.quality_array: np.ndarray = np.zeros(self.length)  # this will have the correct memory length
        self.mapping_quality: int = 0  # initialize at 0
        self.read_quality_string: str = ""  # This will hold the read quality string

    def __repr__(self):
        return f"{self.reference_id}: {self.position}-{self.end_point}"

    def __str__(self):
        return f"{self.reference_id}: {self.position}-{self.end_point}"

    def __gt__(self, other):
        if self.reference_id == other.reference_id:
            return self.position > other.position
        else:
            return False

    def __ge__(self, other):
        if self.reference_id == other.reference_id:
            return self.position >= other.position
        else:
            return False

    def __lt__(self, other):
        if self.reference_id == other.reference_id:
            return self.position < other.position
        else:
            return False

    def __le__(self, other):
        if self.reference_id == other.reference_id:
            return self.position <= other.position
        else:
            return False

    def __ne__(self, other):
        if self.reference_id == other.reference_id:
            return self.position != other.position or self.end_point != other.end_point
        else:
            return True

    def __eq__(self, other):
        if self.reference_id == other.reference_id:
            return self.position == other.position and self.end_point == other.end_point
        else:
            return False

    def __len__(self):
        return self.length

    def update_quality_array(
            self,
            reference: Seq,
            alternate: Seq,
            location: int,
            variant_type: str,
            quality_scores: list,
            quality_score: int = 0
    ):
        """
        This updates the quality score based on the error model. Uniform for mutations, random (but low) for errors
        :param reference: The reference sequence for the variant
        :param alternate: The alternate sequence for the variant
        :param location: The first position of the mutation, in 0-based coordinates
        :param variant_type: Either "mutations" or "errors"
        :param quality_scores: The possible quality scores, used to adjust for mutations and errors
        :param rng: The random number generator for the run
        :param quality_score: The quality score to use, since this has already been calculated for mutations.
            This must be included if type is 'mutation'
        :return: None, updates the quality array in place
        """
        if variant_type == "mutation":
            new_quality_score = [quality_score] * len(alternate)
        else:
            # Since we have an error here, we'll choose a min score
            original_quality_score = self.quality_array[location: location+len(reference)]
            new_quality_score = original_quality_score.copy()
            low_score = min(quality_scores)
            # If is insertion
            if len(alternate) > 1:
                # Generate extra scores for insertions
                # Original ref is unaffected, so it's quality score remains the same
                new_quality_score = [low_score] * (len(alternate) - 1)
            # If is deletion
            elif len(reference) > 1 and len(alternate) == 1:
                new_quality_score = []
            # SNP
            else:
                new_quality_score = [low_score]

        # Replace the given quality score with the new one
        self.quality_array = \
            np.concatenate((self.quality_array[:location],
                            np.array(new_quality_score),
                            self.quality_array[location+len(reference):]))

    def apply_errors(self, mutated_sequence: MutableSeq, err_model: SequencingErrorModel):
        """
        This function applies errors to a sequence and calls the update_quality_array function after

        :param mutated_sequence: The sequence to add errors to.
        :param err_model: The error model for this run,
        :return: None, The sequence, with errors applied
        """
        for error in self.errors:
            # Replace the entire ref sequence with the entire alt sequence
            mutated_sequence = \
                mutated_sequence[:error.location] + error.alt + mutated_sequence[error.location+len(error.ref):]
            # update quality score for error
            self.update_quality_array(
                error.ref,
                error.alt,
                error.location,
                "error",
                err_model.quality_scores,
            )

        return mutated_sequence

    def apply_mutations(self, mutated_sequence: MutableSeq, quality_scores: list, mut_model: MutationModel):
        """
        Applying mutations involves one extra step, because of polyploidism. There may be more than one mutation
        at a given location, so it is formulated as a list. We then pick one at random for this read.

        :param mutated_sequence: The sequence we are applying the variant to
        :param quality_scores: The possible quality scores for this run (to update quality scores)
        :param mut_model: The mutation model for the run, used for the zygosity and for the rng
        :return: mutated sequence, with mutations applied
        """

        for location in self.mutations:
            variant_to_apply = mut_model.rng.choice(self.mutations[location])

            # Fetch parameters
            qual_score = variant_to_apply.get_qual_score()
            # Find the position within the read for this variant, cast as a python int, instead of a numpy int
            position = int(variant_to_apply.get_0_location() - self.position)
            # Figure out if the variant is in this read. If 'to_mutate' selects any 1, then it is mutated.
            to_mutate = mut_model.rng.choice(variant_to_apply.genotype)
            # If a 1 was selected, then apply the variant, else return the original sequence
            if to_mutate:
                if type(variant_to_apply) == Insertion or type(variant_to_apply) == SingleNucleotideVariant:
                    reference_length = 1
                    alternate = variant_to_apply.get_alt()
                elif type(variant_to_apply) == Deletion:
                    reference_length = variant_to_apply.length
                    try:
                        alternate = mutated_sequence[position]
                    except:
                        print(f"sequence: {mutated_sequence}")
                        print(f"position: {position}")
                        raise
                else:
                    reference_length = variant_to_apply.get_ref_len()
                    alternate = mutated_sequence.get_alt()

                # Replace the entire ref with the entire alt
                mutated_sequence = mutated_sequence[:position] + alternate + mutated_sequence[position + reference_length:]

                self.update_quality_array(
                    self.reference_segment[location: location+reference_length],
                    alternate,
                    location,
                    "mutation",
                    quality_scores,
                    qual_score
                )

        return mutated_sequence

    def apply_variants_for_final_output(
            self,
            error_model: SequencingErrorModel,
            mutation_model: MutationModel
    ) -> Seq:
        """
        Gets mutated sequence to output for fastq/bam

        :param error_model: The error model for the run
        :param mutation_model: The mutation model for the run
        :return: the mutated sequence
        """

        rng = mutation_model.rng
        # The read sequence should have some padding still at this point
        read_sequence = MutableSeq(self.read_sequence.seq)
        # Deal with N's
        for i in range(len(read_sequence)):
            if read_sequence[i] == 'N':
                # for now replace 'N' masked regions with pure noise. We may refine this in time
                read_sequence[i] = rng.choice(['A', 'C', 'G', 'T'])
                self.quality_array[i] = min(error_model.quality_scores)

        if self.mutations:
            read_sequence = self.apply_mutations(read_sequence, list(error_model.quality_scores), mutation_model)
        if self.errors:
            read_sequence = self.apply_errors(read_sequence, error_model)

        read_sequence = read_sequence[:self.length]
        self.quality_array = self.quality_array[:self.length]

        return Seq(read_sequence)

    def contains(self, test_pos: int):
        return self.position <= test_pos < self.end_point

    def calculate_flags(self, paired_ended_run):
        """
        Calculates the flags for the read

        :param paired_ended_run: Whether the entire run was done in paired-ended mode
        """
        flag = 0
        if paired_ended_run:
            flag += 1
            if self.is_paired:
                # Whether this read has a mapped pair (proper pair)
                flag += 2
            # Flag 4 is if this read is unmapped, which isn't a situation we'll encounter in the simulation
            if not self.is_paired:
                # If this read's mate is unmapped (usually because it was off the end)
                flag += 8
            if self.is_reverse:
                # Flag to indicate this is the reverse strand
                flag += 16
            elif self.is_paired:
                # Flag that indicates that the mate is reversed
                # (which is always the case, if it exists, in this simulation)
                flag += 32
            if not self.is_reverse and self.is_paired:
                flag += 64
            if self.is_reverse and self.is_paired:
                flag += 128
            # None of the other potential samflags are relevant to this simulation
        return flag

    def finalize_read_and_write(
            self,
            err_model: SequencingErrorModel,
            mut_model: MutationModel,
            fastq_handle: TextIO,
            produce_fastq: bool,
    ):
        """
        Writes the record to the temporary fastq file

        :param err_model: The error model for the run
        :param mut_model: The mutation model for the run
        :param fastq_handle: the path to the fastq model to write the read
        :param produce_fastq: If true, this will write out the temp fastqs. If false, this will only write out the tsams
            to create the bam files.
        """

        # Generate quality scores to cover the extended segment. We'll trim later
        self.quality_array = err_model.get_quality_scores(len(self.reference_segment))

        if self.is_reverse:
            read_quality_array = self.quality_array[::-1]
            read_sequence = self.reference_segment.reverse_complement()
        else:
            read_quality_array = self.quality_array
            read_sequence = self.reference_segment

        # Get errors for the read and update the quality score
        self.errors = err_model.get_sequencing_errors(self.length, read_sequence,
                                                      read_quality_array)

        # Update the read as needed:
        self.quality_array = read_quality_array
        self.read_sequence = read_sequence

        self.read_sequence = self.apply_variants_for_final_output(err_model, mut_model)

        self.read_quality_string = "".join([chr(x + self.quality_offset) for x in self.quality_array])
        self.mapping_quality = err_model.rng.poisson(70)

        if produce_fastq:
            fastq_handle.write(f'@{self.name}\n')
            fastq_handle.write(f'{str(self.read_sequence)}\n')
            fastq_handle.write('+\n')
            fastq_handle.write(f'{self.read_quality_string}\n')

    def make_cigar(self):
        """
        Aligns the reference and mutated sequences.
        """

        # These parameters were set to minimize breaks in the mutated sequence and find the best
        # alignment from there.
        raw_alignment = pairwise2.align.localms(
            self.reference_segment, self.read_sequence, match=1, mismatch=-1, open=-0.5, extend=-0.1,
            penalize_extend_when_opening=True
        )

        alignment = format_alignment(*raw_alignment[0], full_sequences=True).split()
        aligned_template_seq = alignment[0]
        aligned_mut_seq = alignment[-2]
        cig_count = 0
        curr_char = ''
        cig_string = ''
        # Find first match
        start = min([aligned_mut_seq.find(x) for x in ALLOWED_NUCL])
        for char in range(start, start + len(self.read_sequence)):
            if aligned_template_seq[char] == '-':  # insertion
                if curr_char == 'I':  # more insertions
                    cig_count = cig_count + 1
                else:  # new insertion
                    cig_string = cig_string + str(cig_count) + curr_char
                    curr_char = 'I'
                    cig_count = 1
            elif aligned_mut_seq[char] == '-':  # deletion
                if curr_char == 'D':  # more deletions
                    cig_count = cig_count + 1
                else:  # new deletion
                    cig_string = cig_string + str(cig_count) + curr_char
                    curr_char = 'D'
                    cig_count = 1
            else:  # match
                if curr_char == 'M':  # more matches
                    cig_count = cig_count + 1
                else:  # new match
                    # If there is anything before this, add it to the string and increment,
                    # else, just increment
                    if not cig_count == 0:
                        cig_string = cig_string + str(cig_count) + curr_char
                    curr_char = 'M'
                    cig_count = 1
        return cig_string + str(cig_count) + curr_char

    def get_mpos(self):
        """
        Get the mate position of the read
        """
        if self.is_paired:
            if self.is_reverse:
                return self.raw_read[0]
            else:
                return self.raw_read[2]
        else:
            return 0

    def get_tlen(self):
        """
        Get the template length for the read
        """
        if self.is_paired:
            length = self.raw_read[3] - self.raw_read[0] + 1
            if length < 0:
                return 0
            else:
                if self.is_reverse:
                    return -length
                else:
                    return length
        else:
            return 0
