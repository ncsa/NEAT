"""
Class for reads. Each read is described by a name, start position, end position, quality array, mutations, errors,
and whether it is on the reverse strand.

In addition, we attach the reference sequence for later retrieval.

Methods allow comparisons between reads, based on chromosome, start and end. Also, there are methods to retrieve
both the reference sequence and the read and the actual read sequence.
"""
import copy
import logging
import numpy as np

from pathlib import Path
from bisect import bisect_right
from Bio import SeqRecord
from typing import TextIO, Iterator
from Bio.Seq import Seq, MutableSeq

from ...common import open_output
from ...models import SequencingErrorModel, ErrorContainer
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
    :param position: First position of the read
    :param end_point: End point of the read
    :param is_reverse: Whether the read is reversed
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
        self.read_sequence: Seq | None = None
        self.errors: list[ErrorContainer] | None = None
        self.mutations: dict[int: list] | None = None
        self.quality_array: np.ndarray | None = None

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
            self, alternate, location, variant_type: str, model: SequencingErrorModel, quality_score: int = 0
    ):
        """
        This updates the quality score based on the error model. Uniform for mutations, random (but low) for errors
        :param alternate: The alternate sequence for the variant introduced
        :param location: The first position of the mutation, in 0-based coordinates
        :param variant_type: Either "mutations" or "errors"
        :param model: The error model, used to generate any quality scores needed to fill out array
        :param quality_score: The quality score to use, since this has already been calculated for mutations.
            This must be included if type is 'mutation'
        :return: None, updates the quality array in place
        """
        if variant_type == "mutation":
            new_qual = [quality_score] * len(alternate)
        else:
            # Since we have an error here, we'll retain the original quality score and add adjustments as needed
            new_qual = [self.quality_array[location]]
            if len(alternate) > 1:
                # Generate extra scores for insertions
                scores = model.quality_scores
                low_scores = scores[:bisect_right(scores, max(scores)/2)]
                new_qual.extend(model.rng.choice(low_scores, size=len(alternate)-1))

        # Replace the given quality score with the new one
        self.quality_array = \
            self.quality_array[:location] + \
            new_qual + \
            self.quality_array[location+1:]

        if len(self.quality_array) < self.length:
            # Just in case we need to fill out the quality score array
            self.quality_array.append(model.rng.choice(model.quality_scores, size=self.length-len(self.quality_array)))

    def apply_errors(self, mutated_sequence, err_model):
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
            self.update_quality_array(error.alt, error.location, "error", err_model)

        return mutated_sequence

    def apply_mutations(self, mutated_sequence, err_model):
        """
        Applying mutations involves one extra step, because of polyploidism. There may be more than one mutation
        at a given location, so it is formulated as a list. We then pick one at random for this read.

        :param mutated_sequence: The sequence we are applying the variant to
        :param err_model: The error_model for the run, used for the quality score and for the rng
        :return: mutated sequence, with mutations applied
        """
        # Start at the right position based on if this is the forward or reverse read.
        segment_start = int(self.raw_read[2]) if self.is_reverse else int(self.raw_read[0])

        for location in self.mutations:
            # There may be more than one variant to apply, because of polyploidism.
            # So we'll pick one at random for this read, each with equal probability.
            # We may want to add a model to this.
            variant_to_apply = err_model.rng.choice(self.mutations[location])
            # Fetch parameters
            qual_score = variant_to_apply.get_qual_score()
            position = variant_to_apply.get_0_location() - segment_start
            if type(variant_to_apply) == Insertion or type(variant_to_apply) == SingleNucleotideVariant:
                reference_length = 1
                alternate = variant_to_apply.get_alt()
            elif type(variant_to_apply) == Deletion:
                reference_length = variant_to_apply.length
                alternate = mutated_sequence[position]
            else:
                reference_length = variant_to_apply.get_ref_len()
                alternate = mutated_sequence.get_alt()

            # Replace the entire ref with the entire alt
            mutated_sequence = \
                mutated_sequence[:position] + alternate + mutated_sequence[position + reference_length:]

            self.update_quality_array(
                alternate, location, "mutation", err_model, qual_score
            )

        return mutated_sequence

    def get_read(self, error_model) -> Seq:
        """
        Gets mutated sequence to output for fastq/bam

        :param error_model: The error model for the run
        :return: the mutated sequence
        """

        mutated_sequence = self.reference_segment
        if self.mutations:
            mutated_sequence = self.apply_mutations(mutated_sequence, error_model)
        if self.errors:
            mutated_sequence = self.apply_errors(mutated_sequence, error_model)

        mutated_sequence = mutated_sequence[:self.length]
        self.quality_array = self.quality_array[:self.length]

        return Seq(mutated_sequence)

    def contains(self, test_pos: int):
        return self.position <= test_pos < self.end_point

    def calculate_flags(self):
        return 1

    def write_record(
            self, err_model: SequencingErrorModel,
            fastq_handle: TextIO,
            tsam_handle: TextIO,
            produce_fastq: bool,
            produce_tsam: bool
    ):
        """
        Writes the record to the temporary fastq file

        :param err_model: The error model for the run
        :param fastq_handle: the path to the fastq model to write the read
        :param tsam_handle: This will write the corresponding tsam record, assuming the produce_bam option is turned on.
        :param produce_fastq: If true, this will write out the temp fastqs. If false, this will only write out the tsams
            to create the bam files.
        :param produce_tsam: If true, then this will produce the corresponding tsam record, for later bam construction,
            using the same data as in the fastq.
        """

        # Generate quality scores to cover the extended segment. We'll trim later
        self.quality_array = err_model.get_quality_scores(len(self.reference_segment))

        # Get errors for the read
        self.errors = err_model.get_sequencing_errors(self.length, self.reference_segment,
                                                      self.quality_array)

        read_sequence = self.get_read(err_model)

        if self.is_reverse:
            self.quality_array = self.quality_array[::-1]
            read_sequence = read_sequence.reverse_complement()
            temp_position = self.position
            self.position = len(self.reference_segment) - temp_position - self.length
            self.end_point = len(self.reference_segment) - temp_position

        read_quality_score = "".join([chr(x + self.quality_offset) for x in self.quality_array])

        if produce_tsam:
            qname = self.name
            flag = self.calculate_flags()
            rname = self.reference_id
            pos = self.position + 1
            mapq = "".join(self.quality_array)
            fake_cigar = self.reference_segment
            mrnm = "="
            mpos = self.raw_read[2]
            isize = self.raw_read[2] - self.raw_read[1]
            seq = read_sequence
            qual = read_quality_score
            tag = ""
            vtype = ""
            value = ""

            tsam_handle.write('Hello world\n')

        if produce_fastq:
            fastq_handle.write(f'@{self.name}\n')
            fastq_handle.write(f'{str(read_sequence)}\n')
            fastq_handle.write('+\n')
            fastq_handle.write(f'{read_quality_score}\n')
