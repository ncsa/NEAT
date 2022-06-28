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

from Bio import SeqRecord
from Bio.Seq import Seq

_LOG = logging.getLogger(__name__)


class Read:
    """
    A class representing a particular read

    :param name: name for the read
    :param reference: The reference this read is drawn from
    :param position: First position of the read
    :param end_point: End point of the read
    :param is_reverse: Whether the read is reversed
    """
    def __init__(self,
                 name: str,
                 reference: SeqRecord,
                 position: int,
                 end_point: int,
                 quality_array: np.ndarray[int, ...] = None,
                 mutations: list = None,
                 errors: list = None,
                 is_reverse: bool = False):

        self.name = name
        self.reference = reference
        self.position = position
        self.end = end_point
        self.length = end_point - position
        self.errors = errors
        self.mutations = mutations
        self.quality_array = quality_array
        if is_reverse:
            self.quality_array = self.quality_array[::-1]

    def __repr__(self):
        return f"{self.reference.id}: {self.position}-{self.end}"

    def __str__(self):
        return f"{self.reference.id}: {self.position}-{self.end}"

    def __gt__(self, other):
        if self.reference.id == other.reference.id:
            return self.position > other.position
        else:
            return False

    def __ge__(self, other):
        if self.reference.id == other.reference.id:
            return self.position >= other.position
        else:
            return False

    def __lt__(self, other):
        if self.reference.id == other.reference.id:
            return self.position < other.position
        else:
            return False

    def __le__(self, other):
        if self.reference.id == other.reference.id:
            return self.position <= other.position
        else:
            return False

    def __ne__(self, other):
        if self.reference.id == other.reference.id:
            return self.position != other.position or self.end != other.end
        else:
            return True

    def __eq__(self, other):
        if self.reference.id == other.reference.id:
            return self.position == other.position and self.end == other.end
        else:
            return False

    def __len__(self):
        return self.length

    def get_reference(self):
        return Seq(self.reference[self.position: self.end])

    def get_read(self) -> Seq:
        """
        TODO I don't actually know what needs to happen here yet. It's more complicated than this, for sure. I think
           I need to keep track of how many bases are added, fill in the quality score array for insertions,
           and remove quality score items for deletions.

        :return: the mutated sequence
        """
        ref_seq = self.get_reference()
        mut_seq = copy.deepcopy(str(ref_seq))
        for mut in self.mutations:
            # apply mutations to read
            mut_seq[mut[0]] = mut[2]
        for error in self.errors:
            mut_seq[error[0]] = error[2]

        return Seq(mut_seq)

    def contains(self, test_pos: int):
        return self.position <= test_pos < self.end

    def get_reverse_complement(self):
        return self.get_read().reverse_complement()
