from error_handling import log_mssg, premature_exit
from Bio import SeqRecord


class Fragment:
    def __init__(self, reference: SeqRecord, position, length):
        self.reference = reference
        self.position = position
        self.length = length

    def get_fragment(self):
        yield self.reference[self.position: self.position + self.length]
