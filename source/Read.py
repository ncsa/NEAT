from Bio import SeqRecord
from source.Models import Models


class Read:
    def __init__(self, name: str, reference: SeqRecord, position: int, end_point: int):

        self.name = name
        self.reference = reference
        self.position = position
        self.end = end_point
        self.length = end_point - position
        self.errors = None
        self.quality_scores = None

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

    def get_fragment(self):
        return self.reference[self.position: self.end]

    def contains(self, test_pos: int):
        return self.position <= test_pos < self.end

    def generate_errors(self):
        return self

    def generate_quals(self, quality_score_model: Models):
        for i in range(self.length):
