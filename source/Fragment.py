from Bio import SeqRecord


class Fragment:
    def __init__(self, reference: SeqRecord, position: int, length: int):

        self.reference = reference
        self.position = position
        self.length = length
        self.errors = None

    def __repr__(self):
        return f"{self.reference.id}: {self.position}-{self.position+self.length}"

    def __str__(self):
        return f"{self.reference.id}: {self.position}-{self.position + self.length}"

    def __gt__(self, other):
        return self.position > other.position

    def __ge__(self, other):
        return self.position >= other.position

    def __lt__(self, other):
        return self.position < other.position

    def __le__(self, other):
        return self.position <= other.position

    def __ne__(self, other):
        return self.position != other.position or self.length != other.length

    def __eq__(self, other):
        return self.position == other.position and self.length == other.length

    def __len__(self):
        return self.length

    def get_fragment(self):
        return self.reference[self.position: self.position + self.length]

    def contains(self, test_pos: int):
        return self.position <= test_pos < self.position + self.length


