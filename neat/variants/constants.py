from .insertion import Insertion
from .deletion import Deletion
from .single_nucleotide_variant import SingleNucleotideVariant
from .unknown_variant import UnknownVariant

VariantTypes = [Insertion, Deletion, SingleNucleotideVariant, UnknownVariant]
