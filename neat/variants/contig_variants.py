"""
The object storing variants for the current contig.
"""

__all__ = [
    "ContigVariants"
]

import numpy as np
import logging

from bisect import bisect
from Bio.SeqRecord import SeqRecord

from .base_variant import BaseVariant
from ..common import get_genotype_string
from .constants import Insertion, Deletion, SingleNucleotideVariant

_LOG = logging.getLogger(__name__)


class ContigVariants:
    """
    This class will keep track of all the variants in the contig, both those added by the vcf
    and those added by simulation. We output each contig after processing, so a larger scope is not required.

    The dictionary variants_dict will be keyed by locations, so that variants can be output in sorted order,
    as determined by the sorted_locations variable. The genotype is the overall genotype for each key location.

    Using init here resets all parameters for each contig, even though there is no input.
    """
    def __init__(self):
        self.variant_locations: list[int] = []
        self.contig_variants: dict[int, list] = {}
        self.default_dict = {
            "ID": '.',
            "FILTER": "PASS",
            "INFO": '.',
            "FORMAT": 'GT',
        }
        self.all_dels = []
        self.all_ins = []

    def check_if_del(self, other):
        for deletion in self.all_dels:
            if np.array_equal(other.genotype, deletion.genotype) and deletion.contains(other):
                return deletion
        return None

    def check_if_ins(self, other):
        for insert in self.all_ins:
            if np.array_equal(other.genotype, insert.genotype) and insert.contains(other):
                return insert
        return None

    def add_location(self, input_location):
        if input_location not in self.variant_locations:
            self.variant_locations.insert(bisect(self.variant_locations, input_location), input_location)
            self.contig_variants[input_location] = []

    def add_variant(self, input_variant):
        if input_variant.position1 not in self.variant_locations:
            self.add_location(input_variant.position1)
        else:
            if self.find_dups(input_variant):
                return 1
        if type(input_variant) == Insertion:
            self.all_ins.append(input_variant)
        elif type(input_variant) == Deletion:
            self.all_dels.append(input_variant)

        self.contig_variants[input_variant.position1].append(input_variant)
        return 0

    def compile_genotypes_for_location(self, location) -> np.ndarray:
        variants_of_interest = self.contig_variants[location]
        ret_genotype = np.zeros(len(variants_of_interest[0].genotype))
        for item in variants_of_interest:
            ret_genotype[np.where(item.genotype == 1)] = 1
        return ret_genotype

    def generate_field(self, variant, field):
        """
        Generates a metadata field based of defaults or inputs.

        :param variant: The variant to generate the field for
        :param field: The field to try to get
        :return: the field from the metadata or the default value for this run.
        """
        try:
            return variant.metadata[field]
        except KeyError:
            return self.default_dict[field]

    def find_dups(self, variant):
        """
        Checks if the given genotype is already present in a list of variants.

        :param variant: A variant to check for duplicates
        :return: True or False if found or not
        """
        for existing_var in self.contig_variants[variant.position1]:
            if np.array_equal(variant.genotype, existing_var.genotype):
                return True

        return False

    @staticmethod
    def get_ref_alt(variant: BaseVariant, reference: SeqRecord, block_start: int) -> tuple[str, str]:
        """
        Turns variants at a location into a vcf output. Only outputs items it has data for,
        but in the correct vcf format and order. Note that for a VCF, Insertion and Deletion
        type variant positions start at the first unaltered base, so that the alternate or ref is
        not empty. Also, VCF numbering is 1-based rather than 0-based. For an insertion, for example
        this means the location is location - 1 + 1, to account for both of the previous facts.

        We may improve this later, but for now, we will output one variant per line.

        :param variant: The variant to retrieve ref and alternate for
        :param reference: The reference SeqRecord for this contig
        :param block_start: Where, in relation to the overall contig, this SeqRecord starts
        :return: A list of vcf lines based on the variants for that location
        """
        ref: str
        alt: str

        position1 = variant.position1 - block_start
        _LOG.debug(f"Getting ref and alt from {variant.position1}")
        if type(variant) == Insertion:
            ref = reference[position1]
            alt = variant.get_alt()
        elif type(variant) == Deletion:
            ref = str(reference[position1: position1 + variant.length].seq)
            alt = reference[position1]
        elif type(variant) == SingleNucleotideVariant:
            ref = reference[position1]
            alt = variant.get_alt()
        else:
            # Unknown types must have an explicit ref and alternate in metadata
            ref = variant.metadata['REF']
            alt = variant.get_alt()

        return ref, alt

    @staticmethod
    def get_sample_info(variant):
        try:
            return variant.metadata["NEAT_sample"]
        except KeyError:
            return get_genotype_string(variant.genotype)

    def remove_variant(self, variant):
        if variant.position in self.variant_locations:
            if variant in self.contig_variants[variant.position]:
                self.contig_variants[variant.position].remove(variant)
            if not self.contig_variants[variant.position]:
                self.variant_locations.remove(variant.position)

    def __getitem__(self, input_location: int) -> list:
        """
        Fetches all variants at input location
        :param input_location: The location to retrieve variants from
        :return: list of variants at that location
        """
        return self.contig_variants[input_location]

    def __contains__(self, item):
        return item in self.contig_variants
