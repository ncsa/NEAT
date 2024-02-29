"""
Handles local files for the larger loop.
"""
__all__ = [
    "write_local_file"
]

import logging

from Bio import SeqRecord
from Bio.Seq import MutableSeq
from pathlib import Path
from copy import deepcopy

from .options import Options
from ...variants import ContigVariants
from ...common import open_output

_LOG = logging.getLogger(__name__)


def write_local_file(vcf_filename: str or Path,
                     variant_data: ContigVariants,
                     reference: SeqRecord,
                     targeted_regions: list,
                     discarded_regions: list,
                     options: Options,
                     fasta_filename: str or Path = None):
    """

    :param vcf_filename: The path to write the local variant file
    :param variant_data: The variant data for this config
    :param reference: The reference for this contig
    :param targeted_regions: Regions to target in this contig
    :param discarded_regions: Regions to discard in this contig
    :param options: The options for this run
    :param fasta_filename: The filename for the optional fasta file
    """

    with open_output(vcf_filename) as tmp:
        filtered_by_target = 0
        filtered_by_discard = 0
        n_added = 0

        local_fasta = LocalFasta(fasta_filename, reference, options)
        for loc in variant_data.variant_locations:
            for variant in variant_data[loc]:
                if not variant.is_input:
                    target_region = find_region(loc, targeted_regions)
                    # For runs without a targeted bed, target_region[2] is 1, so the following will always fail
                    if options.rng.random() >= target_region[2]:
                        _LOG.debug(f'Variant filtered out by target regions bed: {reference.id}: '
                                   f'{variant}')
                        filtered_by_target += 1
                        continue

                    # Now we check the discard regions. All regions are set to a factor of 1 if there is not discard bed
                    discard_region = find_region(loc, discarded_regions)
                    if discard_region[2] == 0:
                        _LOG.debug(f'Variant filtered out by discard regions bed: {reference.id}: '
                                   f'{variant}')
                        filtered_by_discard += 1
                        continue

                ref, alt = variant_data.get_ref_alt(variant, reference)
                variant.metadata["REF"] = ref
                variant.metadata["ALT"] = alt
                sample = variant_data.get_sample_info(variant)

                if options.produce_fasta:
                    local_fasta.add_variant(variant)
                # +1 to position because the VCF uses 1-based coordinates
                line = f"{reference.id}\t" \
                       f"{variant.position1 + 1}\t" \
                       f"{variant_data.generate_field(variant, 'ID')}\t" \
                       f"{ref}\t" \
                       f"{alt}\t" \
                       f"{variant.qual_score}\t" \
                       f"{variant_data.generate_field(variant, 'FILTER')}\t" \
                       f"{variant_data.generate_field(variant, 'INFO')}\t" \
                       f"{variant_data.generate_field(variant, 'FORMAT')}\t" \
                       f"{sample}\n"
                tmp.write(line)
                n_added += 1

    if options.produce_fasta:
        local_fasta.write_fasta()

    if filtered_by_target:
        _LOG.info(f'{filtered_by_target} variants excluded because '
                  f'of target regions with discard off-target enabled')

    if filtered_by_discard:
        _LOG.info(f'{filtered_by_discard} variants excluded because '
                  f'of target regions with discard off-target enabled')

    _LOG.info(f'Finished outputting temp vcf/fasta')

    _LOG.debug(f"Added {n_added} mutations to the reference.")

    return local_fasta.filenames


class LocalFasta:
    """
    Stores info about local fasta file and includes a method to write records.

    :param filename:
    :param reference:
    :param options:
    """
    def __init__(self, filename: Path, reference: SeqRecord, options: Options):
        self.options = options
        mutable_ref_seq = MutableSeq(reference.seq)

        self.filenames = [filename.parent / f"{filename.stem}_{x+1}{filename.suffix}" for x in range(options.ploidy)]
        self.names = [f"{reference.id}_{k+1}" for k in range(options.ploidy)]
        self.mutated_references = [deepcopy(mutable_ref_seq) for _ in range(options.ploidy)]
        self.offset = [0] * options.ploidy

    def add_variant(self, variant):
        """
        Adds a particular variant to the fasta file(s).

        :param variant: The variant to add
        """
        genotype = variant.genotype
        for k in range(len(self.mutated_references)):
            if genotype[k]:
                position = variant.position1 + self.offset[k]
                self.mutated_references[k][position: position+len(variant.metadata['REF'])] = variant.metadata['ALT']
                # offset is a running total of the position modification caused by insertions and
                # deletions. We update it after inserting the variant,
                # so that the next one is in the correct position.
                self.offset[k] += len(variant.metadata['ALT']) - len(variant.metadata['REF'])

    def write_fasta(self):
        """
        Needs to take the input given and make a fasta record. We'll assume if a read is input, it is
        a full read, not a fragment or something. This can either write a list of reads or a chromosome or
        both. If it gets no data, nothing will be written.
        """
        for i in range(len(self.filenames)):
            with open_output(self.filenames[i]) as out_fasta:
                out_fasta.write(f'>{self.names[i]}\n')

                for j in range(0, len(self.mutated_references[i]), 80):
                    print(*self.mutated_references[i][j: j+80], sep="", file=out_fasta, end='\n')


def find_region(location: int, regions_dict: list[tuple[int, int, int | float]]) -> tuple | None:
    """
    This function finds the region that this variant is in.

    :param location: The location of the variant
    :param regions_dict: The dictionary of all regions of interest in tuples as (start, end, factor) where factor
        gives information about retaining the region
    :return: The region containing this location
    """
    for region in regions_dict:
        # Check that this location is valid.
        # Note that bed coordinates are half-open
        if region[0] <= location < region[1]:
            return region

    return None
