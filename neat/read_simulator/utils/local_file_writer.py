"""
Handles local files for the larger loop.
"""
__all__ = [
    "write_local_file"
]

import logging
import time

from Bio import SeqRecord
from Bio.Seq import MutableSeq
from pathlib import Path
from copy import deepcopy

from .options import Options
from ...variants import ContigVariants
from ...common import open_output

_LOG = logging.getLogger(__name__)


def write_local_file(
        vcf_filename: str or Path,
        variant_data: ContigVariants,
        reference: SeqRecord,
        targeted_regions: list,
        discarded_regions: list
):
    """

    :param vcf_filename: The path to write the local variant file
    :param variant_data: The variant data for this config
    :param reference: The reference for this contig
    :param targeted_regions: Regions to target in this contig
    :param discarded_regions: Regions to discard in this contig
    """

    with open_output(vcf_filename) as tmp:
        t = time.time()
        _LOG.info("Writing temp files")
        filtered_by_target = 0
        filtered_by_discard = 0
        n_added = 0

        for loc in variant_data.variant_locations:
            for variant in variant_data[loc]:
                if not variant.is_input:
                    # Todo make sure this section still works
                    target_region = find_region(loc, targeted_regions)
                    # This will be True in targeted regions, if a bed is present, or everywhere if not bed is present.
                    # Anything outside defined target regions will be marked false and this `if` will activate.
                    if not target_region[2]:
                        _LOG.debug(f'Variant filtered out by target regions bed: {reference.id}: '
                                   f'{variant}')
                        filtered_by_target += 1
                        continue

                    discard_region = find_region(loc, discarded_regions)
                    # This will be True if the discard bed was present and this region was within a discard bed region.
                    if discard_region[2]:
                        _LOG.debug(f'Variant filtered out by discard regions bed: {reference.id}: '
                                   f'{variant}')
                        filtered_by_discard += 1
                        continue

                ref, alt = variant_data.get_ref_alt(variant, reference)
                variant.metadata["REF"] = ref
                variant.metadata["ALT"] = alt
                sample = variant_data.get_sample_info(variant)

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

    if filtered_by_target > 0:
        _LOG.info(f'{filtered_by_target} variants excluded because '
                  f'of target regions with discard off-target enabled')

    if filtered_by_discard > 0:
        _LOG.info(f'{filtered_by_discard} variants excluded because '
                  f'of target regions with discard off-target enabled')

    _LOG.info(f'Finished outputting temp files')

    _LOG.info(f"Added {n_added} mutations in {(time.time() - t) / 60:.2f} m.")


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
