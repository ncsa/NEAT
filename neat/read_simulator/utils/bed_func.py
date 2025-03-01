"""
This function parses the various bed files that can be input into NEAT. It has a special subtask if the input
is a mutation rate dictionary, which will parse out the input mutation rates associated with each region.
"""
import logging
import pathlib
import sys

from Bio.File import _IndexedSeqFileDict

from .options import Options
from ...common import open_input

__all__ = [
    "parse_beds",
    "fill_out_bed_dict"
]

_LOG = logging.getLogger(__name__)


def parse_beds(options: Options, ref_keys_counts: dict, average_mutation_rate: float) -> list:
    """
    This single function parses the three possible bed file types for NEAT.

    :param options: The options object for this run
    :param ref_keys_counts: A dictionary containing the reference keys and the associated lengths.
    :param average_mutation_rate: The average mutation rate from the model or user input for this run
    :return: target_dict, discard_dict, mutation_rate_dict
    """

    return_list = []
    bed_files = (options.target_bed, options.discard_bed, options.mutation_bed)
    # These numbers indicate the factors for target_bed, discard_bed, and mutation_bed, respectively
    processing_structure = {
        0: ("target",
            # For target dict, True will indicate region is to be included (either within a target bed or
            # if there was no target bed entered.
            # For example (0, 10, False), (10, 100, True) would indicade that 0-10 fell outside a target region and will
            # be discarded, whereas 10-100 is within a targeted region and will be included. A single tuple of
            # (0, chromosome_end, True) would indicate either the entire chromosome was targeted or there was no bed
            # file. Default is to target all regions off all contigs.
            True,
            bool(options.target_bed)),
        1: ("discard",
            # For the discard dict, false indicates the region is not within a discard bed (i.e., it should be included
            # in the final results), whereas True indicates a region to be discarded, according to the bed.
            # For example (0, 10, False), (10, 100, True) would mean we keep reads from 0-10 and discard the
            # ones from 10-100. Default is that there are no discard regions and we keep everything.
            False,
            bool(options.discard_bed)),
        2: ("mutation",
            average_mutation_rate,  # Default value will be this for all sections
            bool(options.mutation_bed))
    }

    for i in range(len(bed_files)):
        temp_bed_dict = parse_single_bed(
            bed_files[i],
            ref_keys_counts,
            processing_structure[i]
        )
        # If there was a bed this function will fill in any gaps the bed might have missed, so we have a complete map
        # of each chromosome. The uniform case in parse_single_bed handles this in cases of no input bed.
        if processing_structure[i][2]:
            final_dict = fill_out_bed_dict(ref_keys_counts, temp_bed_dict, processing_structure[i])
        else:
            final_dict = temp_bed_dict

        return_list.append(final_dict)

    return return_list


def parse_single_bed(input_bed: str,
                     ref_keys_counts: dict,
                     process_key: tuple[str, float, bool],
                     ) -> dict:
    """
    Will parse a bed file, returning a dataframe of chromosomes that are found in the reference,
    and the corresponding positions. Some beds may have mutation in the fourth column. If specified
    then we will also import the mutation data.

    :param input_bed: A bed file containing the regions of interest
    :param ref_keys_counts: A dict of chromosomes to check (key), along with their lengths (value)
    :param process_key: Indicates which bed we are processing and the factor we need to process it.
    :return: a dictionary of chromosomes: [(pos1, pos2, factor), etc...]
    """
    ret_dict = {x: [] for x in ref_keys_counts}
    in_bed_only = []
    key_pool = list(ref_keys_counts.keys())

    # process_key[2] is True when there is an input bed, False otherwise.
    if process_key[2]:
        # Pathlib will help us stay machine-agnostic to the degree possible
        input_bed = pathlib.Path(input_bed)
        with open_input(input_bed) as f:
            for line in f:
                if not line.startswith(('@', '#', "\n")):
                    # Note: on targeted and discard regions, we really only need chrom, pos1, and pos2
                    # But for the mutation rate regions, we need a fourth column of meta_data,
                    # So we have to check them all, though it won't be used for targeted and discard regions
                    line_list = line.strip().split('\t')[:4]
                    try:
                        [my_chr, pos1, pos2] = line_list[:3]
                    except ValueError:
                        _LOG.error(f"Improperly formatted bed file line {line}")
                        sys.exit(1)
                    # Bed file chromosome names must match the reference.
                    try:
                        assert my_chr in key_pool
                    except AssertionError:
                        in_bed_only.append(my_chr)
                        _LOG.warning(
                            f"Found chromosome in BED file that isn't in Reference file, skipping {my_chr}"
                        )
                        continue
                    # If this is the mutation bed, we need to process an extra column.
                    if process_key[0] == "mutation":
                        # here we append the metadata, if present
                        index = line_list[3].find('mut_rate=')
                        # Improperly formatted mutation rate bed file
                        if index == -1:
                            _LOG.error(f"Invalid mutation rate: {my_chr}: ({pos1}, {pos2})")
                            _LOG.error('4th column of mutation rate bed must be a semicolon list of key, value '
                                       'pairs, with one key being mut_rate, e.g., "foo=bar;mut_rate=0.001;do=re".')
                            sys.exit(1)

                        # +9 because that's len('mut_rate='). Whatever is that should be our mutation rate.
                        mut_rate = line_list[3][index + 9:]
                        # We'll trim anything after the mutation rate and call it good. These should be ; separated
                        try:
                            mut_rate = float(mut_rate.split(';')[0])
                        except ValueError:
                            _LOG.error(f"Invalid mutation rate: {my_chr}: ({pos1}, {pos2})")
                            _LOG.error('4th column of mutation rate bed must be a semicolon list of key, value '
                                       'pairs, with one key being mut_rate, e.g., "foo=bar;mut_rate=0.001;do=re".')
                            sys.exit(1)

                        if mut_rate > 0.3:
                            _LOG.warning("Found a mutation rate > 0.3. This is unusual.")

                        ret_dict[my_chr].append((int(pos1), int(pos2), mut_rate))
                    elif process_key[0] == "discard":
                        # True here means a discard bed was present and this region was marked for discard.
                        ret_dict[my_chr].append((int(pos1), int(pos2), True))
                    else:  # target dict
                        # True here means that this section is to be included (either within a target bed or because
                        # no targeted regions were defined)
                        ret_dict[my_chr].append((int(pos1), int(pos2), True))

        # some validation
        in_ref_only = [k for k in ref_keys_counts if k not in ret_dict]
        if in_ref_only:
            _LOG.warning(f'Warning: Reference contains sequences not found in BED file {input_bed}. '
                         'These chromosomes will be omitted from the outputs.')
            _LOG.debug(f"In reference only regions: {in_ref_only}")

        if in_bed_only:
            _LOG.warning(f'BED file {input_bed} contains sequence names '
                         'not found in reference. These regions will be ignored.')
            _LOG.debug(f'Regions ignored: {list(set(in_bed_only))}')

    else:
        for my_chr in ref_keys_counts:
            ret_dict[my_chr].append((0, ref_keys_counts[my_chr], process_key[1]))

    return ret_dict


def fill_out_bed_dict(ref_keys_counts: dict,
                      region_dict: dict | None,
                      default_factor: tuple[str, bool, bool],
                      ) -> dict:
    """
    This parses the dict derived from the bed file and fills in any gaps, so it can be more easily cycled through
    later.

    The input to this function is the dict for a single chromosome.

    :param ref_keys_counts: A dictionary with the keys as the contigs from the reference dictionary, and the values as
        the length of the contig. Needed to calculate max values.
    :param region_dict: A dict based on the input from the bed file, with keys being all the chromosomes
        present in the reference.
    :param default_factor: This is a tuple representing the type and any extra information. The extra info is for
        mutation rates. If beds were input, then these values come from the bed, else they are set to one value
        across the contig
    :return: A tuple with (start, end, some_factor) for each region in the genome.
    """

    ret_dict = {}

    for contig in region_dict:
        regions = []
        max_value = ref_keys_counts[contig]
        start = 0

        # The trivial case of no data yet for this contig the region gets filled out as 0 to len(contig_sequence)
        # with the default value applied to the whole region
        if not region_dict[contig]:
            regions.append((start, max_value, default_factor[1]))
        else:
            for region in region_dict[contig]:
                if region[0] > start:
                    if default_factor[0] == "discard" or default_factor[0] == "mutation":
                        regions.append((start, region[0], default_factor[1]))
                    else:  # target regions get assigned "False" values outside the targets, if bed is present
                        regions.append((start, region[0], False))
                    start = region[1]
                    regions.append(region)
                elif region[0] == start:
                    regions.append(region)
                    start = region[1]

            # If the region ends short of the end, this fills in to the end
            if regions[-1][1] != max_value:
                if default_factor[0] == "discard" or default_factor[0] == "mutation":
                    regions.append((start, max_value, default_factor[1]))
                else:  # target regions get assigned "False" values outside the targets, if bed is present
                    regions.append((start, max_value, False))

        ret_dict[contig] = regions

    return ret_dict
