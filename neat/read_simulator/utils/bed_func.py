"""
This function parses the various bed files that can be input into NEAT. It has a special subtask if the input
is a mutation rate dictionary, which will parse out the input mutation rates associated with each region.
"""
import pathlib
import logging
from ...common import open_input

__all__ = [
    "parse_bed"
]

_LOG = logging.getLogger(__name__)


def parse_bed(input_bed: str,
              chromosomes: list,
              mut_rate_present: bool) -> dict:
    """
    Will parse a bed file, returning a dataframe of chromosomes that are found in the reference,
    and the corresponding positions. Some beds may have mutation in the fourth column. If specified
    then we will also import the mutation data as a column in the dataframe. These dataframes will be
    indexed on the chromosome by row, so 1 row per chromosome.

    :param input_bed: A bed file containing the regions of interest
    :param chromosomes: A list of chromosomes to check
    :param mut_rate_present: A true or false if this bed has mut rate regions included.
    :return: a dictionary of chromosomes: [(pos1, pos2, mutation_rate), etc...]
    """
    ret_dict = {x: {'regions': [], 'rates': []} for x in chromosomes}
    in_bed_only = []

    if input_bed:
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
                        raise
                    # Trying not to 'fix' bed files, but an easy case is if the vcf uses 'chr' and the bed doesn't, in
                    # which we frequently see in human data, we'll just put chr on there. Anything more complicated
                    # will be on the user to correct their beds first.
                    try:
                        assert my_chr in chromosomes
                    except AssertionError:
                        my_chr = f'chr{my_chr}'

                    try:
                        assert my_chr in chromosomes
                    except AssertionError:
                        _LOG.warning("Found chromosome in BED file that isn't in Reference file, skipping")
                        continue
                    # Check to see if this chromosome from the bed is even in the ref.
                    # If not, there's no point in including it
                    in_ref = my_chr in chromosomes
                    if not in_ref:
                        in_bed_only.append(my_chr)
                        continue

                    # If it's in ref, then it's in the dict, so add the positions.
                    if mut_rate_present:
                        # here we append the metadata, if present
                        index = line_list[3].find('mut_rate=')
                        if index == -1:
                            _LOG.warning(f"Skipping mutation rate region: {my_chr}: ({pos1}, {pos2})")
                            _LOG.debug(f'4th column of mutation rate bed must be a semicolon list of key, value '
                                       f'pairs, with one key being mut_rate, e.g., "foo=bar;mut_rate=0.001;do=re".')
                            continue

                        # +9 because that's len('mut_rate='). Whatever is that should be our mutation rate.
                        mut_rate = line_list[3][index + 9:]
                        # We'll trim anything after the mutation rate and call it good. These should be ; separated
                        try:
                            mut_rate = float(mut_rate.split(';')[0])
                        except ValueError:
                            _LOG.error(f"Invalid mutation rate: {my_chr}: ({pos1}, {pos2})")
                            _LOG.debug(f'4th column of mutation rate bed must be a semicolon list of key, value '
                                       f'pairs, with one key being mut_rate, e.g., "foo=bar;mut_rate=0.001;do=re".')
                            raise

                        if mut_rate > 0.3:
                            _LOG.warning("Found a mutation rate > 0.3. This is unusual.")

                        ret_dict[my_chr]['regions'].append((int(pos1), int(pos2)))
                        ret_dict[my_chr]['rates'].append(mut_rate)
                    else:
                        ret_dict[my_chr].append((int(pos1), int(pos2)))

        # some validation
        in_ref_only = [k for k in chromosomes if k not in ret_dict]
        if in_ref_only:
            _LOG.warning(f'Warning: Reference contains sequences not found in BED file {input_bed}. '
                         f'These chromosomes will be omitted from the outputs.')
            _LOG.debug(f"In reference only regions: {in_ref_only}")

        if in_bed_only:
            _LOG.warning(f'BED file {input_bed} contains sequence names '
                         f'not found in reference. These regions will be ignored.')
            _LOG.debug(f'Regions ignored: {in_bed_only}')

    # Returns an empty dict if there is no file to process.
    return ret_dict
