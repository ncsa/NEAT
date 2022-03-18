import gzip
import pathlib
import re

from source.error_handling import premature_exit, log_mssg


def parse_bed(input_bed: str, chromosomes: list,
              begins_with_chr: bool,
              mut_rate_present: bool) -> dict:
    """
    Will parse a bed file, returning a dataframe of chromosomes that are found in the reference,
    and the corresponding positions. Some beds may have mutation in the fourth column. If specified
    then we will also import the mutation data as a column in the dataframe. These dataframes will be
    indexed on the chromosome by row, so 1 row per chromosome.
    :param input_bed: A bed file containing the regions of interest
    :param chromosomes: A list of chromosomes to check
    :param begins_with_chr: A true or false if this reference has human like chromosomes (chr1, chr2, etc), used because
    it's a common thing with human datasets to either use or not use "chr" before the number and
    no program seems to agree on what to do with this fact. We'll check for both cases, if chroms in the ref
    begin with 'chr'
    :param mut_rate_present: A true or false if this bed has mut rate regions included.
    :return: a dictionary of chromosomes: [(pos1, pos2, mutation_rate), etc...]
    """
    ret_dict = {x: [] for x in chromosomes}
    in_bed_only = []

    if input_bed:
        # Pathlib will help us stay machine-agnostic to the degree possible
        input_bed = pathlib.Path(input_bed)
        if input_bed.suffix == '.gz':
            f = gzip.open(input_bed, 'r')
        else:
            f = open(input_bed, 'r')
        try:
            for line in f:
                if not line.startswith(('@', '#', "\n")):
                    # Note: on targeted and discard regions, we really only need chrom, pos1, and pos2
                    # But for the mutation rate regions, we need a fourth column of meta_data,
                    # So we have to check them all, though it won't be used for targeted and discard regions
                    line_list = line.strip().split('\t')[:4]
                    try:
                        [my_chr, pos1, pos2] = line_list[:3]
                    except ValueError:
                        log_mssg(f"Improperly formatted bed file line {line}", "error")
                        premature_exit(1)
                    # Trying not to 'fix' bed files, but an easy case is if the vcf uses 'chr' and the bed doesn't, in
                    # which we frequently see in human data, we'll just put chr on there. Anything more complicated
                    # will be on the user to correct their beds first.
                    if begins_with_chr and not my_chr.startswith('chr'):
                        my_chr = 'chr' + my_chr
                    # Check to see if this chromosome from the bed is even in the ref.
                    # If not, there's not point in including it
                    in_ref = my_chr in chromosomes
                    if not in_ref:
                        in_bed_only.append(my_chr)
                        continue

                    # If it's in ref, then its in the dict, so add the positions.
                    if mut_rate_present:
                        # here we append the metadata, if present
                        index = line_list[3].find('mut_rate=')
                        if index == -1:
                            log_mssg(f"Skipping mutation rate region: {my_chr}: ({pos1}, {pos2})", 'warning')
                            log_mssg(f'4th column of mutation rate bed must be a semicolon list of key, value '
                                     f'pairs, with one key being mut_rate, e.g., "foo=bar;mut_rate=0.001;do=re".',
                                     'debug')
                            continue

                        # +9 because that's len('mut_rate='). Whatever is that should be our mutation rate.
                        mut_rate = line_list[3][index + 9:]
                        # We'll trim anything after the mutation rate and call it good. These should be ; separated
                        try:
                            mut_rate = float(mut_rate.split(';')[0])
                        except ValueError:
                            log_mssg(f"Invalid mutation rate: {my_chr}: ({pos1}, {pos2})", 'error')
                            log_mssg(f'4th column of mutation rate bed must be a semicolon list of key, value '
                                     f'pairs, with one key being mut_rate, e.g., "foo=bar;mut_rate=0.001;do=re".',
                                     'debug')
                            premature_exit(1)

                        if mut_rate > 0.3:
                            log_mssg("Found a mutation rate > 0.3. This is unusual.", 'Warning')

                        ret_dict[my_chr].append((int(pos1), int(pos2), mut_rate))
                    else:
                        ret_dict[my_chr].append((int(pos1), int(pos2)))

        except UnicodeDecodeError:
            # this error is to specifically detect files that are zipped using something other than gzip
            # Or if the file is really gzipped but has been misnamed.
            # The standard error here might be ambiguous, so we're using an except.
            log_mssg("Input bed files must be a valid bed files or a valid gzipped bed file. "
                     "Gzipped files must end in extension '.gz'", 'error')
            premature_exit(1)

        # some validation
        in_ref_only = [k for k in chromosomes if k not in ret_dict]
        if in_ref_only:
            log_mssg(f'Warning: Reference contains sequences not found in BED file {input_bed}. '
                          f'These chromosomes will be omitted from the outputs.', 'warning')
            log_mssg(f"In reference only regions: {in_ref_only}", 'debug')

        if in_bed_only:
            log_mssg(f'BED file {input_bed} contains sequence names '
                          f'not found in reference. These regions will be ignored.', 'warning')
            log_mssg(f'Regions ignored: {in_bed_only}', 'debug')

    # Returns an empty dict if there is no file to process.
    return ret_dict
