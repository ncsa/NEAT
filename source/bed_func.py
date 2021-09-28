import pathlib
import gzip
import re
import pandas as pd
from source.error_handling import premature_exit, print_and_log


def parse_bed(input_bed: str, chromosomes: list,
              begins_with_chr: bool,
              mut_rate_present: bool,
              debug: bool = False) -> dict:
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
    :param debug: True or false if we are in debugging mode and need more info
    :return: a dictionary of chromosomes, -1, pos1, and pos2
    """
    ret_dict = {}
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
                    except ValueError as e:
                        print_and_log(f"Improperly formatted bed file line {line}", "error")
                        premature_exit(1)
                    # Trying not to 'fix' bed files, but an easy case is if the vcf uses 'chr' and the bed doesn't, in
                    # which we frequently see in human data, we'll just put chr on there. Anything more complicated
                    # will be on the user to correct their beds first.
                    if begins_with_chr and not my_chr.startswith('chr'):
                        my_chr = 'chr' + my_chr
                    # Check to see if this chromosome from the bed is even in the ref.
                    # If not, there's not point in including it
                    in_ref = False
                    for k in chromosomes:
                        if my_chr == k:
                            in_ref = True
                            break
                    if not in_ref:
                        in_bed_only.append(my_chr)
                        continue
                    # If it's in ref, add chrom to dict if needed, then add the positions.
                    if my_chr not in ret_dict:
                        if not mut_rate_present:
                            ret_dict[my_chr] = [-1]
                        else:
                            ret_dict[my_chr] = []
                    # here we append the metadata, if present
                    if mut_rate_present:
                        ret_dict[my_chr].append([int(pos1), int(pos2), str(line_list[3])])
                    else:
                        ret_dict[my_chr].extend([int(pos1), int(pos2)])
        except UnicodeDecodeError:
            # this error is to specifically detect files that are zipped using something other than gzip
            # Or if the file is really gzipped but has been misnamed.
            # The standard error here might be ambiguous, so we're using an except.
            print_and_log("Input bed files must be a valid bed files or a valid gzipped bed file. "
                          "Gzipped files must end in extension '.gz'", 'error')
            premature_exit(1)

        # some validation
        in_ref_only = [k for k in chromosomes if k not in ret_dict]
        if in_ref_only:
            print_and_log(f'Warning: Reference contains sequences not found in BED file {input_bed}. '
                          f'These chromosomes will be omitted from the outputs.', 'warning')
            if debug:
                print_and_log(f"In reference only regions: {in_ref_only}", 'debug')

        if in_bed_only:
            print_and_log(f'BED file {input_bed} contains sequence names '
                          f'not found in reference. These regions will be ignored.', 'warning')
            if debug:
                print_and_log(f'Regions ignored: {in_bed_only}', 'debug')

        if mut_rate_present:
            # Assuming there are mutation rate regions, this will require one extra value for the return
            mut_rate_dict = {}
            printed_warning = False
            if ret_dict:
                for key in ret_dict.keys():
                    # We have to find each of the regions in each chromosome.
                    for region in ret_dict[key]:
                        try:
                            meta_data = region[2]
                            # We allow for more than one mutation rate, but let's put in a warning so the user knows
                            # something is up with their mutation rate file
                            mut_str = re.findall(r"mut_rate=.*?(?=;)", meta_data + ';')
                            if len(mut_str) > 1:
                                print_and_log("Found mutation rate record with more than one mut_rate value. "
                                              "Using the smallest number", 'warning')
                                if debug:
                                    print_and_log(f"Record with multiple mut_rates: {region}", 'debug')

                            # Skip if the mutation string was empty
                            if not mut_str:
                                print_and_log("Mutation bed record with no mutation rate, skipping", "warning")
                                continue

                            (pos1, pos2) = region[0:2]
                            if pos2 - pos1 >= 0:
                                # mut_rate = #_mutations / length_of_region
                                temp_rate = min([float(k.split('=')[1]) for k in mut_str])
                                if temp_rate > 0.3:
                                    print_and_log("Found a mutation rate > 0.3. This is unusual", 'Warning')
                                if key not in mut_rate_dict:
                                    mut_rate_dict[key] = {'regions': [-1], 'values': [0.0]}
                                mut_rate_dict[key]['values'].extend([temp_rate * (pos2 - pos1)] * 2)
                                mut_rate_dict[key]['regions'].extend([pos1, pos2])
                            else:
                                print_and_log("Skipping invalid region (end - start < 0)", 'warning')
                                if debug:
                                    print_and_log(f'region skipped: {region}', 'debug')

                            if not mut_str and (not printed_warning or debug):
                                print_and_log(f"Warning: Mutation rate record(s) in bed {input_bed} "
                                              f"with no mutation rate in the file. Mutation rates must be in the fourth "
                                              f"column and of the form 'mut_rate=X.XXX'. Skipping records with no "
                                              f"mut_rate.", 'warning')
                                if debug:
                                    print_and_log(f"Record with a problem: {key}, {ret_dict[key]}", 'debug')
                                printed_warning = True
                                continue

                        except IndexError:
                            print_and_log("Malformed mutation bed file. Must be of the format "
                                          "'chromosome\tpos1\tpos2\tmut_rate=X.XX'.", 'error')

            return pd.DataFrame.from_dict(mut_rate_dict, orient='index')
    # Return empty dicts if there is no file to process.
    if mut_rate_present:
        # if we reached this point, then there was no mut_rates, so we're good to just return an empty df.
        return pd.DataFrame()
    else:
        # If there was no file, it returns an empty dict
        return pd.DataFrame.from_dict(ret_dict, orient='index')