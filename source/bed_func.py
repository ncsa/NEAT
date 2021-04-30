import pathlib
import gzip
import sys


def parse_bed(input_bed: str, chromosomes: list,
              begins_with_chr: bool, debug: bool) -> dict:
    """
    Will parse a bed file, returning a dictionary of chromosomes that are found in the reference,
    and the corresponding positions. Some beds may have metadata in the fourth column, in which case,
    we'll throw that in too, in case it's needed.
    :param input_bed: A bed file containing the regions of interest
    :param chromosomes: A list of chromosomes to check
    :param begins_with_chr: A true or false if this reference has human like chromosomes (chr1, chr2, etc), used because
    it's a common thing with human datasets to either use or not use "chr" before the number and no one seems to agree.
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
                    meta_data_present = False
                    line_list = line.strip().split('\t')[:4]
                    if len(line_list) > 3:
                        meta_data_present = True
                    try:
                        [my_chr, pos1, pos2] = line_list[:3]
                    except ValueError as e:
                        print(sys.exc_info()[2])
                        print(e)
                        print("ERROR: possibly due to an error in the format of the BED file record:")
                        if debug:
                            print(f'Record that cause the problem: {line}')
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
                        ret_dict[my_chr] = [-1]
                    # here we append the metadata, if present
                    if meta_data_present:
                        ret_dict[my_chr].append([int(pos1), int(pos2), str(line_list[3])])
                    else:
                        ret_dict[my_chr].extend([int(pos1), int(pos2)])
        except UnicodeDecodeError as e:
            # this error is to specifically detect files that are zipped using something other than gzip
            # Or if the file is really gzipped but has been misnamed.
            # The standard error here might be ambiguous, so we're using an except.
            print(e)
            print("Error: input bed files must be a valid bed files or a valid gzipped bed file. "
                  "Gzipped files must end in extension '.gz'")
            sys.exit(1)


        # some validation
        in_ref_only = [k for k in chromosomes if k not in ret_dict]
        if in_ref_only:
            print(f'Warning: Reference contains sequences not found in BED file {input_bed}. '
                  f'These chromosomes will be omitted from the outputs.')
            if debug:
                print(f'\nFound in reference only: {in_ref_only}')

        if in_bed_only:
            print(f'Warning: BED file {input_bed} contains sequence names '
                  f'not found in reference. These regions will be ignored.\n\tIf this is unexpected, check that '
                  f'the names in the BED file exactly match the names in the reference file.')
            if debug:
                print(f'Regions ignored: {in_bed_only}')

    return ret_dict
