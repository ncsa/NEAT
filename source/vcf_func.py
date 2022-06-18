import gzip

from Bio import SeqIO
from numpy.random import Generator

from source.error_handling import premature_exit, log_mssg
from source.ploid_functions import pick_ploids
from source.constants_and_defaults import ALLOWED_NUCL
from source.Options import Options


def retrieve_attribute_from_format(my_record, value, which_sample):
    # the first sample is index 9, but if there is a tumor, we'll add one to look in that column
    which_sample_index = 9 + which_sample
    index = my_record[8].split(':').index(value)
    # Apply the index corresponding to GT to the sample column to get the genotype, then split into ploids.
    ret = my_record[which_sample_index].split(':')[index].replace('/', '|').split('|')
    return [int(x) for x in ret]


def parse_input_vcf(vcf_path: str, ploidy: int, homozygous_frequency: float,
                    reference: SeqIO, options: Options, tumor_normal: bool = False) -> dict:

    log_mssg(f"Parsing input vcf {vcf_path}", 'debug')
    # Read in the raw vcf using pandas' csv reader.
    if vcf_path.endswith('.gz'):
        f = gzip.open(vcf_path, 'rt')
    else:
        f = open(vcf_path, 'r')

    n_skipped = 0
    mismatched = 0
    ret_dict = {}
    # maximum number of columns we are interested in. Used for trimming unwanted samples.
    max_col = 7
    random_ploid = False
    for line in f:
        # skip headers
        if line.startswith('##'):
            continue
        # Process the columns row
        elif line.startswith('#CHROM'):
            columns = line.strip().strip('#').split('\t')

            # Anything after FORMAT is a sample column
            sample_columns = []
            has_format = False
            if 'FORMAT' in columns:
                has_format = True
                max_col += 1
                sample_columns = columns[columns.index('FORMAT') + 1:]
                if not sample_columns:
                    log_mssg('Sample column missing in in input vcf.', 'error')
                    premature_exit(1)
            else:
                log_mssg('Missing format column, using WP for genotype if present, '
                         'otherwise genotype will be generated randomly', 'info')

            # Recode the sample columns to match the index of the dictionary we are generating
            # We only output 1 sample column for normal runs, 2 for tumor_normal. Those will be indices 7 and 8
            # in the output dictionary, se we hardcode those indices now for later retrieval

            if sample_columns:
                if not tumor_normal:
                    sample_columns = {sample_columns[0]: 7}
                    max_col += 1
                # If the code got here, we're dealing with a cancer sample
                elif len(sample_columns) == 1:
                    log_mssg(f'Tumor-Normal samples require '
                             f'both a tumor and normal sample column in the VCF. {list(sample_columns)}', 'error')
                    premature_exit(1)
                else:
                    normals = [label for label in sample_columns if 'normal' in label.lower()]
                    tumors = [label for label in sample_columns if 'tumor' in label.lower()]
                    if not (tumors and normals):
                        log_mssg("Input VCF for cancer must contain a column with a label containing 'tumor' "
                                 "and 'normal' (case-insensitive).", 'error')
                        premature_exit(1)
                    sample_columns = {normals[0]: 7, tumors[0]: 8}
                    max_col += 2

        # Process the records rows
        else:
            record = line.strip().split('\t')
            # Decrement the position to get 0-based coordinates
            record[1] = int(record[1]) - 1
            # We'll index these by chromosome and position
            """
            For reference, the columns in a VCF, and their indices:
                CHROM [0]
                POS [1]
                ID [2]
                REF [3]
                ALT [4]
                QUAL [5]
                FILTER [6]
                INFO [7]
                FORMAT [8, optional]
                SAMPLE1 [9, optional]
                SAMPLE2 [10, optional, cancer only]
            """
            # First, let's check if the chromosome for this record is even in the reference:
            in_ref = record[0] in reference
            if not in_ref:
                log_mssg(f'Skipping variant because the chromosome is not in the reference:\n{line}', 'warning')
                continue
            # Before we get too deep into parsing, let's make sure some basics are covered:
            # Check if there are any not allowed nucleotides
            if any([x for x in record[3] if x not in ALLOWED_NUCL]) or \
                    any([x for x in record[4] if x not in ALLOWED_NUCL]):
                mismatched += 1
                log_mssg(f'Skipping variant because the ref or alt field contains disallowed nucelotides:'
                         f'{record[0]}: {record[1]}, {record[3]}, {record[4]}', 'warning')
                continue
            # We already accounted for shifting to 0-based coordinates, so this should work.
            if record[3] != reference[record[0]][int(record[1]): int(record[1]) + len(record[3])].seq:
                mismatched += 1
                log_mssg(f'Skipping variant because the ref field did not match the reference:'
                         f'{record[0]}: {record[1]}, {record[3]} v '
                         f'{reference[record[0]][int(record[1]): int(record[1]) + len(record[3])].seq}', 'warning')
                continue

            # We'll need the genotype when we generate reads, and output the records, if applicable
            genotype = None
            genotype_tumor = None
            format_column = None
            normal_sample_column = None
            tumor_sample_column = None

            if has_format:
                if "GT" in record[8].split(':'):
                    # the format column will need no update.
                    format_column = record[8]
                    normal_sample_column = record[9]
                    # Retrieve the GT from the first sample in the record
                    genotype = retrieve_attribute_from_format(record, 'GT', 0)
                    if tumor_normal:
                        # Same procedure as above, but with the tumor sample
                        tumor_sample_column = record[10]
                        genotype_tumor = retrieve_attribute_from_format(record, 'GT', 1)
                else:
                    format_column = 'GT:' + record[8]
                    alt_count = len(record[4].split(';'))
                    genotype = pick_ploids(ploidy, homozygous_frequency, alt_count)
                    gt_field = "/".join([str(x) for x in genotype])
                    normal_sample_column = f'{gt_field}:{record[9]}'
                    if tumor_normal:
                        genotype_tumor = pick_ploids(ploidy, homozygous_frequency, alt_count)
                        # Since this is random, if we accidentally pick the same ploid,
                        # let's just shuffle until they are different
                        # But we'll cap it at 10 tries
                        i = 10
                        while genotype_tumor == genotype or i > 0:
                            options.rng.shuffle(genotype_tumor)
                            i -= 1
                        if genotype == genotype_tumor:
                            log_mssg(f'Skipping variant that already had a '
                                     f'variant at that location {record[0]}: {record[1]}')
                            continue
                        gt_field = '/'.join([str(x) for x in genotype_tumor])
                        tumor_sample_column = f'{gt_field}:{record[10]}'
            # "WP" is the legacy code NEAT used for genotype it added. It was found in the INFO field.
            # We're just going to make a sample column in this version of NEAT
            # The logic of the statement is split the info field on ';' which is used as a divider in that field.
            # Most but not all fields also have an '=', so split there too, then look for "WP"
            elif "WP" in [x.split('=') for x in record[7].split(';')]:
                format_column = "GT"
                info_split = record[7].split(';')
                for record in info_split:
                    if record.startswith('WP'):
                        genotype = record.split('=')[1].replace('/', '|').split('|')
                        genotype = [int(x) for x in genotype]
                        # If there was no format column, there's no sample column, so we'll generate one just
                        # with the genotype matching the old WP notation in INFO.
                        normal_sample_column = '/'.join([str(x) for x in genotype])
            else:
                format_column = "GT"
                genotype = pick_ploids(ploidy, homozygous_frequency, alt_count)
                normal_sample_column = '/'.join([str(x) for x in genotype])

            key = (record[0], int(record[1]))
            if key not in ret_dict:
                """
                key to ret_dict:
                    - key = (CHROM, POS)
                    - 2D array value
                        - 0: The actual info from the record
                        - 1: Genotype data for the first sample
                        - 2: Genotype data for the second (tumor) sample (optional)
                """
                ret_dict[key] = [record[2:8] + [format_column, normal_sample_column], genotype]
                if tumor_normal:
                    ret_dict[key][0].extend([tumor_sample_column])
                    ret_dict[key].append(genotype_tumor)
            else:
                # Potential duplicate, so we check the genotype
                if ret_dict[key][1] == genotype or ret_dict[key][2] == genotype_tumor:
                    # This is a duplicate, two alts in the same location on the samp ploid
                    log_mssg(f'Skipping variant because multiple variants found at this location:'
                             f'{key}', 'warning')
                    n_skipped += 1
                else:
                    # It's on a different ploid, so let's add it.
                    ret_dict[key] = [record[2:8] + [format_column, normal_sample_column], genotype]
                    if tumor_normal:
                        ret_dict[key][0].extend([tumor_sample_column])
                        ret_dict[key].append(genotype_tumor)

    log_mssg(f'Found {len(ret_dict)} variants in input VCF.', 'info')
    log_mssg(f'Skipped {n_skipped} variants because of multiples at the same location', 'info')

    return sample_columns, ret_dict
