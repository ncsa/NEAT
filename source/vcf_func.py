import gzip
import io
import logging
import random

import pandas as pd

from source.error_handling import premature_exit, log_mssg


def retrieve_attribute_from_format(my_record, value, which_sample):
    # the first sample is index 9, but if there is a tumor, we'll add one to look in that column
    which_sample_index = 9 + which_sample
    index = my_record[8].split(':').index(value)
    # Apply the index corresponding to GT to the sample column to get the genotype, then split into ploids.
    ret = my_record[which_sample_index].split(':')[index].replace('/', '|').split('|')
    return [int(x) for x in ret]


def parse_input_vcf(vcf_path: str, tumor_normal: bool = False) -> (list, pd.DataFrame):

    log_mssg(f"Parsing input vcf {vcf_path}", 'debug')
    # Read in the raw vcf using pandas' csv reader.
    if vcf_path.endswith('.gz'):
        f = gzip.open(vcf_path, 'rt')
    else:
        f = open(vcf_path, 'r')

    n_skipped = 0
    ret_dict = {}
    # maximum number of columns we are interested in. Used for trimming unwanted samples.
    max_col = 7
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            columns = line.strip().strip('#').split('\t')

            # Anything after FORMAT is a sample column
            sample_columns = []
            format_col = False
            if 'FORMAT' in columns:
                format_col = True
                max_col += 1
                sample_columns = columns[columns.index('FORMAT') + 1:]
                if not sample_columns:
                    log_mssg('Sample column missing in in input vcf, using WP,'
                             'if present or generating genotypes randomly.', 'info')
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

        else:
            # What follows in this block processes the variant line
            record = line.strip().split('\t')
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
                SAMPLE2 [10, optional]
            """
            # We'll need the genotype when we generate reads, though it's less important for the vcf stage.
            genotype = None
            genotype_tumor = None
            if format_col:
                if "GT" in line:
                    # Retrieve the GT from the first sample in the record
                    genotype = retrieve_attribute_from_format(record, 'GT', 0)
                    if tumor_normal:
                        # Same procedure as above, but with the tumor sample
                        genotype_tumor = retrieve_attribute_from_format(record, 'GT', 1)
            elif "WP" in line:
                rec_split = record[7].split(';')
                for rec in rec_split:
                    if rec.startswith('WP'):
                        genotype = rec.split('=')[1].replace('/', '|').split('|')
                        genotype = [int(x) for x in genotype]
            if genotype_tumor:
                genotype = [genotype, genotype_tumor]
            else:
                genotype = [genotype]
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
                # +1 here so we grab the last column of interest in the slice
                ret_dict[key] = [record[2:max_col + 1]] + genotype
            else:
                log_mssg(f'Skipping variant because multiple variants found at this location:'
                         f'{key}', 'warning')
                n_skipped += 1

    log_mssg(f'Found {len(ret_dict)} variants in input VCF.', 'info')
    log_mssg(f'Skipped {n_skipped} variants because of multiples at the same location', 'info')

    return sample_columns, ret_dict
