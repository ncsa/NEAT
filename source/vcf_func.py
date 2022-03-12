import gzip
import io
import logging
import random

import pandas as pd

from source.error_handling import premature_exit, log_mssg


def parse_input_vcf(vcf_path: str, tumor_normal: bool = False, ploidy: int = 2,
                    choose_random_ploid_if_no_gt_found: bool = True) -> (list, pd.DataFrame):

    log_mssg(f"Parsing input vcf {vcf_path}", 'debug')
    # Read in the raw vcf using pandas' csv reader.
    if vcf_path.endswith('.gz'):
        f = gzip.open(vcf_path, 'rt')
    else:
        f = open(vcf_path, 'r')

    columns = None
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            columns = line.strip().strip('#').split('\t')
            break
        else:
            # If we got to this point, we didn't find a header row.
            log_mssg(f'No header found for vcf: {vcf_path}', 'error')
            premature_exit(1)
    f.close()

    # Anything after FORMAT is a sample column
    sample_columns = columns[columns.index('FORMAT') + 1:]

    if sample_columns:
        if len(sample_columns) == 1 and not tumor_normal:
            sample_columns = [sample_columns[0]]
        elif len(sample_columns) >= 1 and not tumor_normal:
            log_mssg('More than one sample column present, only first sample column used.', 'warning')
            sample_columns = [sample_columns[0]]
        elif len(sample_columns) == 1 and tumor_normal:
            log_mssg(f'Tumor-Normal samples require '
                          f'both a tumor and normal sample column in the VCF. {list(sample_columns)}', 'error')
            premature_exit(1)
        elif len(sample_columns) >= 1 and tumor_normal:
            normals = [label for label in sample_columns if 'normal' in label.lower()]
            tumors = [label for label in sample_columns if 'tumor' in label.lower()]
            if not (tumors and normals):
                log_mssg("Input VCF for cancer must contain a column with a label containing 'tumor' "
                              "and 'normal' (case-insensitive).", 'error')
                premature_exit(1)
            if len(normals) > 1 or len(tumors) > 1:
                log_mssg("If more than one tumor or normal column is present, "
                              "only the first of each is used.", 'warning')
            sample_columns = [normals[0], tumors[0]]

        else:
            log_mssg('Unconsidered case encountered while parsing input vcf!', 'critical')
            log_mssg("Reality: Broken", 'error')
            premature_exit(1)
    else:
        log_mssg('Input VCF must have least one sample column', 'error')
        premature_exit(1)

    # We'll use FORMAT to check for genotype at some point.
    use_columns = ['CHROM', 'POS', 'REF', 'ALT', 'FORMAT'] + sample_columns
    data_types = {'CHROM': str, 'POS': int, 'REF': str, 'ALT': str, 'FORMAT': str}
    for col in sample_columns:
        data_types[col] = str

    variants = pd.read_csv(vcf_path, comment="#", sep="\t", header=None,
                           names=columns,
                           usecols=use_columns,
                           dtype=data_types)

    # Make the sample columns easier to recall later:
    if len(sample_columns) == 1:
        variants.rename(columns={sample_columns[0]: "input_sample"}, inplace=True)
    if len(sample_columns) > 1:
        variants.rename(columns={sample_columns[0]: "normal_sample", sample_columns[1]: "tumor_sample"})

    # Convert vcf coordinates (1-based) to reference coordinates (0-based)
    variants.POS = variants.POS - 1
    """
    Note: there used to be a large section here that figured out the GT and AF from the
    input vcf, but it was slow and we never use that info, so I removed it to try to keep
    the memory signature of this vcf lower. There also used to be checks for duplicates,
    but that's built into the code later, so I dropped it.
    """

    # Nan's give errors, so let's fill them out quick.
    variants = variants.fillna('.')

    log_mssg(f'Found {len(variants)} valid variants in input VCF.', 'info')

    return sample_columns, variants
