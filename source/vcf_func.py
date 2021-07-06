import io
import sys
import time
import gzip
import random
import pandas as pd


def parse_vcf(vcf_path: str, tumor_normal: bool = False, ploidy: int = 2,
              include_homs: bool = False, include_fail: bool = False, debug: bool = False,
              choose_random_ploid_if_no_gt_found: bool = True):

    tt = time.time()
    # Read in the raw vcf using pandas' csv reader.
    if vcf_path.endswith('.gz'):
        f = gzip.open(vcf_path)
    else:
        f = open(vcf_path, 'r')

    # quickest way I've found to read in the file:
    lines = [line for line in f if not line.startswith('##')]
    f.close()

    # Check to make sure header row is included
    if not lines[0].startswith('#CHROM'):
        print(f"ERROR: Improper vcf header row for {vcf_path}. Check and re-run.")
        sys.exit(1)
    else:
        lines[0] = lines[0].strip('#')
    # NOTE: if the vcf that is read in does not match the proper format, this read_csv command
    # will throw an error. This means you can't have data with no column header.
    variants = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str},
        sep='\t'
    )

    # the following section is just some sanity checking.
    min_headers = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL']
    for i in range(len(min_headers)):
        if min_headers[i] != variants.columns[i]:
            print(f"ERROR: VCF must contain the following headers, in order: {min_headers}")
            sys.exit(1)
    if debug:
        optional_headers = ['FILTER', 'INFO', 'FORMAT']
        for j in range(len(optional_headers)):
            if optional_headers[j] != variants.columns[j]:
                print(f'Warning: missing optional header: {optional_headers[j]}.'
                      f'Though not required, a full VCF with complete fields will be helpful.')

    # Check for homs and fails, and drop those rows unless otherwise specified
    if not include_homs:
        variants = variants.drop(variants[(variants['ALT'] == '.') |
                                          (variants['ALT'] == '') |
                                          (variants.apply(lambda row: all(j == row.ALT for j in row.REF), axis=1))
                                          ].index)
    if not include_fail:
        variants = variants.drop(variants[(variants['FILTER'] != 'PASS') &
                                          (variants['FILTER'] != '.')
                                          ].index)

    # If FORMAT is present in the vcf, there must be corresponding Sample columns.
    samp_cols = []
    if 'FORMAT' in variants.columns:
        # VCF spec says that all columns after FORMAT are sample columns.
        samp_cols = variants.columns[list(variants.columns).index('FORMAT') + 1:]
        if len(samp_cols):
            if len(samp_cols) == 1 and not tumor_normal:
                variants['sample_split'] = variants[samp_cols[0]].str.split(':')
                samp_cols = ['sample_split']
            elif len(samp_cols) >= 1 and not tumor_normal:
                print('More than one sample column present, only first sample column used.')
                variants['sample_split'] = variants[samp_cols[0]].str.split(':')
                samp_cols = ['sample_split']
            elif len(samp_cols) == 1 and tumor_normal:
                print(f'Tumor-Normal samples require both a tumor and normal column in the VCF. \n'
                      f'Supplied samples = {list(samp_cols)}')
                sys.exit(1)
            elif len(samp_cols) >= 1 and tumor_normal:
                normals = [label for label in samp_cols if 'normal' in label.lower()]
                tumors = [label for label in samp_cols if 'tumor' in label.lower()]
                if not (tumors and normals):
                    print("ERROR: Input VCF for cancer must contain a column with a label containing 'tumor' "
                          "and 'normal' (case-insensitive).")
                    sys.exit(1)
                if len(normals) > 1 or len(tumors) > 1:
                    print("WARNING: If more than one tumor or normal column is present, "
                          "only the first of each is used.")
                samp_cols = [normals[0], tumors[0]]
                variants['normal_sample_split'] = variants[samp_cols[0]].str.split(':')
                variants['tumor_sample_split'] = variants[samp_cols[1]].str.split(':')
                samp_cols = ['normal_sample_split', 'tumor_sample_split']
            else:
                print('ERROR: Unconsidered case: you may have broken reality. Check your VCF for the proper number'
                      'of sample columns.')
                sys.exit(1)
        else:
            print('ERROR: If FORMAT column is present in VCF, there must be at least one sample column.')
            sys.exit(1)
    else:
        print("Warning: Input VCF files must have a FORMAT and SAMPLE column for variant insert to work.")

    # Split fields with multiple datapoints, if present, into lists
    variants['alt_split'] = variants['ALT'].str.split(',')
    variants = variants.explode('alt_split')
    if 'INFO' in variants.columns:
        variants['info_split'] = variants['INFO'].str.split(';')
    if 'FORMAT' in variants.columns:
        variants['format_split'] = variants['FORMAT'].str.split(':')

    # The following block of code looks for allele frequencies in the VCF.
    # There may be a more clever way to look through these subfields, but I loop (just once) over all the rows.
    new_column = []
    printed_warning = False
    rows_to_delete = []
    n_skipped = 0
    # TODO find a more efficient way than looping over the rows
    # Nan's give errors, so let's fill them out quick.
    variants = variants.fillna('.')
    for index, row in variants.iterrows():
        af_numbers = []
        gt_numbers = []
        # Looking for AF in INFO field to get the allele frequency.
        if 'INFO' in variants.columns:
            found_af = False
            found_wp = False
            for info in row['info_split']:
                # breaks the loop once everthing is found
                if found_af and found_wp:
                    break

                # Looking for allele frequency (AF) in the info field
                if 'AF' in info and not found_af:
                    # In case they try to do something like "AF=0.5;AF=0.25" instead of "AF=0.5,0.25"
                    if not printed_warning and debug:
                        print('Note: NEAT only uses the first AF in the info field of input VCF.')
                        print_message = True
                    found_af = True
                    for i in info.split('=')[1].split(','):
                        # If they didn't supply a number, but instead a '.' or a missing value or something, then
                        # the try/except block should catch it.
                        try:
                            af_numbers.append(float(i))
                        except ValueError:
                            print(f"Warning: format is off for AF on this row: {list(row)} \n"
                                  f"Proceeding without AF for this record.")
                            af_numbers.append(None)
                # WP is NEAT's ploidy indicator. Check if it's there.
                elif 'WP' in info and not found_wp:
                    found_wp =True
                    for j in info.split('=')[1].split(','):
                        # If they didn't supply a number, but instead a '.' or a missing value or something, then
                        # the try/except block should catch it.
                        gt_numbers.append(j)

        # If WP is present, we'll use that, otherwise look for GT in the FORMAT field.
        if 'FORMAT' in variants.columns and not gt_numbers:
            # Assuming there was data in the FORMAT column, this will find the first GT.
            if row['format_split'] != ['.']:
                for format_item in row['format_split']:
                    # GT is the usual genotype indicator. will also need to search for the FORMAT field for this
                    if 'GT' in format_item:
                        if tumor_normal:
                            # Look for the normal and tumor sample GTs by finding the index of the item in question
                            to_add = [row['normal_sample_split']
                                      [row['format_split'].index(format_item)].replace('|', '/'),
                                      row['tumor_sample_split']
                                      [row['format_split'].index(format_item)].replace('|', '/')]
                            gt_numbers.append(to_add)
                        else:
                            # Append the corresponding item from the sample.
                            gt_numbers.append(row['sample_split'][row['format_split'].index(format_item)])
                            # We've found it, so we can quit looking.
                        break

            # If there is not GT or WP present, then we either choose a random ploid or skip
            else:
                if choose_random_ploid_if_no_gt_found:
                    if not printed_warning:
                        print('Warning: Found variants without a GT field, assuming heterozygous...')
                        printed_warning = True
                    tmp_list = []
                    for i in range(len(row['alt_split'])):
                        tmp = ['0'] * ploidy
                        tmp[random.randint(0, ploidy-1)] = '1'
                        tmp_list.append('/'.join(tmp))
                    gt_numbers.extend(tmp_list)
                else:
                    if not printed_warning:
                        print('Warning: Found variants without a GT field, ignoring variants...')
                        printed_warning = True
                    rows_to_delete.append(index)
        # Trim unnecessary sequences from alleles
        while (len(row['REF']) > 1) and \
                (all([n[-1] == row['REF'] for n in row['alt_split']])) and \
                (all([len(n) > 1 for n in row['alt_split']])):
            variants.loc[index, 'REF'] = variants.loc[index, 'REF'][:-1]
            variants.loc[index, 'alt_split'] = [[n[:-1] for n in variants.loc[index, 'alt_split']]]
        # Check that the data are consistent
        if af_numbers:
            if len(af_numbers) < len(row['alt_split']):
                print(f"ERROR: allele frequency (AF) field in INFO must match number of alternate alleles: "
                      f"{list(row)}")
                sys.exit(1)
        else:
            # Used None value if no AF was supplied
            af_numbers.extend([None] * max([len(row['alt_split']), 1]))
        if not gt_numbers:
            rows_to_delete.append(index)
            if debug:
                print(f'Skipping row because no genotype found:\n{row}')
        else:
            # drop variants that aren't actually used
            for gt in gt_numbers:
                if gt == '0/0':
                    rows_to_delete.append(index)
                    if debug:
                        print(f'Skipping row because of 0/0 genotype:\n{row}')
        # Append column to form new AF and GT columns of the dataframe
        new_column.append([af_numbers, gt_numbers])
    # Add the new data to the table
    variants['AF'] = pd.DataFrame(new_column)[0]
    variants['GT'] = pd.DataFrame(new_column)[1]
    # drop rows with no genotype numbers
    variants = variants.drop(rows_to_delete)
    n_skipped += len(rows_to_delete)
    # drop rows with position <= 0
    variants = variants.drop(variants[variants["POS"] <= 0].index)
    n_skipped += len(variants[variants["POS"] <= 0].index)
    # Delete rows where they try to insert more than one variant
    n_skipped_because_hash = 0
    variants = variants.loc[~variants.duplicated(subset=['CHROM', 'POS'])]
    n_skipped_because_hash += len(variants.loc[variants.duplicated(subset=['CHROM', 'POS'])].index)

    variants = variants.sort_values(by=['CHROM', 'POS'])

    print(f'Found {len(variants.index)} valid variants in input vcf.')
    print(f' * {n_skipped} variants skipped: (qual filtered / ref genotypes / invalid syntax)')
    print(f' * {n_skipped_because_hash} variants skipped due to multiple variants found per position')
    print(f'vcf reading took: {int(time.time() - tt)} (sec)')
    print('--------------------------------')
    return list(samp_cols), variants
