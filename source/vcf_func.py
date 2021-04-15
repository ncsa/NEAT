import io
import sys
import time
import re
import random
import pandas as pd


def parse_line(vcf_line, col_dict, col_samp):
    # these were in the original. Not sure the point other than debugging.
    include_homs = False
    include_fail = False

    # check if we want to proceed...
    reference_allele = vcf_line[col_dict['REF']]
    alternate_allele = vcf_line[col_dict['ALT']]
    # enough columns?
    if len(vcf_line) != len(col_dict):
        return None
    # exclude homs / filtered?
    if not include_homs and alternate_allele == '.' or alternate_allele == '' or alternate_allele == reference_allele:
        return None
    if not include_fail and vcf_line[col_dict['FILTER']] != 'PASS' and vcf_line[col_dict['FILTER']] != '.':
        return None

    #	default vals
    alt_alleles = [alternate_allele]
    alt_freqs = []

    gt_per_samp = []

    #	any alt alleles?
    alt_split = alternate_allele.split(',')
    if len(alt_split) > 1:
        alt_alleles = alt_split

    #	check INFO for AF
    af = None
    if 'INFO' in col_dict and ';AF=' in ';' + vcf_line[col_dict['INFO']]:
        info = vcf_line[col_dict['INFO']] + ';'
        af = re.findall(r"AF=.*?(?=;)", info)[0][3:]
    if af is not None:
        af_splt = af.split(',')
        while (len(af_splt) < len(alt_alleles)):  # are we lacking enough AF values for some reason?
            af_splt.append(af_splt[-1])  # phone it in.
        if len(af_splt) != 0 and af_splt[0] != '.' and af_splt[0] != '':  # missing data, yay
            alt_freqs = [float(n) for n in af_splt]
    else:
        alt_freqs = [None] * max([len(alt_alleles), 1])

    gt_per_samp = None
    #	if available (i.e. we simulated it) look for WP in info
    if len(col_samp) == 0 and 'INFO' in col_dict and 'WP=' in vcf_line[col_dict['INFO']]:
        info = vcf_line[col_dict['INFO']] + ';'
        gt_per_samp = [re.findall(r"WP=.*?(?=;)", info)[0][3:]]
    else:
        #	if no sample columns, check info for GT
        if len(col_samp) == 0 and 'INFO' in col_dict and 'GT=' in vcf_line[col_dict['INFO']]:
            info = vcf_line[col_dict['INFO']] + ';'
            gt_per_samp = [re.findall(r"GT=.*?(?=;)", info)[0][3:]]
        elif len(col_samp):
            fmt = ':' + vcf_line[col_dict['FORMAT']] + ':'
            if ':GT:' in fmt:
                gt_ind = fmt.split(':').index('GT')
                gt_per_samp = [vcf_line[col_samp[iii]].split(':')[gt_ind - 1] for iii in range(len(col_samp))]
                for i in range(len(gt_per_samp)):
                    gt_per_samp[i] = gt_per_samp[i].replace('.', '0')
        if gt_per_samp is None:
            gt_per_samp = [None] * max([len(col_samp), 1])

    return alt_alleles, alt_freqs, gt_per_samp


def parse_vcf(vcf_path: str, tumor_normal: bool = False, ploidy: int = 2,
              include_homs: bool = False, include_fail: bool = False, debug: bool = False,
              choose_random_ploid_if_no_gt_found: bool = True):

    if debug:
        print(f"Choosing random ploid if no GT found? {'yes' if choose_random_ploid_if_no_gt_found else 'no'}")

    tt = time.time()
    print('--------------------------------')
    print('reading input VCF...\n', flush=True)

    # Don't probably need any of these
    col_dict = {}
    col_samp = []
    n_skipped = 0
    n_skipped_because_hash = 0
    all_vars = {}  # [ref][pos]
    samp_names = []

    # These I need
    samp_cols = []
    printed_warning = False

    # Read in the raw vcf using pandas' csv reader.
    with open(vcf_path, 'r') as f:
        # quickest way I've found to read in the file:
        lines = [line for line in f if not line.startswith('##')]
        # Check to make sure header row is included
        if not lines[0].startswith('#CHROM'):
            print(f"ERROR: Improper vcf header row for {vcf_path}. Check and re-run.")
            sys.exit(1)
        # NOTE: if the vcf that is read in does not match the proper format, this read_csv command
        # will throw an error. This means you can't have data with no column header.
        variants = pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})

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
    if 'FORMAT' in variants.columns:
        samp_cols = variants.columns[list(variants.columns).index('FORMAT') + 1:]
        if len(samp_cols):
            if len(samp_cols) == 1 and not tumor_normal:
                variants['sample_split'] = variants[samp_cols[0]].str.split(':')
            elif len(samp_cols) >= 1 and not tumor_normal:
                print('More than one sample column present, only first sample column used.')
                variants['sample_split'] = variants[samp_cols[0]].str.split(':')
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
            else:
                print('ERROR: Unconsidered case: you may have broken reality. Check your VCF for the proper number'
                      'of sample columns.')
                sys.exit(1)
        else:
            print('ERROR: If FORMAT column is present in VCF, there must be at least one sample column.')
            sys.exit(1)

    # Split fields with multiple datapoints, if present, into lists
    variants['alt_split'] = variants['ALT'].str.split(',')
    if 'INFO' in variants.columns:
        variants['info_split'] = variants['INFO'].str.split(';')
    if 'FORMAT' in variants.columns:
        variants['format_split'] = variants['FORMAT'].str.split(':')

    # The following block of code looks for allele frequencies in the VCF.
    # There may be a more clever way to look through these subfields, but I loop (just once) over all theh rows.
    new_column = []
    print_message = False
    # TODO find a more efficient way than looping over the rows
    for row in variants.iterrows():
        af_numbers = []
        gt_numbers = []
        if 'INFO' in variants.columns:
            for info in row[1]['info_split']:
                # Looking for allele frequency (AF) in the info field
                if 'AF' in info:
                    # In case they try to do something like "AF=0.5;AF=0.25" instead of "AF=0.5,0.25"
                    if not print_message and debug:
                        print('Note: NEAT only uses the first AF in the info field of input VCF.')
                        print_message = True
                    for i in info.split('=')[1].split(','):
                        # If they didn't supply a number, but instead a '.' or a missing value or something, then
                        # the try/except block should catch it.
                        try:
                            af_numbers.append(float(i))
                        except ValueError:
                            print(f"Warning: format is off for AF on this row: {list(row)} \n"
                                  f"Proceeding without AF for this record.")
                            af_numbers.append(None)
                # WP is NEAT's ploidy indicator.
                elif 'WP' in info:
                    for j in info.split('=')[1].split(','):
                        # If they didn't supply a number, but instead a '.' or a missing value or something, then
                        # the try/except block should catch it.
                        gt_numbers.append(j)
        if 'FORMAT' in variants.columns and not gt_numbers:
            for format_item in row[1]['format_split']:
                # GT is the usual genotype indicator. will also need to search for the FORMAT field for this
                if 'GT' in format_item:
                    if tumor_normal:
                        pass
                    else:
                        # Append the corresponding item from the sample.
                        gt_numbers.append(row[1]['sample_split'][row[1]['format_split'].index(format_item)])
                        # We've found it, so we can quit looking.
                        break

        # Check that the data are consistent
        if af_numbers:
            if len(af_numbers) < len(row[1]['alt_split']):
                print(f"ERROR: allele frequency (AF) field in INFO must match number of alternate alleles: "
                      f"{list(row[1])}")
                sys.exit(1)
        else:
            # Used None value if no AF was supplied
            af_numbers.extend([None] * max([len(row[1]['alt_split']), 1]))
        if not gt_numbers:
            if choose_random_ploid_if_no_gt_found:
                pass
            else:
                gt_numbers.extend([None] * max([len(samp_cols), 1]))

        new_column.append([af_numbers, gt_numbers])
    # Add the new data to the table
    variants['AF'] = pd.DataFrame(new_column)[0]
    variants['GT'] = pd.DataFrame(new_column)[1]





    # TODO figure out what the pl_out stuff is for and convert to pandas.
    with open(vcf_path, 'r') as f:
        for line in f:
            if line[0] != '#':
                splt = line.strip().split('\t')
                pl_out = parse_line(splt, col_dict, col_samp)
                if pl_out is None:
                    n_skipped += 1
                else:
                    (aa, af, gt) = pl_out

                    # make sure at least one allele somewhere contains the variant
                    if tumor_normal:
                        gt_eval = gt[:2]
                    else:
                        gt_eval = gt[:1]
                    # For some reason this had an additional "if True" inserted. I guess it was supposed to be an option
                    # the user could set but was never implemented.
                    if None in gt_eval:
                        if choose_random_ploid_if_no_gt_found:
                            if not printed_warning:
                                print('Warning: Found variants without a GT field, assuming heterozygous...')
                                printed_warning = True
                            for i in range(len(gt_eval)):
                                tmp = ['0'] * ploidy
                                tmp[random.randint(0, ploidy - 1)] = '1'
                                gt_eval[i] = '/'.join(tmp)
                        else:
                            # skip because no GT field was found
                            n_skipped += 1
                            continue
                    non_reference = False
                    for gtVal in gt_eval:
                        if gtVal is not None:
                            if '1' in gtVal:
                                non_reference = True
                    if not non_reference:
                        # skip if no genotype actually contains this variant
                        n_skipped += 1
                        continue

                    chrom = splt[0]
                    pos = int(splt[1])
                    ref = splt[3]
                    # skip if position is <= 0
                    if pos <= 0:
                        n_skipped += 1
                        continue

                    # hash variants to avoid inserting duplicates (there are some messy VCFs out there...)
                    if chrom not in all_vars:
                        all_vars[chrom] = {}
                    if pos not in all_vars[chrom]:
                        all_vars[chrom][pos] = (pos, ref, aa, af, gt_eval)
                    else:
                        n_skipped_because_hash += 1
            else:
                if line[1] != '#':
                    cols = line[1:-1].split('\t')
                    for i in range(len(cols)):
                        if 'FORMAT' in col_dict:
                            col_samp.append(i)
                        col_dict[cols[i]] = i
                    if len(col_samp):
                        samp_names = cols[-len(col_samp):]
                        if len(col_samp) == 1:
                            pass
                        elif len(col_samp) == 2 and tumor_normal:
                            print('Detected 2 sample columns in input VCF, assuming tumor/normal.')
                        else:
                            print(
                                'Warning: Multiple sample columns present in input VCF. By default genReads uses '
                                'only the first column.')
                    else:
                        samp_names = ['Unknown']
                    if tumor_normal:
                        # tumorInd  = samp_names.index('TUMOR')
                        # normalInd = samp_names.index('NORMAL')
                        if 'NORMAL' not in samp_names or 'TUMOR' not in samp_names:
                            print('\n\nERROR: Input VCF must have a "NORMAL" and "TUMOR" column.\n')


    vars_out = {}
    for r in all_vars.keys():
        vars_out[r] = [list(all_vars[r][k]) for k in sorted(all_vars[r].keys())]
        # prune unnecessary sequence from ref/alt alleles
        for i in range(len(vars_out[r])):
            while len(vars_out[r][i][1]) > 1 and all([n[-1] == vars_out[r][i][1][-1] for n in vars_out[r][i][2]]) \
                    and all([len(n) > 1 for n in vars_out[r][i][2]]):
                vars_out[r][i][1] = vars_out[r][i][1][:-1]
                vars_out[r][i][2] = [n[:-1] for n in vars_out[r][i][2]]
            vars_out[r][i] = tuple(vars_out[r][i])

    print('found', sum([len(n) for n in all_vars.values()]), 'valid variants in input vcf.')
    print(' *', n_skipped, 'variants skipped: (qual filtered / ref genotypes / invalid syntax)')
    print(' *', n_skipped_because_hash, 'variants skipped due to multiple variants found per position')
    print(f'vcf reading took: {int(time.time() - tt)} (sec)')
    print('--------------------------------')
    return samp_names, vars_out
