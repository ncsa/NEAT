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


def parse_vcf(vcf_path: str, tumor_normal: bool = False, ploidy: int = 2):
    # this var was in the orig. May have just been a debugging thing.
    # I think this is trying to implement a check on if the genotype
    choose_random_ploid_if_no_gt_found = True
    # Pseudocode for above line
    tt = time.time()
    print('--------------------------------')
    print('reading input VCF...\n', flush=True)

    col_dict = {}
    col_samp = []
    n_skipped = 0
    n_skipped_because_hash = 0
    all_vars = {}  # [ref][pos]
    samp_names = []
    printed_warning = False
    with open(vcf_path, 'r') as f:
        lines = [line for line in f if not line.startswith('##')]
        # sanity check
        if not lines[0].startswith('#CHROM'):
            print("WARNING: Improper vcf header row. Check and re-run.")

        # NOTE: if the vcf that is read in does not match the proper format
        variants = pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})

    # the following section is just some sanity checking.
    valid_headers = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL']
    for i in range(len(valid_headers)):
        if valid_headers[i] != variants.columns[i]:
            print(f"VCF must contain the following headers, in order: {valid_headers}")
            sys.exit(1)

    if 'FORMAT' in variants.columns:
        samp_cols = variants.columns[list(variants.columns).index('FORMAT') + 1:]

        if len(samp_cols):
            if len(samp_cols) == 1 and not tumor_normal:
                pass
            elif len(samp_cols) >= 1 and not tumor_normal:
                print('More than one sample column present, only first sample column used.')
                samp_cols = [samp_cols[0]]
            elif len(samp_cols) == 1 and tumor_normal:
                print('Tumor-Normal samples require both a tumor and normal column in the VCF.')
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
            else:
                print('ERROR: Unconsidered case: you may have broken reality.')
                sys.exit(1)
        else:
            print('ERROR: If FORMAT column is present in VCF, there must be at least one sample column.')
            sys.exit(1)


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
