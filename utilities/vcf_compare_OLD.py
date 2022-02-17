#!/usr/bin/env python
# encoding: utf-8

# Python 3 ready

""" **************************************************

vcf_compare_OLD.py

- compare vcf file produced by workflow to golden vcf produced by simulator

Written by:		Zach Stephens
Date:			January 20, 2015
Contact:		zstephe2@illinois.edu

Bug Fixes:      Jenna Kalleberg
Date:           January 26, 2022
Contac:         jenna.kalleberg@mail.missouri.edu

************************************************** """

import sys, copy, time, bisect, re, argparse
import numpy as np
import argparse
## DEBUG NOTE: updated libraries
import gzip, shutil, csv
from pathlib import Path
from Bio.Seq import MutableSeq

EV_BPRANGE = 50  # how far to either side of a particular variant location do we want to check for equivalents?

DEFAULT_QUAL = -666  # if we can't find a qual score, use this instead so we know it's missing

MAX_VAL = 9999999999999  # an unreasonably large value that no reference fasta could concievably be longer than

## DEBUG NOTE: argparse module can not accept a single '%' in the help string
##             only a double '%%' will work
# DESC = """%(prog): vcf comparison script."""
VERS = 0.1

parser = argparse.ArgumentParser(usage='python {} [options] -r <ref.fa> -g <golden.vcf> -w <workflow.vcf>'.format(sys.argv[0]),
                                 description='{0} v{1}: vcf comparison script.'.format(sys.argv[0], str(VERS)),
                                 # version="%prog v" + str(VERS),
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

## DEBUG NOTE: argparse module returns defaults automatically now
##             and the '%' character was causing a bug
## DEBUG NOTE: 'version' is no longer a valid object in argparse
##             but I made the description behave as intended
parser.add_argument('-r', help='* Reference Fasta', dest='reff', action='store', metavar='<ref.fa>')
parser.add_argument('-g', help='* Golden VCF', dest='golden_vcf', action='store', metavar='<golden.vcf>')
parser.add_argument('-w', help='* Workflow VCF', dest='workflow_vcf', action='store', metavar='<workflow.vcf>')
parser.add_argument('-o', help='* Output Prefix', dest='outfile', action='store', metavar='<prefix>')
parser.add_argument('-m', help='Mappability Track', dest='map_track', action='store', metavar='<track.bed>')
parser.add_argument('-M', help='Maptrack Min Len', dest='map_track_min_len', action='store', metavar='<int>')
parser.add_argument('-t', help='Targetted Regions', dest='target_reg', action='store', metavar='<regions.bed>')
parser.add_argument('-T', help='Min Region Len', dest='min_reg_len', action='store', metavar='<int>')
parser.add_argument('-c', help='Coverage Filter Threshold', dest='dp_thresh', default=15, action='store',
                    metavar='<int>')
parser.add_argument('-a', help='Allele Freq Filter Threshold', dest='af_thresh', default=0.3, action='store',
                    metavar='<float>')

parser.add_argument('--vcf-out', help="Output Match/FN/FP variants", dest='vcf_out', default=False, action='store_true')
parser.add_argument('--no-plot', help="No plotting", dest='no_plot', default=False, action='store_true')
parser.add_argument('--incl-homs', help="Include homozygous ref calls", dest='include_homs', default=False,
                    action='store_true')
parser.add_argument('--incl-fail', help="Include calls that failed filters", dest='include_fail',
                    default=False,
                    action='store_true')
parser.add_argument('--fast', help="No equivalent variant detection", dest='fast', default=False,
                    action='store_true')

## DEBUG NOTE: (opts,args) is from deprecated "optparse"
##             updated to current syntak
# (opts, args) = parser.parse_args()
args = parser.parse_args()

reference = args.reff
golden_vcf = args.golden_vcf
workflow_vcf = args.workflow_vcf
out_prefix = args.outfile
maptrack = args.map_track
min_read_len = args.map_track_min_len
bedfile = args.target_reg
dp_thresh = int(args.dp_thresh)
af_thresh = float(args.af_thresh)
vcf_out = args.vcf_out
no_plot = args.no_plot
include_homs = args.include_homs
include_fail = args.include_fail
fast = args.fast

## DEBUG NOTE: missing a sys call here
if len(sys.argv[1:]) == 0:
    parser.print_help()
    sys.exit(1)

## DEBUG NOTE: variable name bug here
if args.target_reg is not None:
    min_region_len = int(args.min_reg_len)
else:
    min_region_len = None

if min_read_len is None:
    min_read_len = 0
else:
    min_read_len = int(min_read_len)

if reference is None:
    print('Error: No reference provided.')
    sys.exit(1)
if golden_vcf is None:
    print('Error: No golden VCF provided.')
    sys.exit(1)
if workflow_vcf is None:
    print('Error: No workflow VCF provided.')
    sys.exit(1)
if out_prefix is None:
    print('Error: No output prefix provided.')
    sys.exit(1)
if (bedfile is not None and min_region_len is None) or (bedfile is None and min_region_len is not None):
    print('Error: Both -t and -T must be specified')
    sys.exit(1)

## DEBUG NOTE: enabled decompression of gzip vcf files automatically
def gunzip_shutil(source_filepath, dest_filepath, block_size=65536):
    """
      Original Source: https://stackoverflow.com/questions/52332897/how-to-extract-a-gz-file-in-python
    """
    with gzip.open(source_filepath, 'rb') as s_file, \
        open(dest_filepath, 'wb') as d_file:
            shutil.copyfileobj(s_file, d_file, block_size)

## DEBUG NOTE: enabled compression of vcf files automatically
def gzip_shutil(source_filepath, dest_filepath, block_size=65536):
    with open(source_filepath, 'rb') as s_file, \
        gzip.open(dest_filepath, 'wb') as d_file:
            shutil.copyfileobj(s_file, d_file, block_size)

## DEBUG NOTE: decompress golden_vcf, if necessary
if '.gz' in golden_vcf:
    print("decompressing .gz input file {}".format(golden_vcf))
    golden_vcf_path = Path(golden_vcf)
    outpath = golden_vcf_path.parent
    outfile = golden_vcf_path.stem
    outpath_file = str(outpath / outfile)
    # print(outpath_file)
    gunzip_shutil(golden_vcf, outpath_file)
    golden_vcf = outpath_file
    print("Success: Input file decompressed to {}".format(golden_vcf))
else:
    workflow_vcf_path = False

## DEBUG NOTE: decompress workflow_vcf, if necessary
if '.gz' in workflow_vcf:
    print("decompressing .gz input file {}".format(workflow_vcf))
    workflow_vcf_path = Path(workflow_vcf)
    outpath = workflow_vcf_path.parent
    outfile = workflow_vcf_path.stem
    outpath_file = str(outpath / outfile)
    # print(outpath_file)
    gunzip_shutil(workflow_vcf, outpath_file)
    workflow_vcf = outpath_file
    print("Success: Input file decompressed to {}.".format(workflow_vcf))
else:
    workflow_vcf_path = False

if no_plot is False:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as mpl
    from matplotlib_venn import venn2, venn3
    import warnings

    warnings.filterwarnings("ignore", category=UserWarning, module='matplotlib_venn')

AF_STEPS = 20
AF_KEYS = np.linspace(0.0, 1.0, AF_STEPS + 1)

def quantize_AF(af):
    if af >= 1.0:
        return AF_STEPS
    elif af <= 0.0:
        return 0
    else:
        return int(af * AF_STEPS)

VCF_HEADER = '##fileformat=VCFv4.1\n##reference=' + reference + '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n'

DP_TOKENS = ['DP', 'DPU', 'DPI']  # in the order that we'll look for them

def parse_line(splt, col_dict, col_samp):
    #	check if we want to proceed..
    ra = splt[col_dict['REF']]
    aa = splt[col_dict['ALT']]
    if not (include_homs) and (aa == '.' or aa == '' or aa == ra):
        return None
    if not (include_fail) and (splt[col_dict['FILTER']] != 'PASS' and splt[col_dict['FILTER']] != '.'):
        return None
    
    #	default vals
    cov = None
    qual = DEFAULT_QUAL
    alt_alleles = []
    alt_freqs = [None]
    
    #	any alt alleles?
    alt_split = aa.split(',')
    if len(alt_split) > 1:
        alt_alleles = alt_split

    #	cov
    for dp_tok in DP_TOKENS:
        #	check INFO for DP first
        if 'INFO' in col_dict and dp_tok + '=' in splt[col_dict['INFO']]:
            cov = int(re.findall(re.escape(dp_tok) + r"=[0-9]+", splt[col_dict['INFO']])[0][3:])
        #	check FORMAT/SAMPLE for DP second:
        elif 'FORMAT' in col_dict and len(col_samp):
            format = splt[col_dict['FORMAT']] + ':'
            if ':' + dp_tok + ':' in format:
                dp_ind = format.split(':').index(dp_tok)
                cov = int(splt[col_samp[0]].split(':')[dp_ind])
        if cov is not None:
            break
    
    #	check INFO for AF first
    af = None
    if 'INFO' in col_dict and ';AF=' in ';' + splt[col_dict['INFO']]:
        info = splt[col_dict['INFO']] + ';'
        af = re.findall(r"AF=.*?(?=;)", info)[0][3:]
    #	check FORMAT/SAMPLE for AF second:
    elif 'FORMAT' in col_dict and len(col_samp):
        format = splt[col_dict['FORMAT']] + ':'
        if ':AF:' in format:
            af_ind = splt[col_dict['FORMAT']].split(':').index('AF')
            af = splt[col_samp[0]].split(':')[af_ind]
    
    if af is not None:
        af_splt = af.split(',')
        while (len(af_splt) < len(alt_alleles)):  # are we lacking enough AF values for some reason?
            af_splt.append(af_splt[-1])  # phone it in.
        if len(af_splt) != 0 and af_splt[0] != '.' and af_splt[0] != '':  # missing data, yay
            alt_freqs = [float(n) for n in af_splt]
    else:
        alt_freqs = [None] * max([len(alt_alleles), 1])

    #	get QUAL if it's interesting
    if 'QUAL' in col_dict and splt[col_dict['QUAL']] != '.':
        qual = float(splt[col_dict['QUAL']])
    
    return (cov, qual, alt_alleles, alt_freqs)

def parse_vcf(vcf_filename, ref_name, targ_regions_FL, out_file, out_bool):
    v_hashed = {}
    v_pos_hash = {}
    v_alts = {}
    v_cov = {}
    v_af = {}
    v_qual = {}
    v_targ_len = {}
    n_below_min_r_len = 0
    line_unique = 0  # number of lines in vcf file containing unique variant
    hash_coll = 0  # number of times we saw a hash collision ("per line" so non-unique alt alleles don't get counted multiple times)
    var_filtered = 0  # number of variants excluded due to filters (e.g. hom-refs, qual)
    var_merged = 0  # number of variants we merged into another due to having the same position specified
    col_dict = {}
    col_samp = []
    for line in open(vcf_filename, 'r'):
        if line[0] != '#':
            if len(col_dict) == 0:
                print('\n\nError: VCF has no header?\n' + vcf_filename + '\n\n')
                sys.exit(1)
            
            splt = line[:-1].split('\t')
            if splt[0] == ref_name:
                var = (int(splt[1]), splt[3], splt[4])
                targ_ind = bisect.bisect(targ_regions_FL, var[0])
                if targ_ind % 2 == 1:
                    targ_Len = targ_regions_FL[targ_ind] - targ_regions_FL[targ_ind - 1]
                    if (bedfile is not None and targ_Len >= min_region_len) or bedfile is None:
                        pl_out = parse_line(splt, col_dict, col_samp)
                        if pl_out is None:
                            var_filtered += 1
                            continue
                        (cov, qual, aa, af) = pl_out
                        
                        if var not in v_hashed:
                            v_pos = var[0]
                            if v_pos in v_pos_hash:
                                if len(aa) == 0:
                                    aa = [var[2]]
                                aa.extend([n[2] for n in v_hashed.keys() if n[0] == v_pos])
                                var_merged += 1
                            v_pos_hash[v_pos] = 1
                            
                            if len(aa):
                                all_vars = [(var[0], var[1], n) for n in aa]
                                for i in range(len(all_vars)):
                                    v_hashed[all_vars[i]] = 1
                                    # if all_vars[i] not in v_alts:
                                    #	v_alts[all_vars[i]] = []
                                    # v_alts[all_vars[i]].extend(all_vars)
                                    v_alts[all_vars[i]] = all_vars
                            
                            else:
                                v_hashed[var] = 1
                            
                            if cov is not None:
                                v_cov[var] = cov
                            v_af[var] = af[0]  # only use first AF, even if multiple. fix this later?
                            v_qual[var] = qual
                            v_targ_len[var] = targ_Len
                            line_unique += 1
                        
                        else:
                            hash_coll += 1
                    else:
                        n_below_min_r_len += 1
        
        elif line[1] != '#': 
            cols = line[1:-1].split('\t')
            for i in range(len(cols)):
                if 'FORMAT' in col_dict:
                    col_samp.append(i)
                col_dict[cols[i]] = i
            ## DEBUG NOTE: This code ends up writing the field headers 
            ##             for every single contig in the ref, yikes!
            # if vcf_out and out_bool:
            #     print(line)
                # out_bool = False
                # out_file.write(line)
    return (
        v_hashed, v_alts, v_cov, v_af, v_qual, v_targ_len, n_below_min_r_len, line_unique, var_filtered, var_merged,
        hash_coll)

def condense_by_pos(list_in):
    var_list_of_interest = [n for n in list_in]
    ind_count = {}
    for n in var_list_of_interest:
        c = n[0]
        if c not in ind_count:
            ind_count[c] = 0
        ind_count[c] += 1
    # non_unique_dict = {n:[] for n in sorted(ind_count.keys()) if ind_count[n] > 1}		# the source 2.7 way
    non_unique_dict = {}
    for n in sorted(ind_count.keys()):
        if ind_count[n] > 1:
            non_unique_dict[n] = []
    del_list = []
    for i in range(len(var_list_of_interest)):
        if var_list_of_interest[i][0] in non_unique_dict:
            non_unique_dict[var_list_of_interest[i][0]].append(var_list_of_interest[i])
            del_list.append(i)
    del_list = sorted(del_list, reverse=True)
    for di in del_list:
        del var_list_of_interest[di]
    for v in non_unique_dict.values():
        var = (v[0][0], v[0][1], ','.join([n[2] for n in v[::-1]]))
        var_list_of_interest.append(var)
    return var_list_of_interest

def main():
    global bedfile
    ref = []
    f = open(reference, 'r')
    n_lines = 0
    prev_r = None
    prev_p = None
    ref_inds = []
    sys.stdout.write('\nindexing reference fasta... ')
    sys.stdout.flush()
    tt = time.time()
    while 1:
        n_lines += 1
        data = f.readline()
        if not data:
            ref_inds.append((prev_r, prev_p, f.tell() - len(data)))
            break
        if data[0] == '>':
            if prev_p is not None:
                ref_inds.append((prev_r, prev_p, f.tell() - len(data)))
            prev_p = f.tell()
            prev_r = data[1:-1]
    print('{0:.3f} (sec)'.format(time.time() - tt))
    # ref_inds = [('chrM', 6, 16909), ('chr1', 16915, 254252549), ('chr2', 254252555, 502315916), ('chr3', 502315922, 704298801), ('chr4', 704298807, 899276169), ('chr5', 899276175, 1083809741), ('chr6', 1083809747, 1258347116), ('chr7', 1258347122, 1420668559), ('chr8', 1420668565, 1569959868), ('chr9', 1569959874, 1713997574), ('chr10', 1713997581, 1852243023), ('chr11', 1852243030, 1989949677), ('chr12', 1989949684, 2126478617), ('chr13', 2126478624, 2243951900), ('chr14', 2243951907, 2353448438), ('chr15', 2353448445, 2458030465), ('chr16', 2458030472, 2550192321), ('chr17', 2550192328, 2633011443), ('chr18', 2633011450, 2712650243), ('chr19', 2712650250, 2772961813), ('chr20', 2772961820, 2837247851), ('chr21', 2837247858, 2886340351), ('chr22', 2886340358, 2938671016), ('chrX', 2938671022, 3097046994), ('chrY', 3097047000, 3157608038)]

    zt_v = 0  # total golden variants
    zt_w = 0  # total workflow variants
    zn_p = 0  # total perfect matches
    zf_p = 0  # total false positives
    zn_f = 0  # total false negatives
    zn_e = 0  # total equivalent variants detected
    zgF = 0  # total golden variants that were filtered and excluded
    zgR = 0  # total golden variants that were excluded for being redundant
    zgM = 0  # total golden variants that were merged into a single position
    zwF = 0  # total workflow variants that were filtered and excluded
    zwR = 0  # total workflow variants that were excluded for being redundant
    zwM = 0  # total workflow variants that were merged into a single position
    if bedfile is not None:
        zb_m = 0

    mappability_vs_FN = {0: 0,
                         1: 0}  # [0] = # of FNs that were in mappable regions, [1] = # of FNs that were in unmappable regions
    coverage_vs_FN = {}  # [C] = # of FNs that were covered by C reads
    allele_bal_vs_FN = {}  # [AF] = # of FNs that were heterozygous genotypes with allele freq AF (quantized to multiples of 1/AF_STEPS)
    for n in AF_KEYS:
        allele_bal_vs_FN[n] = 0

    #
    #	read in mappability track
    #
    mappability_tracks = {}  # indexed by chr string (e.g. 'chr1'), has boolean array
    prev_Ref = ''
    relevant_regions = []
    if maptrack is not None:
        mtf = open(maptrack, 'r')
        for line in mtf:
            splt = line.strip().split('\t')
            if prev_Ref != '' and splt[0] != prev_Ref:
                # do stuff
                if len(relevant_regions):
                    my_track = [0] * (relevant_regions[-1][1] + 100)
                    for r in relevant_regions:
                        for ri in range(r[0], r[1]):
                            my_track[ri] = 1
                    mappability_tracks[prev_Ref] = [n for n in my_track]
                #
                relevant_regions = []
            if int(splt[3]) >= min_read_len:
                relevant_regions.append((int(splt[1]), int(splt[2])))
            prev_Ref = splt[0]
        mtf.close()
        # do stuff
        if len(relevant_regions):
            my_track = [0] * (relevant_regions[-1][1] + 100)
            for r in relevant_regions:
                for ri in range(r[0], r[1]):
                    my_track[ri] = 1
            mappability_tracks[prev_Ref] = [n for n in my_track]
    #

    #
    #	init vcf output, if desired
    #
    vcfo2 = None
    vcfo3 = None
    global vcfo2_first_time
    global vcfo3_first_time
    vcfo2_first_time = False
    vcfo3_first_time = False
    if vcf_out:
        vcfo2 = open(out_prefix + '_FN.vcf', 'w')
        vcfo3 = open(out_prefix + '_FP.vcf', 'w')
        vcfo2_first_time = True
        vcfo3_first_time = True

    #
    #	data for plotting FN analysis
    #
    set1 = []
    set2 = []
    set3 = []
    var_adj = 0

    #
    #
    #	For each sequence in reference fasta...
    #
    #
    for n_RI in ref_inds:

        ref_name = n_RI[0]
        if not fast:
            f.seek(n_RI[1])
            print('reading ' + ref_name + '...', end=' ')
            my_dat = f.read(n_RI[2] - n_RI[1]).split('\n')
            my_len = sum([len(m) for m in my_dat])
            if sys.version_info >= (2, 7):
                print('{:,} bp'.format(my_len))
            else:
                print('{0:} bp'.format(my_len))
            in_width = len(my_dat[0])
            if len(my_dat[-1]) == 0:  # if last line is empty, remove it.
                del my_dat[-1]
            if in_width * (len(my_dat) - 1) + len(my_dat[-1]) != my_len:
                print('fasta column-width not consistent.')
                print(my_len, in_width * (len(my_dat) - 1) + len(my_dat[-1]))
                for i in range(len(my_dat)):
                    if len(my_dat[i]) != in_width:
                        print(i, len(my_dat[i]), in_width)
                sys.exit(1)
            
            ## DEBUG NOTE: .tomutable method depreciated, updated to current
            # my_dat = Seq(''.join(my_dat)).upper().tomutable()
            my_dat = MutableSeq(''.join(my_dat)).upper()
            my_len = len(my_dat)

        #
        #	Parse relevant targeted regions
        #
        targ_regions_fl = []
        if bedfile is not None:
            bedfile = open(bedfile, 'r')
            for line in bedfile:
                splt = line.split('\t')
                if splt[0] == ref_name:
                    targ_regions_fl.extend((int(splt[1]), int(splt[2])))
            bedfile.close()
        else:
            targ_regions_fl = [-1, MAX_VAL + 1]

        #
        #	Parse vcf files
        #
        sys.stdout.write('comparing variation in ' + ref_name + '... ')
        sys.stdout.flush()
        tt = time.time()

        (correct_hashed, correct_alts, correct_cov, correct_AF, correct_qual, correct_targ_len, correct_below_min_R_len,
         correct_unique, g_filtered, g_merged, g_redundant) = parse_vcf(golden_vcf, ref_name, targ_regions_fl, vcfo2,
                                                                        vcfo2_first_time)
        (workflow_hashed, workflow_alts, workflow_COV, workflow_AF, workflow_qual, workflow_tar_len,
         workflow_below_min_R_len,
         workflow_unique, w_filtered, w_merged, w_redundant) = parse_vcf(workflow_vcf, ref_name, targ_regions_fl, vcfo3,
                                                                         vcfo3_first_time)
        zgF += g_filtered
        zgR += g_redundant
        zgM += g_merged
        zwF += w_filtered
        zwR += w_redundant
        zwM += w_merged

        #
        #	Deduce which variants are FP / FN
        #
        solved_inds = {}
        for var in correct_hashed.keys():
            if var in workflow_hashed or var[0] in solved_inds:
                correct_hashed[var] = 2
                workflow_hashed[var] = 20
                solved_inds[var[0]] = True
        for var in list(correct_hashed.keys()) + list(workflow_hashed.keys()):
            if var[0] in solved_inds:
                correct_hashed[var] = 2
                workflow_hashed[var] = 2
        n_perfect = len(solved_inds)

        # correct_hashed[var] = 1: were not found
        #                    = 2: should be discluded because we were found
        #                    = 3: should be discluded because an alt was found
        not_found = [n for n in sorted(correct_hashed.keys()) if correct_hashed[n] == 1]
        fp_variants = [n for n in sorted(workflow_hashed.keys()) if workflow_hashed[n] == 1]

        #
        #	condense all variants who have alternate alleles and were *not* found to have perfect matches
        #	into a single variant again. These will not be included in the candidates for equivalency checking. Sorry!
        #
        not_found = condense_by_pos(not_found)
        fp_variants = condense_by_pos(fp_variants)

        #
        #	tally up some values, if there are no golden variants lets save some CPU cycles and move to the next ref
        #
        tot_golden_variants = n_perfect + len(not_found)
        total_workflow_variants = n_perfect + len(fp_variants)
        if tot_golden_variants == 0:
            zf_p += len(fp_variants)
            zt_w += total_workflow_variants
            print('{0:.3f} (sec)'.format(time.time() - tt))
            continue

        #
        #	let's check for equivalent variants
        #
        if fast == False:
            del_list_i = []
            del_list_j = []
            regions_to_check = []
            for i in range(len(fp_variants)):
                pos = fp_variants[i][0]
                regions_to_check.append((max([pos - EV_BPRANGE - 1, 0]), min([pos + EV_BPRANGE, len(my_dat) - 1])))
            
            for n in regions_to_check:
                ref_section = my_dat[n[0]:n[1]]

                FP_within = []
                for i in range(len(fp_variants)):
                    m = fp_variants[i]
                    if n[0] < m[0] < n[1]:
                        FP_within.append((m, i))
                FP_within = sorted(FP_within)
                adj = 0
                alt_section = copy.deepcopy(ref_section)
                for (m, i) in FP_within:
                    lr = len(m[1])
                    la = len(m[2])
                    d_pos = m[0] - n[0] + adj
                    alt_section = alt_section[:d_pos - 1] + m[2] + alt_section[d_pos - 1 + lr:]
                    adj += la - lr
                
                nf_within = []
                for j in range(len(not_found)):
                    m = not_found[j]
                    if n[0] < m[0] < n[1]:
                        nf_within.append((m, j))
                nf_within = sorted(nf_within)
                adj = 0
                alt_section2 = copy.deepcopy(ref_section)
                for (m, j) in nf_within:
                    lr = len(m[1])
                    la = len(m[2])
                    d_pos = m[0] - n[0] + adj
                    alt_section2 = alt_section2[:d_pos - 1] + m[2] + alt_section2[d_pos - 1 + lr:]
                    adj += la - lr
                
                if alt_section == alt_section2:
                    for (m, i) in FP_within:
                        if i not in del_list_i:
                            del_list_i.append(i)
                    for (m, j) in nf_within:
                        if j not in del_list_j:
                            del_list_j.append(j)
            
            n_equiv = 0
            for i in sorted(list(set(del_list_i)), reverse=True):
                del fp_variants[i]
            for j in sorted(list(set(del_list_j)), reverse=True):
                del not_found[j]
                n_equiv += 1
            n_perfect += n_equiv
        
        #
        #	Tally up errors and whatnot
        #
        zt_v += tot_golden_variants
        zt_w += total_workflow_variants
        zn_p += n_perfect
        zf_p += len(fp_variants)
        zn_f += len(not_found)
        if fast is False:
            zn_e += n_equiv
        if bedfile is not None:
            zb_m += correct_below_min_R_len
        
        #
        #	try to identify a reason for FN variants:
        #
        
        venn_data = [[0, 0, 0] for n in not_found]  # [i] = (unmappable, low cov, low het)
        for i in range(len(not_found)):
            var = not_found[i]
            no_reason = True

            #	mappability?
            if maptrack is not None:
                if ref_name in mappability_tracks and var[0] < len(mappability_tracks[ref_name]):
                    if mappability_tracks[ref_name][var[0]]:
                        mappability_vs_FN[1] += 1
                        venn_data[i][0] = 1
                        no_reason = False
                    else:
                        mappability_vs_FN[0] += 1

            #	coverage?
            if var in correct_cov:
                c = correct_cov[var]
                if c is not None:
                    if c not in coverage_vs_FN:
                        coverage_vs_FN[c] = 0
                    coverage_vs_FN[c] += 1
                    if c < dp_thresh:
                        venn_data[i][1] = 1
                        no_reason = False
            
            #	heterozygous genotype messing things up?
            # if var in correct_AF:
            #	a = correct_AF[var]
            #	if a != None:
            #		a = AF_KEYS[quantize_AF(a)]
            #		if a not in allele_bal_vs_FN:
            #			allele_bal_vs_FN[a] = 0
            #		allele_bal_vs_FN[a] += 1
            #		if a < AF_THRESH:
            #			venn_data[i][2] = 1

            #	no reason?
            if no_reason:
                venn_data[i][2] += 1
        
        for i in range(len(not_found)):
            if venn_data[i][0]:
                set1.append(i + var_adj)
            if venn_data[i][1]:
                set2.append(i + var_adj)
            if venn_data[i][2]:
                set3.append(i + var_adj)
        var_adj += len(not_found)

        #
        #	if desired, write out vcf files.
        #
        not_found = sorted(not_found)
        fp_variants = sorted(fp_variants)
        if vcf_out:
            for line in open(golden_vcf, 'r'):
                ## DEBUG NOTE: Without the ## lines, unable to open output with bcftools
                if line[0] == '#' and line[1] == '#':
                    vcfo2.write(line)
                ## DEBUG NOTE: Ensure that field header line is only written once
                elif line[0] == '#' and line[1] != '#':
                    vcfo2.write(line)
                # Writing the CHROM POSITION row data 
                elif line[0] != '#':
                    splt = line.split('\t')
                    if splt[0] == ref_name:
                        var = (int(splt[1]), splt[3], splt[4])
                        if var in not_found:
                            vcfo2.write(line)
            for line in open(workflow_vcf, 'r'):
                ## DEBUG NOTE: Without the ## lines, unable to open output with bcftools
                if line[0] == '#' and line[1] == '#':
                    vcfo2.write(line)
                ## DEBUG NOTE: Ensure that field header line is only written once
                elif line[0] == '#' and line[1] != '#':
                    vcfo3.write(line) 
                elif line[0] != '#':
                    splt = line.split('\t')
                    if splt[0] == ref_name:
                        var = (int(splt[1]), splt[3], splt[4])
                        if var in fp_variants:
                            vcfo3.write(line)
        print('{0:.3f} (sec)'.format(time.time() - tt))
    
    #
    #	close vcf output
    #
    print('')
    if vcf_out:
        FN_vcf = out_prefix + '_FN.vcf'
        FP_vcf = out_prefix + '_FP.vcf'
        print(f'SUCCESS: saving false negative positions here ({FN_vcf})')
        print(f'SUCCESS: saving false positive positions here ({FP_vcf})')
        vcfo2.close()
        vcfo3.close()
        ## Compress the FN output vcfs
        zipFN_vcf = FN_vcf + '.gz'
        gzip_shutil(FN_vcf, zipFN_vcf)
        ## Delete the old uncompressed vcf
        FN_vcf_path = Path(FN_vcf)
        if FN_vcf_path.exists():
            print(f"SUCCESS: deleting duplicate, unzipped vcf ({str(FN_vcf_path.name)})")
            FN_vcf_path.unlink()
        ## Compress the FP output vcf
        zipFP_vcf = FP_vcf + '.gz'
        gzip_shutil(FP_vcf, zipFP_vcf)
        ## Delete the old uncompressed vcf
        FP_vcf_path = Path(FP_vcf)
        if FP_vcf_path.exists():
            print(f"SUCCESS: deleting duplicate, unzipped vcf ({str(FP_vcf_path.name)})")
            FP_vcf_path.unlink()
    ## NOTE: Do not uncomment if running multiple comparisions
    #        with the same golden file
    # if '.gz' in str(golden_vcf_path):
    #     outfile = golden_vcf_path.stem
    #     outpath = golden_vcf_path.parent / outfile
    #     print("deleting duplicate file created {}".format(str(outfile)))
    #     outpath.unlink()
    
    ##
    ## DEBUG NOTE: If file was uncompressed, delete duplicate file 
    ## 
    if workflow_vcf_path and '.gz' in str(workflow_vcf_path):
        outfile = workflow_vcf_path.stem
        outpath = workflow_vcf_path.parent / outfile
        print(f"SUCCESS: deleting duplicate, unzipped vcf file ({str(outfile)})")
        outpath.unlink() 
    #
    #	plot some FN stuff
    #
    if no_plot == False:
        n_detected = len(set(set1 + set2 + set3))
        set1 = set(set1)
        set2 = set(set2)
        set3 = set(set3)
        if len(set1):
            s1 = 'Unmappable'
        else:
            s1 = ''
        if len(set2):
            s2 = 'DP < ' + str(dp_thresh)
        else:
            s2 = ''
        # if len(set3): s3 = 'AF < '+str(AF_THRESH)
        if len(set3):
            s3 = 'Unknown'
        else:
            s3 = ''
        mpl.figure(0)
        tstr1 = 'False Negative Variants (Missed Detections)'
        # tstr2 = str(n_detected)+' / '+str(zn_f)+' FN variants categorized'
        tstr2 = ''
        if maptrack is not None:
            v = venn3([set1, set2, set3], (s1, s2, s3))
        else:
            v = venn2([set2, set3], (s2, s3))
        mpl.figtext(0.5, 0.95, tstr1, fontdict={'size': 14, 'weight': 'bold'}, horizontalalignment='center')
        mpl.figtext(0.5, 0.03, tstr2, fontdict={'size': 14, 'weight': 'bold'}, horizontalalignment='center')

        ouf = out_prefix + '_FNvenn.pdf'
        print(f'SUCCESS: saving false negative venn plot here ({ouf})')
        mpl.savefig(ouf)
    
    ##
    ##  Write results to a csv file
    ##
    results_dict = {}
    results_dict['TotalGoldenVariants'] = zt_v
    results_dict['FilteredGoldenVariants'] = zgF
    results_dict['MergedGoldenVariants'] = zgM
    results_dict['RedundantGoldenVariants'] = zgR
    results_dict['TotalWorkflowVariants'] = zt_w
    results_dict['FilteredWorkflowVariants'] = zwF
    results_dict['MergedWorkflowVariants'] = zwM
    results_dict['RedundantWorkflowVariants'] = zwR
    
    if zt_v > 0 and zt_w > 0:
        results_dict['PerfectMatches'] = zn_p
        results_dict['FN_variants'] = zn_f
        results_dict['FP_variants'] = zf_p
    if not fast:
        results_dict['EquivalentVariants'] = zn_e
    if no_plot:
        results_dict['#Unmappable'] = len(set1)
        results_dict['#LowCoverage'] = len(set2)
        results_dict['#Unknown'] = len(set3)
    
    if len(results_dict) > 0:
        results_csv = out_prefix + '_summary_metrics.csv'
        with open(results_csv, 'w') as file:
            writer = csv.writer(file)
            for key, value in results_dict.items():
                writer.writerow([key, value])
    
    #
    #	spit out results to console
    #
    print('\n**********************************\n')
    if bedfile is not None:
        print('ONLY CONSIDERING VARIANTS FOUND WITHIN TARGETED REGIONS\n\n')
    print('Total Golden Variants:  ', zt_v, '\t[', zgF, 'filtered,', zgM, 'merged,', zgR, 'redundant ]')
    print('Total Workflow Variants:', zt_w, '\t[', zwF, 'filtered,', zwM, 'merged,', zwR, 'redundant ]')
    print('')
    if zt_v > 0 and zt_w > 0:
        print('Perfect Matches:', zn_p, '({0:.2f}%)'.format(100. * float(zn_p) / zt_v))
        print('FN variants:    ', zn_f, '({0:.2f}%)'.format(100. * float(zn_f) / zt_v))
        print('FP variants:    ', zf_p)  # ,'({0:.2f}%)'.format(100.*float(zf_p)/zt_w)
    if not fast:
        print('\nNumber of equivalent variants denoted differently between the two vcfs:', zn_e)
    if bedfile is not None:
        print('\nNumber of golden variants located in targeted regions that were too small to be sampled from:', zb_m)
    if fast:
        print(
            "\nWarning! Running with '--fast' means that identical variants denoted differently between the two vcfs will not be detected! The values above may be lower than the true accuracy.")
    # if NO_PLOT:
    if no_plot:
        print('\n#unmappable:  ', len(set1))
        print('#low_coverage:', len(set2))
        print('#unknown:     ', len(set3))
    print('\n**********************************\n')


if __name__ == '__main__':
    main()
