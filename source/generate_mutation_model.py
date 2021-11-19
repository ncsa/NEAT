#! /usr/bin/env python

import argparse
import os.path
import sys
import copy
import gzip
import pathlib
import re
import Bio
from Bio import SeqIO
import numpy as np
import pickle
import pybedtools
import tempfile
import json
import io
import pandas as pd
from memory_profiler import profile

from constants_and_models import HUMAN_WHITELIST, ALL_TRI, ALL_IND, ALLOWED_NUCL


def read_fasta(fasta_file):
    return SeqIO.index(fasta_file, 'fasta')


def extract_header(vcf_file):
    if vcf_file.suffix == '.gz':
        f = gzip.open(vcf_file, 'r')
    else:
        f = open(vcf_file, 'r')

    ret = []
    for line in f:
        if line.startswith('##'):
            ret.append(line.strip())
        elif line.startswith('#CHROM'):
            temp = line.strip().strip("#").split('\t')
            ret.append(temp)
            break

    if not ret:
        print(f'{PROG} - No header found: invalid VCF file.')
        sys.exit(1)

    return ret


def read_and_filter_variants(vcf_file, column_names: list, reference_idx, input_bed: str):

    variants = pd.read_csv(vcf_file, comment="#", sep="\t", header=None)

    if input_bed:
        variants = variants.iloc[:, 3:]
        # In case there were any mismatches (i.e., items in the bed that are not in the vcf), pybedtools
        # would already have pointed them out, so here we drop them
        variants = variants.drop(variants[variants[3] == '.'].index)

    # If there is a problem with the columns, pandas will throw an error here.
    variants = variants.set_axis(column_names, axis=1)

    variant_chroms = variants['CHROM'].unique()
    matching_chroms = []
    for ref_name in reference_idx.keys():
        if ref_name in variant_chroms:
            matching_chroms.append(ref_name)

    ret = variants[variants['CHROM'].isin(matching_chroms)]

    # We'll go ahead and filter out multiple alts, and variants where both REF and ALT are
    # more than one base. This was done to make the trinucelotide context counts make more sense.
    multi_alts = ret[ret['ALT'].str.contains(',')].index
    complex_vars = ret[(ret['REF'].apply(len) > 1) &
                       (ret['ALT'].apply(len) > 1)].index
    ret = ret.drop(multi_alts)
    ret = ret.drop(complex_vars)
    return ret, matching_chroms


def cluster_list(list_to_cluster: list, delta: float) -> list:
    """
    Clusters a sorted list
    :param list_to_cluster: a sorted list
    :param delta: the value to compare list items to
    :return: a clustered list of values
    """
    out_list = [[list_to_cluster[0]]]
    previous_value = list_to_cluster[0]
    current_index = 0
    for item in list_to_cluster[1:]:
        if item - previous_value <= delta:
            out_list[current_index].append(item)
        else:
            current_index += 1
            out_list.append([])
            out_list[current_index].append(item)
        previous_value = item
    return out_list


def count_trinucleotides(reference_idx, input_bed, trinuc_counts, matching_chroms: list, save_trinuc_file: bool):

    # how many times do we observe each trinucleotide in the reference (and input bed region, if present)?
    trinuc_ref_count = {}

    # Count Trinucleotides in reference, based on bed or not
    print(f'{PROG} - Counting trinucleotides in reference...')
    # Count the total number of bases spanned
    track_len = 0

    if input_bed:
        print(f"{PROG} - since you're using a bed input, we have to count trinucs in bed region even if "
              f"you already have a trinuc count file for the reference...")
        if pathlib.Path(input_bed).suffix == ".gz":
            f = gzip.open(input_bed, 'r')
        else:
            f = open(input_bed, 'r')
        for line in f:
            if line.startswith('#'):
                continue
            record = line.strip().split('\t')
            track_len += int(record[2]) - int(record[1]) + 1
            if record[0] in reference_idx.keys():
                for i in range(int(record[1]), int(record[2]) - 1):
                    trinuc = reference_idx[record[0]][i:i + 3].seq
                    if trinuc not in ALL_TRI:
                        continue
                    if trinuc not in trinuc_ref_count:
                        trinuc_ref_count[trinuc] = 0
                    trinuc_ref_count[trinuc] += 1
        if save_trinuc_file:
            print(f'{PROG} - Warning: since we are using bed input, no trinuc file will be saved.')

    # Solution to attribute error (needs to be checked)

    elif not trinuc_counts.is_file():
        for ref_name in matching_chroms:
            sub_seq = reference_idx[ref_name].seq
            for trinuc in ALL_TRI:
                if ref_name not in trinuc_ref_count:
                    trinuc_ref_count[ref_name] = {}
                if trinuc not in trinuc_ref_count:
                    trinuc_ref_count[ref_name][trinuc] = 0
                trinuc_ref_count[ref_name][trinuc] += sub_seq.count_overlap(trinuc)
        if save_trinuc_file:
            with open(trinuc_counts, 'w') as countfile:
                print(f'{PROG} - Saving trinuc counts to file...')
                countfile.write(json.dumps(trinuc_ref_count))

    else:
        print(f'{PROG} - Found counts file, {trinuc_counts}, using that.')
        trinuc_ref_count = json.load(open(trinuc_counts))
        if save_trinuc_file:
            print(f'{PROG} - Warning: existing trinucelotide file will not be changed or overwritten.')

    return trinuc_ref_count, track_len


def find_caf(candidate_field: str) -> float:
    # if parsing a dbsnp vcf, and no CAF= is found in info tag, use this as default val for population freq
    vcf_default_pop_freq = 0.00001
    info_split = [a.split('=') for a in candidate_field.split(';')]
    for item in info_split:
        if item[0].upper() == 'CAF':
            return float(item[1])
    return vcf_default_pop_freq


def main(reference_idx, vcf_file, columns: list, trinuc_count_file, display_counts: bool,
         out_file, input_bed: str, human_flag: bool, common_variants: bool):
    """
    This function generates the mutation model suitable for use in gen_reads. At the moment it must be run as a
    separate utility.
    """
    # [(trinuc_a, trinuc_b)] = # of times we observed a mutation from trinuc_a into trinuc_b
    trinuc_transition_count = {}

    # total count of SNPs
    snp_count = 0
    # overall SNP transition probabilities
    snp_transition_count = {}
    # total count of indels, indexed by length
    indel_count = {}
    # tabulate how much non-N reference sequence we've eaten through
    total_reflen = 0
    # detect variants that occur in a significant percentage of the input samples (pos,ref,alt,pop_fraction)
    common_variants = []
    # tabulate how many unique donors we've encountered (this is useful for identifying common variants)
    total_donors = {}
    # identify regions that have significantly higher local mutation rates than the average
    high_mut_regions = []

    # Clean up and simplify reference index
    # simplify naming and filter out actual human genomes from scaffolding
    if human_flag:
        for key in reference_idx.keys():
            if ' ' in key:
                key_to_check = key.split(' ')[0]
            elif '|' in key:
                key_to_check = key.split('|')[0]
            else:
                key_to_check = key
            if key_to_check not in HUMAN_WHITELIST:
                del reference_idx[key]

    if not reference_idx:
        print(f"{PROG} - No valid contigs detected. If using whitelist, all contigs may have been filtered out.")
        print(f"Note: Names must follow convention of starting with the contig name, "
              "followed by a space or pipe, then supplemental information (>chr1 dna:chromosome). "
              "Check contig names and try again.")
        sys.exit(1)

    # Pre-parsing to find all the matching chromosomes between ref and vcf
    print(f'{PROG} - Processing VCF file...')
    matching_variants, matching_chromosomes = read_and_filter_variants(vcf_file, columns, reference_idx, input_bed)

    # Check to make sure there are some matches, that not everything got filtered out.
    if matching_variants.empty:
        print(f"{PROG} - Found no chromosomes in common between VCF, Fasta, and/or BED. "
              f'Check that files use the same naming convention for chromosomes.')
        sys.exit(1)

    # Starting position of the actual reference, since vcf is 1-based.
    matching_variants['chr_start'] = matching_variants['POS'] - 1

    trinuc_ref_count, bed_track_length = count_trinucleotides(reference_idx, input_bed, trinuc_count_file,
                                                              matching_chromosomes, save_trinuc)

    print(f'{PROG} - Creating mutational model...')
    total_reflen = 0
    for ref_name in matching_chromosomes:
        # Running total of how many non-N bases there are in the reference
        total_reflen += len(reference_idx[ref_name].seq) - reference_idx[ref_name].seq.count('N')

        # list to be used for counting variants that occur multiple times in file (i.e. in multiple samples)
        command_variants = []

        # Create a view that narrows variants list to current ref
        variants_to_process = matching_variants[matching_variants['CHROM'] == ref_name]
        ref_sequence = reference_idx[ref_name].seq

        # Process the variant table
        indel_variants = variants_to_process[variants_to_process['ALT'].apply(len) !=
                                             variants_to_process['REF'].apply(len)]

        snp_variants = variants_to_process[(variants_to_process['REF'].apply(len) == 1) &
                                           (variants_to_process['ALT'].apply(len) == 1)]

        if not snp_variants.empty:
            # only consider positions where ref allele in vcf matches the nucleotide in our reference
            for index, row in snp_variants.iterrows():
                # -2 for POS because the coordinates in a vcf are 1-based,
                # and we want the nucleotide just before the REF base
                trinuc_to_analyze = str(ref_sequence[row['chr_start'] - 1: row['chr_start'] + 2])
                if trinuc_to_analyze not in ALL_TRI:
                    continue
                if row.REF == trinuc_to_analyze[1]:
                    trinuc_ref = trinuc_to_analyze
                    trinuc_alt = trinuc_to_analyze[0] + snp_variants.loc[index, 'ALT'] + trinuc_to_analyze[2]
                    if trinuc_alt not in ALL_TRI:
                        continue

                    key = (trinuc_ref, trinuc_alt)
                    if key not in trinuc_transition_count:
                        trinuc_transition_count[key] = 0
                    trinuc_transition_count[key] += 1
                    snp_count += 1

                    key2 = (str(row.REF), str(row.ALT))
                    if key2 not in snp_transition_count:
                        snp_transition_count[key2] = 0
                    snp_transition_count[key2] += 1

                    my_pop_freq = find_caf(row['INFO'])
                    command_variants.append((row.chr_start, row.REF, row.REF, row.ALT, my_pop_freq))
                else:
                    print(f'{PROG} - Error: ref allele in variant call does not match reference.\n')
                    sys.exit(1)

        # now let's look for indels...
        if not indel_variants.empty:
            for index, row in indel_variants.iterrows():
                if len(row['REF']) != len(row['ALT']):
                    indel_len = len(row['REF']) - len(row['ALT'])
                    if indel_len not in indel_count:
                        indel_count[indel_len] = 0
                    indel_count[indel_len] += 1

                    my_pop_freq = find_caf(row['INFO'])
                    command_variants.append((row.chr_start, row.REF, row.REF, row.ALT, my_pop_freq))

        # if we didn't find anything, skip ahead along to the next reference sequence
        if not len(command_variants):
            print(f'{PROG} - Found no variants for this reference {ref_name}.')
            continue

        # identify common mutations
        percentile_var = 95
        min_value = np.percentile([n[4] for n in command_variants], percentile_var)
        for k in sorted(command_variants):
            if k[4] >= min_value:
                common_variants.append((ref_name, k[0], k[1], k[3], k[4]))
        command_variants = {(n[0], n[1], n[2], n[3]): n[4] for n in command_variants}

        # identify areas that have contained significantly higher random mutation rates
        dist_thresh = 2000
        percentile_clust = 97
        scaler = 1000
        # identify regions with disproportionately more variants in them
        variant_pos = sorted([n[0] for n in command_variants.keys()])
        clustered_pos = cluster_list(variant_pos, dist_thresh)
        by_len = [(len(clustered_pos[i]), min(clustered_pos[i]), max(clustered_pos[i]), i) for i in
                  range(len(clustered_pos))]

        candidate_regions = []
        for n in by_len:
            ref_scalar = int((n[1] - dist_thresh) / float(scaler)) * scaler
            alt_scalar = int((n[2] + dist_thresh) / float(scaler)) * scaler
            candidate_regions.append((n[0] / float(alt_scalar - ref_scalar), max([0, ref_scalar]),
                                      min([len(ref_sequence), alt_scalar])))
        minimum_value = np.percentile([n[0] for n in candidate_regions], percentile_clust)
        for n in candidate_regions:
            if n[0] >= minimum_value:
                high_mut_regions.append((ref_name, n[1], n[2], n[0]))
        # collapse overlapping regions
        for i in range(len(high_mut_regions) - 1, 0, -1):
            if high_mut_regions[i - 1][2] >= high_mut_regions[i][1] and high_mut_regions[i - 1][0] == \
                    high_mut_regions[i][0]:
                # Might need to research a more accurate way to get the mutation rate for this region
                avg_mut_rate = 0.5 * high_mut_regions[i - 1][3] + 0.5 * high_mut_regions[i][
                    3]
                high_mut_regions[i - 1] = (
                    high_mut_regions[i - 1][0], high_mut_regions[i - 1][1], high_mut_regions[i][2], avg_mut_rate)
                del high_mut_regions[i]

    # if for some reason we didn't find any valid input variants, exit gracefully...
    total_var = snp_count + sum(indel_count.values())
    if total_var == 0:
        print(f'{PROG} - Error: No valid variants were found, model could not be created. '
              f'Check that names are compatible.')
        sys.exit(1)

    # COMPUTE PROBABILITIES

    # frequency that each trinuc mutated into anything else
    trinuc_mut_prob = {}
    # frequency that a trinuc mutates into another trinuc, given that it mutated
    trinuc_trans_probs = {}
    # frequency of snp transitions, given a snp occurs.
    snp_trans_freq = {}

    for trinuc in sorted(trinuc_ref_count.keys()):
        my_count = 0
        for k in sorted(trinuc_transition_count.keys()):
            if k[0] == trinuc:
                my_count += trinuc_transition_count[k]

        # Solution to type error (needs to be checked)

        trinuc_ref_list = [v for k, v in trinuc_ref_count[trinuc].items()]

        for i in range(len(trinuc_ref_list)):
            trinuc_mut_prob[trinuc] = my_count / float(trinuc_ref_list[i])

        # Previous code (commented)

        # trinuc_mut_prob[trinuc] = my_count / float(trinuc_ref_count[trinuc])

        for k in sorted(trinuc_transition_count.keys()):
            if k[0] == trinuc:
                trinuc_trans_probs[k] = trinuc_transition_count[k] / float(my_count)

    for n1 in ALLOWED_NUCL:
        rolling_tot = sum([snp_transition_count[(n1, n2)] for n2 in ALLOWED_NUCL if (n1, n2) in snp_transition_count])
        for n2 in ALLOWED_NUCL:
            key2 = (n1, n2)
            if key2 in snp_transition_count:
                snp_trans_freq[key2] = snp_transition_count[key2] / float(rolling_tot)

    # compute average snp and indel frequencies
    snp_freq = snp_count / float(total_var)
    avg_indel_freq = 1. - snp_freq
    indel_freq = {k: (indel_count[k] / float(total_var)) / avg_indel_freq for k in indel_count.keys()}

    if input_bed:
        avg_mut_rate = total_var / bed_track_length
    else:
        avg_mut_rate = total_var / float(total_reflen)

    #	if values weren't found in data, appropriately append null entries
    print_trinuc_warning = False
    for trinuc in ALL_TRI:
        trinuc_mut = [trinuc[0] + n + trinuc[2] for n in ALLOWED_NUCL if n != trinuc[1]]
        if trinuc not in trinuc_mut_prob:
            trinuc_mut_prob[trinuc] = 0.
            print_trinuc_warning = True
        for trinuc2 in trinuc_mut:
            if (trinuc, trinuc2) not in trinuc_trans_probs:
                trinuc_trans_probs[(trinuc, trinuc2)] = 0.
                print_trinuc_warning = True
    if print_trinuc_warning:
        print(f'{PROG} - Warning: Some trinucleotides transitions were not encountered in the input dataset, '
              f'probabilities of 0.0 have been assigned to these events.')

    #	print some stuff

    if display_counts:
        for k in sorted(trinuc_mut_prob.keys()):
          print('p(' + k + ' mutates) =', trinuc_mut_prob[k])

        for k in sorted(trinuc_trans_probs.keys()):
            print('p(' + k[0] + ' --> ' + k[1] + ' | ' + k[0] + ' mutates) =', trinuc_trans_probs[k])

        for k in sorted(indel_freq.keys()):
            if k > 0:
                print('p(ins length = ' + str(abs(k)) + ' | indel occurs) =', indel_freq[k])
            else:
                print('p(del length = ' + str(abs(k)) + ' | indel occurs) =', indel_freq[k])

        for k in sorted(snp_trans_freq.keys()):
            print('p(' + k[0] + ' --> ' + k[1] + ' | SNP occurs) =', snp_trans_freq[k])

    print(f'p(snp)   = {snp_freq}')
    print(f'p(indel) = {avg_indel_freq}')
    print(f'overall average mut rate: {avg_mut_rate}')
    print(f'total variants processed: {total_var}')

    # Convert calculations from what was done to what is needed in gen_reads.
    # This basically does the job of parse_mutation_model
    indel_freq = 1-snp_freq

    # save variables to file
    if common_variants:
        out = {'AVG_MUT_RATE': avg_mut_rate,
               'SNP_FREQ': snp_freq,
               'SNP_TRANS_FREQ': snp_trans_freq,
               'INDEL_FREQ': indel_freq,
               'TRINUC_MUT_PROB': trinuc_mut_prob,
               'TRINUC_TRANS_PROBS': trinuc_trans_probs,
               'COMMON_VARIANTS': common_variants,
               'HIGH_MUT_REGIONS': high_mut_regions}
    else:
        out = {'AVG_MUT_RATE': avg_mut_rate,
               'SNP_FREQ': snp_freq,
               'SNP_TRANS_FREQ': snp_trans_freq,
               'INDEL_FREQ': indel_freq,
               'TRINUC_MUT_PROB': trinuc_mut_prob,
               'TRINUC_TRANS_PROBS': trinuc_trans_probs}

    pickle.dump(out, open(out_file, "wb"))


if __name__ == '__main__':
    PROG = 'generate_mutation_model.py'
    parser = argparse.ArgumentParser(prog=PROG)
    parser.add_argument('reference', type=str, metavar='reference.fa',
                        help="Reference file for organism in fasta format")
    parser.add_argument('mutations', type=str, metavar='mutation.vcf',
                        help="Mutation file for organism in VCF format")
    parser.add_argument('out', type=str, metavar='output_prefix',
                        help="Name of output file (final model will append \'.p\')")
    parser.add_argument('-b', '--bed', type=str, help="Bed file with regions to use in the model")
    parser.add_argument('--outcounts', type=str, help="Path to trinucleotide counts file for reference. Note, this"
                                                      "file will not be used if a bed file is also used.")
    parser.add_argument('--show-trinuc', action='store_true', help='Shows trinucleotide counts, for reference')
    parser.add_argument('--save-trinuc', action='store_true',
                        help='Saves trinucleotide counts to a file called <out>.counts')
    parser.add_argument('--human-sample', action='store_true',
                        help='Only use numbered chroms, X, Y, and MT. '
                             'Omit this flag to include all chroms in reference.')
    parser.add_argument('--skip-common', action='store_true',
                        help="Includes a list of common variants, "
                             "if you want to visualize common variants with plot_mut_model.py.")
    args = parser.parse_args()

    reference = args.reference
    vcf = args.mutations
    out_pickle = args.out
    skip_common = args.skip_common

    # Set bed to None by default. This is important for the main function.
    bed = None
    if args.bed:
        bed = args.bed

    outcounts = None
    if args.outcounts:
        outcounts = pathlib.Path(args.outcounts)
        if not outcounts.is_file():
            print(f'{PROG} - Trinucleotide counts file {str(outcounts)} does not exist.')
            sys.exit(1)

    show_trinuc = args.show_trinuc
    save_trinuc = args.save_trinuc
    is_human = args.human_sample

    if not pathlib.Path(reference).is_file():
        print(f'{PROG} - Input reference is not a file: {reference}')
        sys.exit(1)

    print('Processing reference...')
    reference_index = read_fasta(reference)

    if not pathlib.Path(vcf).is_file():
        print(f'{PROG} - Input VCF is not a file: {vcf}')
        sys.exit(1)

    if bed:
        if not pathlib.Path(bed).is_file():
            print(f'{PROG} - Input BED is not a file: {bed}')

    vcf_header = extract_header(pathlib.Path(vcf))
    vcf_columns = vcf_header[-1]

    vcf_to_process = pathlib.Path(vcf)
    if bed:
        bed_file = pybedtools.BedTool(bed)
        # used bedtools to intersect the bed and vcf. This will require further processing.
        # The fn at the end extracts the filename, which is what the function expects.
        # Also converts to pathlib path.
        print(f'{PROG} - Intersecting bed and vcf.')
        vcf_to_process = pathlib.Path(bed_file.intersect(vcf, wb=True).moveto('temp.vcf').fn)
        print(f'{PROG} - Created temp vcf for processing.')

    outfile = pathlib.Path(f'{out_pickle}.p').resolve()

    if not outfile.parent.is_dir():
        print(f'{PROG} - Unknown parent directory for output: {outfile.resolve().parent}')
        sys.exit(1)

    outcounts_file = pathlib.Path(f'{out_pickle}.counts').resolve()

    main(reference_index, vcf_to_process, vcf_columns, outcounts_file, show_trinuc, outfile,
         bed, is_human, skip_common)

    if os.path.exists('temp.vcf'):
        os.remove('temp.vcf')

    print(f'{PROG} complete! Use {outfile} as input into gen_reads.py.')
