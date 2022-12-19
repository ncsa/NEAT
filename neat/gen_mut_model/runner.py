"""
Creates a mutation model
"""

import json
import os.path
import pathlib
import pickle
import sys

import numpy as np
from numpy import genfromtxt
import pybedtools
from Bio import SeqIO


from pathlib import Path
import logging

from ..models import MutationModel
from .constants import REF_WHITELIST, VALID_NUCL, VALID_TRINUC, VCF_DEFAULT_POP_FREQ, DEF_HOMOZYGOUS_FRQ, DEF_MUT_SUB_MATRIX
from ..common import validate_output_path, validate_input_path, open_input, open_output
from .utils import read_fasta, extract_header, read_and_filter_variants, cluster_list, count_trinucleotides, find_caf

__all__ = [
    "compute_mut_runner"
]

_LOG = logging.getLogger(__name__)


def runner(reference_index, vcf_to_process, vcf_columns: list, outcounts_file, show_trinuc: bool, save_trinuc: bool,
           output, bed: str, human_sample: bool, skip_common: bool):
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
    # total count of insertions, indexed by length
    insert_count = {}
    # total count of insertions, indexed by length
    delete_count = {}
    # tabulate how much non-N reference sequence we've eaten through
    total_reflen = 0
    # detect variants that occur in a significant percentage of the input samples (pos,ref,alt,pop_fraction)
    common_variants = []
    # identify regions that have significantly higher local mutation rates than the average
    high_mut_regions = []

    # Clean up and simplify reference index
    # simplify naming and filter out actual human genomes from scaffolding
    ignore = []
    if human_sample:
        for key in reference_index:
            if key.startswith('chr'):
                key_to_check = int(key[3:])
                if key_to_check not in REF_WHITELIST:
                    ignore.append(key)

    if len(ignore) == len(reference_index):
        _LOG.info("No valid contigs detected. If using whitelist, all contigs may have been filtered out.")
        sys.exit(1)

    # Pre-parsing to find all the matching chromosomes between ref and vcf
    _LOG.info('Processing VCF file...')
    matching_variants, matching_chromosomes = read_and_filter_variants(
        vcf_to_process, vcf_columns, reference_index, bed
    )

    if len(matching_variants) == 0:
        _LOG.info('Found no chromosomes in common between VCF, Fasta, and/or BED. '
                  'Check that files use the same naming convention for chromosomes.')
        sys.exit(1)

    trinuc_ref_count, bed_track_length = count_trinucleotides(reference_index, bed, outcounts_file,
                                                              matching_chromosomes, save_trinuc)
############################
    _LOG.info('Creating mutational model...')
    total_reflen = 0
    for contig in matching_chromosomes:
        # Running total of how many non-N bases there are in the reference
        total_reflen += len(reference_index[contig].seq) - reference_index[contig].seq.count('N')

        # list to be used for counting variants that occur multiple times in file (i.e. in multiple samples)
        vcf_common = []
        vcf_c = []
        variant_to_process = []
        # Create a view that narrows variants list to current ref
        # variants_to_process = matching_variants[matching_variants['CHROM'] == contig]
        ## code for variants to process ##
        for i in range(len(matching_variants)):
            if matching_variants[i][0] == contig:
                variant_to_process.append(matching_variants[i])

        ref_sequence = reference_index[contig].seq

        insert_variants = []
        delete_variants = []
        snp_variants = []

        # Sequential Variants Processing ##
        # new code start
        for variant in range(len(variant_to_process)):
            # ref > Alt then it's a deletion, if the ref < alt, it's an insertion
            if len(variant_to_process[variant][3]) != len(variant_to_process[variant][2]):
                if len(variant_to_process[variant][3]) < len(variant_to_process[variant][2]):
                    delete_variants.append((variant_to_process[variant]))
                else:
                    insert_variants.append(variant_to_process[variant])
            if (len(variant_to_process[variant][2]) == 1) & (len(variant_to_process[variant][3]) == 1):
                snp_variants.append(variant_to_process[variant])
        
        
        if len(snp_variants) != 0:
            for i in range(len(snp_variants)):
                analyze = str(ref_sequence[int(snp_variants[i][1]) - 1: int(snp_variants[i][1]) + 2])
                #analyze is a trinuc from the refrence
                if analyze not in VALID_TRINUC:
                    continue
                if snp_variants[i][2] == analyze[1]:
                    t_ref = analyze
                    #t_alt changes the middle nuc to the ALT
                    t_alt = analyze[0] + snp_variants[i][3] + analyze[2]
                    if t_alt not in VALID_TRINUC:
                        continue

                    key = (t_ref, t_alt)
                    if key not in trinuc_transition_count:
                        trinuc_transition_count[key] = 0
                    trinuc_transition_count[key] += 1
                    snp_count += 1

                    # tracks transition probability
                    key2 = (str(snp_variants[i][3]), str(snp_variants[i][4]))
                    if key2 not in snp_transition_count:
                        snp_transition_count[key2] = 0
                    snp_transition_count[key2] += 1
                    my_pop_freq = find_caf(str(snp_variants[i][4]))
                    vcf_common.append((snp_variants[i][0], snp_variants[i][1], snp_variants[i][2], snp_variants[i][3], my_pop_freq))
                    #vcf_common.append((snp_variants[i][1], snp_variants[i][2], snp_variants[i][2], snp_variants[i][3], my_pop_freq))
                else:
                    _LOG.error('Ref allele in variant call does not match reference.\n')
                    sys.exit(1)

        if len(delete_variants) != 0:
            for i in range(len(delete_variants)):
                del_len = len(delete_variants[i][3]) - len(delete_variants[i][2])
                if del_len not in delete_count:
                    delete_count[del_len] = 0
                delete_count[del_len] += 1
                my_pop_freq = find_caf(delete_variants[i][4])
                vcf_common.append((delete_variants[i][0], delete_variants[i][1], delete_variants[i][2], delete_variants[i][3], my_pop_freq)) 
    ####doube refrence used again???
 
        if len(insert_variants) != 0:
            for i in range(len(insert_variants)):
                insert_len = len(insert_variants[i][3]) - len(insert_variants[i][2])
                if insert_len not in insert_count:
                    insert_count[insert_len] = 0
                insert_count[del_len] += 1
                my_pop_freq = find_caf(insert_variants[i][4])
                vcf_common.append((insert_variants[i][0], insert_variants[i][1], insert_variants[i][2], insert_variants[i][3], my_pop_freq)) 

        # if we didn't find anything, skip ahead along to the next reference sequence
        if not len(vcf_common):
            _LOG.info(f'Found no variants for this reference {contig}.')
            continue

        # identify common mutations
        percentile_var = 95
        min_value = np.percentile([n[4] for n in vcf_common], percentile_var)
        # min = np.percentile([n[4] for n in vcf_c], percentile_var)
        for k in sorted(vcf_common):
            if k[4] >= min_value:
                common_variants.append((contig, k[0], k[1], k[3], k[4]))
        vcf_common = {(n[0], n[1], n[2], n[3]): n[4] for n in vcf_common}

        # identify areas that have contained significantly higher random mutation rates
        dist_thresh = 2000
        percentile_clust = 97
        scaler = 1000
        # identify regions with disproportionately more variants in them
        variant_pos = sorted([n[1] for n in vcf_common.keys()])
        #n0 -> n1
        clustered_pos = cluster_list(variant_pos, dist_thresh)
        by_len = [(len(clustered_pos[i]), min(clustered_pos[i]), max(clustered_pos[i]), i) for i in
                  range(len(clustered_pos))]

        candidate_regions = []
        for n in by_len:
            ref_scalar = int((int(n[1]) - dist_thresh) / float(scaler)) * scaler
            alt_scalar = int((int(n[2]) + dist_thresh) / float(scaler)) * scaler
            candidate_regions.append((n[0] / float(alt_scalar - ref_scalar), max([0, ref_scalar]),
                                      min([len(ref_sequence), alt_scalar])))
        minimum_value = np.percentile([n[0] for n in candidate_regions], percentile_clust)
        for n in candidate_regions:
            if n[0] >= minimum_value:
                high_mut_regions.append((contig, n[1], n[2], n[0]))
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
    total_var = snp_count + sum(insert_count.values()) + sum(delete_count.values())
    if total_var == 0:
        _LOG.error('Error: No valid variants were found, model could not be created. '
                   'Check that names are compatible.')
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

        trinuc_mut_prob[trinuc] = my_count / float(trinuc_ref_count[trinuc])

        for k in sorted(trinuc_transition_count.keys()):
            if k[0] == trinuc:
                trinuc_trans_probs[k] = trinuc_transition_count[k] / float(my_count)

    for n1 in VALID_NUCL:
        rolling_tot = sum([snp_transition_count[(n1, n2)] for n2 in VALID_NUCL if (n1, n2) in snp_transition_count])
        for n2 in VALID_NUCL:
            key2 = (n1, n2)
            if key2 in snp_transition_count:
                snp_trans_freq[key2] = snp_transition_count[key2] / float(rolling_tot)

    # compute average snp and indel frequencies
    snp_freq = snp_count / float(total_var)
    in_freq = sum(insert_count.values()) / float(total_var)
    del_freq = sum(delete_count.values()) / float(total_var)
    # avg_indel_freq = 1. - snp_freq
    # indel_freq = {k: (indel_count[k] / float(total_var)) / avg_indel_freq for k in indel_count.keys()}

    if bed:
        avg_mut_rate = total_var / bed_track_length
    else:
        avg_mut_rate = total_var / float(total_reflen)

    #	if values weren't found in data, appropriately append null entries
    print_trinuc_warning = False
    for trinuc in VALID_TRINUC:
        trinuc_mut = [trinuc[0] + n + trinuc[2] for n in VALID_NUCL if n != trinuc[1]]
        if trinuc not in trinuc_mut_prob:
            trinuc_mut_prob[trinuc] = 0.
            print_trinuc_warning = True
        for trinuc2 in trinuc_mut:
            if (trinuc, trinuc2) not in trinuc_trans_probs:
                trinuc_trans_probs[(trinuc, trinuc2)] = 0.
                print_trinuc_warning = True
    if print_trinuc_warning:
        _LOG.warning('Warning: Some trinucleotides transitions were not encountered in the input dataset, '
                     'probabilities of 0.0 have been assigned to these events.')

    #	print some stuff

    if show_trinuc:
        for k in sorted(trinuc_mut_prob.keys()):
            print('p(' + k + ' mutates) =', trinuc_mut_prob[k])

        for k in sorted(trinuc_trans_probs.keys()):
            print('p(' + k[0] + ' --> ' + k[1] + ' | ' + k[0] + ' mutates) =', trinuc_trans_probs[k])

        for k in sorted(in_freq.keys()):
            print('p(ins length = ' + str(abs(k)) + ' | insertion occurs) =', in_freq[k])

        for k in sorted(del_freq.keys()):
            print('p(del length = ' + str(abs(k)) + ' | deletion occurs) =', del_freq[k])

        for k in sorted(snp_trans_freq.keys()):
            print('p(' + k[0] + ' --> ' + k[1] + ' | SNP occurs) =', snp_trans_freq[k])

    print(f'p(snp)   = {snp_freq}')
    print(f'p(insertion) = {in_freq}')
    print(f'p(deletion) = {del_freq}')
    print(f'overall average mut rate: {avg_mut_rate}')
    print(f'total variants processed: {total_var}')

    # save variables to file
    # if skip_common:
    #     out = {'AVG_MUT_RATE': avg_mut_rate,
    #            'SNP_FREQ': snp_freq,
    #            'SNP_TRANS_FREQ': snp_trans_freq,
    #            'INDEL_FREQ': indel_freq,
    #            'TRINUC_MUT_PROB': trinuc_mut_prob,
    #            'TRINUC_TRANS_PROBS': trinuc_trans_probs}
    # else:
    #     out = {'AVG_MUT_RATE': avg_mut_rate,
    #            'SNP_FREQ': snp_freq,
    #            'SNP_TRANS_FREQ': snp_trans_freq,
    #            'INDEL_FREQ': indel_freq,
    #            'TRINUC_MUT_PROB': trinuc_mut_prob,
    #            'TRINUC_TRANS_PROBS': trinuc_trans_probs,
    #            'COMMON_VARIANTS': common_variants,
    #            'HIGH_MUT_REGIONS': high_mut_regions}

    # # Trying protocol = 4 to maintain backward compatability.
    # pickle.dump(out, gzip.open(out_file, "w"), protocol=4)


    # SNPs = SnvModel(trinuc_trans_matrices=trinuc_trans_probs, trinuc_trans_bias=snp_trans_freq)
    # Insertions = InsertionModel(insert_len_model=in_freq.keys())
    # Deletion = DeletionModel(deletion_len_model=del_freq.keys())

    mut_model = MutationModel()
    mut_model.avg_mut_rate = avg_mut_rate
    mut_model.homozygous_freq = DEF_HOMOZYGOUS_FRQ
    mut_model.variant_probs = {'SNPs': snp_freq, 'Insertions': in_freq, 'Deletions': del_freq}
    mut_model.transition_matrix = DEF_MUT_SUB_MATRIX
    mut_model.trinuc_trans_matrices = trinuc_trans_probs
    mut_model.trinuc_trans_bias = snp_trans_freq
    mut_model.insert_len_model = insert_count
    mut_model.deletion_len_model = delete_count

    print('\nSaving model...')
    with open_output(output, 'w+') as outfile:
        pickle.dump(mut_model, outfile)


def compute_mut_runner(reference: str | Path, mutations: str | Path, bed: str | Path, outcounts: str | Path,
                       show_trinuc: bool, save_trinuc:bool, human_sample: bool, skip_common: bool, output: str | Path,
                       overwrite_output: bool):
    """
    :param reference: (REQ)
    :param mutations: (REQ)
    :param bed: (OPT) 
    :param outcounts: (OPT) 
    :param show_trinuc: (OPT) 
    :param save_trinuc: (OPT) 
    :param human_sample: (OPT) 
    :param skip_common: (OPT) 
    :param output:
    :param overwrite_output: True to overwrite output
    """

    # if not os.path.isfile(reference):
    #     print(f'{PROG} - Input reference is not a file: {reference}')
    #     sys.exit(1)
    validate_input_path(reference)

    # if not os.path.isfile(mutations):
    #     print(f'{PROG} - Input VCF is not a file: {mutations}')
    #     sys.exit(1)
    validate_input_path(mutations)

    if bed:
    #     if not os.path.isfile(bed):
    #         print(f'{PROG} - Input BED is not a file: {bed}')
    #         sys.exit(1)
        validate_input_path(bed)

    if outcounts:
    #     if not os.path.isfile(outcounts):
    #         print(f'{PROG} - Trinucleotide counts file {str(outcounts)} does not exist.')
    #         sys.exit(1)
        validate_input_path(outcounts)

    print('Processing reference...')
    reference_index = read_fasta(reference)

    vcf_header = extract_header(mutations)
    vcf_columns = vcf_header[-1]

    vcf_to_process = pathlib.Path(mutations)

    if bed:
        vcf_columns = ['bed_chr', 'bed_pos1', 'bed_pos2'] + vcf_columns
        bed_file = pybedtools.BedTool(bed)
        # used bedtools to intersect the bed and vcf. This will require further processing.
        # The fn at the end extracts the filename, which is what the function expects.
        # Also converts to pathlib path.
        _LOG.info('Intersecting bed and vcf.')
        # TODO rewrite to remove pybedtools dependency
        vcf_to_process = pathlib.Path(bed_file.intersect(mutations, wb=True).moveto('temp.vcf').fn)
        _LOG.info('Created temp vcf for processing.')

    output_prefix = output
    output = Path(output_prefix + '.pickle.gz')
    validate_output_path(output, overwrite=overwrite_output)

    outcounts_file = pathlib.Path(f'{output}.counts.gz').resolve()
    ##

    runner(reference_index, vcf_to_process, vcf_columns, outcounts_file, show_trinuc, save_trinuc, output,
         bed, human_sample, skip_common)

    if os.path.exists('temp.vcf'): 
        os.remove('temp.vcf')

    _LOG.info(f'Complete! Use {output} as input into gen_reads_runner.py.')
