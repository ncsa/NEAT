#! /usr/bin/env python

import argparse
import sys
import pandas as pd
import pathlib
import re
from Bio import SeqIO
import numpy as np
import pickle

from constants_and_models import HUMAN_WHITELIST, ALL_TRI, ALL_IND, ALLOWED_NUCL


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


def main(reference_idx, vcf_file: str, out_pickle_name: str, bed_file: str, show_trinuc: bool,
         save_trinuc_file: bool, human_flag: bool, is_cancer: bool):
    """
    This function generates the mutation model suitable for use in gen_reads. At the moment it must be run as a
    separate utility.
    """
    # if parsing a dbsnp vcf, and no CAF= is found in info tag, use this as default val for population freq
    vcf_default_pop_freq = 0.00001
    # how many times do we observe each trinucleotide in the reference (and input bed region, if present)?
    trinuc_ref_count = {}
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


    outfile = pathlib.Path(f'{out_pickle_name}.p')
    outcounts = pathlib.Path(f'{out_pickle_name}.counts')

    if not outfile.resolve().parent.is_dir():
        print(f'{PROG} - Unknown parent directory for output: {outfile.resolve().parent}')
        sys.exit(1)

    # Process bed file,
    if bed_file:
        print(f'{PROG} - Processing bed file...')
        col_names = ['chrom', 'start', 'end']
        bed_df = pd.read_csv(bed_file, sep='\t', header=None, usecols=range(3), names=col_names)
        # Adding a couple of columns we'll need for later calculations
        bed_df['coords'] = list(zip(bed_df.start, bed_df.end))
        bed_df['track_len'] = bed_df.end - bed_df.start + 1

    # simplify naming and filter out actual human genomes from scaffolding
    if human_flag:
        for key in reference_idx.keys():
            if key not in HUMAN_WHITELIST:
                del reference_idx[key]

    if not reference_idx:
        print(f"{PROG} - No valid contigs detected. Names must follow convention of starting with the contig name, "
              "followed by a space, then supplemental information (>chr1 dna:chromosome). "
              "Check contig names and try again.")
        sys.exit(1)

    # Pre-parsing to find all the matching chromosomes between ref and vcf
    print(f'{PROG} - Processing VCF file...')
    try:
        variants = pd.read_csv(vcf_file, sep='\t', comment='#', index_col=None, header=None)
    except ValueError:
        print(f'{PROG} - VCF must be in standard VCF format with tab-separated columns')
        sys.exit(1)

    # Narrow chromosomes to those matching the reference
    # This is in part to make sure the names match
    # the set function elminates duplicates, and list returns the set to list format.
    variant_chroms = list(set(variants[0].to_list()))
    matching_chromosomes = []
    for ref_name in reference_idx.keys():
        if ref_name in variant_chroms:
            matching_chromosomes.append(ref_name)
        else:
            continue

    # Check to make sure there are some matches
    if not matching_chromosomes:
        print(f"{PROG} - Found no chromosomes in common between VCF and Fasta. "
              f"Please check chromosome names and try again")
        sys.exit(1)

    # Double check that there are matches
    try:
        matching_variants = variants[variants[0].isin(matching_chromosomes)]
    except ValueError:
        print(f'{PROG} - Problem matching variants with reference.')
        sys.exit(1)

    if matching_variants.empty:
        print(f'{PROG} - There is no overlap between reference and variant file. '
              f'This could be a chromosome naming problem')
        sys.exit(1)

    # If multi-alternate alleles are present, we'll ignore this entry for now.
    # TODO investigate a better way to handle this. Would we need to use allele frequency?
    indices_to_ignore = matching_variants[matching_variants['ALT'].str.contains(',')].index
    matching_variants = matching_variants.drop(indices_to_ignore)

    # Rename header in dataframe for processing
    matching_variants = matching_variants.rename(columns={0: "CHROM", 1: 'chr_start', 2: 'ID', 3: 'REF', 4: 'ALT',
                                                          5: 'QUAL', 6: 'FILTER', 7: 'INFO'})

    # Change the indexing by -1 to match source format indexing (0-based)
    matching_variants['chr_start'] = matching_variants['chr_start'] - 1
    matching_variants['chr_end'] = matching_variants['chr_start'] + len(matching_variants['REF'])

    # Process the variant table
    indices_to_indels = \
        matching_variants.loc[matching_variants.ALT.apply(len) != matching_variants.REF.apply(len)].index

    # indels in vcf don't include the preserved first nucleotide, so lets trim the vcf alleles
    ref_values_to_change = matching_variants.loc[indices_to_indels, 'REF'].copy().str[1:]
    alt_values_to_change = matching_variants.loc[indices_to_indels, 'ALT'].copy().str[1:]
    matching_variants.loc[indices_to_indels, "REF"] = ref_values_to_change
    matching_variants.loc[indices_to_indels, "ALT"] = alt_values_to_change
    matching_variants.replace('', '-', inplace=True)

    # if we encounter a multi-np (i.e. 3 nucl --> 3 different nucl), let's skip it for now...

    # Alt and Ref contain no dashes
    no_dashes = matching_variants[
        ~matching_variants['REF'].str.contains('-') & ~matching_variants['ALT'].str.contains('-')].index
    # Alt and Ref lengths are greater than 1
    long_variants = matching_variants[
        (matching_variants['REF'].apply(len) > 1) & (matching_variants['ALT'].apply(len) > 1)].index
    complex_variants = list(set(no_dashes) & set(long_variants))
    matching_variants = matching_variants.drop(complex_variants)

    # This is solely to make regex easier later, since we can't predict where in the line a string will be
    new_info = ';' + matching_variants['INFO'].copy() + ';'
    matching_variants['INFO'] = new_info

    # Now we check that the bed and vcf have matching regions
    # This also checks that the vcf and bed have the same naming conventions and cuts out scaffolding.
    if bed:
        bed_chroms = list(set(bed_df['chrom']))
        matching_bed_keys = list(set(bed_chroms) & set(variant_chroms))
        try:
            matching_bed = bed_df[bed_df['chrom'].isin(matching_bed_keys)]
        except ValueError:
            print(f'{PROG} - Problem matching bed chromosomes to variant file.')
            sys.exti(1)

        if matching_bed.empty:
            print(f"{PROG} - There is no overlap between bed and variant file. "
                  f"This could be a chromosome naming problem")
            sys.exit(1)

    # Count Trinucleotides in reference, based on bed or not
    print(f'{PROG} - Counting trinucleotides in reference...')

    if bed:
        print(f"{PROG} - since you're using a bed input, we have to count trinucs in bed region even if "
              f"you already have a trinuc count file for the reference...")
        for ref_name in matching_chromosomes:
            sub_bed = matching_bed[matching_bed['chrom'] == ref_name]
            sub_regions = sub_bed['coords'].to_list()
            for sr in sub_regions:
                sub_seq = reference_idx[ref_name][sr[0]: sr[1]].seq
                for trinuc in ALL_TRI:
                    if trinuc not in trinuc_ref_count:
                        trinuc_ref_count[trinuc] = 0
                    trinuc_ref_count[trinuc] += sub_seq.count_overlap(trinuc)

    elif not outcounts.is_file():
        for ref_name in matching_chromosomes:
            sub_seq = reference_idx[ref_name].seq
            for trinuc in ALL_TRI:
                if trinuc not in trinuc_ref_count:
                    trinuc_ref_count[trinuc] = 0
                trinuc_ref_count[trinuc] += sub_seq.count_overlap(trinuc)
    else:
        print(f'{PROG} - Found counts file, {outcounts}, using that.')

    # Load and process variants in each reference sequence individually, for memory reasons...
    print(f'{PROG} - Creating mutational model...')
    for ref_name in matching_chromosomes:
        # Count the number of non-N nucleotides for the reference
        total_reflen += len(reference_idx[ref_name].seq) - reference_idx[ref_name].seq.count('N')

        # list to be used for counting variants that occur multiple times in file (i.e. in multiple samples)
        vdat_common = []

        # Create a view that narrows variants list to current ref
        variants_to_process = matching_variants[matching_variants["CHROM"] == ref_name].copy()
        ref_sequence = str(reference_idx[ref_name].seq)

        # we want only snps
        # so, no '-' characters allowed, and chrStart must be same as chrEnd
        snp_df = variants_to_process[~variants_to_process.index.isin(indices_to_indels)]
        snp_df = snp_df.loc[snp_df['chr_start'] == snp_df['chr_end']]
        if bed:
            bed_to_process = matching_bed[matching_bed['chrom'] == ref_name].copy()
            # TODO fix this line (need the intersection of these two, I think)
            snp_df = bed_to_process.join(snp_df)

        if not snp_df.empty:
            # only consider positions where ref allele in vcf matches the nucleotide in our reference
            for index, row in snp_df.iterrows():
                trinuc_to_analyze = str(ref_sequence[row['chr_start'] - 1: row['chr_start'] + 2])
                if trinuc_to_analyze not in ALL_TRI:
                    continue
                if row.REF == trinuc_to_analyze[1]:
                    trinuc_ref = trinuc_to_analyze
                    trinuc_alt = trinuc_to_analyze[0] + snp_df.loc[index, 'ALT'] + trinuc_to_analyze[2]
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

                    my_pop_freq = vcf_default_pop_freq
                    if ';CAF=' in snp_df.loc[index, 'INFO']:
                        caf_str = re.findall(r";CAF=.*?(?=;)", row.INFO)[0]
                        if ',' in caf_str:
                            my_pop_freq = float(caf_str[5:].split(',')[1])
                    vdat_common.append(
                        (row.chr_start, row.REF, row.REF, row.ALT, my_pop_freq))
                else:
                    print(f'{PROG} - Error: ref allele in variant call does not match reference.\n')
                    sys.exit(1)

        # now let's look for indels...
        indel_df = variants_to_process[variants_to_process.index.isin(indices_to_indels)]
        if not indel_df.empty:
            for index, row in indel_df.iterrows():
                if "-" in row.REF:
                    len_ref = 0
                else:
                    len_ref = len(row.REF)
                if "-" in row.ALT:
                    len_alt = 0
                else:
                    len_alt = len(row.ALT)
                if len_ref != len_alt:
                    indel_len = len_alt - len_ref
                    if indel_len not in indel_count:
                        indel_count[indel_len] = 0
                    indel_count[indel_len] += 1

                    my_pop_freq = vcf_default_pop_freq
                    if ';CAF=' in row.INFO:
                        caf_str = re.findall(r";CAF=.*?(?=;)", row.INFO)[0]
                        if ',' in caf_str:
                            my_pop_freq = float(caf_str[5:].split(',')[1])
                    vdat_common.append((row.chr_start, row.REF, row.REF, row.ALT, my_pop_freq))

        # if we didn't find anything, skip ahead along to the next reference sequence
        if not len(vdat_common):
            print(f'{PROG} - Found no variants for this reference {ref_name}.')
            continue

        # identify common mutations
        percentile_var = 95
        min_value = np.percentile([n[4] for n in vdat_common], percentile_var)
        for k in sorted(vdat_common):
            if k[4] >= min_value:
                common_variants.append((ref_name, k[0], k[1], k[3], k[4]))
        vdat_common = {(n[0], n[1], n[2], n[3]): n[4] for n in vdat_common}

        # identify areas that have contained significantly higher random mutation rates
        dist_thresh = 2000
        percentile_clust = 97
        scaler = 1000
        # identify regions with disproportionately more variants in them
        VARIANT_POS = sorted([n[0] for n in vdat_common.keys()])
        clustered_pos = cluster_list(VARIANT_POS, dist_thresh)
        by_len = [(len(clustered_pos[i]), min(clustered_pos[i]), max(clustered_pos[i]), i) for i in
                  range(len(clustered_pos))]
        # Not sure what this was intended to do or why it is commented out. Leaving it here for now.
        # by_len  = sorted(by_len,reverse=True)
        # minLen = int(np.percentile([n[0] for n in by_len],percentile_clust))
        # by_len  = [n for n in by_len if n[0] >= minLen]
        candidate_regions = []
        for n in by_len:
            bi = int((n[1] - dist_thresh) / float(scaler)) * scaler
            bf = int((n[2] + dist_thresh) / float(scaler)) * scaler
            candidate_regions.append((n[0] / float(bf - bi), max([0, bi]), min([len(reference_idx[ref_name]), bf])))
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

    # if we didn't count ref trinucs because we found file, read in ref counts from file now
    if outcounts.is_file():
        print(f'{PROG} - Reading pre-computed trinuc counts {outcounts}')
        f = open(outcounts, 'r')
        for line in f:
            splt = line.strip().split('\t')
            trinuc_ref_count[splt[0]] = int(splt[1])
        f.close()
    # otherwise, save trinuc counts to file, if desired
    elif save_trinuc_file:
        if bed:
            print(f'{PROG} - Unable to save trinuc counts to file because using input bed region...')
        else:
            print(f'{PROG} - Saving trinuc counts to file...')
            f = open(outcounts, 'w')
            for trinuc in sorted(trinuc_ref_count.keys()):
                f.write(trinuc + '\t' + str(trinuc_ref_count[trinuc]) + '\n')
            f.close()

    # if for some reason we didn't find any valid input variants, exit gracefully...
    total_var = snp_count + sum(indel_count.values())
    if total_var == 0:
        print(f'{PROG} - Error: No valid variants were found, model could not be created. '
              f'Check that names are compatible.')
        sys.exit(1)

    ###	COMPUTE PROBABILITIES

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

    if bed:
        track_sum = float(bed_df['track_len'].sum())
        avg_mut_rate = total_var / track_sum
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

    if show_trinuc:
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

    # In the original code, this value is basically hard-coded in (though it is lost among the
    # the other code and hard to see that). So we will hardcode it in here, too.
    # TODO calculate this value?
    homozygous_freq = 0.01
    if is_cancer:
        homozygous_freq = 0.2

    # Convert calculations from what was done to what is needed in gen_reads.
    # This basically does the job of parse_mutation_model
    indel_freq = 1-snp_freq

    # save variables to file
    out = [avg_mut_rate,
           homozygous_freq,
           snp_freq,
           snp_trans_freq,
           indel_freq,
           trinuc_mut_prob,
           trinuc_trans_probs,
           common_variants,
           high_mut_regions]
    pickle.dump(out, open(out_pickle_name, "wb"))


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
    parser.add_argument('--show-trinuc', action='store_true', help='Shows trinucleotide counts for reference')
    parser.add_argument('--save-trinuc', action='store_true',
                        help='Saves trinucleotide counts to a file called <out>.counts')
    parser.add_argument('--human-sample', action='store_true',
                        help='Only use numbered chroms, X, Y, and MT. '
                             'Omit this flag to include all chroms in reference.')
    parser.add_argument('--is-cancer', action='store_true',
                        help='If this is a cancer sample, this flag will help give the correct output')
    args = parser.parse_args()

    reference = args.reference
    vcf = pathlib.Path(args.mutations)
    out_pickle = args.out
    bed = None
    if args.bed:
        bed = pathlib.Path(args.bed)
    show_trinuc = args.show_trinuc
    save_trinuc = args.save_trinuc
    is_human = args.human_sample
    is_cancer = args.is_cancer

    if not pathlib.Path(reference).is_file():
        print(f'{PROG} - Input reference is not a file: {reference}')
        sys.exit(1)

    print('Processing reference...')
    reference_index = SeqIO.index(reference, "fasta")

    if not vcf.is_file():
        print(f'{PROG} - Input VCF is not a file: {vcf}')
        sys.exit(1)

    if bed:
        if not bed.is_file():
            print(f'{PROG} - Input BED is not a file: {bed}')

    main(reference_index, vcf, out_pickle, bed, show_trinuc, save_trinuc, is_human, is_cancer)
