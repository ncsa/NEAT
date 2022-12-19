"""
Utilities used by the generate mutation model function
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

from .constants import VALID_TRINUC
from ..common import open_input, open_output

__all__ = [
    'read_fasta', 'extract_header', 'read_and_filter_variants',
    'cluster_list', 'count_trinucleotides', 'find_caf'
]

_LOG = logging.getLogger(__name__)


def read_fasta(fasta_file):
    return SeqIO.index(fasta_file, 'fasta')


# In common, a way to validate path/open file
def extract_header(vcf_file):
    ret = []
    with open_input(vcf_file) as f:
        for line in f:
            if line.startswith('##'):
                ret.append(line.strip())
            elif line.startswith('#CHROM'):
                temp = line.strip().strip("#").split('\t')
                ret.append(temp)
                break

    if not ret:
        _LOG.error("No header found: invalid VCF file.")
        sys.exit(1)

    return ret


def read_and_filter_variants(vcf_file, column_names: list, reference_index, bed: str):
    """
    finds all the matching chromosomes between ref and vcf(mutation file)
    """

    variant_chroms = []
    matching_chroms = []

    final_data = []

    header_count = 0

    # TODO Read in the bed file. We'll need a list of chromosomes found in the bed, and a list of regions.b
    #  During this filtering check, any variant must be included in the reference AND
    #  be within a region of the bed file.
    with open_input(vcf_file) as vcf:
        for line in vcf:
            if line.split('\t')[0][0] != '#':
                columns = [x for x in line.split('\t') if line.split('\t')[0][0] != '#']

                """
                columns[0]: CHROM
                columns[1]: POS
                columns[2]: ID
                columns[3]: REF
                columns[4]: ALT
                columns[4]: QUAL
                columns[6]: FILTER
                columns[7]: INFO
                """

                # Add chromosome to variant chroms if it isn't already
                chrom = columns[0]
                if chrom not in variant_chroms:
                    variant_chroms.append(chrom)

                # if a chromosome is present in the reference and the vcf file -> add it to matching_chrom list.
                if chrom not in matching_chroms:
                    if chrom in reference_index.keys():
                        matching_chroms.append(chrom)

                # If CHROM is present in matching_chroms, then the variant is a candidate
                if chrom in matching_chroms:
                    # multi-allelic, we'll just take the first one
                    if ',' in columns[4]:
                        variant_alt = columns[4].split(',')[0]
                    else:
                        variant_alt = columns[4]
                    variant_ref = columns[3]

                    add_variant = False
                    # SNP
                    if len(variant_ref) == 1 and len(variant_alt) == 1:
                        add_variant = True
                    # deletion
                    elif len(variant_ref) > 1 and len(variant_alt) == 1:
                        add_variant = True
                    # insertion
                    elif len(variant_ref) == 1 and len(variant_alt) > 1:
                        add_variant = True
                    # Note that Other variants are ommitted
                    if add_variant:
                        final_data.append([chrom, str(int(columns[1]) - 1), variant_ref, variant_alt, columns[7]])

            else:
                header_count += 1

        print(f"Variant chroms not in ref: {list(set(reference_index.keys()) - set(variant_chroms))}")
        print(f"Matching chroms: {matching_chroms}")

    return final_data, matching_chroms


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
        if int(item) - int(previous_value) <= delta:
            out_list[current_index].append(item)
        else:
            current_index += 1
            out_list.append([])
            out_list[current_index].append(item)
        previous_value = item
    return out_list


def count_trinucleotides(reference_index, bed, trinuc_counts, matching_chroms: list, save_trinuc_file: bool):
    # how many times do we observe each trinucleotide in the reference (and input bed region, if present)?
    trinuc_ref_count = {}

    # Count Trinucleotides in reference, based on bed or not
    # print(f'{PROG} - Counting trinucleotides in reference...')
    # Count the total number of bases spanned
    track_len = 0

    if bed:
        _LOG.info("since you're using a bed input, we have to count trinucs in bed region even if "
                  "you already have a trinuc count file for the reference...")
        with open_input(bed) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                record = line.strip().split('\t')
                track_len += int(record[2]) - int(record[1]) + 1
                if record[0] in reference_index.keys():
                    for i in range(int(record[1]), int(record[2]) - 1):
                        trinuc = reference_index[record[0]][i:i + 3].seq
                        if trinuc not in VALID_TRINUC:
                            continue
                        if trinuc not in trinuc_ref_count:
                            trinuc_ref_count[trinuc] = 0
                        trinuc_ref_count[trinuc] += 1
        if save_trinuc_file:
            _LOG.warning("Since we are using bed input, no trinuc file will be saved.")

    # Solution to attribute error (needs to be checked)
    # TODO remove ref_name from this dict
    elif not trinuc_counts.is_file():
        for ref_name in matching_chroms:
            sub_seq = reference_index[ref_name].seq
            for trinuc in VALID_TRINUC:
                if trinuc not in trinuc_ref_count:
                    trinuc_ref_count[trinuc] = 0
                trinuc_ref_count[trinuc] += sub_seq.count_overlap(trinuc)
        if save_trinuc_file:
            with open_output(trinuc_counts, 'w') as countfile:
                _LOG.info('Saving trinuc counts to file...')
                countfile.write(json.dumps(trinuc_ref_count))

    else:
        _LOG.info(f'Found counts file, {trinuc_counts}, using that.')
        with open_input(trinuc_counts) as counts:
            trinuc_ref_count = json.load(counts)
        if save_trinuc_file:
            _LOG.warning('Existing trinucelotide file will not be changed or overwritten.')

    return trinuc_ref_count, track_len


def find_caf(candidate_field: str) -> float:
    # if parsing a dbsnp vcf, and no CAF= is found in info tag, use this as default val for population freq
    vcf_default_pop_freq = 0.00001
    info_split = [a.split('=') for a in candidate_field.split(';')]
    for item in info_split:
        if item[0].upper() == 'CAF':
            if ',' in item[1]:
                return float(item[1].split(',')[1])
    return vcf_default_pop_freq

