"""
Utilities used by the generate mutation model function
"""

import json
import sys
import pickle
import gzip
import time  # Import time for potential future use here

import numpy as np

from pathlib import Path
import logging

from ..common import (
    open_input,
    open_output,
    ALL_TRINUCS,
    ALLOWED_NUCL,
    ALL_CONTEXTS,
    DINUC_IND,
    NUC_IND,
    TRINUC_IND,
)
from .constants import VCF_DEFAULT_POP_FREQ

__all__ = [
    "extract_header",
    "read_and_filter_variants",
    "convert_trinuc_transition_matrix",
    "cluster_list",
    "count_trinucleotides",
    "find_caf",
    "convert_snp_transition_matrix",
    "convert_trinuc_mutation_dict",
    "check_homozygous",
]

_LOG = logging.getLogger(__name__)


def extract_header(vcf_file):
    """
    This function extracts the header lines into a list from an input VCF file
    :param Path vcf_file: A valid vcf file path
    :return list: A list of header lines
    """
    ret = []
    with open_input(vcf_file) as f:
        for line in f:
            if line.startswith("##"):
                ret.append(line.strip())
            elif line.startswith("#CHROM"):
                temp = line.strip().strip("#").split("\t")
                ret.append(temp)
                break

    if not ret:
        _LOG.error("No header found: invalid VCF file.")
        sys.exit(1)

    return ret


def read_and_filter_variants(vcf_file, reference_index, ignore):
    """
    Finds all the matching chromosomes between ref and vcf(mutation file)

    :param Path vcf_file: Full path to input vcf
    :param dict reference_index: SeqIO dictionary of the index
    :param list ignore: list of contigs to ignore
    :return list, list: list of matching variants and a list of matching chromosomes
    """
    start_time = time.time()  # Start timing this function
    variant_chroms = []
    matching_chroms = []

    final_data = []

    with open_input(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            else:
                columns = line.strip().split("\t")

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
                if chrom in ignore:
                    continue
                if chrom not in variant_chroms:
                    variant_chroms.append(chrom)

                # if a chromosome is present in the reference and the vcf file -> add it to matching_chrom list.
                if chrom not in matching_chroms:
                    if chrom in reference_index:
                        matching_chroms.append(chrom)

                # If CHROM is present in matching_chroms, then the variant is a candidate
                if chrom in matching_chroms:
                    # multi-allelic, we'll just take the first one
                    if "," in columns[4]:
                        variant_alt = columns[4].split(",")[0]
                    else:
                        variant_alt = columns[4]
                    variant_ref = columns[3]

                    # If this variant has invalid characters, we'll skip it as too complex
                    if any(
                        [char for char in variant_ref if char not in ALLOWED_NUCL]
                    ) or any(
                        [char for char in variant_alt if char not in ALLOWED_NUCL]
                    ):
                        continue

                    # initialize flags
                    add_variant = False
                    type_flag = "SNP"

                    # SNP
                    if len(variant_ref) == 1 and len(variant_alt) == 1:
                        add_variant = True
                    # deletion
                    elif len(variant_ref) > 1 and len(variant_alt) == 1:
                        add_variant = True
                        type_flag = "DEL"
                    # insertion
                    elif len(variant_ref) == 1 and len(variant_alt) > 1:
                        add_variant = True
                        type_flag = "INS"
                    # Other variant types are skipped

                    """
                    Structure:
                    final_data[i][0]: Chromosome (CHROM) column
                    final_data[i][1]: Reference-relative position (0-based)
                    final_data[i][2]: Variant REF column
                    final_data[i][3]: Variant ALT column
                    final_data[i][4]: variant type flag ("SNP", "INS", "DEL")
                    final_data[i][5]: Variant INFO field
                    """
                    if add_variant:
                        # Ensure all necessary columns exist before accessing index 7 (INFO)
                        if len(columns) > 7:
                            info_field = columns[7]
                        else:
                            info_field = (
                                "."  # Or handle missing INFO field appropriately
                            )
                        final_data.append(
                            [
                                chrom,
                                str(int(columns[1]) - 1),
                                variant_ref,
                                variant_alt,
                                type_flag,
                                info_field,
                            ]
                        )

        print(
            f"Variant chroms not in ref: {list(set(reference_index) - set(variant_chroms))}"
        )
        print(f"Matching chroms: {matching_chroms}")

    end_time = time.time()  # End timing
    _LOG.debug(
        f"read_and_filter_variants execution time: {end_time - start_time:.4f} seconds"
    )  # Log duration
    return final_data, matching_chroms


def cluster_list(list_to_cluster, delta):
    """
    Clusters a sorted list
    :param list list_to_cluster: a sorted list
    :param float delta: the value to compare list items to
    :return list: a clustered list of values
    """
    start_time = time.time()  # Start timing
    if not list_to_cluster:  # Handle empty list case
        return []
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
    end_time = time.time()  # End timing
    _LOG.debug(
        f"cluster_list execution time: {end_time - start_time:.4f} seconds"
    )  # Log duration
    return out_list


def count_trinucleotides(
    reference_index, bed, trinuc_counts, matching_chroms, save_trinuc_file
):
    """
    Counts the frequency of the various trinucleotide combinations in the dataset

    :param dict reference_index: The index of the fasta, from SeqIO
    :param Path bed: Full path to bed file, if using
    :param Path trinuc_counts: Full path to existing trinculeotide counts file, if input, or the output path
    :param list matching_chroms: List of matching chromosomes
    :param save_trinuc_file: Boolean determines whether to save trinucleotide counts file.
    :return dict, int: A dictionary of trinculeotides and counts, the number of bases spanned
    """
    start_time = time.time()  # Start timing
    # Data format: TRINUC: COUNT (e.g., "AAA": 12), where COUNT is the number of times trinuc was observed
    # in the reference sequences.
    trinuc_ref_count = {}

    # Count the total number of bases spanned
    track_len = 0

    trinuc_counts_exists = False
    if trinuc_counts:
        if trinuc_counts.is_file():
            _LOG.info(
                "Trinucleotide counts file exists, skipping save to avoid overwriting data."
            )
            trinuc_counts_exists = True
            save_trinuc_file = False

    if bed:
        _LOG.info("Counting trinucleotide combinations in bed regions")
        if trinuc_counts_exists:
            _LOG.warning(
                "Ignoring trinucleotide counts file, to restrict counts to bed regions"
            )
        if save_trinuc_file:
            _LOG.warning("Using bed input, no trinuc counts file will be saved.")
        current_line = 0
        with open_input(bed) as f:
            print("line ", current_line, " of ", len(matching_chroms))
            for line in f:
                if line.startswith("#"):
                    continue
                record = line.strip().split("\t")
                track_len += int(record[2]) - int(record[1]) + 1
                if record[0] in reference_index:
                    for i in range(int(record[1]), int(record[2]) - 1):
                        trinuc = reference_index[record[0]][i : i + 3].seq
                        trinuc_ref_count = count_trinuc(trinuc, trinuc_ref_count)
                current_line += 1

    # Solution to attribute error (needs to be checked)
    elif not trinuc_counts_exists:
        _LOG.info("Counting trinucleotide combinations in reference.")
        for ref_name in matching_chroms:
            contig_seq = reference_index[ref_name].seq
            for i in range(len(contig_seq)):
                trinuc = contig_seq[i : i + 3]
                trinuc_ref_count = count_trinuc(trinuc, trinuc_ref_count)

        if save_trinuc_file:
            with open_output(trinuc_counts, "w") as countfile:
                _LOG.info("Saving trinuc counts to file...")
                # Convert all values to writable formats
                pickle.dump(trinuc_ref_count, countfile)

    else:
        _LOG.info(f"Loading file: {trinuc_counts}.")
        with gzip.open(trinuc_counts, "rb") as counts:
            trinuc_ref_count = pickle.load(counts)
        if save_trinuc_file:
            _LOG.warning("Existing trinucelotide file will not be changed.")

    end_time = time.time()  # End timing
    _LOG.debug(
        f"count_trinucleotides execution time: {end_time - start_time:.4f} seconds"
    )  # Log duration
    return trinuc_ref_count, track_len


def count_trinuc(trinuc_seq, ref_count_dict):
    """
    Helper function that adds trinuc to dict, if it's valid, and counts it.
    :param Sequence trinuc_seq: Three letter sequence
    :param dict ref_count_dict:
    :return dict: the updated ref_count_dict
    """
    if trinuc_seq not in ALL_TRINUCS:
        return ref_count_dict
    if trinuc_seq not in ref_count_dict:
        ref_count_dict[trinuc_seq] = 0
    ref_count_dict[trinuc_seq] += 1
    return ref_count_dict


def find_caf(candidate_field):
    """
    Finds the CAF field, if present in INFO, otherwise returns the default
    :param str candidate_field: The vcf info to parse
    :return float: The given allele frequency
    """
    info_split = [a.split("=") for a in candidate_field.split(";")]
    for item in info_split:
        if item[0].upper() == "CAF":
            if "," in item[1]:
                return float(item[1].split(",")[1])


def check_homozygous(info_field):
    """
    Checks info field for genotype to try to detect if this sample is homozygous

    :param str info_field: The info field of an input variant
    :return int: 1 if homozygous, 0 otherwise
    """
    # Minimal overhead, timing might not be very informative unless called extremely often
    # start_time = time.time()
    fields = info_field.strip().split(":")
    genotype = None
    for item in fields:
        if len(item.split("/")) == 1:
            if len(item.split("|")) == 1:
                continue
            else:
                genotype = [int(x) for x in item.split("|")]
                break
        else:
            genotype = [int(x) for x in item.split("/")]
            break

    if genotype:
        gt = genotype[0]
        for allele in genotype[1:]:
            if allele != gt:
                return 0
        return 1

    return 0


def convert_trinuc_transition_matrix(trans_probs):
    """
    Convert the calculated transitions into a transition matrix for the mutation model.
    Each matrix is the probability of one

    :param dict trans_probs: A dictionary of trinuc1 -> trinuc2 probabilities
    :return np.ndarray: The transition matrix for each trinucleotide, indexed by our common scheme.
    """
    start_time = time.time()  # Start timing
    ret_matrix = np.zeros((16, 4, 4))
    for key, value in trans_probs.items():
        # Construct the context to fill out, and fetch the index
        context = DINUC_IND[key[0][0] + "_" + key[0][2]]
        # Get indexes of the ref and alt
        mutation_ref = NUC_IND[key[0][1]]
        mutation_alt = NUC_IND[key[1][1]]
        if ret_matrix[context][mutation_ref][mutation_alt] == 0.0:
            ret_matrix[context][mutation_ref][mutation_alt] = value
        else:
            _LOG.error("Repeat Trinuc detected.")
            _LOG.debug(
                f"Error on {ALL_CONTEXTS[context]}: "
                f"{ALLOWED_NUCL[mutation_ref]} -> {ALLOWED_NUCL[mutation_alt]}"
            )
            sys.exit(1)

    end_time = time.time()  # End timing
    _LOG.debug(
        f"convert_trinuc_transition_matrix execution time: {end_time - start_time:.4f} seconds"
    )  # Log duration
    return ret_matrix


def convert_snp_transition_matrix(snp_matrix):
    """
    Converts the dict form, used for display, to the final indexed form optimized for speed.

    :param dict snp_matrix: A dictionary of snp transition probabilities
    :return np.ndarray: Indexed by our common indexing scheme
    """
    start_time = time.time()  # Start timing
    ret_matrix = np.zeros((4, 4))
    for key, value in snp_matrix.items():
        ret_matrix[NUC_IND[key[0]]][NUC_IND[key[1]]] = value

    end_time = time.time()  # End timing
    _LOG.debug(
        f"convert_snp_transition_matrix execution time: {end_time - start_time:.4f} seconds"
    )  # Log duration
    return ret_matrix


def convert_trinuc_mutation_dict(trinuc_mut_dict):
    """
    Converts the dict form, used for display, to the final indexed form optimized for speed.

    :param dict trinuc_mut_dict: A dictionary of trinuc mutation probabilities
    :return np.ndarray: Indexed by our common indexing scheme
    """
    start_time = time.time()  # Start timing
    ret_array = np.zeros(64)
    for trinuc in ALL_TRINUCS:
        ret_array[TRINUC_IND[trinuc]] = trinuc_mut_dict[trinuc]

    end_time = time.time()  # End timing
    _LOG.debug(
        f"convert_trinuc_mutation_dict execution time: {end_time - start_time:.4f} seconds"
    )  # Log duration
    return ret_array
