"""
Helper functions for VCF files. The primary purpose is to retrieve the data from an input vcf file,
to ensure that those variants get inserted into the reads.
"""

import logging
from pathlib import Path

import numpy as np
import sys

from Bio import SeqIO
from Bio.File import _IndexedSeqFileDict

from .options import Options
from ...common import open_input, pick_ploids, get_genotype_string
from ...variants import ContigVariants, SingleNucleotideVariant, Insertion, Deletion, UnknownVariant

__all__ = [
    "parse_input_vcf"
]

_LOG = logging.getLogger(__name__)


def retrieve_genotype(my_record: list, is_cancer=False):
    # the first sample is index 9, but if there is a tumor, we'll add one to look in that column
    which_sample_index = 9 + (1 if is_cancer else 0)
    index = my_record[8].split(':').index('GT')
    # Apply the index corresponding to GT to the sample column to get the genotype, then split into ploids.
    ret = my_record[which_sample_index].split(':')[index].replace('/', '|').split('|')
    return np.array([int(x) for x in ret])


def variant_genotype(ploidy, full_genotype, which_alt):
    new_genotype = np.zeros(ploidy)
    for i in range(len(full_genotype)):
        if full_genotype[i] == which_alt:
            new_genotype[i] = 1
    return new_genotype


def parse_input_vcf(
        input_dict: dict[str, ContigVariants],
        vcf_path: Path,
        ploidy: int,
        homozygous_frequency: float,
        reference: _IndexedSeqFileDict,
        options: Options
) -> dict:
    """
    key to input_dict:
        - input_dict = {contig_1: ContigVariants_object, contit_2: ContigVariants_object, ...}

        - For each contig in the reference, there is a ContigVariants object that will be updated
          here, for use in the generate_variants module.
        - extra items from input vcf:
            - metadata: The actual info from the record
            - genotype: Genotype data for the first sample
            - genotype_tumor: Genotype data for the second (tumor) sample (optional)

    :param input_dict: Dicitonary of contig variants objects holding data for each contig
    :param vcf_path: Path to input vcf file
    :param ploidy: number of copies of each chromosome in the dataset.
    :param homozygous_frequency: Chance of an allele appearing on more than one ploid
    :param reference: The reference index object for this run.
    :param options: Options for this run
    :return: A dictionary with sample name data read from the vcf
    """

    _LOG.info(f"Parsing input vcf {vcf_path}")

    n_skipped = 0
    mismatched = 0
    records_found = 0
    # maximum number of columns we are interested in. Used for trimming unwanted samples.
    max_col = 7
    with open_input(vcf_path) as f:
        for line in f:
            # skip headers
            if line.startswith('##'):
                continue
            # Process the header row
            elif line.startswith('#CHROM'):
                columns = line.strip().strip('#').split('\t')

                # Anything after FORMAT is a sample column
                sample_columns = []
                has_format = False
                if 'FORMAT' in columns:
                    has_format = True
                    max_col += 1
                    sample_columns = columns[columns.index('FORMAT') + 1:]
                    if not sample_columns:
                        _LOG.error('Input vcf has FORMAT column but no sample columns.')
                        sys.exit(1)
                else:
                    _LOG.warning('Missing format column in vcf, using WP for genotype if present, '
                                 'otherwise genotype will be generated randomly')

                # Recode the sample columns to match the index of the dictionary we are generating
                # We only output 1 sample column for normal runs, 2 for tumor_normal. Those will be indices 7 and 8
                # in the output dictionary, se we hardcode those indices now for later retrieval

                if sample_columns:
                    sample_columns = {sample_columns[0]: 7}
                    max_col += 1

            # Process the records rows
            else:
                # list of variants from this line:
                line_variants = []

                record = line.strip().split('\t')
                # Decrement the position to get 0-based coordinates
                record[1] = int(record[1]) - 1
                # We'll index these by chromosome and position
                """
                For reference, the columns in a VCF, and their indices:
                    CHROM [0]
                    POS [1]
                    ID [2]
                    REF [3]
                    ALT [4]
                    QUAL [5]
                    FILTER [6]
                    INFO [7]
                    FORMAT [8, optional]
                    SAMPLE1 [9, optional]
                    SAMPLE2 [10, optional, cancer only]
                """
                # First, let's check if the chromosome for this record is even in the reference. Since input_dict is
                # constructed from the reference, the keys list is the same.
                in_ref = record[0] in input_dict.keys()
                if not in_ref:
                    _LOG.warning(f'Skipping variant because the chromosome is not in the reference:\n{line}')
                    continue

                reference_string = reference[record[0]][int(record[1]): int(record[1]) + len(record[3])].seq.upper()
                # We already accounted for shifting to 0-based coordinates, so this should work.
                if record[3] != str(reference_string):
                    mismatched += 1
                    _LOG.warning(f'Skipping variant because the ref field did not match the reference:'
                                 f'{record[0]}: {record[1]}, {record[3]} v '
                                 f'{reference_string}')
                    continue

                # Quality score could be missing, in this case, we treat it as ground truth and assign a default score
                default_qual = "42"
                if record[5] == ".":
                    record[5] = default_qual

                # We'll need the genotype when we generate reads, and output the records, if applicable
                genotype = None
                normal_sample_field = None
                tumor_sample_field = None

                if has_format:
                    if "GT" in record[8].split(':'):
                        # the format column will need no update.
                        format_column = record[8]
                        normal_sample_field = record[9]
                        # Retrieve the GT from the first sample in the record
                        genotype = retrieve_genotype(record)

                    elif "WP" in [x.split('=') for x in record[7].split(';')]:
                        """
                        "WP" is the legacy code NEAT used for genotype it added. It was found in the INFO field.
                        We're just going to make a sample column in this version of NEAT
                        The logic of the statement is split the info field on ';' which is used as a divider in that field.
                        Most but not all fields also have an '=', so split there too, then look for "WP"
                        """
                        format_column = f"GT:{record[8]}"
                        for record in record[7].split(';'):
                            if record.startswith('WP'):
                                genotype = record.split('=')[1].replace('/', '|').split('|')
                                genotype = np.array([int(x) for x in genotype])
                                normal_sample_field = f"{get_genotype_string(genotype)}:{record[9]}"

                    else:
                        format_column = 'GT:' + record[8]
                        alt_count = len(record[4].split(';'))
                        genotype = pick_ploids(ploidy, homozygous_frequency, alt_count, options.rng)
                        gt_field = get_genotype_string(genotype)
                        normal_sample_field = f'{gt_field}:{record[9]}'

                elif "WP" in [x.split('=') for x in record[7].split(';')]:
                    """
                    "WP" is the legacy code NEAT used for genotype it added. It was found in the INFO field.
                    We're just going to make a sample column in this version of NEAT
                    The logic of the statement is split the info field on ';' which is used as a divider in that field.
                    Most but not all fields also have an '=', so split there too, then look for "WP"
                    """
                    format_column = "GT"
                    info_split = record[7].split(';')
                    for record in info_split:
                        if record.startswith('WP'):
                            genotype = record.split('=')[1].replace('/', '|').split('|')
                            genotype = np.array([int(x) for x in genotype])
                            normal_sample_field = get_genotype_string(genotype)

                else:
                    # If there was no format column, there's no sample column, so we'll generate one
                    format_column = "GT"
                    alt_count = len(record[4].split(';'))
                    genotype = pick_ploids(ploidy, homozygous_frequency, alt_count, options.rng)
                    normal_sample_field = get_genotype_string(genotype)

                chrom = record[0]
                location = int(record[1])

                # first we'll attempt to classify this variant:
                ref = record[3]
                alts = record[4].split(',')
                data = {"REF": ref,
                        "ID": record[2],
                        "FILTER": record[6],
                        "INFO": record[7],
                        "FORMAT": format_column,
                        "QUAL": record[5],
                        "NEAT_sample": normal_sample_field,
                        "NEAT_cancer_sample": tumor_sample_field}

                count = 0
                for alt in alts:
                    count += 1
                    # This temp genotype teases out only the ploids with this particular variant
                    temp_genotype = variant_genotype(options.ploidy, genotype, count)
                    if len(ref) == len(alt) == 1:
                        # Type = SNV
                        temp_variant = SingleNucleotideVariant(
                            location, alt, temp_genotype, record[5], is_input=True, kwargs=data
                        )
                    elif len(ref) > len(alt) and ref.startswith(alt):
                        # type = deletion
                        temp_variant = Deletion(
                            location, len(alt), temp_genotype, record[5], is_input=True, kwargs=data
                        )
                    elif len(alt) > len(ref) and alt.startswith(ref):
                        # type = insertion
                        temp_variant = Insertion(
                            location, len(alt), alt, temp_genotype, record[5], is_input=True, kwargs=data)
                    else:
                        # Unknown variant type.
                        # We'll need the alternate, so we'll add it to data.
                        data["ALT"] = alt
                        temp_variant = UnknownVariant(location, temp_genotype, record[5], is_input=True, kwargs=data
                                                      )

                    rc = input_dict[chrom].add_variant(temp_variant)
                    if rc == 1:
                        _LOG.warning(f"Input variant skipped because a variant already existed at that location:"
                                     f"{chrom}: {location} ({temp_variant})")
                        n_skipped += 1
                    else:
                        records_found += 1

    _LOG.info(f'Found {records_found} variants in input VCF.')
    _LOG.info(f'Skipped {n_skipped} variants because of multiples at the same location')
    _LOG.info(f'Skipped {mismatched} variants because of a mismatch between Ref and reference.')

    return sample_columns
