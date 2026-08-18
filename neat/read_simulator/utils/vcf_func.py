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
    "parse_input_vcf",
    "collect_header_declarations"
]

_LOG = logging.getLogger(__name__)


def retrieve_genotype(my_record: list, is_cancer=False):
    """
    Reads the GT subfield of a record's sample column into a per-ploid array of allele indices.

    :param my_record: The VCF record, already split on tabs
    :param is_cancer: True to read the tumor sample column instead of the normal one
    :return: An array with one allele index per ploid, or None if the genotype is a no-call.
             VCF 4.2 allows any allele to be given as '.', which several callers and truth-set
             conventions emit for a record that is present but unassigned to a genotype. Such a
             record carries no ploid assignment for us to use, so the caller has to supply one.
    """
    # the first sample is index 9, but if there is a tumor, we'll add one to look in that column
    which_sample_index = 9 + (1 if is_cancer else 0)
    index = my_record[8].split(':').index('GT')
    # Apply the index corresponding to GT to the sample column to get the genotype, then split into ploids.
    ret = my_record[which_sample_index].split(':')[index].replace('/', '|').split('|')
    if any(x == '.' for x in ret):
        return None
    return np.array([int(x) for x in ret])


def replace_genotype_field(sample_field: str, format_column: str, genotype: np.ndarray):
    """
    Substitutes a genotype into the GT subfield of a sample column, leaving the other subfields
    alone. Used when the input GT could not be used and we generated one instead, so that the
    output VCF reports the genotype the reads were actually built from.

    :param sample_field: The sample column from the input record
    :param format_column: The FORMAT column, used to locate GT among the subfields
    :param genotype: The genotype to write in
    :return: The sample column with its GT subfield replaced
    """
    subfields = sample_field.split(':')
    gt_index = format_column.split(':').index('GT')
    # A VCF sample column may drop trailing subfields, so GT's slot might not exist yet.
    while len(subfields) <= gt_index:
        subfields.append('.')
    subfields[gt_index] = get_genotype_string(genotype)
    return ':'.join(subfields)


def prepend_genotype_field(sample_field: str, format_column: str, genotype: np.ndarray):
    """
    Builds the FORMAT and sample columns for a record NEAT had to generate a genotype for, with GT
    first as the spec requires.

    A FORMAT of '.' has no keys to preserve, so GT replaces it rather than being prepended to it:
    'GT:.' names a key called '.', which is not valid VCF and which downstream tools reject.

    :param sample_field: The sample column from the input record
    :param format_column: The FORMAT column from the input record
    :param genotype: The genotype NEAT generated for this record
    :return: The new (FORMAT, sample) column pair
    """
    gt_field = get_genotype_string(genotype)
    if format_column == '.':
        return 'GT', gt_field
    return f'GT:{format_column}', f'{gt_field}:{sample_field}'


def collect_header_declarations(vcf_path):
    """
    Works out the INFO, FILTER and FORMAT declarations the golden VCF needs in order to carry an
    input VCF's fields.

    parse_input_vcf copies each record's ID, FILTER, INFO and FORMAT fields onto the variant so
    they survive into the golden VCF. A record using a key the output header never declares is not
    a valid VCF, and tools reject it — bcftools, which NEAT runs to sort the golden VCF, fails on
    the first such record — so the declarations have to travel with the fields they describe.

    An input VCF is not required to be tidy about this: undeclared INFO and FORMAT keys are common
    (the H1N1 example shipped with NEAT uses AF and PP without declaring either). Anything the
    input's own header leaves out is declared permissively here, as Type=String with an unspecified
    Number, which is the same assumption bcftools makes when it reads such a file.

    :param vcf_path: Path to the input vcf file
    :return: Header lines to add to the golden VCF, without trailing newlines
    """
    declarations = []
    declared = set()
    used = {"##INFO": set(), "##FILTER": set(), "##FORMAT": set()}
    # PASS and the missing value need no declaration, and GT is declared by the output writer.
    never_declare = {"##INFO": {"."}, "##FILTER": {".", "PASS"}, "##FORMAT": {".", "GT"}}

    with open_input(vcf_path) as f:
        for line in f:
            if line.startswith('##'):
                if line.startswith(('##INFO=', '##FILTER=', '##FORMAT=')):
                    key = _declaration_id(line.strip())
                    # The output writer emits its own GT declaration, and a key declared twice is
                    # itself invalid, so keep only the first of each.
                    if key[1] not in never_declare.get(key[0], ()) and key not in declared:
                        declarations.append(line.strip())
                        declared.add(key)
                continue
            if line.startswith('#'):
                # The '#CHROM' line; records follow.
                continue

            record = line.strip().split('\t')
            if len(record) < 8:
                continue
            for entry in record[6].split(';'):
                used["##FILTER"].add(entry)
            if record[7] != '.':
                for entry in record[7].split(';'):
                    # A flag-style key has no '=' and so no value.
                    used["##INFO"].add(entry.split('=')[0])
            if len(record) > 8 and record[8] != '.':
                for entry in record[8].split(':'):
                    used["##FORMAT"].add(entry)

    # Fill in whatever the input used but never declared.
    source = Path(vcf_path).name
    for kind in ("##FILTER", "##INFO", "##FORMAT"):
        for identifier in sorted(used[kind]):
            if (not identifier or identifier in never_declare[kind]
                    or (kind, identifier) in declared):
                continue
            if kind == "##FILTER":
                declarations.append(f'{kind}=<ID={identifier},'
                                    f'Description="Carried over from {source}">')
            else:
                declarations.append(f'{kind}=<ID={identifier},Number=.,Type=String,'
                                    f'Description="Carried over from {source}">')

    return declarations


def _declaration_id(line: str):
    """
    Identifies a VCF header declaration by its kind and ID.

    :param line: A '##INFO=<ID=...>'-style header line
    :return: A (kind, id) tuple, or (kind, None) if the line has no ID field
    """
    kind = line.split('=', 1)[0]
    if 'ID=' in line:
        return kind, line.split('ID=', 1)[1].split(',')[0].rstrip('>')
    return kind, None


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

    :param input_dict: Dictionary of contig variants objects holding data for each contig
    :param vcf_path: Path to input vcf file
    :param ploidy: number of copies of each chromosome in the dataset.
    :param reference: The reference index object for this run.
    :param options: Options for this run
    :return: A dictionary with sample name data read from the vcf
    """

    _LOG.info(f"Parsing input vcf {vcf_path}")

    n_skipped = 0
    mismatched = 0
    records_found = 0
    n_nocall = 0
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
                # in the output dictionary, so we hard code those indices now for later retrieval

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
                        if genotype is None:
                            # A no-call GT is valid VCF but says nothing about which ploids carry
                            # the variant, so we treat it the same as a record whose FORMAT has no
                            # GT at all and generate one. The sample column is rewritten to match,
                            # so the reads and the output VCF agree on the genotype.
                            alt_count = len(record[4].split(','))
                            genotype = pick_ploids(ploidy, 0.001, alt_count, options.rng)
                            normal_sample_field = replace_genotype_field(record[9], format_column, genotype)
                            n_nocall += 1
                            _LOG.warning(f'No-call genotype in input VCF, assigned one at random '
                                         f'({get_genotype_string(genotype)}): '
                                         f'{record[0]}: {record[1] + 1}')

                    elif "WP" in [x.split('=')[0] for x in record[7].split(';') if '=' in x]:
                        """
                        "WP" is the legacy code NEAT used for genotype it added. It was found in the INFO field.
                        We're just going to make a sample column in this version of NEAT
                        The logic of the statement is split the info field on ';' which is used as a divider in that field.
                        Most but not all fields also have an '=', so split there too, then look for "WP"
                        """
                        for info_item in record[7].split(';'):
                            if info_item.startswith('WP') and '=' in info_item:
                                genotype = info_item.split('=')[1].replace('/', '|').split('|')
                                genotype = np.array([int(x) for x in genotype])
                                format_column, normal_sample_field = prepend_genotype_field(
                                    record[9], record[8], genotype
                                )
                            elif info_item.startswith('WP'):
                                _LOG.error(f'Malformed WP field in INFO (missing value): {record[7]}')
                                sys.exit(1)

                    else:
                        alt_count = len(record[4].split(','))
                        genotype = pick_ploids(ploidy, 0.001, alt_count, options.rng)
                        format_column, normal_sample_field = prepend_genotype_field(
                            record[9], record[8], genotype
                        )

                elif "WP" in [x.split('=')[0] for x in record[7].split(';') if '=' in x]:
                    """
                    "WP" is the legacy code NEAT used for genotype it added. It was found in the INFO field.
                    We're just going to make a sample column in this version of NEAT
                    The logic of the statement is split the info field on ';' which is used as a divider in that field.
                    Most but not all fields also have an '=', so split there too, then look for "WP"
                    """
                    format_column = "GT"
                    for info_item in record[7].split(';'):
                        if info_item.startswith('WP') and '=' in info_item:
                            genotype = info_item.split('=')[1].replace('/', '|').split('|')
                            genotype = np.array([int(x) for x in genotype])
                            normal_sample_field = get_genotype_string(genotype)
                        elif info_item.startswith('WP'):
                            _LOG.error(f'Malformed WP field in INFO (missing value): {record[7]}')
                            sys.exit(1)

                else:
                    # If there was no format column, there's no sample column, so we'll generate one
                    format_column = "GT"
                    alt_count = len(record[4].split(','))
                    genotype = pick_ploids(ploidy, 0.001, alt_count, options.rng)
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
                    if ref == alt:
                        _LOG.warning(
                            f"Skipping variant at {chrom}:{location + 1} — REF == ALT ({ref!r}). "
                            f"This is not a valid variant."
                        )
                        n_skipped += 1
                        continue
                    # This temp genotype teases out only the ploids with this particular variant
                    temp_genotype = variant_genotype(options.ploidy, genotype, count)
                    if len(ref) == len(alt) == 1:
                        # Type = SNV
                        temp_variant = SingleNucleotideVariant(
                            location, alt, temp_genotype, record[5], is_input=True, **data
                        )
                    elif len(ref) > len(alt) and ref.startswith(alt):
                        # type = deletion
                        # length is the full reference span (number of bases the deletion
                        # covers, including the shared anchor base), matching how deletions
                        # are reconstructed on output and generated by the mutation model.
                        temp_variant = Deletion(
                            location, len(ref), temp_genotype, record[5], is_input=True, **data
                        )
                    elif len(alt) > len(ref) and alt.startswith(ref):
                        # type = insertion
                        temp_variant = Insertion(
                            location, len(alt), alt, temp_genotype, record[5], is_input=True, **data)
                    else:
                        # Unknown variant type.
                        # We'll need the alternate, so we'll add it to data.
                        data["ALT"] = alt
                        temp_variant = UnknownVariant(location, temp_genotype, record[5], is_input=True, **data
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
    if n_nocall:
        _LOG.info(f'Assigned a random genotype to {n_nocall} variants with a no-call genotype.')

    return sample_columns
