import logging
import pathlib
import time

from Bio import SeqRecord
from Bio.Seq import Seq
from numpy.random import Generator

from ...models import SequencingErrorModel, GcModel, FragmentLengthModel
from .options import Options
from ...common import open_input, ALLOWED_NUCL
from .read import Read

__all__ = [
    'generate_reads'
]

_LOG = logging.getLogger(__name__)


def is_low_coverage(input_list, bias_vector, test_val):
    # We'll give ourselves a padding of 10% to prevent oversampling
    allowance = 0.9
    for i in range(len(input_list)):
        if input_list[i] < test_val * bias_vector[i] * allowance:
            return True
    return False


def is_too_many_n(segment):
    n = segment.count('N')
    return n/len(segment) >= 0.2


def create_windows(sequence: SeqRecord, size: int, overlap: int):
    """
    Create a list of windows
    :param sequence: Sequence to split
    :param size: size of windows
    :param overlap: size of overlap between windows
    :return: list of windows
    """
    windows = []
    for i in range(0, len(sequence), size):
        if i < overlap and i + size + overlap < len(sequence):
            windows.append((i, i+size+overlap))
        if i + size + overlap < len(sequence):
            windows.append((i-overlap, i+size+overlap))
        else:
            windows.append((i, len(sequence)))

    return windows


def process_read(myread: Read,
                 mylist: list[Read, ...],
                 rng: Generator,
                 err_model: SequencingErrorModel,
                 coverage: list):
    """
    Inserts a read into a list, if it is not already there

    :param myread: The read to insert
    :param mylist: The list in which to insert myread
    :param rng: The random number generator for the run
    :param err_model: The error model for the run
    :param coverage: list of read counts per base
    """

    if not any(x for x in mylist if x == myread):
        # Generate quality scores
        myread.quality_array = err_model.get_quality_scores(myread.length, rng)
        # Get errors
        myread.errors = err_model.get_sequencing_errors(myread.quality_array, rng)
        mylist.append(myread)
        # This counts the new coverage of the bases, only if we added the read.
        for n in range(myread.position, myread.end):
            coverage[n] += 1


def generate_reads(reference: SeqRecord,
                   error_model: SequencingErrorModel,
                   gc_bias: GcModel,
                   fraglen_model: FragmentLengthModel,
                   input_vcf: str,
                   temporary_directory: str | pathlib.Path,
                   targeted_regions: list,
                   discarded_regions: list,
                   mutation_rates: list,
                   options: Options,
                   chrom: str):
    # apply errors and mutations and generate reads
    chrom_fastq_r1 = temporary_directory / f'{chrom}_tmp_r1.fq'
    chrom_fastq_r2 = None
    if options.paired_ended:
        chrom_fastq_r2 = temporary_directory / f'{chrom}_tmp_r2.fq'

    read_length = int(options.read_len)
    average_fragment_size = int(options.read_len)
    if options.paired_ended:
        average_fragment_size = int(fraglen_model.fragment_mean)

    vars_to_insert = {}
    frag_mean = None
    frag_std = None
    if options.paired_ended:
        frag_mean = fraglen_model.fragment_mean
        frag_std = fraglen_model.fragment_st_dev
    with open_input(input_vcf) as input_variants:
        for line in input_variants:
            if line.startswith("@") or line.startswith("#"):
                continue
            line_split = line.strip().split('\t')
            # Since these vars are on the same chromosome, we can index by position
            # We don't need the ID for this so let's skip it
            # The remaining fields will be:
            #   - vars_to_insert[position][0] = REF
            #   - vars_to_insert[position][1] = ALT
            #   - vars_to_insert[position][2] = QUAL
            #   - vars_to_insert[position][3] = FILTER
            #   - vars_to_insert[position][4] = INFO
            #   - vars_to_insert[position][5] = FORMAT
            #   - vars_to_insert[position][6] = SAMPLE_1
            #   - vars_to_insert[position][7] = SAMPLE_2 (optional)
            vars_to_insert[line_split[1]] = line_split[3:]

    _LOG.info(f'Sampling reads...')
    start = time.time()

    base_name = f'{pathlib.Path(options.output).name}-{chrom}'
    total_bp_spanned = len(reference)
    coverage_vector = [0.0] * len(reference)
    coverage_bias = gc_bias.create_coverage_bias_vector(reference.seq)
    window_size = min(total_bp_spanned, average_fragment_size * 100)
    overlap = average_fragment_size
    windows = create_windows(reference, window_size, overlap)
    min_length = read_length
    max_length = window_size
    if options.paired_ended:
        max_length = int(min(float(max_length), frag_mean + (6.0 * frag_std)))

    left_reads = []
    right_reads = []
    previous_percent = 0
    i = 0  # Counter number of loops
    for window in windows:
        current_window_size = window[1] - window[0]
        random_positions = list(range(window[0], window[1]))
        options.rng.shuffle(random_positions)
        for position in random_positions:

            # Simple progress tracker
            current_percent = (i * 100)//current_window_size
            if current_percent > previous_percent:
                print(f'{current_percent}%', end='\r')
                previous_percent = current_percent

            # A sample segment to check for N concentration
            segment = reference[i: i + average_fragment_size]

            # Make sure the segments will come out to at around 80% valid bases
            if is_too_many_n(segment.seq):
                i += 1
                continue

            # check if we are in a targeted region
            # check if we are in a discard region

            dist = [read_length] * (int(options.coverage) * 3)
            if options.paired_ended:
                # Generate a pool of possible fragment lengths based on the mean and standard deviation
                dist = options.rng.normal(loc=frag_mean,
                                          scale=frag_std,
                                          size=int(options.coverage) * 3)
                # filter down to lengths between the min and max, then round to the nearest int
                dist = [round(m) for m in dist if min_length <= m <= max_length]

            for fragment_length in dist:
                end_point = position + fragment_length
                low_cov = is_low_coverage(coverage_vector[position: len(reference)], coverage_bias,
                                          options.coverage)
                if not low_cov:
                    break
                # need to make sure we have enough bases to generate reads toward the end
                if end_point >= len(reference):
                    if position + read_length > len(reference):
                        # If we can't squeeze in a minimum read, and we still need to cover the end of the genome,
                        # we'll shift position back to make sure
                        # we don't get caught in a loop trying to finish up the coverage vector.
                        position = len(reference) - read_length

                    # Since we have enough room, we'll just grab to the end of the reference
                    fragment_length = len(reference) - position

                if 'N' in segment:
                    modified_segment = ""
                    for base in segment:
                        if base not in ALLOWED_NUCL:
                            modified_segment += options.rng.choice(ALLOWED_NUCL)
                        else:
                            modified_segment += base

                    segment = Seq(modified_segment)

                read_name = f'{base_name}-{str(i)}'
                read1 = Read(read_name, reference, position, position + read_length)
                process_read(read1, left_reads, options.rng, error_model, coverage_vector)

                if options.paired_ended:
                    read2 = Read(read_name, reference, position + fragment_length - read_length,
                                 position + fragment_length, is_reverse=True)

                    process_read(read2, right_reads, options.rng, error_model, coverage_vector)

            i += 1

    print("100%")

    print([x for x in left_reads + right_reads if x.contains(1)])

    _LOG.info(f"Finished sampling reads in {time.time() - start} seconds")
    return chrom_fastq_r1, chrom_fastq_r2
