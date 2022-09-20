import logging
import time
import numpy as np

from math import ceil, floor
from pathlib import Path
from Bio import SeqRecord
from Bio.Seq import Seq
from numpy.random import Generator
from bisect import bisect_left, bisect_right
from scipy.stats import ks_2samp
from matplotlib import pyplot as plt

from ...models import SequencingErrorModel, GcModel, FragmentLengthModel
from .options import Options
from ...common import open_input, open_output, ALLOWED_NUCL
from ...variants import ContigVariants
from .read import Read

__all__ = [
    'generate_reads',
    'cover_dataset'
]

_LOG = logging.getLogger(__name__)


def is_too_many_n(segment):
    """
    Checks that the segment is valid and there aren't too many invalid characters in it.

    :param segment: the sequence to check
    :return: True if there are too many invalid characters, or no reads; False otherwise
    """
    if not segment:
        return True
    n = segment.count('N')
    return n / len(segment) >= 0.2


def create_windows(sequence: SeqRecord, size: int, overlap: int):
    """
    Create a list of windows. We might need this for parallelization, so I'm leaving the code in place

    :param sequence: Sequence to split
    :param size: size of windows
    :param overlap: size of overlap between windows
    :return: list of windows
    """
    windows = []
    for i in range(0, len(sequence), size):
        if i < overlap and i + size + overlap < len(sequence):
            windows.append((i, i + size + overlap))
        if i + size + overlap < len(sequence):
            windows.append((i - overlap, i + size + overlap))
        else:
            windows.append((i, len(sequence)))

    return windows


def cover_dataset(
        span_length: int,
        target_vector: np.ndarray,
        window_size: int,
        options: Options,
        read_length_model: FragmentLengthModel,
        fragment_model: FragmentLengthModel | None,
) -> list:
    """
    Covers a dataset to the desired depth in the paired ended case. This is the main algorithm for creating the reads
    to the proper coverage depth. It uses an abstract representation of the reads, by end points.

    :param span_length: The total length the cover needs to span
    :param target_vector: The vector of coverage targets by position
    :param window_size: The size of the window used for gc-bias
    :param options: The options for the run
    :param read_length_model: The model used to generate random read lengths
    :param fragment_model: The fragment model used for to generate random fragment lengths
    :return:
    """

    # Variables to track reads and current coverage
    temp_reads = []
    coverage = [0] * span_length
    read_cover = {x: [] for x in range(span_length)}

    read_pool = read_length_model.generate_fragments(span_length, options.read_len, options.coverage)

    if fragment_model:
        # Paired ended case
        fragment_pool = fragment_model.generate_fragments(span_length, options.read_len, options.coverage)
    else:
        # single ended case
        fragment_pool = read_pool
    padding = max(fragment_pool) // 2

    i = 0

    _LOG.debug("Covering dataset")

    target = options.rng.choice([ceil(max(target_vector)), floor(max(target_vector))])
    while i <= target:
        _LOG.debug(f"Covering layer {i+1} of {target}")
        # In order to make the dataset look random, but still achieve proper coverage, we're going to have some loops
        # Space the reads randomly, and some loops create very regular reads. Outside this, we'll shuffle the reads
        # to give them a random appearance.
        if i % 10 == 0:
            random_mode = 1
        else:
            random_mode = 0

        # For either random or nonrandom mode, we wand to start at a random place just before the span, out to padding
        j = -options.rng.choice(range(padding))
        previous_percent = 0
        k = 0  # independent tracking
        while -padding <= j < span_length:
            k += 1
            # Simple progress tracker
            current_percent = (j * 100) // span_length
            if current_percent > previous_percent:
                print(f'{current_percent}%', end='\r')
                previous_percent = current_percent
            # The structure for these reads will be (left_start, left_end, right_start, right_end) where start and end
            # are ints with start > end. Reads can overlap, so right_start < left_end is possible, but the reads cannot
            # extend past each other, so right_start < left_start and left_end > right_end are not possible.
            if options.paired_ended:
                left_read_length, right_read_length = options.rng.choice(read_pool, size=2, replace=False)
                fragment_length = options.rng.choice(fragment_pool)
            else:
                left_read_length = options.rng.choice(read_pool)
                right_read_length = 0
                fragment_length = left_read_length

            if fragment_length < max(left_read_length, right_read_length):
                # Outlier, skip
                j += fragment_length//2 + (random_mode * options.rng.poisson(1))
                continue

            # Booleans telling us how to deal with the various exceptions below.
            skip_left = False
            skip_right = False

            # This trims the left read to the region of interest
            left_start = max(j, 0)
            left_end = min(span_length, j + left_read_length)

            # Make right read, if paired ended, else set it equal to the default
            if options.paired_ended:
                # The uncorrected left_start is j
                # The uncorrected left_end is j + left_read_length
                # The candidate end points of the paired right read, uncorrected
                uncorrected_right_end = j + fragment_length
                uncorrected_right_start = uncorrected_right_end - right_read_length

                # this trims the raw right read, so that right_start >= left_start and left_end <= right_end
                # this keeps us from getting unrealistic paired ended reads
                right_start = max(uncorrected_right_start, left_start)
                right_end = min(uncorrected_right_end, left_end)

                # Both reads are off the map, so we throw this out:
                if right_end <= 0:
                    j += fragment_length//2 + (random_mode * options.rng.poisson(1))
                    continue

            else:
                skip_right = True
                right_start = 0
                right_end = 0

            # Cases that might make us want to skip these reads:
            # Case 1: Left read end completely out of the area of consideration, unusable
            if left_end <= 0:
                # right read is usable if it's not off the map and has at least half a read
                # length of usable bases
                if right_start < span_length and \
                        (right_end - right_start > right_read_length/2):
                    # skip the left read, keep the right read, indicating an unpaired right read
                    skip_left = True
                # else, skip both
                else:
                    skip_left = skip_right = True

            # Case 2: Truncated left read, unusable
            if left_end - left_start < left_read_length/2:
                # right read is usable if it's not off the map and has at least half a read
                # length of usable bases
                if right_start < span_length and \
                        (right_end - right_start >= right_read_length/2):
                    # skip the left read, keep the right read, indicating an unpaired right read
                    skip_left = True
                # else, skip both
                else:
                    skip_left = skip_right = True

            # Case 3: Right read starts after the end of the span, unusable
            # Note: In single ended mode the next two cases are irrelevant, but it shouldn't cost much time to run them
            if right_start >= span_length:
                # left read is usable if it's not off the map and has at least half a
                # read length of usable bases
                if left_start < span_length and \
                        (left_end - left_start > left_read_length/2):
                    # Set the right read to (0, 0), indicating an unpaired left read
                    skip_right = True
                # else, skip both
                else:
                    skip_left = skip_right = True

            # Case 4: Truncated right read, unusable
            if right_end - right_start < right_read_length/2:
                # left read is usable if it's not off the map and has at least half a
                # read length of usable bases
                if left_start < span_length and \
                        (left_end - left_start > left_read_length/2):
                    # Set the right read to (0, 0), indicating an unpaired left read
                    skip_right = True
                # else, skip both
                else:
                    skip_left = skip_right = True

            if skip_left and skip_right:
                # If neither read is usable, we'll just move to the next
                j += fragment_length//2 + (random_mode * options.rng.poisson(1))
                continue

            if not skip_left:
                temp_left_read = parse_coverage(
                    left_start,
                    left_end,
                    left_read_length,
                    coverage,
                    target_vector
                )
            else:
                temp_left_read = (0, 0)

            if not skip_right:
                temp_right_read = parse_coverage(
                    right_start,
                    right_end,
                    right_read_length,
                    coverage,
                    target_vector,
                )
            else:
                temp_right_read = (0, 0)

            temp_reads.append(temp_left_read + temp_right_read)
            j += left_read_length + (random_mode * options.rng.poisson(1))

            if k > span_length:
                delta = span_length - k
                if delta % 100 == 0:
                    _LOG.debug(f"K is getting large: {k} v {span_length}")

        i += run_one_more(coverage, i, target, options.coverage)

    _LOG.debug("Culling reads to final set.")
    final_reads = culling_procedure(
        temp_reads,
        span_length,
        coverage,
        read_cover,
        target_vector,
        window_size,
        options.rng
    )

    return final_reads


def calculate_coverage(coverage_section: list, target_vector: np.ndarray) -> (list, int):
    """
    Calculates the coverage of a given read, but prevents us from going too far over the target coverage

    :param coverage_section: The coverage section to update
    :param target_vector: The vector with the coverage targets,
        should be at least the same length as the coverage section
    :return: The updated coverage section as a list plus the k we reached in the calculation
    """
    for k in range(len(coverage_section)):
        # We'll give ourselves a little buffer
        if coverage_section[k] <= target_vector[k] + 2:
            coverage_section[k] += 1
        else:
            return coverage_section, k

    return coverage_section, len(coverage_section)


def parse_coverage(
        read_start: int,
        read_end: int,
        read_length: int,
        coverage_vector: list,
        target_vector: np.ndarray,
) -> tuple[int, int]:
    """
    We need to be able to calculate the coverage and determine a course of action based on the results multiple
    times. This function handles a call to generate the coverage, then interprets the coverage results.

    :param read_start: Start point of the read
    :param read_end: End point of the read
    :param read_length: The length of the read (end - start)
    :param coverage_vector: The vector showing the current coverage
    :param target_vector: The vector showing the target coverage values at each location
    :return: A read of the maximum length allowed by coverage or a read of (0,0), indicating it isn't to be used
    """

    temp_coverage = coverage_vector[read_start: read_end]
    temp_target = target_vector[read_start: read_end]

    temp_coverage, new_length = calculate_coverage(temp_coverage, temp_target)

    if new_length > read_length / 2:
        # If we have enough bases to make at least half a read, we'll add it as a full read.
        # In most cases left_length = left_read_length
        coverage_vector[read_start: read_end] = temp_coverage
        return read_start, read_start + new_length
    else:
        # Else, we throw this one out and try to use the paired read
        return 0, 0


def culling_procedure(
        candidate_reads: list,
        span_length: int,
        coverage: list,
        cover_dict: dict,
        target_vector: np.ndarray,
        window_size: int,
        rng: Generator
) -> list:
    """
    The culling procedure for reads. It is based on the median values (median as opposed to mean to
    allow for datasets with lots of zeros) of the coverage and target vectors. Basically, we want to
    throw away reads until we get to the right median value. We start by trying to throw out reads
    where the coverage is too high. If that doesn't work, we throw them out randomly one by one until the
    median is in acceptable parameters.

    :param candidate_reads: The list of reads in the following format: (left_start, left_end, right_start, right_end)
    :param span_length: The total length spanned by this cover
    :param coverage: The calculated coverage vector for the left_reads + right_reads dataset
    :param cover_dict: a dictionary of the cover at each location
    :param target_vector: The target coverage values, by location
    :param window_size: The size of the window used to compute GC bias
    :param rng: The random number generator for the run
    :return: Two lists, left and right reads, culled to proper coverage depth.
    """

    reads_culled = 0
    target_med = np.median(target_vector)
    target_std = np.std(target_vector)

    #  First check if any culling is necessary
    if np.median(coverage) <= target_med + target_std:
        _LOG.debug(f"Culled {reads_culled} reads")
        return candidate_reads

    for i in range(0, span_length, window_size):
        window_target = target_vector[i: i+window_size]
        window_target_median = np.median(window_target)

        window_coverage = coverage[i: i+window_size]
        window_coverage_median = np.median(window_coverage)

        # Trying to identify areas of high coverage we can trim from
        while window_coverage_median > window_target_median + (2 * target_std):
            reads_culled += 1
            # toss a read from the coverage window
            position_to_cull = rng.choice(range(i, i+window_size))
            cull_read(position_to_cull, coverage, candidate_reads, rng)
            # Recalculate median for the window
            window_coverage_median = np.median(coverage[i: i+window_size])
            # Test to make sure this doesn't completely mess up our calculations
            # (needed for paired-end in particular)
            if np.median(coverage) <= target_med:
                break

    # If we didn't find enough areas of too high coverage to bring our total under an acceptable parameter,
    # we'll throw out a couple random reads
    if np.median(coverage) > target_med + (2 * target_std):
        # We'll pare slightly more than we need to
        while np.median(coverage) > target_med + (2 * target_std) and len(candidate_reads) > 0:
            reads_culled += 1
            # Pick a random read to toss and update coverage
            cull_read(candidate_reads, coverage, candidate_reads, rng)

    _LOG.debug(f"Culled {reads_culled} reads")
    return candidate_reads


def cull_read(position, coverage_vector, total_reads, rng):
    """
    This function picks a read to toss and updates coverage.

    :param position: A position we want to thin out
    :param coverage_vector: The current coverage vector for the area
    :param total_reads: The entire reads dataset, from which to remove the culled read
    :param rng: The random number generator for the run
    :return: None. Updates total_reads in place
    """
    # Pick a random read to toss and update coverage
    depth = coverage_vector[position]
    random_selection = rng.choice(range(depth))
    current_depth = 0
    # Hoping this will be faster because i don't cycle through every read every time
    read_to_cull = (0, 0, 0, 0)
    for i in range(len(total_reads)):
        found = False
        if total_reads[i][0] <= position < total_reads[i][1]:
            current_depth += 1
            found = True
        elif total_reads[i][2] <= position < total_reads[i][3]:
            current_depth += 1
            found = True
        if current_depth == random_selection and found:
            read_to_cull = total_reads[i]
            break

    if not any(read_to_cull):
        read_to_cull = tuple(rng.choice(total_reads))

    for k in range(read_to_cull[0], read_to_cull[1]):
        coverage_vector[k] -= 1
    for m in range(read_to_cull[2], read_to_cull[3]):
        coverage_vector[m] -= 1
    # toss the read
    total_reads.remove(read_to_cull)


def run_one_more(coverage, index, max_runs, coverage_target):
    """
    This checks if our index as at the max (i.e., we are on the last loop). If we are and we haven't hit our median
    coverage target (to within rounding error), then we'll run one more loop

    :param coverage: The coverage vector for the simulation
    :param index: The current index of the run
    :param max_runs: The maximum number of loops we are running
    :param coverage_target: The median target for the dataset
    :return: -1 if we need to rerun a loop, 1 otherwise
    """
    if index == max_runs:
        if np.median(coverage) < coverage_target + 0.5:
            return -1

    return 1


def find_applicable_mutations(my_read: Read, all_variants: ContigVariants) -> dict:
    """
    Scans the variants dict for appropriate mutations.

    :param my_read: The read object to add the mutations to
    :param all_variants: All the variants for the dataset
    :return: A list of relevant mutations
    """
    return_dict = {}
    left = bisect_left(all_variants.variant_locations, my_read.position)
    right = bisect_right(all_variants.variant_locations, my_read.end_point - 1)
    subset = all_variants.variant_locations[left: right]
    for index in subset:
        return_dict[index] = all_variants[index]
    return return_dict


def replace_n(segment: Seq, rng: Generator) -> Seq:
    """
    Replaces invalid characters with random valid characters

    :param segment: A Seq object containing a DNA sequence
    :param rng: The random number generator for the run
    :return: The modified sequence object
    """
    modified_segment = ""
    for base in segment:
        if base not in ALLOWED_NUCL:
            modified_segment += rng.choice(ALLOWED_NUCL)
        else:
            modified_segment += base

    return Seq(modified_segment)


def modify_target_coverage(included_regions: list, excluded_regions: list, coverage_vector: np.ndarray):
    """
    Modifies the coverage vector by applying the list of regions. For this version, areas
    outside the regions have coverage adjusted by the off_target_percent

    :param included_regions: A list of intervals to target, extracted from a bed file
    :param excluded_regions: A list of regions to throw out, extracted from a bed file
    :param coverage_vector: The target coverage vector, which will be modified
    :return: The updated target coverage vector.
    """

    # this will tabulate values for included regions first, then for excluded regions. Hopefully both are not present.
    # If only one is present, the other should not change it.
    for region in included_regions + excluded_regions:
        coverage_vector[region[0]: region[1]] = coverage_vector[region[0]: region[1]] * region[2]

    return coverage_vector


def generate_reads(reference: SeqRecord,
                   error_model: SequencingErrorModel,
                   gc_bias: GcModel,
                   fraglen_model: FragmentLengthModel,
                   readlen_model: FragmentLengthModel,
                   contig_variants: ContigVariants,
                   temporary_directory: str | Path,
                   targeted_regions: list,
                   discarded_regions: list,
                   options: Options,
                   chrom: str,
                   ref_start: int = 0
                   ) -> tuple:
    """
    This will generate reads given a set of parameters for the run. The reads will output in a fastq.

    :param reference: The reference segment that reads will be drawn from.
    :param error_model: The error model for this run
    :param gc_bias: The GC-Bias model for this run
    :param fraglen_model: The fragment length model for this run
    :param readlen_model: The read length model for this run
    :param contig_variants: An object containing all input and randomly generated variants to be included.
    :param temporary_directory: The directory where to store temporary files for the run
    :param targeted_regions: A list of regions to target for the run (at a rate defined in the options
        file or 2% retained by default)
    :param discarded_regions: A list of regions to discard for the run
    :param options: The options entered for this run by the user
    :param chrom: The chromosome this reference segment originates from
    :param ref_start: The start point for this reference segment. Default is 0 and this is currently not fully
        implemented, to be used for parallelization.
    :return: A tuple of the filenames for the temp files created
    """
    # Set up files for use. May not need r2, but it's there if we do.
    chrom_fastq_r1 = temporary_directory / f'{chrom}_tmp_r1.fq.bgz'
    chrom_fastq_r2 = temporary_directory / f'{chrom}_tmp_r2.fq.bgz'

    # set up a temporary 'sam' file for processing by generate_bam, if the option is set
    tsam = temporary_directory / f'{chrom}.tsam.gz'

    _LOG.info(f'Sampling reads...')
    start = time.time()

    base_name = f'{Path(options.output).name}-{chrom}'
    target_coverage_vector = np.full(shape=len(reference), fill_value=options.coverage)
    if not options.no_coverage_bias:
        target_coverage_vector = gc_bias.create_coverage_bias_vector(reference.seq)

    target_coverage_vector = modify_target_coverage(targeted_regions, discarded_regions, target_coverage_vector)

    reads = cover_dataset(
        len(reference),
        target_coverage_vector,
        gc_bias.window_size,
        options,
        readlen_model,
        fraglen_model,
    )

    paired_reads = sorted([x for x in reads if any(x[0:2]) and any(x[2:4])])
    singletons = sorted([x for x in reads if x not in paired_reads and any(x)])

    # I'm worried about sorting the sam later, so trying to take care of that now
    sam_order = np.concatenate((paired_reads, singletons))
    sam_read_order = [""] * len(sam_order)

    previous_percent = 0

    _LOG.debug("Writing fastq(s) and optional tsam, if indicated")
    with (
        open_output(chrom_fastq_r1) as fq1,
        open_output(chrom_fastq_r2) as fq2,
        open_output(tsam) as temp_sam
    ):

        for i in np.arange(len(reads)):
            # Simple progress tracker
            current_percent = (i * 100) // len(reads)
            if current_percent > previous_percent:
                print(f'{current_percent}%', end='\r')
                previous_percent = current_percent

            # Added some padding after, in case there are deletions
            segments = [reference[reads[i][0]: reads[i][1] + 50].seq, reference[reads[i][2]: reads[i][3] + 50].seq]
            # Check for N concentration
            # Make sure the segments will come out to at around 80% valid bases
            # So as not to mess up the calculations, we will skip the padding we added above.
            if not is_too_many_n(segments[0][:-50] + segments[1][:-50]):
                read_name = f'{base_name}-{str(i)}'
                # reads[i] = (left_start, left_end, right_start, right_end)
                if any(reads[i][0:2]):
                    segment = segments[0]
                    segment = replace_n(segment, options.rng)

                    if any(reads[i][2:4]):
                        is_paired = True
                    else:
                        is_paired = False

                    read1 = Read(name=read_name,
                                 raw_read=reads[i],
                                 reference_segment=segment,
                                 reference_id=reference.id,
                                 position=reads[i][0] + ref_start,
                                 end_point=reads[i][1] + ref_start,
                                 quality_offset=options.quality_offset,
                                 is_paired=is_paired
                                 )

                    # The zero, zero slice of the np.where command is because this array is a 2d numpy array, so
                    # to get the first index of the first full array that matched, we need [0][0]
                    # we only have to do this once since reads[i] contains info on both reads.
                    # The idea here is to have the order established above, find where this fits in the order, and set
                    # The read. This will allow us to find these reads when we write out the full file
                    sam_read_order[(np.where(np.isin(sam_order[:, 1], reads[i])))[0][0]] = read_name
                    # TODO There may be a more efficient way to add mutations
                    read1.mutations = find_applicable_mutations(read1, contig_variants)

                    read1.write_record(error_model, fq1, temp_sam, options.produce_fastq, options.produce_bam)

                if any(reads[i][2:4]):
                    segment = segments[1]
                    segment = replace_n(segment, options.rng)

                    if any(reads[i][0:2]):
                        is_paired = True
                    else:
                        is_paired = False

                    read2 = Read(name=read_name,
                                 raw_read=reads[i],
                                 reference_segment=segment,
                                 reference_id=reference.id,
                                 position=reads[i][2] + ref_start,
                                 end_point=reads[i][3] + ref_start,
                                 quality_offset=options.quality_offset,
                                 is_reverse=True,
                                 is_paired=is_paired
                                 )

                    # TODO There may be a more efficient way to add mutations
                    read2.mutations = find_applicable_mutations(read2, contig_variants)
                    read2.write_record(error_model, fq2, temp_sam, options.produce_fastq, options.produce_bam)

    print("100%")

    _LOG.info(f"Finished sampling reads in {time.time() - start} seconds")
    return chrom_fastq_r1, chrom_fastq_r2, tsam, sam_read_order
