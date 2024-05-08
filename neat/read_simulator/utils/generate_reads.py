import logging
import time
import pickle
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
    'generate_reads'
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
    n = segment.upper().count('N')
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
        options: Options,
        gc_bias_model: GcModel,
        fragment_model: FragmentLengthModel | None,
) -> list:
    """
    Covers a dataset to the desired depth in the paired ended case. This is the main algorithm for creating the reads
    to the proper coverage depth. It uses an abstract representation of the reads, by end points.

    :param span_length: The total length the cover needs to span
    :param target_vector: The vector of coverage targets by position
    :param options: The options for the run
    :param gc_bias_model: The gc-bias model being used for this run.
    :param fragment_model: The fragment model used for to generate random fragment lengths
    """

    # Variables to track reads and current coverage
    temp_reads = []
    final_coverage = [0] * len(target_vector)
    final_reads = []

    # Determine potential fragment lengths
    if fragment_model:
        # Paired ended case
        fragment_pool = fragment_model.generate_fragments(span_length, options.read_len, options.coverage)
    else:
        # single ended case
        fragment_pool = [options.read_len]

    _LOG.debug("Covering dataset")

    # We need a window to cover. If window size of the gc model is greater than read length, then we'll call it 1,
    # meaning 1 section of the target vector.
    read_window_width = max(max(fragment_pool)//gc_bias_model.window_size, 1)

    total_number_reads = 0
    # for each position, the set of reads covering that position.
    for i in range(len(target_vector)-read_window_width):
        section_id = (i, i+read_window_width)
        subsection_target_vector = target_vector[i: i+read_window_width]
        subsection_target = ceil(np.mean(subsection_target_vector))

        # Check if we even need to bother this section
        section_max = max(final_coverage[i: i+read_window_width])
        if section_max >= subsection_target:
            # we have enough reads here, move on.
            i += 1
            continue

        # We need to know where we are along the reference contig. We are i * read_window_widths in,
        # and each read_window is gc_bias_model.window_size long.
        start_coordinate = i * read_window_width * gc_bias_model.window_size
        # read_count counts the reads added to this section
        read_count = 0
        # 1.5x as many as needed to give the final set a bit of randomness.
        temp_pool_size = ceil(subsection_target * 1.5)
        # debug variable
        k = 0
        while read_count < temp_pool_size:
            # The structure for these reads will be (left_start, left_end, right_start, right_end)
            # where start and end are ints with end > start. Reads can overlap, so right_start < left_end
            # is possible, but the reads cannot extend past each other, so right_start < left_start and
            # left_end > right_end are not possible.
            if options.paired_ended:
                left_read_length = right_read_length = options.read_len
                # Select a random fragment length from the pool.
                fragment_length = options.rng.choice(fragment_pool)
            else:
                left_read_length = options.read_len
                right_read_length = 0
                fragment_length = left_read_length

            # Booleans telling us how to deal with the various exceptions below.
            skip_left = False
            skip_right = False

            left_start = start_coordinate
            # This trims the left read to the region of interest
            left_end = min(span_length, start_coordinate + left_read_length)

            # Make right read, if paired ended, else set it equal to the default
            if options.paired_ended:
                # The candidate end points of the paired right read, uncorrected. It starts from the far end, so we
                # know the end based on the fragment length and we count in from there.
                uncorrected_right_end = start_coordinate + fragment_length
                uncorrected_right_start = uncorrected_right_end - right_read_length

                # this trims the raw right read, so that right_start >= left_start and left_end <= right_end
                # this keeps us from getting unrealistic paired ended reads
                right_start = max(uncorrected_right_start, left_start)
                right_end = min(max(uncorrected_right_end, left_end), span_length)

                # Both reads are off the map, so we throw this out:
                if right_end <= 0:
                    # Move the start coordinate slightly to see if this keeps things moving.
                    start_coordinate = min(
                        max(start_coordinate + options.rng.integers(-10, 10), 0),
                        span_length - fragment_length - 1
                    )
                    k += 1
                    if k >= 1000:
                        _LOG.debug(f"skipped {k} times")
                    continue

            else:
                skip_right = True
                right_start = 0
                right_end = 0

            # Cases that might make us want to skip these reads:

            # Case 1: Truncated left read, unusable
            if left_end - left_start < left_read_length:
                # right read is usable if it's not off the map and is at least a read length long
                if right_start < span_length and \
                        (right_end - right_start >= right_read_length):
                    # skip the left read, keep the right read, indicating an unpaired right read
                    skip_left = True
                # else, skip both
                else:
                    skip_left = skip_right = True

            # Case 2: Truncated right read, unusable
            if right_end - right_start < right_read_length:
                # left read is usable if it's not off the map and has at least half a
                # read length of usable bases
                if left_start < span_length and \
                        (left_end - left_start >= left_read_length):
                    # Set the right read to (0, 0), indicating an unpaired left read
                    skip_right = True
                # else, skip both
                else:
                    skip_left = skip_right = True

            if skip_left and skip_right:
                # If neither read is usable, we'll just move to the next
                # Move the start coordinate slightly to see if this keeps things moving.
                start_coordinate = min(
                    max(start_coordinate + options.rng.integers(-10, 10), 0),
                    span_length - fragment_length - 1
                )
                k += 1
                if k >= 1000:
                    _LOG.debug(f"skipped {k} times")
                continue

            if not skip_left:
                temp_left_read = (left_start, left_end)
            else:
                temp_left_read = (0, 0)

            if not skip_right:
                temp_right_read = (right_start, right_end)
            else:
                temp_right_read = (0, 0)

            read_to_add = temp_left_read + temp_right_read
            temp_reads.append(read_to_add)
            read_count += 1
            # bump it a bit to give it a little variety.
            start_coordinate = min(
                max(start_coordinate + options.rng.integers(-10, 10), 0),
                span_length - fragment_length - 1
            )

        # Subset this set of temporary reads to make sure we don't overdo it.
        reads_to_add, final_coverage = final_subsetting(
            candidate_reads=temp_reads,
            target=subsection_target,
            coverage_vector=final_coverage,
            coordinates=section_id,
            rng=options.rng
        )
        final_reads.extend(reads_to_add)

    return final_reads


def final_subsetting(
        candidate_reads: list,
        target: int,
        coverage_vector: list,
        coordinates: tuple,
        rng: Generator
) -> (list, list):
    """
    The final subsetting for reads. We basically just calculate if we need to subset, then how many final reads we
    should have, and we draw them randomly without replacement until we have enough.

    :param candidate_reads: The list of reads in the following format: (left_start, left_end, right_start, right_end)
    :param target: The target coverage value for this section
    :param coverage_vector: The coverage vector, will be updated with new counts.
    :param coordinates: The indices of interest.
    :param rng: The random number generator for the run
    :return: Two lists, left and right reads, culled to proper coverage depth.
    """

    # Throw out fully empty reads.
    filtered_candidates = [x for x in candidate_reads if any(x)]

    # if we're under target, we want to make up the difference.
    # grab however many we need.
    ret_set = rng.choice(filtered_candidates, size=target, replace=False)

    # Update the coverage vector
    for i in range(coordinates[0], coordinates[1]):
        coverage_vector[i] += target

    return ret_set, coverage_vector


def find_applicable_mutations(my_read: Read, all_variants: ContigVariants) -> dict:
    """
    Scans the variants' dict for appropriate mutations.

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
    # This takes care of soft masking
    segment_no_mask = segment.upper()
    for base in segment_no_mask:
        if base in ALLOWED_NUCL:
            modified_segment += base
        else:
            modified_segment += rng.choice(ALLOWED_NUCL)

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


def merge_sort(my_array: np.ndarray):
    """
    This sorts the reads in reverse position, merging as it goes, in order to get the proper order

    :param my_array: the array to be sorted.
    :return: the sorted array
    """
    ret_array = my_array[my_array[:, 3].argsort()]
    ret_array = ret_array[ret_array[:, 2].argsort(kind='mergesort')]
    ret_array = ret_array[ret_array[:, 1].argsort(kind='mergesort')]
    ret_array = ret_array[ret_array[:, 0].argsort(kind='mergesort')]
    return ret_array


def generate_reads(reference: SeqRecord,
                   reads_pickle: str,
                   error_model_1: SequencingErrorModel,
                   error_model_2: SequencingErrorModel | None,
                   gc_bias: GcModel,
                   fraglen_model: FragmentLengthModel,
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
    :param reads_pickle: The file to put the reads generated into, for bam creation.
    :param error_model_1: The error model for this run, the forward strand
    :param error_model_2: The error model for this run, reverse strand
    :param gc_bias: The GC-Bias model for this run
    :param fraglen_model: The fragment length model for this run
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

    _LOG.info(f'Sampling reads...')
    start = time.time()

    base_name = f'{Path(options.output).name}-{chrom}'
    # We want to bin the coverage into window-sized segments to speed up calculations.
    # This divides the segment into len(reference) // window_size rounded up to the nearest int.
    target_shape = ceil(len(reference) // gc_bias.window_size)
    # Assume uniform coverage
    target_coverage_vector = np.full(shape=target_shape, fill_value=options.coverage)
    if not options.no_coverage_bias:
        pass
        # I'm trying to move this section into cover_dataset.
        # target_coverage_vector = gc_bias.create_coverage_bias_vector(target_coverage_vector, reference.seq.upper())

    # Apply the targeting/discarded rates.
    target_coverage_vector = modify_target_coverage(targeted_regions, discarded_regions, target_coverage_vector)

    reads = cover_dataset(
        len(reference),
        target_coverage_vector,
        options,
        gc_bias,
        fraglen_model,
    )

    # Reads that are paired
    paired_reads = np.asarray([tuple(x) for x in reads if any(x[0:2]) and any(x[2:4])])
    if paired_reads.size:
        paired_reads = merge_sort(paired_reads)

    # singletons
    if paired_reads.any():
        singletons = np.asarray([tuple(x) for x in reads if x not in paired_reads and any(x)])
    else:
        singletons = np.asarray([tuple(x) for x in reads if any(x)])
    if singletons.size:
        singletons = merge_sort(singletons)

    # determine sam read order. It should be paired reads, then singletons, unles one or the other is missing.
    if paired_reads.size and singletons.size:
        sam_read_order = np.concatenate((paired_reads, singletons))
    elif paired_reads.size:
        # if singletons is empty
        sam_read_order = paired_reads
    else:
        # if there are no paired reads
        sam_read_order = singletons

    final_sam_dict = {}

    if len(sam_read_order) > 0:
        _LOG.debug(f"Paired percentage = {len(paired_reads) / len(sam_read_order)}")
    else:
        _LOG.debug("Paired percentage = 0 (no sam read orders)")

    _LOG.debug("Writing fastq(s) and optional tsam, if indicated")
    with (
        open_output(chrom_fastq_r1) as fq1,
        open_output(chrom_fastq_r2) as fq2
    ):

        for i in range(len(reads)):
            print(f'{i/len(reads):.2%}', end='\r')
            # Added some padding after, in case there are deletions
            segments = [reference[reads[i][0]: reads[i][1] + 50].seq.upper(),
                        reference[reads[i][2]: reads[i][3] + 50].seq.upper()]
            # Check for N concentration
            # Make sure the segments will come out to at around 80% valid bases
            # So as not to mess up the calculations, we will skip the padding we added above.
            if not is_too_many_n(segments[0][:-50] + segments[1][:-50]):
                read_name = f'{base_name}-{str(i)}'

                if options.produce_bam:
                    # This index gives the position of the read for the final bam file
                    # [0][0] gives us the first index where this occurs (which should be the only one)
                    final_order_read_index = np.where(sam_read_order == reads[i])[0][0]
                    if final_order_read_index not in final_sam_dict:
                        final_sam_dict[final_order_read_index] = []

                # reads[i] = (left_start, left_end, right_start, right_end)
                # If there is a read one:
                if any(reads[i][0:2]):
                    segment = segments[0]
                    segment = replace_n(segment, options.rng)

                    if any(reads[i][2:4]):
                        is_paired = True
                    else:
                        is_paired = False

                    read1 = Read(name=read_name + "/1",
                                 raw_read=reads[i],
                                 reference_segment=segment,
                                 reference_id=reference.id,
                                 position=reads[i][0] + ref_start,
                                 end_point=reads[i][1] + ref_start,
                                 quality_offset=options.quality_offset,
                                 is_paired=is_paired
                                 )

                    read1.mutations = find_applicable_mutations(read1, contig_variants)

                    read1.finalize_read_and_write(error_model_1, fq1, options.produce_fastq)

                    if options.produce_bam:
                        # Save this for later
                        final_sam_dict[final_order_read_index].append(read1)

                elif options.produce_bam:
                    # If there's no read1, append a 0 placeholder
                    final_sam_dict[final_order_read_index].append(0)

                if any(reads[i][2:4]):
                    segment = segments[1]
                    segment = replace_n(segment, options.rng)

                    if any(reads[i][0:2]):
                        is_paired = True
                    else:
                        is_paired = False

                    read2 = Read(name=read_name + "/2",
                                 raw_read=reads[i],
                                 reference_segment=segment,
                                 reference_id=reference.id,
                                 position=reads[i][2] + ref_start,
                                 end_point=reads[i][3] + ref_start,
                                 quality_offset=options.quality_offset,
                                 is_reverse=True,
                                 is_paired=is_paired
                                 )

                    read2.mutations = find_applicable_mutations(read2, contig_variants)
                    read2.finalize_read_and_write(error_model_2, fq2, options.produce_fastq)

                    if options.produce_bam:
                        # Save this for later
                        final_sam_dict[final_order_read_index].append(read2)

                elif options.produce_bam:
                    # If there's no read2, append a 0 placeholder
                    final_sam_dict[final_order_read_index].append(0)

    if options.produce_bam:
        with open_output(reads_pickle) as reads:
            pickle.dump(final_sam_dict, reads)

    _LOG.info(f"Finished sampling reads in {time.time() - start} seconds")
    return chrom_fastq_r1, chrom_fastq_r2
