import logging
import pickle
import time

from math import ceil
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from bisect import bisect_left, bisect_right

from .output_file_writer import OutputFileWriter
from ...common import open_output
from ...models import SequencingErrorModel, FragmentLengthModel, TraditionalQualityModel
from .options import Options
from ...variants import ContigVariants
from .read import Read

__all__ = [
    'generate_reads',
    'cover_dataset',
    'overlaps',
]

_LOG = logging.getLogger(__name__)


def cover_dataset(
        span_length: int,
        options: Options,
        fragment_model: FragmentLengthModel | None,
) -> list:
    """
    Covers a dataset to the desired depth in the paired ended case. This is the main algorithm for creating the reads
    to the proper coverage depth. It uses an abstract representation of the reads, by end points.

    :param span_length: The total length the cover needs to span
    :param options: The options for the run
    :param fragment_model: The fragment model used for to generate random fragment lengths
    """

    final_reads = []
    # sanity check
    if span_length/fragment_model.fragment_mean < 5:
        _LOG.warning("The fragment mean is relatively large compared to the chromosome size. You may need to increase "
                     "standard deviation, or decrease fragment mean, if NEAT cannot complete successfully.")
    # precompute how many reads we want
    # The numerator is the total number of base pair calls needed.
    # Divide that by read length gives the number of reads needed
    number_reads_per_layer = ceil(span_length / fragment_model.fragment_mean)
    number_reads = number_reads_per_layer * options.coverage

    # step 1: Divide the span up into segments drawn from the fragment pool. Assign reads based on that.
    # step 2: repeat above until number of reads exceeds number_reads * 1.5
    # step 3: shuffle pool, then draw number_reads (or number_reads/2 for paired ended) reads to be our reads
    read_count = 0
    loop_count = 0

    while read_count <= number_reads + 10:
        start = 0
        loop_count += 1
        # We use fragments to model the DNA
        fragment_pool = fragment_model.generate_fragments(number_reads_per_layer + 11, options.rng)

        temp_fragments = []
        # Mapping random fragments onto genome
        i = 0
        while start < span_length:
            # We take the first element and put it back on the end to create an endless pool of fragments to draw from
            fragment = fragment_pool[i]
            i = (i + 1) % len(fragment_pool)
            end = min(start + fragment, span_length)

            # Ensure the read is long enough to form a read, else we will not use it.
            if (end > options.read_len) and (end - start > options.read_len + 10):
                temp_fragments.append((start, end))
            start = end

        # Generating reads from fragments
        for fragment in temp_fragments:
            read_start = fragment[0]
            read_end = read_start + options.read_len
            read1 = (read_start, read_end)
            if options.paired_ended:
                # This will be valid because of the check above
                read2 = (fragment[1] - options.read_len, fragment[1])
            else:
                read2 = (0, 0)
            # The structure for these reads will be (left_start, left_end, right_start, right_end)
            # where start and end are ints with end > start. Reads can overlap, so right_start < left_end
            # is possible, but the reads cannot extend past each other, so right_start < left_start and
            # left_end > right_end are not possible.

            # sanity check that we haven't created an unrealistic read:
            read = read1 + read2
            final_reads.append(read)
            read_count += 1

    _LOG.debug(f"Coverage required {loop_count} loops")
    # Now we shuffle them to add some randomness
    options.rng.shuffle(final_reads)
    # And only return the number we needed
    return final_reads[:number_reads]


def find_applicable_mutations(my_read: Read, all_variants: ContigVariants) -> dict:
    """
    Scans the variants' dict for appropriate mutations.

    :param my_read: The read object to add the mutations to
    :param all_variants: All the variants for the contig
    :return: A list of relevant mutations
    """
    return_dict = {}
    left = bisect_left(all_variants.variant_locations, my_read.position)
    right = bisect_right(all_variants.variant_locations, my_read.end_point - 1)
    subset = all_variants.variant_locations[left: right]
    for index in subset:
        return_dict[index] = all_variants[index]
    return return_dict


def overlaps(test_interval: tuple[int, int], comparison_interval: tuple[int, int]) -> bool:
    """
    This function checks if the read overlaps with an input interval.
    :param test_interval: the interval to test, expressed as a tuple of end points
        (understood to be a half-open interval)
    :param comparison_interval: the interval to check against, expressed as a tuple of end points

    Four situations where we can say there is an overlap:
       1. The comparison interval contains the test interval start point
       2. The comparison interval contains the test interval end point
       3. The comparison interval contains both start and end points of the test interval
       4. The comparison interval is within the test interval
    Although 3 is really just a special case of 1 and 2, so we don't need a separate check

    If the read is equal to the interval, then all of these will be trivially true,
    and we don't need a separate check.
    """
    return (comparison_interval[0] < test_interval[1] < comparison_interval[1]) or \
           (comparison_interval[0] <= test_interval[0] < comparison_interval[1]) or \
           (test_interval[0] <= comparison_interval[0] and test_interval[1] >= comparison_interval[1])


def generate_reads(
        thread_index: int,
        reference: SeqRecord,
        error_model: SequencingErrorModel,
        qual_model: TraditionalQualityModel,
        fraglen_model: FragmentLengthModel,
        contig_variants: ContigVariants,
        targeted_regions: list,
        discarded_regions: list,
        options: Options,
        ofw: OutputFileWriter,
        contig_name: str,
        contig_index: int,
        ref_start: int,
):
    """
    This will generate reads given a set of parameters for the run. The reads will output in a fastq.

    :param thread_index: Index of current thread
    :param reference: The reference segment that reads will be drawn from.
    :param error_model: The error model for this run, the forward strand
    :param qual_model: The quality score model for this run, forward strand
    :param fraglen_model: The fragment length model for this run
    :param contig_variants: An object containing all input and randomly generated variants to be included.
    :param targeted_regions: A list of regions to target for the run (at a rate defined in the options
        file or 2% retained by default)
    :param discarded_regions: A list of regions to discard for the run
    :param options: The options entered for this run by the user
    :param ofw: the output file writer for the run
    :param contig_name: The name of the chromosome this ref segment originates from
    :param contig_index: The index of the above chromosome within the overall bam header
    :param ref_start: The start point for this reference segment. Default is 0 and this is currently not fully
        implemented, to be used for parallelization.
    :return: A tuple of the filenames for the temp files created
    """
    # _LOG.info(f'Sampling reads for thread {thread_index}...')
    start_time = time.time()

    # _LOG.debug("Covering dataset.")
    t = time.time()
    reads = cover_dataset(
        len(reference),
        options,
        fraglen_model,
    )
    # _LOG.debug(f"Dataset coverage took: {(time.time() - t)/60:.2f} m")

    # _LOG.debug("Writing fastq(s) and optional bam, if indicated")
    t = time.time()

    reads_to_write = []

    for i in range(len(reads)):
        # First thing we'll do is check to see if this read is filtered out by a bed file
        read1, read2 = (reads[i][0], reads[i][1]), (reads[i][2], reads[i][3])
        found_read = False
        # For no target bed, there wil only be one region to check and this will complete very quickly
        for region in targeted_regions:
            # If this is a false region, we can skip it
            if not region[2]:
                continue
            # We need to make sure this hasn't been filtered already, so if any coordinate is nonzero (falsey)
            if any(read1):
                # Check if read1 is in this targeted region (any part of it overlaps)
                if overlaps(read1, (region[0], region[1])):
                    found_read = True
            # Again, make sure it hasn't been filtered, or this is a single-ended read
            if any(read2):
                if overlaps(read2, (region[0], region[1])):
                    found_read = True
        # This read was outside targeted regions
        if not found_read:
            # Filter out this read
            continue

        # If there was no discard bed, this will complete very quickly
        discard_read = False
        for region in discarded_regions:
            # If region[2] is False then this region is not being discarded and we can skip it
            if not region[2]:
                continue
            # Check to make sure the read isn't already filtered out
            if any(read1):
                if overlaps(read1, (region[0], region[1])):
                    discard_read = True

            if any(read2):
                if overlaps(read2, (region[0], region[1])):
                    discard_read = True
        if discard_read:
            # toss the whole fragment
            continue

        if not any(read2) and options.paired_ended:
            # Marked paired, but no read 2 so we toss this one.
            continue
        block_read = read1 + read2
        raw_read = tuple(x + ref_start for x in block_read)
        # +1 to account for sam indexing
        read_name = f'NEAT_generated_{contig_name}_{thread_index}_{raw_read[0]+1:010d}_{raw_read[3]+1:010d}'

        # add a small amount of padding to the end to account for deletions.
        # Trying out this method of using the read-length, which for the default neat run gives ~30.
        padding = options.read_len//5
        segment = reference[read1[0]: read1[1] + padding].seq

        # if we're at the end of the contig, this may not pick up the full padding
        actual_padding = len(segment) - options.read_len

        read_1 = Read(
            name=read_name + "/1",
            raw_read=raw_read,
            reference_segment=segment,
            reference_id=contig_name,
            ref_id_index=contig_index,
            position=read1[0] + ref_start,
            end_point=read1[1] + ref_start,
            padding=actual_padding,
            run_read_len=options.read_len,
            is_paired=options.paired_ended,
        )

        read_1.mutations = find_applicable_mutations(read_1, contig_variants)
        read_1.finalize_read_and_write(
            error_model,
            qual_model,
            ofw.files_to_write[ofw.fq1],
            options.quality_offset,
            options.produce_fastq,
            options.rng
        )

        # skip over read 2 for single ended reads.
        if options.paired_ended:
            # Padding, as above
            padding = options.read_len//5
            start_coordinate = max((read2[0] - padding), 0)
            # this ensures that we get a segment with NEAT-recognized bases
            segment = reference[start_coordinate: read2[1]].seq
            # See note above
            actual_padding = len(segment) - options.read_len

            read_2 = Read(
                name=read_name + "/2",
                raw_read=raw_read,
                reference_segment=segment,
                reference_id=contig_name,
                ref_id_index=contig_index,
                position=read2[0] + ref_start,
                end_point=read2[1] + ref_start,
                padding=actual_padding,
                run_read_len=options.read_len,
                is_reverse=True,
                is_paired=options.paired_ended
            )

            read_2.mutations = find_applicable_mutations(read_2, contig_variants)

            read_2.finalize_read_and_write(
                error_model,
                qual_model,
                ofw.files_to_write[ofw.fq2],
                options.quality_offset,
                options.produce_fastq,
                options.rng
            )
            reads_to_write.append((read_1, read_2))
        else:
            reads_to_write.append((read_1, None))

    _LOG.info(f"Finished sampling reads for thread {thread_index} in {(time.time() - start_time)/60:.2f} m")
    return reads_to_write
