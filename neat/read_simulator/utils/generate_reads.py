import logging
import time
import pickle
import numpy as np

from math import ceil
from pathlib import Path
from Bio import SeqRecord
from Bio.Seq import Seq
from numpy.random import Generator
from bisect import bisect_left, bisect_right

from ...models import SequencingErrorModel, GcModel, FragmentLengthModel, MutationModel
from .options import Options
from ...common import open_input, open_output, ALLOWED_NUCL
from ...variants import ContigVariants
from .read import Read

__all__ = [
    'generate_reads',
    'cover_dataset',
    'overlaps',
]

_LOG = logging.getLogger(__name__)


"""
Currently this code block is unused, but I think it might be useful for parallelization, so I wanted to leave it in.

def create_windows(sequence: SeqRecord, size: int, overlap: int):
    \"""
    Create a list of windows. Currently Unused. We might need this for parallelization, so I'm leaving the code in place

    :param sequence: Sequence to split
    :param size: size of windows
    :param overlap: size of overlap between windows
    :return: list of windows
    \"""
    windows = []
    for i in range(0, len(sequence), size):
        if i < overlap and i + size + overlap < len(sequence):
            windows.append((i, i + size + overlap))
        if i + size + overlap < len(sequence):
            windows.append((i - overlap, i + size + overlap))
        else:
            windows.append((i, len(sequence)))

    return windows
"""


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

    final_reads = set()
    # precompute how many reads we want
    # The numerator is the total number of base pair calls needed.
    # Divide that by read length gives the number of reads needed
    number_reads = ceil((span_length * options.coverage) / options.read_len)

    # We use fragments to model the DNA
    fragment_pool = fragment_model.generate_fragments(span_length, number_reads * 3)

    # step 1: Divide the span up into segments drawn froam the fragment pool. Assign reads based on that.
    # step 2: repeat above until number of reads exceeds number_reads * 1.5
    # step 3: shuffle pool, then draw number_reads (or number_reads/2 for paired ended) reads to be our reads
    read_count = 0
    while read_count <= number_reads:
        start = 0
        temp_fragments = []
        # Breaking the gename into fragments
        while start < span_length:
            # We take the first element and put it back on the end to create an endless pool of fragments to draw from
            fragment = fragment_pool.pop(0)
            end = min(start + fragment, span_length)
            # these are equivalent of reads we expect the machine to filter out, but we won't actually use it
            if end - start < options.read_len:
                # add some random flavor to try to keep it to falling into a loop
                if options.rng.normal() < 0.5:
                    fragment_pool.insert(fragment, len(fragment_pool)//2)
                else:
                    fragment_pool.insert(fragment, len(fragment_pool) - 3)
            else:
                fragment_pool.append(fragment)
                temp_fragments.append((start, end))
            start = end

        # Generating reads from fragments
        for fragment in temp_fragments:
            read_start = fragment[0]
            read_end = read_start + options.read_len
            # This filters out those small fragments, to give the dataset some realistic variety
            if read_end > fragment[1]:
                continue
            else:
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
                read = read1 + read2
                if read not in final_reads:
                    final_reads.add(read)
                    read_count += 1

    # Convert set to final list
    final_reads = list(final_reads)
    # Now we shuffle them to add some randomness
    options.rng.shuffle(final_reads)
    # And only return the number we needed
    if options.paired_ended:
        # Since each read is actually 2 reads, we only need to return half as many. But to cheat a few extra, we scale
        # that down slightly to 1.85 reads per read. This factor is arbitrary and may even be a function. But let's see
        # how well this estimate works
        return final_reads[:ceil(number_reads/1.85)]
    else:
        # Each read lacks a pair, so we need the full number of single ended reads
        return final_reads[:number_reads]


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


def generate_reads(reference: SeqRecord,
                   reads_pickle: str,
                   error_model_1: SequencingErrorModel,
                   error_model_2: SequencingErrorModel | None,
                   mutation_model: MutationModel,
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
    :param mutation_model: The mutation model for this run
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

    _LOG.debug("Covering dataset.")
    t = time.process_time()
    reads = cover_dataset(
        len(reference),
        options,
        fraglen_model,
    )
    _LOG.debug(f"Dataset coverage took: {(time.process_time() - t)/60} m")

    # Filters left to apply: N filter... not sure yet what to do with those (output random low quality reads, maybe)
    # Target regions
    # discard regions
    # GC-bias figure out some way to throw out 75% of reads with high GC or AT content (see original gc bias model).
    # heterozygous variants. Okay, so I've been applying the variant to all but for heterozygous we need to apply the
    # variant to only some

    # These will hold the values as inserted.
    properly_paired_reads = []
    singletons = []

    _LOG.debug("Writing fastq(s) and optional tsam, if indicated")
    t = time.process_time()
    with (
        open_output(chrom_fastq_r1) as fq1,
        open_output(chrom_fastq_r2) as fq2
    ):

        for i in range(len(reads)):
            print(f'{i/len(reads):.2%}', end='\r')
            # First thing we'll do is check to see if this read is filtered out by a bed file
            read1, read2 = (reads[i][0], reads[i][1]), (reads[i][2], reads[i][3])
            found_read1, found_read2 = False, False
            # For no target bed, there wil only be one region to check and this will complete very quickly
            for region in targeted_regions:
                # If this is a false region, we can skip it
                if not region[2]:
                    continue
                # We need to make sure this hasn't been filtered already, so if any coordinate is nonzero (falsey)
                if any(read1):
                    # Check if read1 is in this targeted region (any part of it overlaps)
                    if overlaps(read1, (region[0], region[1])):
                        found_read1 = True
                # Again, make sure it hasn't been filtered, or this is a single-ended read
                elif any(read2):
                    if overlaps(read2, (region[0], region[1])):
                        found_read2 = True
            # This read was outside of targeted regions
            if not found_read1:
                # Filter out this read
                read1 = (0, 0)
            if not found_read2:
                # Note that for single ended reads, it will never find read2 and this does nothing (it's already (0,0))
                read2 = (0, 0)

            # If there was no discard bed, this will complete very quickly
            discard_read1, discard_read2 = False, False
            for region in discarded_regions:
                # If region[2] is False then this region is not being discarded and we can skip it
                if not region[2]:
                    continue
                # Check to make sure the read isn't already filtered out
                if any(read1):
                    if overlaps(read1, (region[0], region[1])):
                        discard_read1 = True
                # No need to worry about already filtered reads
                if any(read2):
                    if overlaps(read2, (region[0], region[1])):
                        discard_read2 = True
            if discard_read1:
                read1 = (0, 0)
            if discard_read2:
                read2 = (0, 0)

            # This raw read will replace the original reads[i], and is the raw read with filters applied.
            raw_read = read1 + read2

            # If both reads were filtered out, we can move along
            if not any(raw_read):
                continue
            else:
                # We must have at least 1 read that has data
                # If only read1 or read2 is absent, this is a singleton
                read1_is_singleton = False
                read2_is_singleton = False
                properly_paired = False
                if not any(read2):
                    # Note that this includes all single ended reads that passed filter
                    read1_is_singleton = True
                elif not any(read1):
                    # read1 got filtered out
                    read2_is_singleton = True
                else:
                    properly_paired = True

            read_name = f'{base_name}-{str(i)}'

            # If the other read is marked as a singleton, then this one was filtered out, or these are single-ended
            if not read2_is_singleton:
                # It's properly paired if it's not a singleton
                is_paired = not read1_is_singleton
                # add a small amount of padding to the end to account for deletions
                # 30 is basically arbitrary. This may be a function of read length or determinable?
                # Though we don't figure out the variants until later
                segment = reference[read1[0]: read1[1] + 30]

                read1 = Read(name=read_name + "/1",
                             raw_read=raw_read,
                             reference_segment=segment,
                             reference_id=reference.id,
                             position=read1[0] + ref_start,
                             end_point=read1[1] + ref_start,
                             quality_offset=options.quality_offset,
                             is_paired=is_paired
                             )

                read1.mutations = find_applicable_mutations(read1, contig_variants)
                read1.finalize_read_and_write(error_model_1, mutation_model, fq1, options.produce_fastq)

            # if read1 is a sinleton then these are single-ended reads or this one was filtered out, se we skip
            if not read1_is_singleton:
                is_paired = not read2_is_singleton
                # Because this segment reads back to front, we need padding at the beginning.
                # 30 is arbitrary. This may be a function of read length or determinable?
                # Though we don't figure out the variants until later
                segment = reference[read2[0] - 30: read2[1]]

                read2 = Read(name=read_name + "/2",
                             raw_read=reads[i],
                             reference_segment=segment,
                             reference_id=reference.id,
                             position=read2[0] + ref_start,
                             end_point=read2[1] + ref_start,
                             quality_offset=options.quality_offset,
                             is_reverse=True,
                             is_paired=is_paired
                             )

                read2.mutations = find_applicable_mutations(read2, contig_variants)
                read2.finalize_read_and_write(error_model_2, mutation_model, fq2, options.produce_fastq)
            if properly_paired:
                properly_paired_reads.append((read1, read2))
            elif read1_is_singleton:
                # This will be the choice for all single-ended reads
                singletons.append(read1)
            else:
                singletons.append(read2)

    _LOG.debug(f"Fastqs written in: {(time.process_time() - t)/60:.2f} m")

    if options.produce_bam:
        # this will give us the proper read order of the elements, for the sam. They are easier to sort now
        properly_paired_reads = sorted(properly_paired_reads)
        singletons = sorted(singletons)
        sam_order = properly_paired_reads + singletons

        with open_output(reads_pickle) as reads:
            pickle.dump(sam_order, reads)

        if options.paired_ended:
            _LOG.debug(f"Properly paired percentage = {len(properly_paired_reads)/len(sam_order)}")

    _LOG.info(f"Finished sampling reads in {(time.process_time() - start)/60:.2f} m")
    return chrom_fastq_r1, chrom_fastq_r2
