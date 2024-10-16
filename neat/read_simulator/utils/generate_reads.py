import logging
import time
import pickle
import sys

from math import ceil
from pathlib import Path
from Bio import SeqRecord
from bisect import bisect_left, bisect_right

from ...models import SequencingErrorModel, FragmentLengthModel, MutationModel, QualityScoreModel
from .options import Options
from ...common import open_output
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

    final_reads = set()
    # sanity check
    if span_length/fragment_model.fragment_mean < 5:
        _LOG.warning("The fragment mean is relatively large compared to the chromosome size. You may need to increase "
                     "standard deviation, or decrease fragment mean, if NEAT cannot complete successfully.")
    # precompute how many reads we want
    # The numerator is the total number of base pair calls needed.
    # Divide that by read length gives the number of reads needed
    number_reads = ceil((span_length * options.coverage) / options.read_len)

    # We use fragments to model the DNA
    fragment_pool = fragment_model.generate_fragments(number_reads * 3)

    # step 1: Divide the span up into segments drawn from the fragment pool. Assign reads based on that.
    # step 2: repeat above until number of reads exceeds number_reads * 1.5
    # step 3: shuffle pool, then draw number_reads (or number_reads/2 for paired ended) reads to be our reads
    read_count = 0
    loop_count = 0
    while read_count <= number_reads:
        start = 0
        loop_count += 1
        # if loop_count > options.coverage * 100:
        #     _LOG.error("The selected fragment mean and standard deviation are causing NEAT to get stuck.")
        #     _LOG.error("Please try adjusting fragment mean or standard deviation to see if that fixes the issue.")
        #     _LOG.error(f"parameters:\n"
        #                f"chromosome length: {span_length}\n"
        #                f"read length: {options.read_len}\n"
        #                f"fragment mean: {options.fragment_mean}\n"
        #                f"fragment standard deviation: {options.fragment_st_dev}")
        #     sys.exit(1)
        temp_fragments = []
        # trying to get enough variability to harden NEAT against edge cases.
        if loop_count % 10 == 0:
            fragment_model.rng.shuffle(fragment_pool)
        # Mapping random fragments onto genome
        i = 0
        while start < span_length:
            # We take the first element and put it back on the end to create an endless pool of fragments to draw from
            fragment = fragment_pool[i]
            i = (i + 1) % len(fragment_pool)
            end = min(start + fragment, span_length)
            # Ensure the read is long enough to form a read, else we will not use it.
            if end - start > options.read_len:
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

                # sanity check that we haven't created an unrealistic read:
                insert_size = read2[0] - read1[1]
                if insert_size > 2 * options.read_len:
                    # Probably an outlier fragment length. We'll just pitch one of the reads
                    # and consider it lost to the ages.
                    if fragment_model.rng.choice((True, False)):
                        read1 = (0, 0)
                    else:
                        read2 = (0, 0)
                read = read1 + read2
                if read not in final_reads:
                    final_reads.add(read)
                    read_count += 1

    # Convert set to final list
    final_reads = list(final_reads)
    # Now we shuffle them to add some randomness
    fragment_model.rng.shuffle(final_reads)
    # And only return the number we needed
    _LOG.debug(f"Coverage required {loop_count} loops")
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
        reference: SeqRecord,
        reads_pickle: str,
        error_model_1: SequencingErrorModel,
        error_model_2: SequencingErrorModel | None,
        mut_model: MutationModel,
        markov_model: QualityScoreModel,
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

    :param reference: The reference segment that reads will be drawn from
    :param reads_pickle: The file to put the reads generated into, for bam creation
    :param error_model_1: The error model for this run, the forward strand
    :param error_model_2: The error model for this run, reverse strand
    :param mut_model: The mutation model for this run
    :param markov_model: The quality score model for this run
    :param fraglen_model: The fragment length model for this run
    :param contig_variants: An object containing all input and randomly generated variants to be included
    :param temporary_directory: The directory where to store temporary files for the run
    :param targeted_regions: A list of regions to target for the run (at a rate defined in the options
        file or 2% retained by default)
    :param discarded_regions: A list of regions to discard for the run
    :param options: The options entered for this run by the user
    :param chrom: The chromosome this reference segment originates from
    :param ref_start: The start point for this reference segment. Default is 0 and this is currently not fully
        implemented, to be used for parallelization
    :return: A tuple of the filenames for the temp files created
    """

    _LOG.debug(f"Generating reads for contig: {chrom}")
    _LOG.debug(f"Mutation model: {mut_model}, Mutation model RNG: {getattr(mut_model, 'rng', None)}")

    # Check if mut_model is None or doesn't have an RNG
    assert mut_model is not None, "Mutation model is None"
    assert mut_model.rng is not None, "Mutation model RNG is None"

    # Set up files for use. May not need r2, but it's there if we do.
    # We will separate the properly paired and the singletons.
    # For now, we are making an assumption that the chromosome name contains no invalid characters for bash file names
    # such as `*` or `:` even though those are technically allowed.
    # TODO We'll need to add some checks to ensure that this is the case.
    chrom_fastq_r1_paired = temporary_directory / f'{chrom}_r1_paired.fq.bgz'
    chrom_fastq_r1_single = temporary_directory / f'{chrom}_r1_single.fq.bgz'
    chrom_fastq_r2_paired = temporary_directory / f'{chrom}_r2_paired.fq.bgz'
    chrom_fastq_r2_single = temporary_directory / f'{chrom}_r2_single.fq.bgz'

    _LOG.info(f'Sampling reads...')
    start_time = time.time()

    base_name = f'NEAT-generated_{chrom}'

    _LOG.debug("Covering dataset.")
    t = time.time()
    reads = cover_dataset(
        len(reference),
        options,
        fraglen_model,
    )
    _LOG.debug(f"Dataset coverage took: {(time.time() - t)/60:.2f} m")

    # These will hold the values as inserted.
    properly_paired_reads = []
    singletons = []

    _LOG.debug("Writing fastq(s) and optional tsam, if indicated")
    t = time.time()

    # remove the combined model (ONLY use Markov)

    # Generate the quality scores using Markov model or traditional model
    if options.use_markov:
        markov_preds_df = QualityScoreModel(options.qual_score_model)
    else:
        markov_preds_df = None  # This is set to None and will be skipped in traditional mode

    with (
        open_output(chrom_fastq_r1_paired) as fq1_paired,
        open_output(chrom_fastq_r1_single) as fq1_single,
        open_output(chrom_fastq_r2_paired) as fq2_paired,
        open_output(chrom_fastq_r2_single) as fq2_single
    ):

        for i in range(len(reads)):
            print(f'{i/len(reads):.2%}', end='\r')
            # First thing we'll do is check to see if this read is filtered out by a bed file
            read1, read2 = (reads[i][0], reads[i][1]), (reads[i][2], reads[i][3])
            found_read1, found_read2 = False, False
            # For no target bed, there will only be one region to check and this will complete very quickly
            for region in targeted_regions:
                if not region[2]:
                    continue
                if any(read1) and overlaps(read1, (region[0], region[1])):
                    found_read1 = True
                if any(read2) and overlaps(read2, (region[0], region[1])):
                    found_read2 = True
            if not found_read1:
                read1 = (0, 0)
            if not found_read2:
                read2 = (0, 0)

            # Filter out by discarded regions
            discard_read1, discard_read2 = False, False
            for region in discarded_regions:
                if not region[2]:
                    continue
                if any(read1) and overlaps(read1, (region[0], region[1])):
                    discard_read1 = True
                if any(read2) and overlaps(read2, (region[0], region[1])):
                    discard_read2 = True
            if discard_read1:
                read1 = (0, 0)
            if discard_read2:
                read2 = (0, 0)

            raw_read = read1 + read2
            if not any(raw_read):
                continue

            read1_is_singleton = not any(read2)
            read2_is_singleton = not any(read1)
            properly_paired = not (read1_is_singleton or read2_is_singleton)

            read_name = f'{base_name}_{str(i+1)}'

            # Apply Markov chain quality scores if enabled, else use traditional model
            if options.use_markov:
                pred_qualities = markov_preds_df.iloc[i % len(markov_preds_df)].values
            else:
                pred_qualities = None

            # First read (if it's not a singleton)
            if not read2_is_singleton:
                is_paired = not read1_is_singleton
                padding = options.read_len // 5
                segment = reference[read1[0]: read1[1] + padding].seq
                actual_padding = len(segment) - options.read_len

                read_1 = Read(
                    name=read_name + "/1",
                    raw_read=raw_read,
                    reference_segment=segment,
                    reference_id=reference.id,
                    position=read1[0] + ref_start,
                    end_point=read1[1] + ref_start,
                    padding=actual_padding,
                    is_paired=is_paired
                )

                if options.use_markov:
                    read_1.qualities = pred_qualities[:len(segment)]
                read_1.mutations = find_applicable_mutations(read_1, contig_variants)

                handle = fq1_paired if is_paired else fq1_single
                read_1.finalize_read_and_write(
                    error_model_1,
                    mut_model,
                    markov_model,
                    handle,
                    options.quality_offset,
                    options.produce_fastq
                )

            # Second read (if not a singleton)
            if not read1_is_singleton:
                is_paired = not read2_is_singleton
                padding = options.read_len // 5
                start_coordinate = max((read2[0] - padding), 0)
                segment = reference[start_coordinate: read2[1]].seq
                actual_padding = len(segment) - options.read_len

                read_2 = Read(
                    name=read_name + "/2",
                    raw_read=reads[i],
                    reference_segment=segment,
                    reference_id=reference.id,
                    position=read2[0] + ref_start,
                    end_point=read2[1] + ref_start,
                    padding=actual_padding,
                    is_reverse=True,
                    is_paired=is_paired
                )

                if options.use_markov:
                    read_2.qualities = pred_qualities[:len(segment)]
                read_2.mutations = find_applicable_mutations(read_2, contig_variants)

                handle = fq2_paired if is_paired else fq2_single
                read_2.finalize_read_and_write(
                    error_model_2,
                    mut_model,
                    markov_model,
                    handle,
                    options.quality_offset,
                    options.produce_fastq
                )

            if properly_paired:
                properly_paired_reads.append((read_1, read_2))
            elif read1_is_singleton:
                singletons.append((read_1, None))
            else:
                singletons.append((None, read_2))

    _LOG.info(f"Contig fastq(s) written in: {(time.time() - t)/60:.2f} m")

    if options.produce_bam:
        properly_paired_reads = sorted(properly_paired_reads)
        singletons = sorted(singletons)
        sam_order = properly_paired_reads + singletons

        with open_output(reads_pickle) as reads:
            pickle.dump(sam_order, reads)

        if options.paired_ended:
            _LOG.debug(f"Properly paired percentage = {len(properly_paired_reads)/len(sam_order)}")

    _LOG.info(f"Finished sampling reads in {(time.time() - start_time)/60:.2f} m")
    return chrom_fastq_r1_paired, chrom_fastq_r1_single, chrom_fastq_r2_paired, chrom_fastq_r2_single
