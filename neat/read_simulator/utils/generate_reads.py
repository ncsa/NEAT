import logging
import pickle
import time
from math import ceil
from pathlib import Path

import numpy as np
from Bio.SeqRecord import SeqRecord
from bisect import bisect_left, bisect_right

from .output_file_writer import OutputFileWriter
from ...common import open_output
from ...models import SequencingErrorModel, FragmentLengthModel, TraditionalQualityModel, GCBiasModel
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
        reference: SeqRecord,
        options: Options,
        fragment_model: FragmentLengthModel | None,
        gc_model: GCBiasModel | None,
) -> list:
    """
    Covers a dataset to the desired depth in the paired ended case. This is the main algorithm for creating the reads
    to the proper coverage depth. It uses an abstract representation of the reads, by end points.

    :param reference: The reference sequence for this block
    :param options: The options for the run
    :param fragment_model: The fragment model used for to generate random fragment lengths
    :param gc_model: The GC bias model used for fragment selection
    """

    final_reads = []
    span_length = len(reference)
    # sanity check
    if span_length / fragment_model.fragment_mean < 5:
        _LOG.warning("The fragment mean is relatively large compared to the chromosome size. You may need to increase "
                     "standard deviation, or decrease fragment mean, if NEAT cannot complete successfully.")

    # precompute how many reads we want
    # The numerator is the total number of base pair calls needed.
    # Divide that by read length gives the number of reads needed
    if options.paired_ended:
        # TODO use gc bias to skew this number. Calculate at the runner level.
        number_reads = ceil(span_length * options.coverage / (2 * options.read_len))
    else:
        number_reads = ceil(span_length * options.coverage / options.read_len)

    if gc_model and not gc_model.is_uniform:
        # CDF-based sampling for GC bias
        window_size = gc_model.window_size
        if span_length <= window_size:
            # Fallback to uniform if region is too short
            return _uniform_sampling(span_length, number_reads, options, fragment_model)

        # Build prefix sum of weights
        # We need weights for every possible start position that can produce a read.
        # For single-ended, start must be in [0, span_length - read_len]
        max_start = span_length - options.read_len
        if max_start < 0:
            return []
            
        n_positions = max_start + 1
        weights = np.zeros(n_positions, dtype=float)
        
        # Initial window GC count
        seq_str = str(reference.seq).upper()
        effective_window = min(window_size, span_length)
        
        # The window for start position i is [i, i + effective_window]
        # BUT we must ensure i + effective_window does not exceed span_length.
        # Sliding window logic:
        # i = 0: window [0, effective_window]
        # last i = n_positions - 1: window [max_start, max_start + effective_window]
        # max_start + effective_window = span_length - read_len + min(window_size, span_length)
        # If window_size > read_len, then max_start + window_size > span_length!
        
        # We need to cap the window end at span_length.
        
        gc_count = seq_str[:effective_window].count('G') + seq_str[:effective_window].count('C')
        n_count = seq_str[:effective_window].count('N')
        
        def get_weight(gc, n, win_len):
            called = win_len - n
            if called <= 0:
                return 1.0
            return gc_model.get_weight(gc / called)

        weights[0] = get_weight(gc_count, n_count, effective_window)
        
        # Sliding window
        for i in range(1, n_positions):
            # Start position i means window [i, min(i + window_size, span_length)]
            # Previous window was [i-1, min(i-1 + window_size, span_length)]
            
            # Outgoing base is always i-1
            outgoing = seq_str[i-1]
            if outgoing in 'GC':
                gc_count -= 1
            elif outgoing == 'N':
                n_count -= 1
                
            # Incoming base depends on whether the window is still expanding or sliding
            # current_window_end = min(i + window_size, span_length)
            # prev_window_end = min(i - 1 + window_size, span_length)
            
            current_window_len = effective_window
            if i + window_size <= span_length:
                # Still sliding full window
                incoming = seq_str[i + window_size - 1]
                if incoming in 'GC':
                    gc_count += 1
                elif incoming == 'N':
                    n_count += 1
            else:
                # Window is shrinking at the end of the span
                current_window_len = span_length - i
                # No new base incoming, outgoing already removed
                
            weights[i] = get_weight(gc_count, n_count, current_window_len)
            
        prefix_sum = np.cumsum(weights)
        total_weight = prefix_sum[-1]
        
        if total_weight == 0:
            _LOG.debug("All positions in region have zero GC bias weight; no fragments generated.")
            return []

        mean_weight = total_weight / n_positions

        # Scale total reads by the mean/max ratio so that GC-rich regions receive
        # proportionally more reads than AT-rich regions across the genome.
        # For a uniform model mean_weight == max_weight, so this is a no-op.
        number_reads = ceil(number_reads * mean_weight / gc_model.max_weight)

        if number_reads == 0:
            return []

        # Batch CDF sampling: oversample so fragment-length filtering still leaves
        # enough reads. All operations are vectorized.
        n_candidates = number_reads * 10
        uniform_vals = options.rng.random(n_candidates) * total_weight
        starts = np.searchsorted(prefix_sum, uniform_vals).astype(int)
        starts = np.clip(starts, 0, max_start)

        frag_lengths = np.array(fragment_model.generate_fragments(n_candidates, options.rng))
        ends = np.minimum(starts + frag_lengths, span_length)

        valid = ends - starts >= options.read_len
        if options.paired_ended:
            valid &= ends - starts >= options.read_len + 10

        valid_starts = starts[valid][:number_reads]
        valid_ends = ends[valid][:number_reads]

        placed = len(valid_starts)
        if placed < number_reads:
            _LOG.warning(
                f"GC-weighted fragment generation placed {placed}/{number_reads} fragments "
                f"(short fragments filtered)."
            )

        for s, e in zip(valid_starts.tolist(), valid_ends.tolist()):
            read1 = (s, s + options.read_len)
            read2 = (e - options.read_len, e) if options.paired_ended else (0, 0)
            final_reads.append(read1 + read2)
             
    else:
        # Uniform sampling
        final_reads = _uniform_sampling(span_length, number_reads, options, fragment_model)

    # Now we shuffle them to add some randomness
    options.rng.shuffle(final_reads)
    return final_reads


def _uniform_sampling(span_length, number_reads, options, fragment_model):
    if span_length <= options.read_len:
        return []

    max_start = span_length - options.read_len
    # Oversample by 100x so fragment-length filtering (e.g. when frag mean < read_len)
    # still yields enough valid reads. All operations are vectorized.
    n_candidates = number_reads * 100
    starts = options.rng.integers(0, max_start + 1, size=n_candidates)
    frag_lengths = np.array(fragment_model.generate_fragments(n_candidates, options.rng))
    ends = np.minimum(starts + frag_lengths, span_length)

    valid = ends - starts >= options.read_len
    if options.paired_ended:
        valid &= ends - starts >= options.read_len + 10

    valid_starts = starts[valid][:number_reads]
    valid_ends = ends[valid][:number_reads]

    final_reads = []
    for s, e in zip(valid_starts.tolist(), valid_ends.tolist()):
        read1 = (s, s + options.read_len)
        read2 = (e - options.read_len, e) if options.paired_ended else (0, 0)
        final_reads.append(read1 + read2)
    return final_reads


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
        errors_per_read: int,
        qual_model: TraditionalQualityModel,
        fraglen_model: FragmentLengthModel,
        gc_model: GCBiasModel,
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
    :param errors_per_read: Total number of errors to add to contig
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

    if len(reference) < options.read_len:
        _LOG.warning(
            f"Contig '{contig_name}' (length {len(reference)}) is shorter than read_len "
            f"({options.read_len}). Skipping contig."
        )
        return []

    # _LOG.debug("Covering dataset.")
    t = time.time()
    reads = cover_dataset(
        reference,
        options,
        fraglen_model,
        gc_model,
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
            segment_start=read1[0] + ref_start,
            is_paired=options.paired_ended,
        )

        read_1.mutations = find_applicable_mutations(read_1, contig_variants)
        if options.produce_fastq:
            fastq_handle = ofw.files_to_write[ofw.fq1]
        else:
            fastq_handle = None
        read_1.finalize_read_and_write(
            error_model,
            qual_model,
            fastq_handle,
            options.quality_offset,
            options.produce_fastq,
            errors_per_read,
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
                segment_start=start_coordinate + ref_start,
                is_reverse=True,
                is_paired=options.paired_ended
            )

            read_2.mutations = find_applicable_mutations(read_2, contig_variants)
            if options.produce_fastq:
                fastq_handle = ofw.files_to_write[ofw.fq2]
            else:
                fastq_handle = None
            read_2.finalize_read_and_write(
                error_model,
                qual_model,
                fastq_handle,
                options.quality_offset,
                options.produce_fastq,
                errors_per_read,
                options.rng
            )
            reads_to_write.append((read_1, read_2))
        else:
            reads_to_write.append((read_1, None))

    _LOG.info(f"Finished sampling reads for thread {thread_index} in {(time.time() - start_time)/60:.2f} m")
    return reads_to_write
