import heapq
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
        *,
        responsibility_length: int | None = None,
) -> list:
    """
    Covers a dataset to the desired depth in the paired ended case. This is the main algorithm for creating the reads
    to the proper coverage depth. It uses an abstract representation of the reads, by end points.

    :param reference: The reference sequence for this block
    :param options: The options for the run
    :param fragment_model: The fragment model used for to generate random fragment lengths
    :param gc_model: The GC bias model used for fragment selection
    :param responsibility_length: The number of bases at the start of `reference` that this
        chunk is responsible for placing read1 starts in. Reads still extend into the
        trailing overlap region for context. Defaults to len(reference) (the chunk owns
        its full reference). For non-final chunks under sub-contig parallelism, this is
        the chunk step (chunk_size - overlap) — restricting read1 positions to the
        non-overlapping range so that per-chunk BAMs can be byte-concatenated into a
        coordinate-sorted whole without re-sorting.
    """

    final_reads = []
    span_length = len(reference)
    # Number of bases this chunk owns for read1 placement. Defaults to the full chunk.
    if responsibility_length is None:
        responsibility_length = span_length
    # Last valid read1 start position. Bounded by both the chunk's responsibility (so
    # reads don't appear in the next chunk's range) and by the requirement that the read
    # fit in available reference (so a read of length read_len has its tail in-bounds).
    max_start = min(responsibility_length - 1, span_length - options.read_len)
    # sanity check
    if span_length / fragment_model.fragment_mean < 5:
        _LOG.warning("The fragment mean is relatively large compared to the chromosome size. You may need to increase "
                     "standard deviation, or decrease fragment mean, if NEAT cannot complete successfully.")

    # precompute how many reads we want
    # The numerator is the total number of base pair calls needed. Coverage is scaled by
    # responsibility_length, not span_length, so chunks don't over-sample their overlap
    # region (which is also covered by the next chunk).
    if options.paired_ended:
        # TODO use gc bias to skew this number. Calculate at the runner level.
        number_reads = ceil(responsibility_length * options.coverage / (2 * options.read_len))
    else:
        number_reads = ceil(responsibility_length * options.coverage / options.read_len)

    if gc_model and not gc_model.is_uniform:
        # CDF-based sampling for GC bias
        window_size = gc_model.window_size
        if span_length <= window_size:
            # Fallback to uniform if region is too short
            return _uniform_sampling(span_length, number_reads, options, fragment_model,
                                     max_start=max_start)

        # Build prefix sum of weights only over the positions this chunk owns.
        if max_start < 0:
            return []
            
        n_positions = max_start + 1

        # Vectorized sliding window via cumulative sums — O(n) numpy, no Python loop.
        seq_arr = np.frombuffer(str(reference.seq).upper().encode(), dtype=np.uint8)
        gc_mask = (seq_arr == ord('G')) | (seq_arr == ord('C'))
        n_mask  = seq_arr == ord('N')
        gc_cumsum = np.empty(span_length + 1, dtype=np.int32)
        n_cumsum  = np.empty(span_length + 1, dtype=np.int32)
        gc_cumsum[0] = 0;  np.cumsum(gc_mask, out=gc_cumsum[1:])
        n_cumsum[0]  = 0;  np.cumsum(n_mask,  out=n_cumsum[1:])

        positions   = np.arange(n_positions, dtype=np.int32)
        window_ends = np.minimum(positions + window_size, span_length)
        gc_counts_arr = (gc_cumsum[window_ends] - gc_cumsum[positions]).astype(np.float32)
        n_counts_arr  = (n_cumsum[window_ends]  - n_cumsum[positions]).astype(np.float32)
        window_lens   = (window_ends - positions).astype(np.float32)
        called        = np.maximum(window_lens - n_counts_arr, 1.0)
        gc_indices    = np.clip(np.rint(gc_counts_arr / called * 100).astype(np.int16), 0, 100)
        weights       = np.asarray(gc_model.weights, dtype=np.float32)[gc_indices]

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

        # Batch CDF sampling with adaptive retry (same pattern as _uniform_sampling).
        min_frag = options.read_len + (10 if options.paired_ended else 0)
        # For PE, the read2 record must stay within this chunk's responsibility so the
        # cat-stitched output remains coordinate-sorted (read2.position = e - read_len).
        # Cap e at responsibility_length + read_len so read2.position <= responsibility_length,
        # which lies at-or-before the next chunk's first read1 position. Falls back to
        # span_length for SE and for the final chunk (where responsibility = span_length).
        if options.paired_ended:
            e_limit = min(responsibility_length + options.read_len, span_length)
        else:
            e_limit = span_length
        acc_starts: list[np.ndarray] = []
        acc_ends:   list[np.ndarray] = []
        collected = 0
        n_batch = number_reads * 2

        while collected < number_reads:
            uv = options.rng.random(n_batch) * total_weight
            s = np.clip(np.searchsorted(prefix_sum, uv).astype(int), 0, max_start)
            fl = np.array(fragment_model.generate_fragments(n_batch, options.rng))
            e = np.minimum(s + fl, e_limit)
            mask = e - s >= min_frag
            acc_starts.append(s[mask])
            acc_ends.append(e[mask])
            collected += int(mask.sum())
            n_batch = max(10, (number_reads - collected) * 5)

        valid_starts = np.concatenate(acc_starts)[:number_reads]
        valid_ends   = np.concatenate(acc_ends)[:number_reads]

        for s, e in zip(valid_starts.tolist(), valid_ends.tolist()):
            read1 = (s, s + options.read_len)
            read2 = (e - options.read_len, e) if options.paired_ended else (0, 0)
            final_reads.append(read1 + read2)
             
    else:
        # Uniform sampling
        final_reads = _uniform_sampling(
            span_length, number_reads, options, fragment_model,
            max_start=max_start, responsibility_length=responsibility_length,
        )

    # FASTQ is written in the natural fragment-sampling order (no explicit shuffle).
    # The BAM is sorted to coordinate order at the BAM-write boundary in single_runner,
    # which interleaves PE mate reads correctly. Users who want randomized FASTQ
    # ordering can pipe the output through `seqkit shuffle` — see README "FASTQ
    # output order".
    return final_reads


def _uniform_sampling(span_length, number_reads, options, fragment_model, *,
                      max_start=None, responsibility_length=None):
    if span_length <= options.read_len:
        return []

    # If caller didn't restrict the read1 placement range, default to the full
    # in-bounds span (last valid r1.position = span_length - read_len).
    if max_start is None:
        max_start = span_length - options.read_len
    if max_start < 0:
        return []
    if responsibility_length is None:
        responsibility_length = span_length
    min_frag = options.read_len + (10 if options.paired_ended else 0)
    # For PE, cap e so read2.position stays within this chunk's responsibility (see GC path
    # for the full rationale).
    if options.paired_ended:
        e_limit = min(responsibility_length + options.read_len, span_length)
    else:
        e_limit = span_length

    # First batch: 2× candidates covers >99 % of cases when frag_mean >> read_len.
    # Retry in small increments only when fragment_mean < read_len (rare).
    acc_starts: list[np.ndarray] = []
    acc_ends:   list[np.ndarray] = []
    collected = 0
    n_batch = number_reads * 2

    while collected < number_reads:
        s = options.rng.integers(0, max_start + 1, size=n_batch)
        fl = np.array(fragment_model.generate_fragments(n_batch, options.rng))
        e = np.minimum(s + fl, e_limit)
        mask = e - s >= min_frag
        acc_starts.append(s[mask])
        acc_ends.append(e[mask])
        collected += int(mask.sum())
        n_batch = max(10, (number_reads - collected) * 5)

    all_starts = np.concatenate(acc_starts)[:number_reads]
    all_ends   = np.concatenate(acc_ends)[:number_reads]

    final_reads = []
    for s, e in zip(all_starts.tolist(), all_ends.tolist()):
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
        responsibility_length: int | None = None,
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
    :return: None. FASTQ and BAM records are streamed directly to the output handles
        on `ofw` as reads are generated; no per-chunk accumulation of Read objects.
    """
    # _LOG.info(f'Sampling reads for thread {thread_index}...')
    start_time = time.time()

    if len(reference) < options.read_len:
        _LOG.warning(
            f"Contig '{contig_name}' (length {len(reference)}) is shorter than read_len "
            f"({options.read_len}). Skipping contig."
        )
        return

    # _LOG.debug("Covering dataset.")
    t = time.time()
    reads = cover_dataset(
        reference,
        options,
        fraglen_model,
        gc_model,
        responsibility_length=responsibility_length,
    )
    # _LOG.debug(f"Dataset coverage took: {(time.time() - t)/60:.2f} m")

    # Process fragments in read1-start order so the BAM emerges coordinate-sorted by
    # construction. In paired-end mode the read2 of each fragment is at a different
    # (typically later) position, so we hold those in a min-heap keyed by position and
    # flush each one before writing the next read1 that would precede it. The heap is
    # bounded by ~(fragment_length / read_length) entries — single-digit reads in
    # practice — so per-worker memory stays constant in chunk size and coverage.
    reads.sort(key=lambda r: r[0])

    # _LOG.debug("Writing fastq(s) and optional bam, if indicated")
    t = time.time()

    bam_handle = ofw.files_to_write[ofw.bam] if options.produce_bam else None
    r2_buffer: list[tuple[int, int, "Read"]] = []  # (position, counter, read)
    r2_counter = 0

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

        # Stream BAM in coordinate order: flush any buffered read2 records whose
        # position lies before this read1, then write read1 itself. Since fragments
        # are sorted by read1.position, read1 positions arrive monotonically.
        if bam_handle is not None:
            while r2_buffer and r2_buffer[0][0] < read_1.position:
                _, _, buffered_r2 = heapq.heappop(r2_buffer)
                ofw.write_bam_record(buffered_r2, contig_index, bam_handle, options.read_len)
            ofw.write_bam_record(read_1, contig_index, bam_handle, options.read_len)

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
            if bam_handle is not None:
                heapq.heappush(r2_buffer, (read_2.position, r2_counter, read_2))
                r2_counter += 1

    # Flush any read2 records still in the buffer — these all have positions at or
    # after the last read1 we wrote, so popping them in heap order gives the correct
    # coordinate-sorted tail.
    if bam_handle is not None:
        while r2_buffer:
            _, _, buffered_r2 = heapq.heappop(r2_buffer)
            ofw.write_bam_record(buffered_r2, contig_index, bam_handle, options.read_len)

    _LOG.info(f"Finished sampling reads for thread {thread_index} in {(time.time() - start_time)/60:.2f} m")
