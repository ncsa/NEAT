import heapq
import logging
import pickle
import time
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

# Shortest insert accepted when short fragments are being kept (adapter readthrough or the
# adapter-free short-insert control). Real libraries size-select shorter molecules out, and a
# read that is almost entirely adapter carries no usable signal anyway.
MIN_SHORT_INSERT = 25

# Ceiling on rejection-sampling rounds in _sample_fragments. Each round draws at least ten
# candidates, so a workable fragment model converges in one or two; hitting the ceiling means
# the model cannot produce fragments this reference and read length can use, which is a
# configuration problem to report rather than one to spin on. This bound replaces the fixed
# spacer lengths FragmentLengthModel used to splice into every batch for the same purpose.
_MAX_SAMPLING_ROUNDS = 50


def _sample_fragments(
        number_reads: int,
        draw_starts,
        fragment_model: FragmentLengthModel,
        options: Options,
        min_frag: int,
        e_limit: int,
) -> tuple:
    """
    Rejection-sample (start, end) fragment windows until number_reads of them clear the floor.

    Shared by both placement strategies; they differ only in how a batch of start positions is
    drawn, which is what ``draw_starts`` supplies. Fragments are drawn in batches and the ones
    whose clamped span falls below ``min_frag`` are discarded, with each subsequent batch sized
    from the shortfall.

    :param number_reads: How many fragment windows are wanted.
    :param draw_starts: Callable taking a batch size and returning that many start positions.
    :param fragment_model: The fragment length model to sample lengths from.
    :param options: Run options (supplies the rng).
    :param min_frag: Shortest acceptable fragment span, from _min_fragment.
    :param e_limit: Upper bound for a fragment end, so mates stay inside this chunk.
    :return: (starts, ends), each a numpy array of at most number_reads entries. Short of that
        only when the fragment model cannot satisfy the floor, which is logged.
    """
    acc_starts: list[np.ndarray] = []
    acc_ends:   list[np.ndarray] = []
    collected = 0
    n_batch = number_reads * 2

    for _ in range(_MAX_SAMPLING_ROUNDS):
        if collected >= number_reads:
            break
        s = draw_starts(n_batch)
        fl = np.array(fragment_model.generate_fragments(n_batch, options.rng))
        e = np.minimum(s + fl, e_limit)
        mask = e - s >= min_frag
        acc_starts.append(s[mask])
        acc_ends.append(e[mask])
        collected += int(mask.sum())
        n_batch = max(10, (number_reads - collected) * 5)

    if collected < number_reads:
        _LOG.warning(
            f"Fragment model (mean {fragment_model.fragment_mean}, "
            f"st dev {fragment_model.fragment_st_dev}) produced only {collected} of "
            f"{number_reads} fragments at or above the {min_frag} bp floor. Coverage in this "
            f"region will be below the requested depth; check that the fragment mean suits the "
            f"read length and reference."
        )

    return (
        np.concatenate(acc_starts)[:number_reads],
        np.concatenate(acc_ends)[:number_reads],
    )


def _min_fragment(options: Options) -> int:
    """
    The shortest fragment accepted by either sampler.

    Ordinarily a fragment must be at least a full read long (plus 10 in paired mode, so mates do
    not entirely overlap); anything shorter is resampled. When short inserts are being kept
    deliberately, that floor drops to MIN_SHORT_INSERT so the left tail of the fragment
    distribution survives instead of being silently truncated away.
    """
    if options.adapters or options.keep_short_fragments:
        return MIN_SHORT_INSERT
    return options.read_len + (10 if options.paired_ended else 0)


def _read_windows(start: int, end: int, options: Options) -> tuple:
    """
    Map a sampled fragment onto its read1 and read2 reference windows.

    Both windows are clamped to the fragment: a read can never extend past the molecule it was
    sequenced from. When the insert is at least read_len long — always the case unless short
    fragments are being kept — the clamps are inert and this returns exactly the windows NEAT has
    always produced. For a shorter insert both mates collapse onto the whole fragment, which is
    what a real sequencer produces: R1 and R2 fully overlap, and each runs on into its adapter.

    :param start: Fragment start, relative to the reference chunk.
    :param end: Fragment end, relative to the reference chunk.
    :param options: Run options (reads read_len and paired_ended).
    :return: (r1_start, r1_end, r2_start, r2_end), with (0, 0) for read2 in single-ended mode.
    """
    read1 = (start, min(start + options.read_len, end))
    read2 = (max(end - options.read_len, start), end) if options.paired_ended else (0, 0)
    return read1 + read2


def _stochastic_round(expected_count: float, rng) -> int:
    """
    Round a fractional read count to an integer without biasing total coverage.

    Each chunk of the genome computes its own fractional expected read count. Rounding
    every chunk up (math.ceil) systematically over-covers: the bias is negligible at
    high coverage but dominates at low or fractional coverage, where a chunk's fair
    share can be well under one read yet still gets rounded up to a whole read. Rounding
    the fractional part probabilistically instead makes E[count] == expected_count, so
    the average coverage matches the requested target at any depth. The run's seeded RNG
    is used, so output remains reproducible for a given seed.

    :param expected_count: The fractional number of reads this chunk should produce.
    :param rng: The run's seeded numpy Generator.
    :return: A non-negative integer read count whose expectation is expected_count.
    """
    if expected_count <= 0:
        return 0
    floor = int(expected_count)
    remainder = expected_count - floor
    return floor + int(rng.random() < remainder)


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
    # region (which is also covered by the next chunk). expected_reads is kept as a float
    # and rounded stochastically (see _stochastic_round) so that summing across many
    # chunks does not over-cover — critical for low/fractional coverage.
    if options.paired_ended:
        expected_reads = responsibility_length * options.coverage / (2 * options.read_len)
    else:
        expected_reads = responsibility_length * options.coverage / options.read_len

    number_reads = _stochastic_round(expected_reads, options.rng)

    if gc_model and not gc_model.is_uniform:
        # CDF-based sampling for GC bias
        window_size = gc_model.window_size
        if span_length <= window_size:
            # Fallback to uniform if region is too short
            return _filter_n_regions(
                _uniform_sampling(span_length, number_reads, options, fragment_model,
                                  max_start=max_start),
                reference, options,
            )

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

        # Scale this chunk's read count by chunk_mean / global_mean. With CDF-biased
        # placement, expected genome-wide average coverage = options.coverage × (chunk_mean /
        # denominator). Using the genome-wide mean as denominator makes that ratio average
        # to 1.0 across chunks, so average coverage matches options.coverage (the documented
        # contract). The global mean is precomputed once at the runner level. When unset
        # (e.g., cover_dataset is called directly outside the runner), we fall back to the
        # chunk's own mean — i.e., no scaling, which is the correct single-chunk behavior.
        denominator = options.gc_global_mean_weight if getattr(
            options, "gc_global_mean_weight", None
        ) else mean_weight
        # Re-round from the unscaled float expectation so the GC scaling doesn't compound
        # a prior rounding step.
        number_reads = _stochastic_round(expected_reads * mean_weight / denominator, options.rng)

        if number_reads == 0:
            return []

        # Batch CDF sampling with adaptive retry (same pattern as _uniform_sampling).
        min_frag = _min_fragment(options)
        # For PE, the read2 record must stay within this chunk's responsibility so the
        # cat-stitched output remains coordinate-sorted (read2.position = e - read_len).
        # Cap e at responsibility_length + read_len so read2.position <= responsibility_length,
        # which lies at-or-before the next chunk's first read1 position. Falls back to
        # span_length for SE and for the final chunk (where responsibility = span_length).
        if options.paired_ended:
            e_limit = min(responsibility_length + options.read_len, span_length)
        else:
            e_limit = span_length
        def draw_cdf_starts(batch_size):
            uv = options.rng.random(batch_size) * total_weight
            return np.clip(np.searchsorted(prefix_sum, uv).astype(int), 0, max_start)

        valid_starts, valid_ends = _sample_fragments(
            number_reads, draw_cdf_starts, fragment_model, options, min_frag, e_limit,
        )

        for s, e in zip(valid_starts.tolist(), valid_ends.tolist()):
            final_reads.append(_read_windows(s, e, options))

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
    return _filter_n_regions(final_reads, reference, options)


def _filter_n_regions(final_reads: list, reference: SeqRecord, options: Options) -> list:
    """
    Drop fragments that fall mostly inside reference 'N' runs (assembly gaps).

    Under the default ``n_handling="exclude"`` policy, a read whose window is at least
    ``options.n_max_fraction`` 'N' is removed so that true gaps get ~zero coverage, mirroring
    real WGS. Reads that merely clip the edge of a gap survive and have their residual N's
    emitted as literal low-quality 'N' base calls in ``Read.convert_masking``. Under the legacy
    ``n_handling="telomere"`` policy this is a no-op (N's are filled with a telomere repeat
    downstream instead).

    The common case — a reference chunk with no 'N' at all — returns the input list unchanged
    after a single vectorized scan, so N-free references are neither reordered nor slowed.

    :param final_reads: List of (r1_start, r1_end, r2_start, r2_end) fragment tuples.
    :param reference: The reference chunk these positions index into.
    :param options: Run options (reads n_handling and n_max_fraction).
    :return: The filtered fragment list (a new list), or the input unchanged when nothing applies.
    """
    if options.n_handling == "telomere" or not final_reads:
        return final_reads

    seq_arr = np.frombuffer(str(reference.seq).upper().encode(), dtype=np.uint8)
    span = seq_arr.size
    n_mask = seq_arr == ord('N')
    if not n_mask.any():
        # No unknown bases in this chunk: nothing to exclude, preserve order exactly.
        return final_reads

    n_cumsum = np.empty(span + 1, dtype=np.int64)
    n_cumsum[0] = 0
    np.cumsum(n_mask, out=n_cumsum[1:])

    arr = np.asarray(final_reads, dtype=np.int64)
    r1s, r1e, r2s, r2e = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]

    def n_fraction(start, end):
        # Windows are half-open and always in-bounds (start >= 0, end <= span).
        length = np.maximum(end - start, 1)
        return (n_cumsum[end] - n_cumsum[start]) / length

    # Keep a fragment only if its read1 window is below the N threshold...
    keep = n_fraction(r1s, r1e) < options.n_max_fraction
    if options.paired_ended:
        # ...and, for paired reads, its read2 window too. Paired mates are always real,
        # in-bounds windows here (read2 = (e - read_len, e) with e >= read_len + 10), so no
        # degenerate-mate guard is needed — the (0, 0) placeholder only occurs single-ended.
        keep &= n_fraction(r2s, r2e) < options.n_max_fraction

    return [tuple(row) for row in arr[keep].tolist()]


def _uniform_sampling(span_length, number_reads, options, fragment_model, *,
                      max_start=None, responsibility_length=None):
    # A chunk's stochastically-rounded read count can legitimately be 0 at low/fractional
    # coverage. Bail early so we don't np.concatenate an empty accumulator list below.
    if number_reads <= 0:
        return []
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
    min_frag = _min_fragment(options)
    # For PE, cap e so read2.position stays within this chunk's responsibility (see GC path
    # for the full rationale).
    if options.paired_ended:
        e_limit = min(responsibility_length + options.read_len, span_length)
    else:
        e_limit = span_length

    # First batch: 2× candidates covers >99 % of cases when frag_mean >> read_len.
    # Retry in small increments only when fragment_mean < read_len (rare).
    def draw_uniform_starts(batch_size):
        return options.rng.integers(0, max_start + 1, size=batch_size)

    all_starts, all_ends = _sample_fragments(
        number_reads, draw_uniform_starts, fragment_model, options, min_frag, e_limit,
    )

    final_reads = []
    for s, e in zip(all_starts.tolist(), all_ends.tolist()):
        final_reads.append(_read_windows(s, e, options))
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
    # (reference start, counter, read), ordered by the position each record will actually carry.
    # For a reverse read that is not read.position: an indel moves where its alignment starts. Both
    # reads go through the buffer, not just read 2 — once the insert is short enough that the mates
    # cover the same window, read 2 can start before read 1 of its own fragment, so writing read 1
    # straight out leaves the BAM unsorted and unindexable.
    bam_buffer: list[tuple[int, int, "Read"]] = []
    bam_counter = 0
    # Fragments arrive in read-1 order, and no read can start earlier than its own fragment start
    # less the deletion headroom (read_len // 5). A whole read length is a generous bound on that,
    # so anything below the watermark can be written out without a later read undercutting it.
    bam_slack = options.read_len

    # Resolved once per chunk. Empty strings whenever readthrough is off, which makes every
    # adapter branch in Read a no-op and keeps output identical to a run without the feature.
    r1_adapter = options.adapter_r1 if options.adapters else ""
    r2_adapter = options.adapter_r2 if options.adapters else ""

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

        # Genomic bases available to read 1. Equals read_len unless the insert was shorter, in
        # which case the rest of the read is adapter (or, under keep_short_fragments alone, the
        # read simply ends early).
        insert_1 = read1[1] - read1[0]
        # add a small amount of padding to the end to account for deletions.
        # Trying out this method of using the read-length, which for the default neat run gives ~30.
        padding = options.read_len//5
        # Deletion headroom is drawn from the reference beyond the read, but for a short insert
        # that reference is past the end of the molecule and was never sequenced, so there is no
        # headroom to take.
        short_insert = insert_1 != options.read_len
        segment_end = read1[1] if short_insert else read1[1] + padding
        segment = reference[read1[0]: segment_end].seq

        # For a full-length insert the headroom is literal — extra reference bases the read pulls
        # in to stay read_len long after a deletion — and at the end of a contig there may be
        # fewer of them than asked for. A short insert has no reference to pull from, so its
        # headroom is a budget instead: a deletion shortens the genomic portion and the adapter
        # readthrough (or, with no adapter, the read simply ends earlier) absorbs the difference.
        # Passing 0 here would make every guard keyed on padding drop short-insert deletions
        # outright. See Read._resync_short_insert_lengths.
        actual_padding = padding if short_insert else len(segment) - insert_1

        read_1 = Read(
            name=read_name + "/1",
            raw_read=raw_read,
            reference_segment=segment,
            reference_id=contig_name,
            ref_id_index=contig_index,
            position=read1[0] + ref_start,
            end_point=read1[1] + ref_start,
            padding=actual_padding,
            run_read_len=options.read_len if r1_adapter else insert_1,
            segment_start=read1[0] + ref_start,
            is_paired=options.paired_ended,
            genomic_len=insert_1,
            adapter_seq=r1_adapter,
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
            options.rng,
            options.n_handling,
        )

        # Stream BAM in coordinate order: write out everything the watermark has cleared, then
        # buffer this read. Fragments arrive sorted by read1.position, so the watermark advances
        # monotonically.
        if bam_handle is not None:
            watermark = read_1.position - bam_slack
            while bam_buffer and bam_buffer[0][0] < watermark:
                _, _, buffered = heapq.heappop(bam_buffer)
                ofw.write_bam_record(buffered, contig_index, bam_handle, buffered.run_read_length)
            heapq.heappush(bam_buffer, (read_1.make_alignment()[1], bam_counter, read_1))
            bam_counter += 1

        # skip over read 2 for single ended reads.
        if options.paired_ended:
            insert_2 = read2[1] - read2[0]
            # Padding, as above
            padding = options.read_len//5
            # Read 2 is reverse, so its deletion headroom sits before the window — and, as for
            # read 1, it does not exist once the window already spans the whole fragment.
            short_insert = insert_2 != options.read_len
            if short_insert:
                start_coordinate = read2[0]
            else:
                start_coordinate = max((read2[0] - padding), 0)
            # this ensures that we get a segment with NEAT-recognized bases
            segment = reference[start_coordinate: read2[1]].seq
            # See note above
            actual_padding = padding if short_insert else len(segment) - insert_2

            read_2 = Read(
                name=read_name + "/2",
                raw_read=raw_read,
                reference_segment=segment,
                reference_id=contig_name,
                ref_id_index=contig_index,
                position=read2[0] + ref_start,
                end_point=read2[1] + ref_start,
                padding=actual_padding,
                run_read_len=options.read_len if r2_adapter else insert_2,
                segment_start=start_coordinate + ref_start,
                is_reverse=True,
                is_paired=options.paired_ended,
                genomic_len=insert_2,
                adapter_seq=r2_adapter,
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
                options.rng,
                options.n_handling,
            )
            if bam_handle is not None:
                heapq.heappush(bam_buffer, (read_2.make_alignment()[1], bam_counter, read_2))
                bam_counter += 1

    # Drain whatever the watermark never cleared; popping in heap order gives the sorted tail.
    if bam_handle is not None:
        while bam_buffer:
            _, _, buffered = heapq.heappop(bam_buffer)
            ofw.write_bam_record(buffered, contig_index, bam_handle, buffered.run_read_length)

    _LOG.info(f"Finished sampling reads for thread {thread_index} in {(time.time() - start_time)/60:.2f} m")
