"""
Utilities for modeling GC bias
"""

import numpy as np
import pysam
from Bio import SeqIO
from pathlib import Path


def compute_genome_wide_gc_mean_weight(reference_path: str | Path, gc_model) -> float:
    """
    Compute the genome-wide mean GC bias weight across every length-`window_size` window
    in the reference. Used at the runner level so that per-chunk read-count scaling can
    divide by the global mean (preserving the user-requested average coverage) rather
    than by the model's max weight (which would under-deliver coverage).

    For a uniform model this returns weights[0] without scanning the reference.
    All-N or sub-window contigs are skipped.
    """
    if gc_model.is_uniform:
        return float(gc_model.weights[0])

    window_size = gc_model.window_size
    weights_arr = np.asarray(gc_model.weights, dtype=np.float32)

    total_weight = 0.0
    total_positions = 0

    with pysam.FastaFile(str(reference_path)) as fa:
        for contig in fa.references:
            span_length = fa.get_reference_length(contig)
            if span_length <= window_size:
                continue
            seq = fa.fetch(contig).upper()
            n_positions = span_length - window_size + 1
            seq_arr = np.frombuffer(seq.encode(), dtype=np.uint8)
            gc_mask = (seq_arr == ord('G')) | (seq_arr == ord('C'))
            n_mask = seq_arr == ord('N')
            gc_cumsum = np.empty(span_length + 1, dtype=np.int64)
            n_cumsum = np.empty(span_length + 1, dtype=np.int64)
            gc_cumsum[0] = 0
            np.cumsum(gc_mask, out=gc_cumsum[1:])
            n_cumsum[0] = 0
            np.cumsum(n_mask, out=n_cumsum[1:])
            positions = np.arange(n_positions, dtype=np.int64)
            window_ends = positions + window_size
            gc_counts = (gc_cumsum[window_ends] - gc_cumsum[positions]).astype(np.float32)
            n_counts = (n_cumsum[window_ends] - n_cumsum[positions]).astype(np.float32)
            called = np.maximum(window_size - n_counts, 1.0)
            gc_indices = np.clip(np.rint(gc_counts / called * 100).astype(np.int16), 0, 100)
            weights = weights_arr[gc_indices]
            total_weight += float(weights.sum())
            total_positions += n_positions

    if total_positions == 0:
        return float(gc_model.weights[0])
    return total_weight / total_positions


def calculate_gc_content(sequence: str) -> float:
    """
    Calculates GC content of a sequence.
    """
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = sum(1 for base in sequence if base in 'ATGC')
    if total_count == 0:
        return 0.0
    return gc_count / total_count


def get_gc_bias_weights(bam_file: str | Path, reference_file: str | Path, window_size: int = 100) -> list[float]:
    """
    Estimates GC bias weights from a BAM file and a reference.

    For each contig, reads are fetched once and accumulated into a per-base
    coverage array.  Windows are then sampled from that array rather than
    issuing one bam.count() call per window, which eliminates O(n_windows)
    index lookups and is significantly faster on large genomes.
    """
    gc_bins = [0] * 101
    gc_counts = [0.0] * 101

    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        ref_index = SeqIO.index(str(reference_file), "fasta")

        for contig in bam.references:
            if contig not in ref_index:
                continue

            contig_seq = str(ref_index[contig].seq).upper()
            contig_len = len(contig_seq)

            # Single pass: accumulate per-base coverage
            coverage = np.zeros(contig_len, dtype=np.int32)
            for read in bam.fetch(contig):
                if read.is_unmapped or read.reference_end is None:
                    continue
                rs = read.reference_start
                re = min(read.reference_end, contig_len)
                coverage[rs:re] += 1

            # Sample windows using the pre-built coverage array
            step = max(window_size, contig_len // 1000)
            for start in range(0, contig_len - window_size, step):
                end = start + window_size
                gc_fraction = calculate_gc_content(contig_seq[start:end])
                gc_percent = int(round(gc_fraction * 100))

                gc_bins[gc_percent] += 1
                gc_counts[gc_percent] += int(coverage[start:end].sum())

    # Normalize: mean reads per window at each GC bin
    weights = [0.0] * 101
    for i in range(101):
        if gc_bins[i] > 0:
            weights[i] = gc_counts[i] / gc_bins[i]

    # Rescale so max weight is 1.0
    max_w = max(weights) if any(weights) else 1.0
    if max_w > 0:
        weights = [w / max_w for w in weights]
    else:
        weights = [1.0] * 101

    # Fill unobserved GC bins with the global mean so every bin has a positive weight
    positive = [w for w in weights if w > 0]
    avg_w = sum(positive) / len(positive) if positive else 1.0
    weights = [w if w > 0 else avg_w for w in weights]

    # Rescale again after fill
    max_w = max(weights)
    weights = [w / max_w for w in weights]

    return weights
