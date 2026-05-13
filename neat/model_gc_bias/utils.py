"""
Utilities for modeling GC bias
"""

import numpy as np
import pysam
from Bio import SeqIO
from pathlib import Path

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
    """
    # 1. Sample windows across the genome
    # 2. For each window, calculate GC content and count reads
    # 3. Aggregate counts by GC percentage
    # 4. Normalize by the number of windows at each GC percentage
    
    gc_bins = [0] * 101
    gc_counts = [0.0] * 101
    
    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        ref_index = SeqIO.index(str(reference_file), "fasta")
        
        for contig in bam.references:
            if contig not in ref_index:
                continue
            
            contig_seq = ref_index[contig].seq
            contig_len = len(contig_seq)
            
            # Sample windows (not all windows to be faster)
            step = max(window_size, contig_len // 1000)
            for start in range(0, contig_len - window_size, step):
                end = start + window_size
                window_seq = contig_seq[start:end]
                gc_fraction = calculate_gc_content(str(window_seq))
                gc_percent = int(round(gc_fraction * 100))
                
                # Count reads in this window
                # Note: count() returns number of alignments overlapping the region
                read_count = bam.count(contig, start, end)
                
                gc_bins[gc_percent] += 1
                gc_counts[gc_percent] += read_count
                
    # Normalize
    weights = [0.0] * 101
    for i in range(101):
        if gc_bins[i] > 0:
            weights[i] = gc_counts[i] / gc_bins[i]
        else:
            weights[i] = 0.0
            
    # Rescale so max weight is 1.0
    max_w = max(weights) if any(weights) else 1.0
    if max_w > 0:
        weights = [w / max_w for w in weights]
    else:
        weights = [1.0] * 101
        
    # Smoothing might be good here, but keeping it simple for now
    # Fill in zeros with neighbors or global average if needed
    avg_w = sum(weights) / sum(1 for w in weights if w > 0) if any(w > 0 for w in weights) else 1.0
    for i in range(101):
        if weights[i] == 0:
            weights[i] = avg_w
            
    # Rescale again
    max_w = max(weights)
    weights = [w / max_w for w in weights]
            
    return weights
