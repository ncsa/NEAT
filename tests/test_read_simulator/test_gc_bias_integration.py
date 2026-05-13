import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from neat.models.gc_bias_model import GCBiasModel
from neat.models.fragment_model import FragmentLengthModel
from neat.read_simulator.utils import Options
from neat.read_simulator.utils.generate_reads import cover_dataset

def test_gc_biased_sampling():
    # Create a reference with two halves: one 100% AT, one 100% GC
    span_length = 10000
    half = span_length // 2
    seq = "A" * half + "G" * half
    reference = SeqRecord(Seq(seq), id="chr1")
    
    # Create a GC bias model that favors GC-rich regions (100% weight for 100% GC, 0% for 0% GC)
    weights = [0.0] * 101
    weights[100] = 1.0
    weights[0] = 0.01 # Allow a tiny bit of AT to avoid infinite loops if possible
    gc_model = GCBiasModel(weights, window_size=100)
    
    options = Options(rng_seed=42)
    options.read_len = 50
    options.paired_ended = False
    options.coverage = 10
    options.overwrite_output = True
    
    fragment_model = FragmentLengthModel(150, 20)
    
    # Generate reads
    reads = cover_dataset(reference, options, fragment_model, gc_model)
    
    # Count how many reads are in the AT half vs GC half
    at_half_count = 0
    gc_half_count = 0
    
    for start, end, _, _ in reads:
        if start < half:
            at_half_count += 1
        else:
            gc_half_count += 1
            
    print(f"AT half reads: {at_half_count}, GC half reads: {gc_half_count}")
    
    # We expect significantly more reads in the GC half
    assert gc_half_count > at_half_count * 4
