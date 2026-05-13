import gzip
import pickle
import textwrap
import pytest
import pysam
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from neat.read_simulator.runner import read_simulator_runner
from neat.model_gc_bias.runner import compute_gc_bias_runner

def _write_ref(path: Path, seq: str) -> Path:
    """Write a FASTA reference."""
    path.write_text(f">chr1\n{seq}\n", encoding="utf-8")
    return path

def _write_config(path: Path, ref_path: Path, **overrides) -> Path:
    defaults = {
        "reference": str(ref_path),
        "produce_fastq": "true",
        "produce_vcf": "false",
        "produce_bam": "false",
        "read_len": 50,
        "coverage": 2,
        "rng_seed": 42,
        "overwrite_output": "true",
        "cleanup_splits": "true",
    }
    defaults.update(overrides)
    lines = "\n".join(f"{k}: {v}" for k, v in defaults.items())
    path.write_text(lines + "\n", encoding="utf-8")
    return path

def create_fake_bam(bam_path, ref_path, read_count_at=10, read_count_gc=1000):
    """Creates a fake BAM file with more reads in GC-rich region."""
    # Reference is 10000bp: 5000bp A, 5000bp G
    header = {'HD': {'VN': '1.0'},
              'SQ': [{'LN': 10000, 'SN': 'chr1'}]}
    
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out_bam:
        # Reads in AT region (0-5000)
        for i in range(read_count_at):
            a = pysam.AlignedSegment()
            a.query_name = f"read_at_{i}"
            a.query_sequence = "A" * 50
            a.reference_id = 0
            a.reference_start = 1000 + i
            a.cigar = ((0, 50),)
            a.qual = "B" * 50
            out_bam.write(a)
            
        # Reads in GC region (5000-10000)
        for i in range(read_count_gc):
            a = pysam.AlignedSegment()
            a.query_name = f"read_gc_{i}"
            a.query_sequence = "G" * 50
            a.reference_id = 0
            a.reference_start = 6000 + i
            a.cigar = ((0, 50),)
            a.qual = "B" * 50
            out_bam.write(a)
    
    pysam.index(str(bam_path))

def test_gc_bias_full_pipeline(tmp_path):
    """
    Test the full GC bias pipeline:
    1. Create a reference and a BAM with GC bias.
    2. Run model-gc-bias to create a model.
    3. Run read-simulator using that model.
    4. Verify the output reads reflect the bias.
    """
    # 1. Setup data
    ref_seq = "A" * 5000 + "G" * 5000
    ref_path = _write_ref(tmp_path / "ref.fa", ref_seq)
    bam_path = tmp_path / "input.bam"
    # Create BAM where GC region has much more reads
    create_fake_bam(bam_path, ref_path, read_count_at=10, read_count_gc=1000)
    
    # 2. Run model-gc-bias
    model_prefix = "gc_model"
    compute_gc_bias_runner(
        bam_file=bam_path,
        reference_file=ref_path,
        output_dir=tmp_path,
        output_prefix=model_prefix,
        window_size=100,
        overwrite=True
    )
    
    model_path = tmp_path / (model_prefix + ".pickle.gz")
    assert model_path.exists()
    
    # 3. Run read-simulator with the generated model
    out_dir = tmp_path / "sim_out"
    cfg_path = _write_config(
        tmp_path / "sim.yml", 
        ref_path, 
        gc_model=str(model_path),
        coverage=1,
        read_len=50
    )
    
    read_simulator_runner(str(cfg_path), str(out_dir), "sim")
    
    # 4. Verify bias in simulated reads
    fq_files = list(out_dir.glob("*.fastq.gz"))
    assert len(fq_files) >= 1
    
    at_reads = 0
    gc_reads = 0
    
    # NEAT_generated_chr1_0_0000000101_0000000151
    # The read name contains the coordinates.
    with gzip.open(fq_files[0], "rt") as f:
        for line in f:
            if line.startswith("@NEAT_generated"):
                parts = line.split("_")
                start = int(parts[4])
                if start < 5000:
                    at_reads += 1
                else:
                    gc_reads += 1
                    
    print(f"Simulated AT reads: {at_reads}, GC reads: {gc_reads}")
    # Since our model was trained on 10x more GC reads, we expect more GC reads here too.
    # It won't be exactly 10x because of how cover_dataset works and smoothing in the model, 
    # but it should be significantly more.
    assert gc_reads > at_reads

def test_gc_bias_all_n_reference(tmp_path):
    """Test GC bias simulation with an all-N reference."""
    ref_seq = "N" * 1000
    ref_path = _write_ref(tmp_path / "ref_n.fa", ref_seq)
    
    # Uniform model
    model_path = tmp_path / "uniform_model.pickle.gz"
    from neat.models.gc_bias_model import get_uniform_gc_model
    get_uniform_gc_model().save(model_path)
    
    out_dir = tmp_path / "sim_n_out"
    cfg_path = _write_config(
        tmp_path / "sim_n.yml", 
        ref_path, 
        gc_model=str(model_path),
        coverage=1,
        read_len=50
    )
    
    # Should not crash and should produce reads (treated as neutral)
    read_simulator_runner(str(cfg_path), str(out_dir), "sim_n")
    
    fq_files = list(out_dir.glob("*.fastq.gz"))
    assert len(fq_files) >= 1
    with gzip.open(fq_files[0], "rt") as f:
        content = f.read()
    assert "@NEAT_generated" in content

def test_gc_bias_multi_contig(tmp_path):
    """Test GC bias simulation with multiple contigs."""
    ref_path = tmp_path / "multi.fa"
    ref_path.write_text(">chr1\n" + "A"*1000 + "\n>chr2\n" + "G"*1000 + "\n")
    
    # Model that favors G
    weights = [0.1] * 101
    weights[100] = 1.0
    from neat.models.gc_bias_model import GCBiasModel
    model = GCBiasModel(weights, window_size=100)
    model_path = tmp_path / "skewed_model.pickle.gz"
    model.save(model_path)
    
    out_dir = tmp_path / "multi_out"
    cfg_path = _write_config(
        tmp_path / "multi.yml", 
        ref_path, 
        gc_model=str(model_path),
        coverage=10,
        read_len=50
    )
    
    read_simulator_runner(str(cfg_path), str(out_dir), "sim_multi")
    
    fq_files = list(out_dir.glob("*.fastq.gz"))
    chr1_reads = 0
    chr2_reads = 0
    
    with gzip.open(fq_files[0], "rt") as f:
        for line in f:
            if line.startswith("@NEAT_generated_chr1"):
                chr1_reads += 1
            elif line.startswith("@NEAT_generated_chr2"):
                chr2_reads += 1
                
    # chr2 is 100% G, chr1 is 100% A. With our model, chr2 should have significantly more reads.
    assert chr2_reads > chr1_reads
