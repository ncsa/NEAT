import gzip
import pickle
import pytest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

from neat.read_simulator.runner import read_simulator_runner
from neat.model_fragment_lengths.runner import compute_fraglen_runner
from neat.model_sequencing_error.runner import model_seq_err_runner

def _write_ref(path: Path, seq: str = "ACGT" * 2500) -> Path:
    """Write a FASTA reference (10kb default)."""
    path.write_text(f">chr1\n{seq}\n", encoding="utf-8")
    return path

def _write_config(path: Path, ref_path: Path, **overrides) -> Path:
    defaults = {
        "reference": str(ref_path),
        "produce_fastq": "true",
        "produce_vcf": "false",
        "produce_bam": "true",
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

def test_fraglen_model_pipeline(tmp_path):
    """Test training a fragment length model and using it in simulation."""
    ref_path = _write_ref(tmp_path / "ref.fa")
    
    # 1. Run simulation to get a BAM
    out_dir1 = tmp_path / "sim1"
    cfg_path1 = _write_config(tmp_path / "sim1.yml", ref_path, paired_ended="true", fragment_mean=300, fragment_st_dev=30)
    read_simulator_runner(str(cfg_path1), str(out_dir1), "sim1")
    
    # NEAT might append .bam to what we thought was the filename if not careful, 
    # but Options.from_cli usually uses the prefix.
    # Looking at runner.py: output_files = list(out_dir.glob("sim1*"))
    print(f"Files in {out_dir1}: {list(out_dir1.glob('*'))}")
    bam_path = out_dir1 / "sim1_golden.bam"
    assert bam_path.exists()
    
    # 2. Train fraglen model from BAM
    model_prefix = "fraglen_model"
    compute_fraglen_runner(
        file=bam_path,
        filter_minreads=1,
        output_dir=tmp_path,
        output_prefix=model_prefix,
        overwrite=True
    )
    
    model_path = tmp_path / (model_prefix + ".pickle.gz")
    assert model_path.exists()
    
    # Verify model content
    with gzip.open(model_path, 'rb') as f:
        model = pickle.load(f)
        assert hasattr(model, 'fragment_mean')
        # Should be around 300
        assert 250 < model.fragment_mean < 350
        
    # 3. Run simulation with this model
    out_dir2 = tmp_path / "sim2"
    cfg_path2 = _write_config(tmp_path / "sim2.yml", ref_path, fragment_model=str(model_path), paired_ended="true")
    read_simulator_runner(str(cfg_path2), str(out_dir2), "sim2")
    assert (out_dir2 / "sim2_r1.fastq.gz").exists()

def test_seq_err_model_pipeline(tmp_path):
    """Test training a sequencing error model and using it in simulation."""
    ref_path = _write_ref(tmp_path / "ref.fa")
    
    # 1. Run simulation to get FASTQs
    out_dir1 = tmp_path / "sim1"
    cfg_path1 = _write_config(tmp_path / "sim1.yml", ref_path, coverage=1)
    read_simulator_runner(str(cfg_path1), str(out_dir1), "sim1")
    
    # Check what files were actually produced
    print(f"Files in {out_dir1}: {list(out_dir1.glob('*'))}")
    # Based on previous failure, maybe it's sim1.fastq.gz if single ended?
    fq_path = out_dir1 / "sim1.fastq.gz"
    if not fq_path.exists():
        fq_path = out_dir1 / "sim1_r1.fastq.gz"
    
    assert fq_path.exists()
    
    # 2. Train sequencing error model from FASTQ
    model_prefix = "seq_err_model"
    model_seq_err_runner(
        files=[str(fq_path)],
        offset=33,
        qual_scores=42,
        max_reads=100,
        overwrite=True,
        output_dir=str(tmp_path),
        output_prefix=model_prefix
    )
    
    model_path = tmp_path / (model_prefix + ".p.gz")
    assert model_path.exists()
    
    # 3. Run simulation with this model
    out_dir2 = tmp_path / "sim2"
    cfg_path2 = _write_config(tmp_path / "sim2.yml", ref_path, error_model=str(model_path))
    read_simulator_runner(str(cfg_path2), str(out_dir2), "sim2")
    assert (out_dir2 / "sim2.fastq.gz").exists()
