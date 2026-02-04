"""
Run the 'read-simulator' end-to-end on a tiny synthetic reference
"""

import sys
import yaml
import subprocess
import tempfile
from pathlib import Path

FA = ">chr1\n" + ("A" * 200) + "\n"

def test_basic_cli():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        ref = td / "ref.fa"
        cfg = td / "neat.yaml"
        out = td / "out"
        ref.write_text(FA)
        cfg_dict = {
            "reference": str(ref),
            "read_len": 50,
            "coverage": 1,
            "produce_fastq": True,
            "produce_bam": False,
            "produce_vcf": False,
            "threads": 1,
            "rng_seed": 1,
            "avg_seq_error": 0.0,
            "mutation_rate": 0.0,
        }
        cfg.write_text(yaml.safe_dump(cfg_dict, sort_keys=False), encoding="utf-8")
        # Use module invocation so Poetry install works in CI
        proc = subprocess.run(
            [sys.executable, "-m", "neat", "read-simulator", "-c", str(cfg), "-o", str(out)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert proc.returncode == 0, f"STDERR:\n{proc.stderr}"
        assert out.exists()
