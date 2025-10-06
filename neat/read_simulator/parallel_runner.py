"""
Split a reference genome, run NEAT simulations in
parallel, and stitch their outputs back together.
"""

import subprocess
import sys
import time
from pathlib import Path
from typing import List, Tuple
from multiprocessing import Pool

from .utils import Options, Options
from .utils.split_inputs import main as split_main
from .utils.stitch_outputs import main as stitch_main
from .single_runner import read_simulator_single

import logging

__all__ = ["main"]

EXTENSIONS = ["gz", "fastq", "bam", "vcf"]

_LOG = logging.getLogger(__name__)


# Parallelization worker
def worker(args: Tuple[List[str], Path]):
    """
    Run one NEAT simulation in its own working directory.

    Each worker is executed in a separate process by the
    :class:`concurrent.futures.ProcessPoolExecutor`.
    """
    cmd, workdir = args
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.check_call(cmd, cwd=workdir)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"NEAT job failed (exit {e.returncode}) in {workdir.name}: {' '.join(cmd)}"
        ) from e
    return workdir


def main(options: Options) -> None:
    """
    This function executes the three high-level steps of parallelization:

    1. Split the reference sequence into chunks
    2. Launch NEAT's read simulator on each chunk
    3. Stitch the resulting outputs back into a single dataset
    """
    # 1) Split (or reuse)
    splits_files = []
    if options.reuse_splits:
        if list(options.splits_dir.glob("*.y*ml")):
            _LOG.info("[parallel] Recycling existing split FASTA/YAML files – skipping split step")
            split_sec = 0.0
    else:
        t0 = time.time()
        _LOG.info("[parallel] Splitting reference...")
        splits_files = split_main(options)
        split_sec = time.time() - t0

    _LOG.setLevel(logging.ERROR)
    # 2) Submit NEAT sims in a deterministic order using the manifest
    _LOG.info(f"[parallel] Launching {len(splits_files)} NEAT job(s) (max {options.threads} in parallel)...")
    t1 = time.time()
    output_opts = []
    for splits_file in splits_files:
        reference = splits_file
        current_output_dir = options.temp_dir_path / splits_file.stem
        current_output_dir.mkdir(parents=True, exist_ok=True)
        # Create local filenames based on fasta indexing scheme.
        fq1 = None
        fq2 = None
        bam = None
        vcf = None
        if options.produce_fastq:
            fq1 = current_output_dir / options.fq1.name
            if options.paired_ended:
                # if paired ended, there must be a fq2
                fq2 = current_output_dir / options.fq2.name
        if options.produce_bam:
            bam = current_output_dir / options.bam.name
        if options.produce_vcf:
            vcf = current_output_dir / options.vcf.name

        current_options = options.copy_with_changes(reference, current_output_dir, fq1, fq2, vcf, bam)
        output_opts.append(current_options)

    pool = Pool()
    pool.map(read_simulator_single, output_opts)

    sim_sec = time.time() - t1

    # 3) Stitch

    _LOG.info("[parallel] Stitching per-chunk outputs...")
    t2 = time.time()

    stitch_main(options, output_opts)
    stitch_sec = time.time() - t2

    if not options.cleanup_splits:
        _LOG.info(f"[parallel] File splits preserved for future runs: {options.splits_dir}")

    _LOG.info(
        f"[parallel] Pipeline complete — stitched files under {options.output_dir}\n"
        f"Timings:\n"
        f"  split : {split_sec:6.1f} s\n"
        f"  sim   : {sim_sec:6.1f} s\n"
        f"  stitch: {stitch_sec:6.1f} s\n"
        f"  total : {split_sec + sim_sec + stitch_sec:6.1f} s"
    )