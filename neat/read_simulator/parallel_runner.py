"""
Split a reference genome, run NEAT simulations in
parallel, and stitch their outputs back together.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple

from .utils import Options
from .utils.split_inputs import main as split_main
from .utils.stitch_outputs import main as stitch_main
from .single_runner import read_simulator_single

import logging

__all__ = ["main"]

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
    need_split = not options.reuse_splits
    splits_files = []
    if not need_split and list(options.existing_splits_dir.glob("*.y*ml")):
        _LOG.info("[parallel] Recycling existing split FASTA/YAML files – skipping split step")
    else:
        t0 = time.time()
        _LOG.info("[parallel] Splitting reference...")
        splits_files = split_main(options)
        split_sec = time.time() - t0

    # 2) Submit NEAT sims in a deterministic order using the manifest
    _LOG.info(f"[parallel] Launching {len(splits_files)} NEAT job(s) (max {options.threads} in parallel)...")
    t1 = time.time()
    futures = []
    output_opts = []
    with ProcessPoolExecutor(max_workers=options.threads) as pool:
        for splits_file in splits_files:
            current_options = options.deep_copy()
            current_options.reference = splits_file
            current_output_dir = options.temp_dir / splits_file.stem
            current_output_dir.mkdir(parents=True, exist_ok=True)
            # Create local filenames based on fasta indexing scheme.
            if options.produce_fastq:
                current_options.fq1 = make_new_name(options.fq1, splits_file, current_output_dir)
                if options.paired_ended:
                    # if paired ended, there must be an fq2
                    current_options.fq2 = make_new_name(options.fq2, splits_file, current_output_dir)
            if options.produce_bam:
                current_options.bam = make_new_name(options.bam, splits_file, current_output_dir)
            if options.produce_vcf:
                current_options.vcf = make_new_name(options.vcf, splits_file, current_output_dir)
            output_opts.append(current_options)
            current_options.output_dir = current_output_dir
            futures.append(pool.submit(worker, read_simulator_single(current_options)))

        for fut in as_completed(futures):
            try:
                wd = fut.result()
                _LOG.info(f" Finished {wd.name}")
            except Exception as exc:
                _LOG.error(f" Simulation failed: {exc}")
                sys.exit(1)
    sim_sec = time.time() - t1

    # 3) Stitch
    _LOG.info("[parallel] Stitching per-chunk outputs...")
    t2 = time.time()

    stitch_main(options, output_opts)
    stitch_sec = time.time() - t2

    if args.cleanup_splits:
        _LOG.info("[parallel] Cleaning up split FASTA/configs...")
        shutil.rmtree(splits_dir, ignore_errors=True)

    _LOG.info(
        f"[parallel] Pipeline complete — stitched files under {args.final_prefix.parent}\n"
        f"Timings:\n"
        f"  split : {split_sec:6.1f} s\n"
        f"  sim   : {sim_sec:6.1f} s\n"
        f"  stitch: {stitch_sec:6.1f} s\n"
        f"  total : {split_sec + sim_sec + stitch_sec:6.1f} s"
    )

EXTENSIONS = ["fastq", "bam", "vcf", "gz"]

def make_new_name(file_path: Path, ref_path: Path, output_dir) -> Path | None:
    name_parts = file_path.name.split('.')
    for index in range(len(name_parts), 0, -1):
        if name_parts[index] in EXTENSIONS:
            continue
        else:
            last_name = name_parts[index]
            last_name += f"_{ref_path.stem}"
            name_parts[index] = last_name
            return Path(output_dir / ".".join(name_parts))
    # If we reach this point, then our file name was empty
    _LOG.error("Error renaming file, bad file name")
    sys.exit(1)
