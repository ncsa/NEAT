"""
Stitch NEAT split‑run outputs into one dataset.
"""

import argparse
import gzip
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Tuple
import yaml

import logging

__all__ = ["main"]

from Bio import SeqIO, bgzf

from neat.common import open_output, open_input
from neat.read_simulator.utils import Options, OutputFileWriter

_LOG = logging.getLogger(__name__)

def concat(files_to_join: List[Path], ofw: OutputFileWriter, file: Path) -> None:
    if not files_to_join:
        # Nothing to do, and no error to throw
        return

    out_handle = ofw.files_to_write[file]
    for f in files_to_join:
        with bgzf.BgzfReader(f) as in_f:
            shutil.copyfileobj(in_f, out_handle)

def merge_bam(bams: List[Path], bam_handle) -> None:
    if not bams:
        return

    for bam in bams:
        with bgzf.BgzfReader(bam) as in_bam:
            for line in in_bam:
                if line.startswith('@'):
                    continue
                else:
                    bam_handle.write(line)

def main(ofw: OutputFileWriter, output_files: list[tuple[int, str, dict[str, Path]]]) -> None:

    fq1_list = []
    fq2_list = []
    bam_list = []
    # Gather all output files from the ops objects
    for (thread_idx, contig_name, file_dict) in output_files:
        if file_dict["fq1"]:
            fq1_list.append(file_dict["fq1"])
        if file_dict["fq2"]:
            fq2_list.append(file_dict["fq2"])
        if file_dict["bam"]:
            bam_list.append(file_dict["bam"])
    # concatenate all files of each type. An empty list will result in no action
    concat(fq1_list, ofw, ofw.fq1)
    concat(fq2_list, ofw, ofw.fq2)
    merge_bam(bam_list, ofw.files_to_write[ofw.bam])

    # Final success message via logging
    _LOG.info("Stitching complete!")
