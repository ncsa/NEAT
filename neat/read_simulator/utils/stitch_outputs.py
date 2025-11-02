"""
Stitch NEAT split‑run outputs into one dataset.
"""

import argparse
import gzip
import pickle
import re
import shutil
import pysam
import subprocess
import sys
from multiprocessing import Pool, Process
from pathlib import Path
from struct import pack
from typing import Iterable, List, Tuple
import yaml

import logging

__all__ = ["main"]

from Bio import SeqIO, bgzf
from Bio.bgzf import BgzfWriter

from neat.common import open_output, open_input
from neat.read_simulator.utils import Options, OutputFileWriter

_LOG = logging.getLogger(__name__)

def concat(files_to_join: List[Path], ofw: OutputFileWriter, read_num: str) -> None:
    if not files_to_join:
        # Nothing to do, and no error to throw
        return

    read_file = ofw.fq1 if read_num == "1" else ofw.fq2
    dest_file = ofw.files_to_write[read_file]

    for f in files_to_join:
        with bgzf.BgzfReader(f) as in_f:
            shutil.copyfileobj(in_f, dest_file)

def merge_vcfs(vcfs: List[Path], dest: BgzfWriter) -> None:
    if not vcfs:
        return
    for vcf in vcfs:
        with bgzf.BgzfReader(vcf) as fh:
            for line in fh:
                if not line.startswith("#"):
                    dest.write(line)

def merge_bam(bam_files: List[Path], ofw: OutputFileWriter, threads: int) -> None:
    if not bam_files:
        return

    # ofw.files_to_write[ofw.bam].close()
    unsorted = ofw.bam.with_suffix(".unsorted.bam")
    pysam.merge("--no-PG", "-@", str(threads), "-f", str(unsorted), *map(str, bam_files))
    pysam.sort("-@", str(threads), "-o", str(ofw.bam), str(unsorted))
    unsorted.unlink(missing_ok=True)

def main(
        ofw: OutputFileWriter,
        output_files: list[tuple[int, str, dict[str, Path]]],
        threads: int | None = None
) -> None:

    fq1_list = []
    fq2_list = []
    bam_list = []
    vcf_list = []
    # Gather all output files from the ops objects
    for (thread_idx,file_dict) in output_files:
        if file_dict["fq1"]:
            fq1_list.append(file_dict["fq1"])
        if file_dict["fq2"]:
            fq2_list.append(file_dict["fq2"])
        if file_dict["vcf"]:
            vcf_list.append(file_dict["vcf"])
        if file_dict["bam"]:
            bam_list.append(file_dict["bam"])
    # concatenate all files of each type. An empty list will result in no action
    if fq1_list:
        concat(fq1_list, ofw, "1")
    if fq2_list:
        concat(fq2_list, ofw, "2")
    if vcf_list:
        merge_vcfs(vcf_list, ofw.files_to_write[ofw.vcf])
    if bam_list:
        merge_bam(bam_list, ofw, threads)
    # Final success message via logging
    _LOG.info("Stitching complete!")
