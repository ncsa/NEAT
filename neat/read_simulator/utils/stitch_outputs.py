"""
Stitch NEAT split‑run outputs into one dataset.
"""

import argparse
import gzip
import pickle
import re
import shutil
import subprocess
import sys
from pathlib import Path
from struct import pack
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

def merge_bam(reads_pickles: List[Path], ofw: OutputFileWriter, contig_dict: dict, read_length: int) -> None:
    if not reads_pickles:
        return
    bam_out = ofw.files_to_write[ofw.bam]
    bam_out.write("BAM\1")
    header = "@HD\tVN:1.4\tSO:coordinate\n"
    for item in ofw.bam_header:
        header += f'@SQ\tSN:{item}\tLN:{str(ofw.bam_header[item])}\n'
    header += "@RG\tID:NEAT\tSM:NEAT\tLB:NEAT\tPL:NEAT\n"
    header_bytes = len(header)
    num_refs = len(ofw.bam_header)
    bam_out.write(pack('<i', header_bytes))
    bam_out.write(header)
    bam_out.write(pack('<i', num_refs))

    for item in ofw.bam_header:
        name_length = len(item) + 1
        bam_out.write(pack('<i', name_length))
        bam_out.write(f'{item}\0')
        bam_out.write(pack('<i', ofw.bam_header[item]))

    for file in reads_pickles:
        contig_reads_data = pickle.load(gzip.open(file))
        for read_data in contig_reads_data:
            read1 = read_data[0]
            read2 = read_data[1]
            if read1:
                ofw.write_bam_record(
                    read1,
                    contig_dict[read1.reference_id],
                    bam_out,
                    read_length
                )
            if read2:
                ofw.write_bam_record(
                    read2,
                    contig_dict[read2.reference_id],
                    bam_out,
                    read_length
                )
    bam_out.close()

def main(
        ofw: OutputFileWriter,
        output_files: list[tuple[int, str, dict[str, Path]]],
        contig_dict: dict | None = None,
        read_length: int | None = None,
) -> None:

    fq1_list = []
    fq2_list = []
    reads_pickles = []
    # Gather all output files from the ops objects
    for (thread_idx,file_dict) in output_files:
        if file_dict["fq1"]:
            fq1_list.append(file_dict["fq1"])
        if file_dict["fq2"]:
            fq2_list.append(file_dict["fq2"])
        if file_dict["reads"]:
            reads_pickles.append(file_dict["reads"])
    # concatenate all files of each type. An empty list will result in no action
    concat(fq1_list, ofw, ofw.fq1)
    concat(fq2_list, ofw, ofw.fq2)
    merge_bam(reads_pickles, ofw, contig_dict, read_length)

    # Final success message via logging
    _LOG.info("Stitching complete!")
