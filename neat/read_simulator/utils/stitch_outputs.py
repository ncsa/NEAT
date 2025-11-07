"""
Stitch NEAT split‑run outputs into one dataset.
"""

import shutil
import pysam
from pathlib import Path
from typing import List

import logging

__all__ = ["main"]

from Bio import bgzf
import gzip
from neat.read_simulator.utils import OutputFileWriter

_LOG = logging.getLogger(__name__)

def concat(files_to_join: List[Path], file_to_write: Path) -> None:
    dest_file = bgzf.BgzfWriter(file_to_write)
    for f in files_to_join:
        with gzip.open(f, 'r') as in_f:
            for line in in_f:
                dest_file.write(line)
    dest_file.close()

def merge_vcfs(vcfs: List[Path], file_to_write: Path) -> None:
    dest = bgzf.BgzfWriter(file_to_write)
    for vcf in vcfs:
        with gzip.open(vcf, 'rt') as fh:
            for line in fh:
                if not line.startswith("#"):
                    dest.write(line)
    dest.close()

def merge_bam(bam_files: List[Path], ofw: OutputFileWriter, threads: int) -> None:
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
        concat(fq1_list, ofw.fq1)
    if fq2_list:
        concat(fq2_list, ofw.fq2)
    if vcf_list:
        merge_vcfs(vcf_list, ofw.vcf)
    if bam_list:
        merge_bam(bam_list, ofw, threads)
    # Final success message via logging
    _LOG.info("Stitching complete!")
