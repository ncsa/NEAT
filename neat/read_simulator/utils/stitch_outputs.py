"""
Stitch NEAT split‑run outputs into one dataset.
"""
import resource
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

def concat(files_to_join: List[Path], dest_file: gzip.GzipFile) -> None:
    for f in files_to_join:
        with gzip.open(f, 'rt') as in_f:
            shutil.copyfileobj(in_f, dest_file)

def merge_vcfs(vcfs: List[Path], ofw: OutputFileWriter) -> None:
    dest = ofw.files_to_write[ofw.vcf]
    for vcf in vcfs:
        with gzip.open(vcf, 'rt') as fh:
            for line in fh:
                if not line.startswith("#"):
                    dest.write(line)

def merge_bam(bam_files: List[Path], ofw: OutputFileWriter, threads: int):
    merged_file = ofw.tmp_dir / "temp_merged.bam"
    intermediate_files = []
    # Note 1000 is abritrary. May need to be a user parameter/adjustable/a function
    for i in range(0, len(bam_files), 500):
        temp_file = str(ofw.tmp_dir / f"temp_merged_{i}.bam")
        pysam.merge("--no-PG", "-f", temp_file, *map(str, bam_files[i:i+500]))
        intermediate_files.append(temp_file)
    pysam.merge("--no-PG", "-f", str(merged_file), *intermediate_files)
    pysam.sort("-@", str(threads), "-m", "1G", "-o", str(ofw.bam), str(merged_file))
    merged_file.unlink(missing_ok=True)

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
        concat(fq1_list, ofw.files_to_write[ofw.fq1])
    if fq2_list:
        concat(fq2_list, ofw.files_to_write[ofw.fq2])
    if vcf_list:
        merge_vcfs(vcf_list, ofw)
    if bam_list:
        merge_bam(bam_list, ofw, threads)
    # Final success message via logging
    _LOG.info("Stitching complete!")
