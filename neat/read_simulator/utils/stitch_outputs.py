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
    if not files_to_join:
        # Nothing to do, and no error to throw
        _LOG.warning(f"Concat called but there are no files to join: {files_to_join}" )
        return

    for f in files_to_join:
        with gzip.open(f, 'rt') as in_f:
            shutil.copyfileobj(in_f, dest_file)

def merge_vcfs(vcfs: List[Path], ofw: OutputFileWriter) -> None:
    dest = ofw.files_to_write[ofw.vcf]
    seen: set[str] = set()
    n_duplicates = 0
    for vcf in vcfs:
        with gzip.open(vcf, 'rt') as fh:
            for line in fh:
                if not line.startswith("#"):
                    normalized = line.rstrip("\r\n")
                    if normalized in seen:
                        n_duplicates += 1
                        continue
                    seen.add(normalized)
                    dest.write(line)
    if n_duplicates:
        _LOG.warning(f"merge_vcfs: removed {n_duplicates} duplicate VCF line(s) during merge.")

def merge_bam(bam_files: List[Path], ofw: OutputFileWriter, threads: int):
    # Per-worker BAMs are coordinate-sorted within themselves, and each chunk owns a
    # non-overlapping reference range for r1 placement (see cover_dataset's
    # responsibility_length), so the BAMs concatenate into a globally coordinate-sorted
    # output without any sort or k-way merge. pysam.cat does a raw BGZF concatenation —
    # no decompression / re-encode — and is typically 10-30x faster than pysam.merge on
    # large outputs. At supercomputer scale this is the difference between the stitch
    # step being I/O-bound (fast) and BGZF-bound (slow).
    pysam.cat("-o", str(ofw.bam), *map(str, bam_files))
    # The .bai index is produced by runner.py after stitching.

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
