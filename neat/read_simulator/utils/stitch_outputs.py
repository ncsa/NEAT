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

from neat.read_simulator.utils import Options

_LOG = logging.getLogger(__name__)

def concat(files_to_join: List[(Path,)], dest: Path) -> None:
    if not files_to_join:
        # Nothing to do, and no error to throw
        return

    with dest.open("wb") as out_f:
        for f in files_to_join:
            with f.open("rb") as in_f:
                shutil.copyfileobj(in_f, out_f)

def merge_bam(bams: List[Path], dest: Path) -> None:
    if not bams:
        return

    import pysam
    unsorted = dest.with_suffix(".unsorted.bam")
    pysam.merge("--no-PG", "-f", str(unsorted), *map(str, bams))
    pysam.sort("-o", str(dest), str(unsorted))
    unsorted.unlink(missing_ok=True)
    """Make BAM header deterministic: drop @PG, sort @SQ/@RG."""
    hdr_text = pysam.view("-H", str(dest), text=True)

    hd: list[str] = []
    sq: list[str] = []
    rg: list[str] = []
    other: list[str] = []
    for line in hdr_text.splitlines():
        if line.startswith("@PG"):
            continue
        if line.startswith("@HD"):
            hd.append(line)
        elif line.startswith("@SQ"):
            sq.append(line)
        elif line.startswith("@RG"):
            parts = [t for t in line.split("\t") if not t.startswith("CL:")]
            line = "\t".join(parts)
            rg.append(line)
        else:
            other.append(line)

    def sn(x: str) -> str:
        for t in x.split("\t"):
            if t.startswith("SN:"):
                return t[3:]
        return x

    def rgid(x: str) -> str:
        for t in x.split("\t"):
            if t.startswith("ID:"):
                return t[3:]
        return x

    sq.sort(key=sn)
    rg.sort(key=rgid)

    new_hdr = "\n".join([
        *(hd or ["@HD\tVN:1.6\tSO:coordinate"]),
        *sq,
        *rg,
        *other,
    ]) + "\n"

    hdr_tmp = dest.with_suffix(".reheader.sam")
    out_tmp = dest.with_suffix(".reheader.bam")
    hdr_tmp.write_text(new_hdr)

    with open(out_tmp, 'wb') as out_f:
        data = pysam.reheader(str(hdr_tmp), str(dest))
        out_f.write(data)

    out_tmp.replace(dest)
    hdr_tmp.unlink(missing_ok=True)

def merge_vcf(vcfs: List[Path], dest: Path) -> None:
    if not vcfs:
        return
    first, *rest = vcfs
    shutil.copy(first, dest)
    with dest.open("ab") as out_f:
        for vcf in rest:
            with vcf.open("rb") as fh:
                for line in fh:
                    if not line.startswith(b"#"):
                        out_f.write(line)

def main(options: Options, thread_options: list[Options]) -> None:
    fq1_list = []
    fq2_list = []
    vcf_list = []
    bam_list = []
    # Gather all output files from the ops objects
    for local_ops in thread_options:
        if local_ops.fq1:
            fq1_list.append(local_ops.fq1)
        if local_ops.fq2:
            fq2_list.append(local_ops.fq2)
        if local_ops.vcf:
            vcf_list.append(local_ops.vcf)
        if local_ops.bam:
            bam_list.append(local_ops.bam)
    # concatenate all files of each type. An empty list will result in no action
    concat(fq1_list, options.fq1)
    concat(fq2_list, options.fq2)
    merge_bam(bam_list, options.bam)
    merge_vcf(vcf_list, options.vcf)

    # Final success message via logging
    _LOG.info("Stitching complete!")
