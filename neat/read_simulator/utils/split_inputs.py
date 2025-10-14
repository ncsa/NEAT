"""
Split a FASTA and NEAT config into per‑contig or fixed‑size chunks ready for simulation.
"""

import shutil
import sys
import gzip
from pathlib import Path
from textwrap import wrap
from typing import Iterator

import logging

__all__ = ["main"]

from Bio import bgzf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from neat.read_simulator.utils import Options

_LOG = logging.getLogger(__name__)


# Utility helpers
def print_stderr(msg: str, *, exit_: bool = False) -> None:
    """Log an error message and optionally exit with status 1."""
    _LOG.error(msg)
    if exit_:
        sys.exit(1)

def disk_bytes_free(path: Path) -> int:
    return shutil.disk_usage(path).free

def write_fasta(contig: str, rec: Seq, out_fa: Path, width: int = 80) -> None:
    """Write a sequence record to a FASTA file, wrapping lines at a given width."""
    with bgzf.BgzfWriter(out_fa) as fh:
        fh.write(f">{contig}\n")
        for chunk in wrap(str(rec), width):
            fh.write(chunk + "\n")

# Chunk generators
def chunk_record(record: Seq, chunk_len: int, overlap: int) -> Iterator[tuple[int, Seq]]:
    """
    Yield overlapping slices of a record.
    """
    seq_len = len(record)
    start: int = 0
    chunkid = 0
    while start < seq_len:
        end = min(start + chunk_len, seq_len)
        sub_seq = record[start:end]
        yield start, sub_seq
        if end == seq_len:
            break
        start = end - overlap
        chunkid += 1

def main(options: Options, reference_index: dict) -> tuple[dict, int]:
    """Perform the splitting of a FASTA and NEAT configuration."""

    overlap = int(options.read_len) * 2

    approx_out_bytes = int(options.reference.stat().st_size * 1.1)
    if disk_bytes_free(options.output_dir) < approx_out_bytes:
        print_stderr(
            f"Not enough free space in {options.output_dir} (need about {approx_out_bytes/1e9:.2f} GB)",
            exit_=True,
        )

    idx = 1
    pad = 10  # zero-pad width for global indices

    written = 0
    # We'll keep track of chunks by contig, to help us out later
    split_fasta_dict: dict[str, dict[tuple[int, int], Path]] = {key: {} for key in reference_index.keys()}
    for contig, seq_record in reference_index.items():
        if options.mode == "contig":
            stem = f"{idx:0{pad}d}__{contig}"
            fa = options.splits_dir / f"{stem}.fa.gz"
            write_fasta(contig, seq_record.seq, fa)
            split_fasta_dict[contig][(0, len(seq_record))] = fa
            idx += 1
            written += 1
        else:
            for start, subseq in chunk_record(seq_record.seq, options.size, overlap):
                stem = f"{idx:0{pad}d}__{contig}"
                fa = options.splits_dir / f"{stem}.fa.gz"
                write_fasta(contig, subseq, fa)
                split_fasta_dict[contig][(start, start+len(subseq))] = fa
                idx += 1
                written += 1

    # Report success via the logger instead of printing to stderr
    _LOG.info(f"Generated {written} FASTAs in {options.splits_dir}")
    return split_fasta_dict, written
