"""
Split a FASTA and NEAT config into per‑contig or fixed‑size chunks ready for simulation.
"""

import shutil
import sys
import resource
from pathlib import Path
from textwrap import wrap
from typing import Iterator

import logging

__all__ = ["main"]

from Bio import bgzf, SeqIO
from Bio.Seq import Seq

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

def write_fasta(contig: str, seq: str, out_fa: Path, width: int = 80) -> None:
    """Write a sequence record to a FASTA file, wrapping lines at a given width."""
    with bgzf.BgzfWriter(out_fa) as fh:
        fh.write(f">{contig}\n")
        # Ensure seq is a string for bgzf.BgzfWriter.write()
        seq_str = str(seq)
        for i in range(0, len(seq_str), width):
            fh.write(seq_str[i:i+width] + "\n")

# Chunk generators
def chunk_record(seq: str, chunk_len: int, overlap: int) -> Iterator[tuple[int, str]]:
    """
    Yield overlapping slices of a record.
    """
    seq_len = len(seq)
    start: int = 0
    while start < seq_len:
        end = min(start + chunk_len, seq_len)
        sub_seq = seq[start:end]
        yield start, sub_seq
        if end == seq_len:
            break
        start = end - overlap

def main(options: Options) -> tuple[dict, int, dict]:
    """Perform the splitting of a FASTA and NEAT configuration."""

    overlap = int(options.read_len) * 2

    idx = 1
    pad = 10  # zero-pad width for global indices

    written = 0
    # We'll keep track of chunks by contig, to help us out later
    split_fasta_dict: dict[str, dict[tuple[int, int], Path]] = {}
    reference_keys_with_lens: dict[str, int] = {}

    for seq_record in SeqIO.parse(str(options.reference), "fasta"):
        contig = seq_record.id
        seq_str = str(seq_record.seq).upper()
        reference_keys_with_lens[contig] = len(seq_str)
        split_fasta_dict[contig] = {}
        
        if options.parallel_mode == "contig":
            stem = f"{idx:0{pad}d}__{contig}"
            fa = options.splits_dir / f"{stem}.fa.gz"
            write_fasta(contig, seq_str, fa)
            split_fasta_dict[contig][(0, len(seq_str))] = fa
            idx += 1
            written += 1
        else:
            for start, subseq in chunk_record(seq_str, options.parallel_block_size, overlap):
                stem = f"{idx:0{pad}d}__{contig}"
                fa = options.splits_dir / f"{stem}.fa.gz"
                write_fasta(contig, subseq, fa)
                split_fasta_dict[contig][(start, start+len(subseq))] = fa
                idx += 1
                written += 1

    # Report success via the logger instead of printing to stderr
    _LOG.info(f"Generated {written} FASTAs in {options.splits_dir}")
    return split_fasta_dict, written, reference_keys_with_lens
