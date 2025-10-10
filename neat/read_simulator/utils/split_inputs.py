"""
Split a FASTA and NEAT config into per‑contig or fixed‑size chunks ready for simulation.
"""

import shutil
import sys
import yaml
from pathlib import Path
from textwrap import wrap
from typing import Iterator, Tuple

import logging

__all__ = ["main"]

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

# FASTA helpers
def parse_fasta(path: Path) -> Iterator[SeqRecord]:
    """Parse a FASTA file into SimpleRecord objects."""
    with path.open() as fh:
        header: str | None = None
        seq_chunks: list[str] = []
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield SeqRecord(Seq("".join(seq_chunks)), header)
                header = line[1:].split(maxsplit=1)[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header:
            yield SeqRecord(Seq("".join(seq_chunks)), header)


def write_fasta(contig: str, rec: Seq, out_fa: Path, width: int = 60) -> None:
    """Write a sequence record to a FASTA file, wrapping lines at a given width."""
    with out_fa.open("w") as fh:
        fh.write(f">{contig}\n")
        for chunk in wrap(str(rec), width):
            fh.write(chunk + "\n")

# Chunk generators
def chunk_record(record: SeqRecord, chunk_len: int, overlap: int) -> Iterator[Seq]:
    """
    Yield overlapping slices of a record.
    """
    seq_len = len(record)
    start = 0
    chunkid = 0
    while start < seq_len:
        end = min(start + chunk_len, seq_len)
        sub_seq = record.seq[start:end]
        yield sub_seq
        if end == seq_len:
            break
        start = end - overlap
        chunkid += 1

def main(options: Options, frag_mean: float, reference_index: dict) -> dict:
    """Perform the splitting of a FASTA and NEAT configuration."""

    overlap = int(frag_mean) * 2

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
    split_fasta_dict = {key: [] for key in reference_index.keys()}
    for contig, seq_record in reference_index.items():
        if options.mode == "contig":
            stem = f"{idx:0{pad}d}__{contig}"
            fa = options.splits_dir / f"{stem}.fa"
            write_fasta(contig, seq_record.seq, fa)
            split_fasta_dict[contig].append(fa)
            idx += 1
            written += 1
        else:
            for subseq in chunk_record(seq_record.seq, options.size, overlap):
                stem = f"{idx:0{pad}d}__{contig}"
                fa = options.splits_dir / f"{stem}.fa"
                write_fasta(contig, subseq, fa)
                split_fasta_dict[contig].append(fa)
                idx += 1
                written += 1

    # Report success via the logger instead of printing to stderr
    _LOG.info(f"Generated {written} FASTAs in {options.splits_dir}")
    return split_fasta_dict
