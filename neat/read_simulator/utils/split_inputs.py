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

from neat.read_simulator.utils import Options

_LOG = logging.getLogger(__name__)


class SimpleRecord:
    """
    Basic stand‑in for Bio.SeqRecord.SeqRecord.
    This lightweight class stores an identifier and a sequence and
    implements ``__len__`` so that sequence length can be computed.
    """

    __slots__ = ("id", "seq")

    def __init__(self, id_: str, seq: str):
        self.id: str = id_
        self.seq: str = seq

    def __len__(self) -> int:
        return len(self.seq)


# Utility helpers
def print_stderr(msg: str, *, exit_: bool = False) -> None:
    """Log an error message and optionally exit with status 1."""
    _LOG.error(msg)
    if exit_:
        sys.exit(1)

def disk_bytes_free(path: Path) -> int:
    return shutil.disk_usage(path).free

# FASTA helpers
def parse_fasta(path: Path) -> Iterator[SimpleRecord]:
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
                    yield SimpleRecord(header, "".join(seq_chunks))
                header = line[1:].split(maxsplit=1)[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header:
            yield SimpleRecord(header, "".join(seq_chunks))


def write_fasta(rec: SimpleRecord, out_fa: Path, width: int = 60) -> None:
    """Write a sequence record to a FASTA file, wrapping lines at a given width."""
    with out_fa.open("w") as fh:
        fh.write(f">{rec.id}\n")
        for chunk in wrap(rec.seq, width):
            fh.write(chunk + "\n")

# Chunk generators
def chunk_record(record: SimpleRecord, chunk_len: int, overlap: int) -> Iterator[Tuple[SimpleRecord, int]]:
    """
    Yield overlapping slices of a record.
    """
    seq_len = len(record)
    start = 0
    chunk_id = 1
    while start < seq_len:
        end = min(start + chunk_len, seq_len)
        sub_seq = record.seq[start:end]
        new_id = f"{record.id}_chunk{chunk_id}"
        yield SimpleRecord(new_id, sub_seq), chunk_id
        if end == seq_len:
            break
        start = end - overlap
        chunk_id += 1

def write_config(template_cfg: dict, out_fa: Path, out_yml: Path) -> None:
    """Write a per‑fragment configuration derived from the template config."""
    cfg_copy = template_cfg.copy()
    cfg_copy["reference"] = str(out_fa.resolve())
    if "name" in cfg_copy:
        cfg_copy["name"] = out_fa.stem
    out_yml.write_text(yaml.safe_dump(cfg_copy, sort_keys=False))

def main(options: Options) -> list:
    """Perform the splitting of a FASTA and NEAT configuration."""
    overlap = options.read_len * 2

    approx_out_bytes = int(options.reference.stat().st_size * 1.1)
    if disk_bytes_free(options.output_dir) < approx_out_bytes:
        print_stderr(
            f"Not enough free space in {options.output_dir} (need about {approx_out_bytes/1e9:.2f} GB)",
            exit_=True,
        )

    idx = 1
    pad = 6  # zero-pad width for global indices

    written = 0
    fasta_files = []
    for rec in parse_fasta(options.reference):
        if options.mode == "contig":
            stem = f"{idx:0{pad}d}__{rec.id}"
            fa = options.splits_dir / f"{stem}.fa"
            write_fasta(rec, fa)
            fasta_files.append(fa)
            idx += 1
            written += 1
        else:
            for subrec, chunk_id in chunk_record(rec, options.size, overlap):
                stem = f"{idx:0{pad}d}__{subrec.id}"
                fa = options.splits_dir / f"{stem}.fa"
                write_fasta(subrec, fa)
                fasta_files.append(fa)
                idx += 1
                written += 1

    # Report success via the logger instead of printing to stderr
    _LOG.info(f"Generated {written} FASTAs in {options.splits_dir}")
    files_written_string = "\n\t".join([str(x) for x in fasta_files])
    _LOG.debug(f'Splits files written: \n\t{files_written_string}')
    return fasta_files
