from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Iterator, List

import yaml
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

__all__ = ["main"]


# Utility helpers
def print_stderr(msg: str, *, exit_: bool = False):
    """Print with an option to exit."""
    print(msg, file=sys.stderr)
    if exit_:
        sys.exit(1)


def disk_bytes_free(path: Path) -> int:
    """Return free bytes on the filesystem where the path lives."""
    return shutil.disk_usage(path).free


# Chunk generators
def chunk_record(
    record: SeqRecord,
    chunk_len: int,
    overlap: int,
) -> Iterator[SeqRecord]:
    """
    Yield chunk_len-bp overlapping slices of record.
    Each chunk overlaps the previous to allow downstream stitching of reads.
    """
    seq_len = len(record)
    start = 0
    chunk_id = 1
    while start < seq_len:
        end = min(start + chunk_len, seq_len)
        chunk_seq = record.seq[start:end]
        new_id = f"{record.id}_chunk{chunk_id}"
        yield SeqRecord(chunk_seq, id=new_id, description="")
        if end == seq_len:
            break
        start = end - overlap
        chunk_id += 1


# File writers
def write_fasta(rec: SeqRecord, out_fa: Path):
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w") as fh:
        SeqIO.write(rec, fh, "fasta")


def write_config(template_cfg: dict, out_fa: Path, out_yml: Path):
    cfg_copy = template_cfg.copy()
    cfg_copy["reference"] = str(out_fa)
    # Optional: set "name" to the FASTA stem if present in template
    if "name" in cfg_copy:
        cfg_copy["name"] = out_fa.stem
    out_yml.write_text(yaml.safe_dump(cfg_copy, sort_keys=False))


# Argument parsing
def parse_args(argv: List[str] | None = None):
    p = argparse.ArgumentParser(
        description="Split a FASTA and NEAT config into per-contig or "
        "fixed-size chunks ready for neat read‑simulator."
    )
    p.add_argument("fasta", type=Path, help="Input multi‑contig FASTA")
    p.add_argument("config", type=Path, help="Original NEAT YAML config")
    p.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help="Directory to write split FASTAs + configs",
    )
    p.add_argument(
        "--by",
        choices=["contig", "size"],
        default="contig",
        help="Mode for splitting is either whole contigs (default) or size‑based chunks",
    )
    p.add_argument(
        "--size",
        type=int,
        default=1_000_000,
        metavar="N",
        help="Target chunk size",
    )
    return p.parse_args(argv)


# Main
def main(argv: List[str] | None = None):
    args = parse_args(argv)

    # Load template config
    template_cfg = yaml.safe_load(args.config.read_text())
    read_len = int(template_cfg.get("read_len", 150))
    overlap = read_len * 2

    # Disk‑space sanity check
    approx_out_bytes = int(args.fasta.stat().st_size * 1.1)
    if disk_bytes_free(args.outdir) < approx_out_bytes:
        print_stderr(
            f"Not enough free space in {args.outdir} "
            f"(need ≈ {approx_out_bytes/1e9:.2f} GB)",
            exit_=True,
        )

    written = 0
    mode = args.by
    outdir = args.outdir

    # Iterate FASTA records
    for rec in SeqIO.parse(args.fasta, "fasta"):
        if mode == "contig":
            fasta_path = outdir / f"{rec.id}.fa"
            yaml_path = fasta_path.with_suffix(".yaml")
            write_fasta(rec, fasta_path)
            write_config(template_cfg, fasta_path, yaml_path)
            written += 1
        else:  # fixed size mode
            for subrec in chunk_record(rec, args.size, overlap):
                fasta_path = outdir / f"{subrec.id}.fa"
                yaml_path = fasta_path.with_suffix(".yaml")
                write_fasta(subrec, fasta_path)
                write_config(template_cfg, fasta_path, yaml_path)
                written += 1

    print_stderr(f"Generated {written} FASTA/YAML pair(s) in {outdir}")


if __name__ == "__main__":
    main()
