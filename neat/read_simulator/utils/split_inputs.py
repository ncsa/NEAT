"""
Description.
"""

import argparse
import shutil
import sys
import yaml
import time
from pathlib import Path
from textwrap import wrap
from typing import Iterator, List

__all__ = ["main"]


class SimpleRecord:
    """Basic stand-in for Bio.SeqRecord.SeqRecord."""
    __slots__ = ("id", "seq")

    def __init__(self, id_: str, seq: str):
        self.id: str = id_
        self.seq: str = seq

    def __len__(self) -> int:
        return len(self.seq)


# Utility helpers
def print_stderr(msg: str, *, exit_: bool = False):
    print(msg, file=sys.stderr)
    if exit_:
        sys.exit(1)


def disk_bytes_free(path: Path) -> int:
    return shutil.disk_usage(path).free


# FASTA helpers
def parse_fasta(path: Path) -> Iterator[SimpleRecord]:
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
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w") as fh:
        fh.write(f">{rec.id}\n")
        for chunk in wrap(rec.seq, width):
            fh.write(chunk + "\n")


# Chunk generators
def chunk_record(record: SimpleRecord, chunk_len: int, overlap: int):
    """Yield overlapping chunk_len slices."""
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


def write_config(template_cfg: dict, out_fa: Path, out_yml: Path):
    cfg_copy = template_cfg.copy()
    cfg_copy["reference"] = str(out_fa.resolve())
    if "name" in cfg_copy:
        cfg_copy["name"] = out_fa.stem
    out_yml.write_text(yaml.safe_dump(cfg_copy, sort_keys=False))


# Argument parsing
def parse_args(argv: List[str] | None = None):
    p = argparse.ArgumentParser(
        description="Split a FASTA and NEAT config into per-contig or fixed-size chunks ready for neat read-simulator.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("fasta", type=Path, help="Input multi-contig FASTA")
    p.add_argument("config", type=Path, help="Original NEAT YAML/YML config")
    p.add_argument("--outdir", type=Path, required=True, help="Directory to write split FASTAs + configs")
    p.add_argument("--by", choices=["contig", "size"], default="contig",
                   help="Split whole contigs (default) or fixed-size chunks")
    p.add_argument("--size", type=int, default=1_000_000,
                   help="Target chunk size when --by size is chosen")
    return p.parse_args(argv)


def main(argv: List[str] | None = None):
    args = parse_args(argv)

    template_cfg = yaml.safe_load(args.config.read_text())
    read_len = int(template_cfg.get("read_len", 150))
    overlap = read_len * 2

    approx_out_bytes = int(args.fasta.stat().st_size * 1.1)
    if disk_bytes_free(args.outdir) < approx_out_bytes:
        print_stderr(
            f"Not enough free space in {args.outdir} (need about {approx_out_bytes/1e9:.2f} GB)",
            exit_=True,
        )

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    manifest_entries: list[dict] = []
    idx = 1
    pad = 6  # zero-pad width for global indices

    mode = args.by
    written = 0

    for rec in parse_fasta(args.fasta):
        if mode == "contig":
            stem = f"{idx:0{pad}d}__{rec.id}"
            fa = outdir / f"{stem}.fa"
            yml = outdir / f"{stem}.yaml"
            write_fasta(rec, fa)
            write_config(template_cfg, fa, yml)
            manifest_entries.append({
                "index": idx,
                "contig": rec.id,
                "chunk": 1,
                "stem": stem,
                "fasta": str(fa),
                "yaml": str(yml),
            })
            idx += 1
            written += 1
        else:
            for subrec, chunk_id in chunk_record(rec, args.size, overlap):
                stem = f"{idx:0{pad}d}__{subrec.id}"
                fa = outdir / f"{stem}.fa"
                yml = outdir / f"{stem}.yaml"
                write_fasta(subrec, fa)
                write_config(template_cfg, fa, yml)
                manifest_entries.append({
                    "index": idx,
                    "contig": rec.id,
                    "chunk": chunk_id,
                    "stem": stem,
                    "fasta": str(fa),
                    "yaml": str(yml),
                })
                idx += 1
                written += 1

    manifest = {
        "version": 1,
        "created": int(time.time()),
        "total": written,
        "entries": manifest_entries,
    }
    (outdir / "manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))

    print_stderr(f"Generated {written} FASTA/YAML pair(s) in {outdir}")


if __name__ == "__main__":
    main()
