from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List
import yaml

__all__ = ["main"]


def print_stderr(msg: str, *, exit_: bool = False):
    """Print to stderr, with an option to exit."""
    print(msg, file=sys.stderr)
    if exit_:
        sys.exit(1)


def gather(paths: Iterable[Path], suffixes: tuple[str, ...]) -> List[Path]:
    """Gather all files in the provided path whose name ends with the given suffix."""
    out: list[Path] = []
    for p in paths:
        if p.is_file() and p.name.endswith(suffixes):
            out.append(p)
        elif p.is_dir():
            for suf in suffixes:
                out.extend(p.rglob(f"*{suf}"))
    return sorted(set(out))


def is_gzipped(file: Path) -> bool:
    return file.suffix in {".gz", ".bgz"}


def concat(files: List[Path], dest: Path, keep_compression: bool):
    if not files:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)

    if keep_compression:
        if not dest.suffix.endswith("gz"):
            dest = dest.with_suffix(dest.suffix + ".gz")
        with dest.open("wb") as out_f:
            for f in files:
                with f.open("rb") as in_f:
                    shutil.copyfileobj(in_f, out_f)
    else:
        # stream out uncompressed lines for speed
        with dest.open("wb") as out_f:
            for f in files:
                opener = gzip.open if is_gzipped(f) else open
                with opener(f, "rb") as in_f:
                    shutil.copyfileobj(in_f, out_f)


def merge_bam(bams: List[Path], dest: Path, samtools: str):
    if not bams:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    samtools = samtools if Path(samtools).exists() else shutil.which(samtools) or samtools
    unsorted = dest.with_suffix(".unsorted.bam")
    cmd_merge = [samtools, "merge", "--no-PG", "-f", str(unsorted), *map(str, bams)]
    subprocess.check_call(cmd_merge)
    subprocess.check_call([samtools, "sort", "-o", str(dest), str(unsorted)])
    subprocess.check_call([samtools, "index", str(dest)])
    unsorted.unlink(missing_ok=True)


def merge_vcf(vcfs: List[Path], dest: Path):
    if not vcfs:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    first, *rest = vcfs
    shutil.copy(first, dest)
    with dest.open("ab") as out_f:
        for vcf in rest:
            with vcf.open("rb") as fh:
                for line in fh:
                    if not line.startswith(b"#"):
                        out_f.write(line)


def parse_args(argv: list[str] | None = None):
    p = argparse.ArgumentParser(description="Stitch NEAT split‑run outputs into one dataset.")
    p.add_argument("-i", "--inputs", nargs="+", type=Path, required=True,
                   help="Input split‑run directories or file prefixes")
    p.add_argument("-o", "--output-prefix", type=Path, required=True,
                   help="Prefix for stitched outputs (no extension)")
    p.add_argument("-c", "--config", type=Path, required=True,
                   help="Original NEAT YAML config to detect paired‑ended, etc.")
    p.add_argument("--samtools", default="samtools", help="Path to samtools binary")
    return p.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)

    # read config to know paired vs. single
    cfg = yaml.safe_load(args.config.read_text())
    paired = bool(cfg.get("paired_ended", False))

    # gather files
    fq_r1 = gather(args.inputs, ("_r1_paired.fq", "_r1_paired.fq.gz", "_r1_paired.fq.bgz"))
    fq_r2 = gather(args.inputs, ("_r2_paired.fq", "_r2_paired.fq.gz", "_r2_paired.fq.bgz")) if paired else []
    fq_single = gather(args.inputs, ("_r1_single.fq", "_r1_single.fq.gz", "_r1_single.fq.bgz"))
    bams = gather(args.inputs, ("_golden.bam",))
    vcfs = gather(args.inputs, ("_golden.vcf",))

    if not any((fq_r1, fq_r2, fq_single, bams, vcfs)):
        print_stderr("No NEAT output files found under supplied inputs", exit_=True)

    # check compression consistency for R1 and R2
    for lane in (fq_r1, fq_r2):
        if lane:
            gz = [is_gzipped(f) for f in lane]
            if any(gz) and not all(gz):
                print_stderr("Mixing compressed & uncompressed FASTQs; aborting stitching", exit_=True)
    keep_comp = fq_r1 and all(is_gzipped(f) for f in fq_r1)

    prefix = args.output_prefix

    # concatenate FASTQs
    concat(fq_r1, prefix.with_name(prefix.name + "_read1.fq"), keep_comp)
    if paired:
        concat(fq_r2, prefix.with_name(prefix.name + "_read2.fq"), keep_comp)
    concat(fq_single, prefix.with_name(prefix.name + "_readSingle.fq"), keep_comp)

    merge_bam(bams, prefix.with_name(prefix.name + "_golden.bam"), args.samtools)
    merge_vcf(vcfs, prefix.with_name(prefix.name + "_golden.vcf"))

    print_stderr("Stitching complete!")


if __name__ == "__main__":
    main()
